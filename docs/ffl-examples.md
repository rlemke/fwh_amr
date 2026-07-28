# FFL Examples — `amr`

Every numbered scenario is a **complete, compilable FFL file**. Copy one into
`my.ffl` and run it against this domain:

```bash
fw ffl run --primary my.ffl \
  --library ~/fw_handlers/fwh_amr/src/amr/ffl/amr.ffl \
  --workflow my.amr.<WorkflowName>
```

A runner serving the `amr` namespace must be up (`fw runner start --domain amr`).
Every block below is compile-checked against `src/amr/ffl/amr.ffl`.

New to the language? Start with the
[FFL grammar](https://github.com/rlemke/facetwork/blob/main/docs/reference/language/grammar.md)
and the [canonical examples](https://github.com/rlemke/facetwork/tree/main/examples/canonical).

---

## The facets at a glance

| Declaration | Signature | Does |
|---|---|---|
| `amr.sources.ListOrganisms` | `(quick: Boolean = false) => (organisms: [String], count: Long)` | The curated pathogen groups — the fan-out unit |
| `amr.sources.DownloadOrganism` | `(group: String, force: Boolean = false, max_rows: Long = 0) => (group, isolate_count, gene_count)` | Ingest ONE NCBI Pathogen Detection group → cached country × year aggregate |
| `amr.maps.BuildResistanceMap` | `(dependency_signal: Long = 0) => (html_path, country_count, organism_count, isolate_count, year_min, year_max)` | Fan-in → year-slider + drug-class world choropleth |
| `amr.maps.BuildGeneSpeciesTable` | `(dependency_signal: Long = 0) => (html_path, row_count, organism_count)` | Fan-in → ranked gene × species table |
| `amr.workflows.BuildOrganismAtlas` | `(group: String, force: Boolean = false)` | One organism, end to end |
| `amr.workflows.BuildAmrFanout` | `(quick, force, max_rows) => (built: [String])` | One ingest per group, in parallel across the fleet |
| `amr.workflows.BuildAmrMaps` | `(force: Boolean = false)` | Fan-in: render map + table over the whole cache |

---

## 1. Run what ships — no FFL to write

```bash
fw ffl seed --include amr

# one pathogen group, end to end (the quick demo)
fw ffl run --primary ~/fw_handlers/fwh_amr/src/amr/ffl/amr.ffl \
  --workflow amr.workflows.BuildOrganismAtlas \
  --inputs '{"group": "Vibrio_cholerae"}'

# the real thing: fan out every group, then render
fw ffl run --primary ~/fw_handlers/fwh_amr/src/amr/ffl/amr.ffl \
  --workflow amr.workflows.BuildAmrFanout
fw ffl run --primary ~/fw_handlers/fwh_amr/src/amr/ffl/amr.ffl \
  --workflow amr.workflows.BuildAmrMaps
```

Write FFL when you want a different *shape* of run — a subset of organisms, your
own error handling, or composition with another domain.

## 2. The smallest workflow you can write

One ingest, one render. Every FFL workflow needs a `namespace`, a `use` per
namespace it calls, and a `yield` back to itself.

```ffl
namespace my.amr {

    use amr.sources
    use amr.maps

    /** Ingest one pathogen group and render the world map over the cache. */
    workflow OneOrganism(group: String = "Vibrio_cholerae") => (html_path: String, isolates: Long) andThen {

        dl = amr.sources.DownloadOrganism(group = $.group)

        map = amr.maps.BuildResistanceMap(dependency_signal = dl.isolate_count)

        yield OneOrganism(html_path = map.html_path, isolates = map.isolate_count)
    }
}
```

Rules visible above: `=>` sits on the **same line** as the closing `)`;
references are always `step.field`; `$.group` reads the container's (here the
workflow's) parameters.

## 3. Fan out over a list — the distributed pattern

`andThen foreach v in <list>` turns one step into N runtime steps that the fleet
claims in parallel. Inside the loop body the loop variable is `$.v`, and `$$`
walks out to the enclosing workflow's parameters.

```ffl
namespace my.amr {

    use amr.sources

    /** One ingest task per curated group — claimed in parallel across the fleet. */
    workflow IngestEverything(quick: Boolean = false) => (built: [String]) andThen {

        groups = amr.sources.ListOrganisms(quick = $.quick) andThen foreach g in $.organisms {

            // $.g   = this iteration's group
            // $$.quick = the workflow parameter, one level out
            dl = amr.sources.DownloadOrganism(group = $.g, force = false)

            yield IngestEverything(built = [dl.group])
        }
    }
}
```

Wall clock is the slowest organism, not the sum — this is the whole point of the
runtime. Watch the fan-out live on the dashboard's execution graph.

## 4. Fan out over *your own* list

The list doesn't have to come from a facet. A `Json` parameter works, so you can
drive the fan-out from `--inputs`:

```ffl
namespace my.amr {

    use amr.sources
    use amr.maps

    /** Ingest exactly the groups you name, then render. */
    workflow IngestChosen(groups: Json) => (built: [String]) andThen foreach g in $.groups {

        dl = amr.sources.DownloadOrganism(group = $.g, max_rows = 50000)

        yield IngestChosen(built = [dl.group])
    }
}
```

```bash
fw ffl run --primary my.ffl --library …/amr.ffl \
  --workflow my.amr.IngestChosen \
  --inputs '{"groups": ["Salmonella_enterica", "Escherichia_coli_Shigella"]}'
```

## 5. Fan-in — sequencing renders after ingest

`BuildResistanceMap` reads the *cache* the ingests wrote, so it needs no value
from them. Steps with no reference between them may run in any order, so make the
dependency explicit by passing any upstream field as `dependency_signal`:

```ffl
namespace my.amr {

    use amr.maps

    /** Render both artifacts, table strictly after the map. */
    workflow RenderBoth() => (map_path: String, table_path: String) andThen {

        map = amr.maps.BuildResistanceMap(dependency_signal = 0)

        // referencing map.country_count is what orders this second
        tbl = amr.maps.BuildGeneSpeciesTable(dependency_signal = map.country_count)

        yield RenderBoth(map_path = map.html_path, table_path = tbl.html_path)
    }
}
```

## 6. Call-time mixins — timeouts and retries

A facet ships its own defaults (`DownloadOrganism` declares
`with Timeout(minutes = 45)`); the **call site** can override them for one
particular use without forking the facet.

```ffl
namespace my.amr {

    use amr.sources
    use amr.maps

    /** Give the biggest organism groups more headroom, and retry NCBI hiccups. */
    workflow PatientIngest(group: String = "Salmonella_enterica") => (isolates: Long) andThen {

        dl = amr.sources.DownloadOrganism(group = $.group, force = true) with Timeout(minutes = 120) with Retry(maxAttempts = 3, backoffSeconds = 120)

        map = amr.maps.BuildResistanceMap(dependency_signal = dl.isolate_count) with Timeout(minutes = 30)

        yield PatientIngest(isolates = map.isolate_count)
    }
}
```

## 7. Don't let one bad group kill the run — `catch`

`catch` fires when a step errors after its retries are exhausted. Inside a
`foreach`, this is per-iteration: one organism failing leaves the others alone.

```ffl
namespace my.amr {

    use amr.sources

    /** Best-effort ingest: record failures instead of failing the run. */
    workflow BestEffortIngest(quick: Boolean = false) => (built: [String]) andThen {

        groups = amr.sources.ListOrganisms(quick = $.quick) andThen foreach g in $.organisms {

            dl = amr.sources.DownloadOrganism(group = $.g) catch {
                yield BestEffortIngest(built = ["failed"])
            }

            yield BestEffortIngest(built = [dl.group])
        }
    }
}
```

## 8. Branch on a result — `when`

A `when` block hangs off the step it inspects: inside a case `$` is that step and
`$$` reaches the workflow. Every `when` needs a default case, last.

```ffl
namespace my.amr {

    use amr.sources
    use amr.maps

    /** Only render if the ingest actually produced isolates. */
    workflow GuardedRender(group: String = "Vibrio_cholerae", min_isolates: Long = 1000) => (status: String, html_path: String) andThen {

        dl = amr.sources.DownloadOrganism(group = $.group) andThen when {
            case $.isolate_count >= $$.min_isolates => {
                map = amr.maps.BuildResistanceMap(dependency_signal = $.isolate_count)
                yield GuardedRender(status = "completed", html_path = map.html_path)
            }
            case _ => {
                yield GuardedRender(status = "too_few_isolates", html_path = "")
            }
        }
    }
}
```

## 9. Compose across domains — publish the map

Facets from different domains compose freely, as long as some runner in the fleet
serves each namespace. `census.Publish` is the generic publisher every map domain
reuses.

```ffl
namespace my.amr {

    use amr.maps
    use census.Publish

    /** Render, then push to the public maps site. */
    workflow RenderAndPublish(repo: String = "rlemke/facetwork-maps") => (pages_url: String) andThen {

        map = amr.maps.BuildResistanceMap(dependency_signal = 0)

        published = census.Publish.PublishWebBundle(
            repo = $.repo,
            prefixes = ["amr/output"],
            dests = ["amr"],
            labels = ["Antibiotic-resistance gene spread"],
            landing_title = "Facetwork maps")

        yield RenderAndPublish(pages_url = published.pages_url)
    }
}
```

Add `--library ~/fw_handlers/fwh_census_us/src/census_us/ffl/census.ffl` when
compiling that one.

---

## Cheat sheet

| You want to… | Write |
|---|---|
| Read a workflow/step parameter | `$.name` (`$$.name` one level out) |
| Read a previous step's result | `stepname.field` |
| Order two independent steps | reference a field of the first from the second |
| Fan out over a list | `step = Facet(…) andThen foreach v in $.list { … }` |
| More time / retries for one call | `… with Timeout(minutes = 120) with Retry(maxAttempts = 3, backoffSeconds = 120)` |
| Handle a step failure | `step = Facet(…) catch { yield … }` |
| Branch | `step = Facet(…) andThen when { case <bool> => { … } case _ => { … } }` |
| Concatenate strings | `a ++ b` |

**Validate before you run:** `afl my.ffl --check` or MCP `fw_validate`. Every error
carries a `rule_id` — fetch `fw://docs/rules/{rule_id}` for a wrong/right pair.

## See also

- [`docs/README.md`](README.md) — per-feature specs for this domain
- [FFL grammar](https://github.com/rlemke/facetwork/blob/main/docs/reference/language/grammar.md) ·
  [canonical examples](https://github.com/rlemke/facetwork/tree/main/examples/canonical) ·
  [relative `$`-scoping](https://github.com/rlemke/facetwork/blob/main/docs/architecture/ffl-relative-scoping.md)
- `src/amr/ffl/amr.ffl` — the source of truth for every signature above
