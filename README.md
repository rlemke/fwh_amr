# fwh_amr — antibiotic-resistance gene spread (geography + time + species)

A Facetwork domain that maps how **antibiotic-resistance genes** are distributed
across **countries**, **years**, and **bacterial species**, from
[NCBI Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens/) — millions of
public bacterial/fungal isolates with precomputed
[AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinderPlus/)
resistance-gene calls. No human-genome privacy concerns; all data is open.

It produces two artifacts:

1. **Resistance map** (`index.html`) — a MapLibre world choropleth of resistance-gene
   *prevalence* (share of sequenced isolates carrying a gene of a chosen drug class)
   with a **year slider** and a **drug-class dropdown** (carbapenemase, colistin `mcr`,
   ESBL, vancomycin, aminoglycoside, …). Critical WHO-priority mechanisms are flagged.
2. **Gene × species table** (`genes.html`) — a ranked, explainable table: each
   resistance gene × pathogen species, with isolate count, country count, and
   first-seen year, so **emerging gene/organism combinations** surface. Critical
   mechanisms sort first; filter by class, "critical only", or free-text search.

## How it works (and what it deliberately does *not* do)

- **Ingests NCBI's precomputed calls, not raw sequence.** Each pathogen group's
  per-isolate `metadata.tsv` already contains the `AMR_genotypes` column (the
  AMRFinderPlus gene calls) alongside `collection_date`, `geo_loc_name`, and
  `scientific_name`. We **stream-parse** that one file per organism (bounded
  memory, even for the ~840 MB *Salmonella* table) — **no assemblies are
  downloaded and AMRFinderPlus is not re-run.**
- **Fan-out / fan-in.** One `DownloadOrganism` per pathogen group runs in parallel
  across the fleet, each caching a compact aggregate into a shared dir; the map and
  table steps fan those in.

### Honesty boundaries (surfaced in every artifact's "About this data" modal)

- Two organisms carrying the same gene is **co-occurrence, not proof of horizontal
  gene transfer** — that needs plasmid/phylogenetic analysis this view does not do.
- Sampling is **heavily surveillance-biased** (US/EU sequence far more isolates); a
  dark country can mean more *surveillance*, not more resistance. Country-years with
  < 5 sequenced isolates are greyed.
- A gene's presence is a **genotype, not a phenotype**.
- Scope is NCBI's **foodborne + healthcare pathogen panel**, not all bacteria.

## Install

```bash
pip install -e ~/fw_handlers/fwh_amr      # or: fw install domain amr --check
```

Discovered by the runner via the `facetwork.domains` entry point (`amr = "amr:domain"`).

## Run

Quick single-organism demo (Vibrio cholerae, ~15 MB metadata, ~1 min):

```bash
fw ffl run --library src/amr/ffl/amr.ffl \
  --workflow amr.workflows.BuildOrganismAtlas \
  --inputs '{"group":"Vibrio_cholerae"}'
```

Full curated panel (14 pathogen groups) — fan-out then fan-in:

```bash
# 1. fan out one DownloadOrganism per pathogen group (parallel across runners)
fw ffl run --library src/amr/ffl/amr.ffl --workflow amr.workflows.BuildAmrFanout
# 2. render the map + table once the aggregates are cached
fw ffl run --library src/amr/ffl/amr.ffl --workflow amr.workflows.BuildAmrMaps
```

Outputs land in the storage backend under `cache/amr/output/{index.html,genes.html}`
(MinIO on the fleet; a local dir in single-box dev).

## Workflows

| Workflow | Does |
|----------|------|
| `BuildOrganismAtlas(group, force)` | one pathogen group end-to-end (download → map → table) |
| `BuildAmrFanout(quick, force, max_rows)` | fan out `DownloadOrganism` over the curated panel (`quick` = demo set only) |
| `BuildAmrMaps(force)` | fan-in: build the resistance map + gene × species table |

## FFL at a glance

The domain is driven from [FFL](https://github.com/rlemke/facetwork/blob/main/docs/reference/language/grammar.md),
Facetwork's workflow language. A step is `name = Facet(args)`; later steps
reference earlier ones as `step.field`; `andThen foreach` fans a step out across
the fleet:

```ffl
namespace my.amr {

    use amr.sources
    use amr.maps

    /** Ingest every curated pathogen group in parallel, then render the map. */
    workflow IngestEverything(quick: Boolean = false) => (built: [String]) andThen {

        groups = amr.sources.ListOrganisms(quick = $.quick) andThen foreach g in $.organisms {

            dl = amr.sources.DownloadOrganism(group = $.g, force = false)

            yield IngestEverything(built = [dl.group])
        }
    }
}
```

```bash
fw ffl run --primary my.ffl --library src/amr/ffl/amr.ffl \
  --workflow my.amr.IngestEverything
```

📖 **[docs/ffl-examples.md](docs/ffl-examples.md)** — the full example gallery:
fan-out/fan-in, call-time mixins (timeout/retry), `catch` so one bad organism
doesn't kill the run, `when` branching, and cross-domain composition (publishing).
Every snippet there is compile-checked.

## Data source

NCBI Pathogen Detection FTP: `ftp.ncbi.nlm.nih.gov/pathogen/Results/<Group>/latest_snps/Metadata/<PDG…>.metadata.tsv`.
The curated pathogen panel and the gene → drug-class classification live in
[`src/amr/organisms.py`](src/amr/organisms.py); the ingest/render logic in
[`src/amr/_lib.py`](src/amr/_lib.py).

Tag topics: `facetwork`, `facetwork-domain`, `amr`, `antimicrobial-resistance`,
`genomics`, `ncbi-pathogen-detection`.
