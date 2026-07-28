"""AMR gene spread — ingest NCBI Pathogen Detection, build map + table.

Pipeline (all backend-aware via :mod:`amr.storage`):

1. ``download_organism(group)`` — stream one organism group's per-isolate
   ``metadata.tsv`` from the NCBI Pathogen Detection FTP (the *precomputed*
   AMRFinderPlus calls in the ``AMR_genotypes`` column — no assemblies fetched,
   AMRFinderPlus not re-run), parse it line-by-line (bounded memory), and cache
   a compact per-organism aggregate: resistance-class prevalence by
   country×year, and a per-gene rollup (first-seen year, isolate count,
   countries). Written to the *shared* organisms dir so a fan-out of one
   download per pathogen writes side-by-side artifacts.
2. ``build_resistance_map()`` — fan-in: merge every cached organism aggregate,
   join country prevalence onto Natural Earth geometry, and render a MapLibre
   world choropleth with a YEAR SLIDER + a drug-class dropdown.
3. ``build_gene_species_table()`` — fan-in: merge the per-gene rollups into a
   ranked, explainable gene × species table (critical mechanisms first), with
   class/critical/search filters and first-seen-year "emerging" highlighting.

Honesty boundaries are stated in each artifact's "About this data" modal:
co-occurrence is not HGT proof; sampling is surveillance-biased; genotype ≠
phenotype; scope is NCBI's foodborne + healthcare panel.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from html import escape

from . import organisms as org
from . import storage as cstore

try:
    import requests
except ImportError:  # pragma: no cover
    requests = None

logger = logging.getLogger("amr")

FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results"
WORLD_GEOJSON_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/"
    "geojson/ne_110m_admin_0_countries.geojson"
)
REPO_URL = "https://github.com/rlemke/fwh_amr"
FFL_URL = f"{REPO_URL}/blob/main/src/amr/ffl/amr.ffl"
USER_AGENT = "facetwork-amr/1.0 (+https://github.com/rlemke/facetwork)"

YEAR_MIN, YEAR_MAX = 1995, 2025  # slider bounds; older isolates are sparse/retrospective
MIN_DENOM = 5  # country-years with fewer sequenced isolates are greyed (noisy prevalence)


@dataclass
class OrganismResult:
    group: str
    label: str
    n_isolates: int
    n_dated: int
    n_geo_matched: int
    n_geo_unmatched: int
    n_genes: int


@dataclass
class MapResult:
    html_path: str
    geojson_path: str
    country_count: int
    organism_count: int
    isolate_count: int
    year_min: int
    year_max: int


@dataclass
class TableResult:
    html_path: str
    row_count: int
    organism_count: int


# ---------------------------------------------------------------------------
# 1. Download + aggregate one organism
# ---------------------------------------------------------------------------


def _resolve_metadata_url(group: str) -> str:
    """Find the ``<PDG…>.metadata.tsv`` under the group's latest_snps/Metadata/."""
    if requests is None:
        raise RuntimeError("requests is required to reach NCBI Pathogen Detection")
    listing = f"{FTP_BASE}/{group}/latest_snps/Metadata/"
    resp = requests.get(listing, headers={"User-Agent": USER_AGENT}, timeout=(30, 90))
    resp.raise_for_status()
    import re
    m = re.search(r'href="(PDG[^"]+\.metadata\.tsv)"', resp.text)
    if not m:
        raise RuntimeError(f"{group}: no metadata.tsv found under {listing}")
    return listing + m.group(1)


def _parse_year(collection_date: str) -> int | None:
    if not collection_date:
        return None
    v = collection_date.strip().strip('"')
    if len(v) < 4 or not v[:4].isdigit():
        return None
    y = int(v[:4])
    return y if 1900 <= y <= YEAR_MAX else None


def _parse_genes(amr_genotypes: str) -> list[str]:
    """Split the AMR_genotypes cell into bare gene symbols."""
    v = (amr_genotypes or "").strip().strip('"').strip()
    if not v or v.upper() == "NULL":
        return []
    out = []
    for tok in v.split(","):
        # values can be "gene" or "gene=COMPLETE" / "gene=PARTIAL" in some versions
        sym = tok.split("=", 1)[0].strip().strip('"').strip()
        if sym and sym.upper() != "NULL":
            out.append(sym)
    return out


def download_organism(group: str, *, force: bool = False, max_rows: int = 0) -> OrganismResult:
    """Stream one organism's metadata.tsv, aggregate, cache. Returns counts."""
    if group not in org.organism_groups():
        logger.warning("amr: %s is not in the curated panel (proceeding anyway)", group)
    out_key = cstore.join(cstore.organisms_root(), f"{group}.json")
    if not force and cstore.exists(out_key):
        with cstore.open_read(out_key) as f:
            blob = json.load(f)
        return OrganismResult(
            group, blob.get("label", group), blob.get("n_isolates", 0),
            blob.get("n_dated", 0), blob.get("n_geo_matched", 0),
            blob.get("n_geo_unmatched", 0), len(blob.get("genes", {})),
        )

    if requests is None:
        raise RuntimeError("requests is required to fetch NCBI Pathogen Detection data")
    url = _resolve_metadata_url(group)
    logger.info("amr: streaming %s", url)

    # prevalence[iso3][year] = [total, {class: count_isolates_with_class}]
    prevalence: dict[str, dict[int, list]] = {}
    # genes[symbol] = {class, first_year, n, countries(set)}
    genes: dict[str, dict] = {}
    n_isolates = n_dated = n_geo_matched = n_geo_unmatched = 0

    resp = requests.get(url, headers={"User-Agent": USER_AGENT}, stream=True, timeout=(30, 300))
    resp.raise_for_status()
    resp.encoding = "utf-8"
    header: list[str] | None = None
    idx: dict[str, int] = {}
    truncated = False
    for raw in resp.iter_lines(decode_unicode=True):
        if raw is None or raw == "":
            continue
        row = raw.split("\t")
        if header is None:
            header = row
            idx = {name: i for i, name in enumerate(header)}
            missing = [c for c in ("collection_date", "geo_loc_name", "scientific_name", "AMR_genotypes") if c not in idx]
            if missing:
                raise RuntimeError(f"{group}: metadata.tsv missing columns {missing}")
            continue
        if len(row) <= idx["AMR_genotypes"]:
            continue
        n_isolates += 1
        if max_rows and n_isolates > max_rows:
            truncated = True
            n_isolates -= 1
            break

        year = _parse_year(row[idx["collection_date"]])
        iso3 = org.country_to_iso3(row[idx["geo_loc_name"]])
        gene_syms = _parse_genes(row[idx["AMR_genotypes"]])
        classes_present: set[str] = set()
        for sym in gene_syms:
            cls = org.classify_gene(sym)
            classes_present.add(cls)
            g = genes.get(sym)
            if g is None:
                g = genes[sym] = {"class": cls, "first_year": year, "n": 0, "countries": set()}
            g["n"] += 1
            if year is not None and (g["first_year"] is None or year < g["first_year"]):
                g["first_year"] = year
            if iso3:
                g["countries"].add(iso3)

        if year is not None:
            n_dated += 1
        if iso3:
            n_geo_matched += 1
        elif row[idx["geo_loc_name"]].strip().strip('"') not in ("", "NULL"):
            n_geo_unmatched += 1

        # prevalence needs a matched country AND a year in-window
        if iso3 and year is not None and YEAR_MIN <= year <= YEAR_MAX:
            slot = prevalence.setdefault(iso3, {}).setdefault(year, [0, {}])
            slot[0] += 1
            for cls in classes_present:
                slot[1][cls] = slot[1].get(cls, 0) + 1

    if truncated:
        logger.warning("amr: %s truncated at max_rows=%d (aggregates are a sample)", group, max_rows)

    genes_out = {
        sym: {
            "class": g["class"],
            "first_year": g["first_year"],
            "n": g["n"],
            "countries": sorted(g["countries"]),
        }
        for sym, g in genes.items()
    }
    blob = {
        "group": group,
        "label": org.organism_label(group),
        "n_isolates": n_isolates,
        "n_dated": n_dated,
        "n_geo_matched": n_geo_matched,
        "n_geo_unmatched": n_geo_unmatched,
        "truncated": truncated,
        "prevalence": {iso: {str(y): v for y, v in yrs.items()} for iso, yrs in prevalence.items()},
        "genes": genes_out,
    }
    with cstore.open_write(out_key, "w") as f:
        json.dump(blob, f, separators=(",", ":"))
    logger.info(
        "amr: %s — %d isolates, %d dated, %d geo-matched (%d unmatched), %d genes",
        group, n_isolates, n_dated, n_geo_matched, n_geo_unmatched, len(genes_out),
    )
    return OrganismResult(group, blob["label"], n_isolates, n_dated,
                          n_geo_matched, n_geo_unmatched, len(genes_out))


# ---------------------------------------------------------------------------
# fan-in helpers
# ---------------------------------------------------------------------------


def _load_organism_aggregates() -> list[dict]:
    root = cstore.organisms_root()
    paths = [p for p in cstore.list_dir(root) if p.endswith(".json")]
    blobs = []
    for p in paths:
        try:
            with cstore.open_read(p) as f:
                blobs.append(json.load(f))
        except Exception as exc:  # pragma: no cover
            logger.warning("amr: skipping unreadable aggregate %s: %s", p, exc)
    if not blobs:
        raise RuntimeError("no organism aggregates cached — run DownloadOrganism first")
    return blobs


def _world_geojson() -> dict:
    cache_key = cstore.join(cstore.cache_root(), "world-countries.geojson")
    if cstore.exists(cache_key):
        with cstore.open_read(cache_key) as f:
            return json.load(f)
    if requests is None:
        raise RuntimeError("requests is required to download world geometry")
    resp = requests.get(WORLD_GEOJSON_URL, headers={"User-Agent": USER_AGENT}, timeout=(30, 120))
    resp.raise_for_status()
    gj = resp.json()
    with cstore.open_write(cache_key, "w") as f:
        json.dump(gj, f, separators=(",", ":"))
    return gj


# ---------------------------------------------------------------------------
# 2. Resistance map (fan-in)
# ---------------------------------------------------------------------------


def build_resistance_map() -> MapResult:
    blobs = _load_organism_aggregates()
    classes = list(org.DRUG_CLASSES.keys())

    # Merge prevalence across organisms: PREV[iso3][year] = [total, {class: count}]
    merged: dict[str, dict[int, list]] = {}
    total_isolates = 0
    ymin, ymax = YEAR_MAX, YEAR_MIN
    for b in blobs:
        total_isolates += b.get("n_isolates", 0)
        for iso, yrs in b.get("prevalence", {}).items():
            for ys, (tot, cls_counts) in yrs.items():
                y = int(ys)
                ymin, ymax = min(ymin, y), max(ymax, y)
                slot = merged.setdefault(iso, {}).setdefault(y, [0, {}])
                slot[0] += tot
                for cls, c in cls_counts.items():
                    slot[1][cls] = slot[1].get(cls, 0) + c
    if ymin > ymax:
        ymin, ymax = YEAR_MIN, YEAR_MAX

    # Compact per-country prevalence object for the client:
    #   PREV[iso3] = {year: {"n": total, cls: pct_int, ...}, ...}
    prev_client: dict[str, dict] = {}
    for iso, yrs in merged.items():
        rec = {}
        for y, (tot, cls_counts) in yrs.items():
            if tot <= 0:
                continue
            entry = {"n": tot}
            for cls, c in cls_counts.items():
                entry[cls] = round(c / tot * 100.0, 1)
            rec[str(y)] = entry
        if rec:
            prev_client[iso] = rec

    world = _world_geojson()
    keep = []
    for ft in world["features"]:
        props = ft["properties"]
        iso = props.get("ISO_A3")
        if not iso or iso == "-99":
            iso = props.get("ADM0_A3")
        name = props.get("NAME", "")
        ft["properties"] = {"NAME": name, "ISO_A3": iso}
        keep.append(ft)
    fc = {"type": "FeatureCollection", "features": keep}

    years = list(range(ymin, ymax + 1))
    html = _render_map_html(fc, prev_client, classes=classes, years=years,
                            organism_count=len(blobs), isolate_count=total_isolates)

    out_dir = cstore.output_root()
    geojson_path = cstore.join(out_dir, "resistance.geojson")
    html_path = cstore.join(out_dir, "index.html")
    with cstore.open_write(geojson_path, "w") as f:
        json.dump(fc, f, separators=(",", ":"))
    with cstore.open_write(html_path, "w") as f:
        f.write(html)
    logger.info("amr map: %d countries with data, %d organisms, %d isolates, %d-%d",
                len(prev_client), len(blobs), total_isolates, ymin, ymax)
    return MapResult(html_path, geojson_path, len(prev_client), len(blobs),
                     total_isolates, ymin, ymax)


# ---------------------------------------------------------------------------
# 3. Gene x species table (fan-in)
# ---------------------------------------------------------------------------


def build_gene_species_table() -> TableResult:
    blobs = _load_organism_aggregates()
    rows = []
    for b in blobs:
        species = b.get("label", b.get("group", "?"))
        for sym, g in b.get("genes", {}).items():
            cls = g.get("class", "other")
            rows.append({
                "gene": sym,
                "class": cls,
                "class_name": org.DRUG_CLASSES.get(cls, {}).get("name", cls),
                "critical": cls in org.CRITICAL_CLASSES,
                "species": species,
                "first_year": g.get("first_year"),
                "n": g.get("n", 0),
                "n_countries": len(g.get("countries", [])),
            })
    # rank: critical first, then by isolate count
    rows.sort(key=lambda r: (not r["critical"], -r["n"]))
    html = _render_table_html(rows, organism_count=len(blobs))
    out_dir = cstore.output_root()
    html_path = cstore.join(out_dir, "genes.html")
    with cstore.open_write(html_path, "w") as f:
        f.write(html)
    logger.info("amr table: %d gene x species rows across %d organisms", len(rows), len(blobs))
    return TableResult(html_path, len(rows), len(blobs))


# ---------------------------------------------------------------------------
# Rendering — shared chrome
# ---------------------------------------------------------------------------


def _ts() -> str:
    from datetime import UTC, datetime
    return datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")


def _attribution(call: str) -> str:
    return (
        '<div style="position:fixed;bottom:10px;left:10px;z-index:9999;'
        "background:rgba(255,255,255,0.92);border-radius:6px;padding:6px 10px;"
        "box-shadow:0 1px 4px rgba(0,0,0,0.2);font:11px system-ui,sans-serif;color:#444;"
        'max-width:480px">Generated by Facetwork workflow '
        '<code style="background:#f0f0f0;padding:0 3px;border-radius:3px">'
        f"{escape(call)}</code> &middot; "
        f'<a href="{escape(FFL_URL)}" target="_blank" rel="noopener" '
        'style="color:#1565c0;text-decoration:none">view FFL</a>'
        f' &middot; <a href="{escape(REPO_URL)}" target="_blank" '
        'rel="noopener" style="color:#1565c0;text-decoration:none">source repo</a>'
        f" &middot; NCBI Pathogen Detection &middot; generated {_ts()}</div>"
    )


def _about_body(kind: str) -> str:
    repo = escape(REPO_URL)
    common = (
        '<p style="margin:0 0 10px"><b>Source:</b> NCBI Pathogen Detection — the open '
        "per-isolate metadata tables (<code>AMR_genotypes</code>, the precomputed "
        "AMRFinderPlus gene calls, plus <code>collection_date</code>, "
        "<code>geo_loc_name</code>, <code>scientific_name</code>). No assemblies are "
        "downloaded and AMRFinderPlus is not re-run — we ingest NCBI's own calls.</p>"
        '<p style="margin:0 0 10px"><b>Read this honestly:</b></p>'
        '<ul style="margin:0 0 10px;padding-left:18px">'
        "<li>Two organisms carrying the <i>same gene</i> is <b>co-occurrence, not proof of "
        "transfer</b> between them — horizontal gene transfer needs plasmid/phylogenetic "
        "analysis this view does not do.</li>"
        "<li>Sampling is <b>heavily surveillance-biased</b>: the US and Europe sequence far "
        "more isolates, so a dark country may mean more <i>surveillance</i>, not more "
        "resistance. Country-years with fewer than "
        f"{MIN_DENOM} sequenced isolates are greyed.</li>"
        "<li>A gene's presence is a <b>genotype, not a phenotype</b> — carriage does not "
        "guarantee clinical resistance.</li>"
        "<li>Scope is NCBI's <b>foodborne + healthcare pathogen</b> panel, not all bacteria.</li>"
        "</ul>"
    )
    if kind == "map":
        lead = (
            '<p style="margin:0 0 10px"><b>Prevalence</b> is the share of sequenced isolates '
            "in a country-year carrying at least one gene of the selected drug class. Pick a "
            "class, drag the year slider or press play; click a country for its numbers.</p>"
        )
    else:
        lead = (
            '<p style="margin:0 0 10px">Each row is one <b>resistance gene × pathogen species</b>: '
            "how many sequenced isolates of that species carry the gene, in how many countries, and "
            "the earliest collection year seen. <b>Critical mechanisms</b> (carbapenemase, colistin "
            "<code>mcr</code>, ESBL, vancomycin) are flagged and sorted first.</p>"
        )
    return (
        lead + common +
        '<p style="margin:0">Generated by Facetwork &middot; '
        f'<a href="{escape(FFL_URL)}" target="_blank" rel="noopener">view FFL</a> &middot; '
        f'<a href="{repo}" target="_blank" rel="noopener">source repo</a>.</p>'
    )


def _modal(kind: str) -> str:
    return (
        '<div id="infomodal" class="modal" style="display:none"><div class="modalcard">'
        '<button id="infoclose" class="modalclose">&times;</button>'
        '<h2>About this data</h2>'
        f'<div class="modalbody">{_about_body(kind)}</div></div></div>'
    )


_MODAL_CSS = (
    ".infobtn{margin-top:8px;width:100%;padding:7px;font-size:13px;cursor:pointer;"
    "background:#fff8e1;border:1px solid #f6c343;border-radius:4px;color:#5d4b00}"
    ".infobtn:hover{background:#fff2c4}"
    ".modal{position:fixed;inset:0;z-index:50;background:rgba(0,0,0,.45);display:flex;"
    "align-items:center;justify-content:center}"
    ".modalcard{background:#fff;max-width:480px;margin:16px;padding:18px 20px 20px;"
    "border-radius:10px;box-shadow:0 4px 24px rgba(0,0,0,.4);position:relative;"
    "font:13px/1.5 system-ui,sans-serif;max-height:82vh;overflow:auto}"
    ".modalcard h2{margin:0 0 8px;font-size:16px} .modalbody b{color:#111}"
    ".modalclose{position:absolute;top:6px;right:10px;border:none;background:none;font-size:24px;"
    "line-height:1;cursor:pointer;color:#999} .modalclose:hover{color:#444}"
)

_MODAL_JS = (
    "const _im=document.getElementById('infomodal');"
    "if(_im){const ib=document.getElementById('infobtn'),ic=document.getElementById('infoclose');"
    "if(ib)ib.onclick=()=>{_im.style.display='flex';};"
    "if(ic)ic.onclick=()=>{_im.style.display='none';};"
    "_im.onclick=e=>{if(e.target===_im)_im.style.display='none';};}"
)


# ---------------------------------------------------------------------------
# Rendering — the map
# ---------------------------------------------------------------------------


def _render_map_html(fc: dict, prev: dict, *, classes: list[str], years: list[int],
                     organism_count: int, isolate_count: int) -> str:
    data_js = json.dumps(fc, separators=(",", ":"))
    prev_js = json.dumps(prev, separators=(",", ":"))
    class_meta = {k: {"name": v["name"], "color": v["color"], "critical": v["critical"]}
                  for k, v in org.DRUG_CLASSES.items() if k in classes}
    classes_js = json.dumps(class_meta, separators=(",", ":"))
    order_js = json.dumps(classes)
    y0, y1 = years[0], years[-1]
    call = "amr.workflows.BuildAmrAtlas"
    default_class = "carbapenem"
    return f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Antibiotic resistance by country, {y0}-{y1}</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<link href="https://unpkg.com/maplibre-gl@4.7.1/dist/maplibre-gl.css" rel="stylesheet">
<script src="https://unpkg.com/maplibre-gl@4.7.1/dist/maplibre-gl.js"></script>
<style>
  html,body,#map{{margin:0;height:100%;width:100%;font-family:system-ui,sans-serif}}
  .panel{{position:absolute;z-index:1;background:rgba(255,255,255,.95);padding:10px 12px;
    border-radius:6px;box-shadow:0 1px 4px rgba(0,0,0,.3);font-size:12px}}
  #ctl{{top:10px;left:10px;max-width:320px}}
  #ctl h3{{margin:0 0 6px;font-size:14px}}
  #ctl select{{width:100%;padding:5px;font-size:13px;margin:4px 0 2px}}
  #legend{{bottom:18px;right:10px}}
  #legend .bar{{height:12px;width:200px;border-radius:2px;margin:5px 0 2px}}
  #legend .tk{{display:flex;justify-content:space-between;font-size:10px;color:#555}}
  #tl{{bottom:14px;left:50%;transform:translateX(-50%);display:flex;align-items:center;gap:10px;
    width:min(560px,80vw)}}
  #tl input[type=range]{{flex:1}} #tl #yr{{font-weight:700;min-width:46px;text-align:center}}
  #tl button{{border:1px solid #aaa;background:#fff;border-radius:5px;padding:4px 10px;cursor:pointer;font-size:14px}}
  .maplibregl-popup-content{{max-width:280px;font-size:12px}}
  .maplibregl-popup-content h4{{margin:0 0 4px;font-size:13px}}
  table.m{{border-collapse:collapse;margin-top:4px}} table.m td{{padding:1px 6px 1px 0}}
  table.m td.v{{text-align:right;font-variant-numeric:tabular-nums}}
  .crit{{color:#b30000;font-weight:700}}
  {_MODAL_CSS}
</style></head>
<body>
<div id="map"></div>
<div id="ctl" class="panel">
  <h3>Antibiotic resistance genes</h3>
  <div style="color:#555">Share of sequenced isolates carrying a gene of the chosen
  drug class, by country and year. Source: NCBI Pathogen Detection
  ({organism_count} pathogen groups, {isolate_count:,} isolates).</div>
  <select id="cls"></select>
  <button id="infobtn" class="infobtn">&#8505;&#65039; About this data (read the caveats)</button>
</div>
<div id="legend" class="panel">
  <b id="legttl">Prevalence</b><div class="bar" id="legbar"></div>
  <div class="tk"><span>0%</span><span>&ge;25%</span></div>
  <div style="font-size:10px;color:#888;margin-top:3px">grey = &lt;{MIN_DENOM} isolates that year</div>
</div>
<div id="tl" class="panel">
  <button id="play" title="Play">&#9654;</button>
  <span id="yr">{y1}</span>
  <input id="year" type="range" min="{y0}" max="{y1}" step="1" value="{y1}">
</div>
{_modal("map")}
{_attribution(call)}
<script>
const DATA={data_js}, PREV={prev_js}, CLASSES={classes_js}, ORDER={order_js};
const Y0={y0}, Y1={y1}, MIN_DENOM={MIN_DENOM}, NODATA='#e6e6e6', LOWN='#cfcfcf';
let year={y1}, cls='{default_class}';
// linear mix white -> class colour by t in [0,1]
function hex2rgb(h){{return [parseInt(h.slice(1,3),16),parseInt(h.slice(3,5),16),parseInt(h.slice(5,7),16)];}}
function mix(c,t){{const a=[255,255,255],b=hex2rgb(c);
  const r=a.map((x,i)=>Math.round(x+(b[i]-x)*t));return 'rgb('+r.join(',')+')';}}
function colorFor(iso){{
  const rec=(PREV[iso]||{{}})[year]; if(!rec) return NODATA;
  if(rec.n<MIN_DENOM) return LOWN;
  const p=rec[cls]; if(p===undefined) return mix(CLASSES[cls].color,0); // 0% -> white
  return mix(CLASSES[cls].color, Math.min(1, p/25));
}}
function matchExpr(){{
  const e=['match',['get','ISO_A3']];
  for(const f of DATA.features){{const iso=f.properties.ISO_A3; if(!iso) continue;
    e.push(iso,colorFor(iso));}}
  e.push(NODATA); return e;
}}
const map=new maplibregl.Map({{container:'map',style:{{version:8,
  sources:{{bm:{{type:'raster',tiles:['https://a.basemaps.cartocdn.com/rastertiles/voyager/{{z}}/{{x}}/{{y}}.png','https://b.basemaps.cartocdn.com/rastertiles/voyager/{{z}}/{{x}}/{{y}}.png'],tileSize:256,attribution:'&copy; OpenStreetMap &copy; CARTO &middot; NCBI Pathogen Detection'}}}},
  layers:[{{id:'bm',type:'raster',source:'bm'}}]}},center:[14,25],zoom:1.4}});
map.addControl(new maplibregl.NavigationControl());
const slider=document.getElementById('year'), ylab=document.getElementById('yr');
const sel=document.getElementById('cls'), legttl=document.getElementById('legttl'), legbar=document.getElementById('legbar');
for(const k of ORDER){{const o=document.createElement('option');o.value=k;
  o.textContent=(CLASSES[k].critical?'\\u26a0 ':'')+CLASSES[k].name;sel.appendChild(o);}}
sel.value=cls;
function legend(){{legttl.innerHTML=(CLASSES[cls].critical?'<span class="crit">&#9888; </span>':'')+CLASSES[cls].name+' prevalence';
  legbar.style.background='linear-gradient(to right,#fff,'+CLASSES[cls].color+')';}}
function popupFor(f,lngLat){{
  const iso=f.properties.ISO_A3, rec=(PREV[iso]||{{}})[year];
  let body;
  if(!rec){{body='<div style="color:#888">No sequenced isolates this year.</div>';}}
  else{{let rowsH='';for(const k of ORDER){{const p=rec[k];if(p===undefined)continue;
      rowsH+=`<tr><td>${{(CLASSES[k].critical?'<span class=crit>&#9888;</span> ':'')+CLASSES[k].name}}</td><td class="v">${{p}}%</td></tr>`;}}
    if(!rowsH)rowsH='<tr><td colspan=2 style="color:#888">no resistance genes detected</td></tr>';
    const note=rec.n<MIN_DENOM?`<div style="color:#b30000;margin-top:3px">only ${{rec.n}} isolate(s) — too few to read prevalence</div>`:'';
    body=`<div style="color:#555">${{rec.n.toLocaleString()}} sequenced isolates in ${{year}}</div>`
      +`<table class="m">${{rowsH}}</table>`+note;}}
  new maplibregl.Popup({{closeButton:true,maxWidth:'280px'}}).setLngLat(lngLat)
    .setHTML(`<h4>${{f.properties.NAME||iso}}</h4>`+body).addTo(map);
}}
map.on('load',()=>{{
  map.addSource('c',{{type:'geojson',data:DATA}});
  map.addLayer({{id:'fill',type:'fill',source:'c',paint:{{'fill-color':matchExpr(),'fill-opacity':0.85}}}});
  map.addLayer({{id:'line',type:'line',source:'c',paint:{{'line-color':'#888','line-width':0.3}}}});
  function apply(){{ylab.textContent=year;legend();map.setPaintProperty('fill','fill-color',matchExpr());}}
  slider.oninput=()=>{{year=+slider.value;apply();}};
  sel.onchange=()=>{{cls=sel.value;apply();}};
  map.on('click','fill',e=>popupFor(e.features[0],e.lngLat));
  map.on('mouseenter','fill',()=>map.getCanvas().style.cursor='pointer');
  map.on('mouseleave','fill',()=>map.getCanvas().style.cursor='');
  let timer=null; const btn=document.getElementById('play');
  btn.onclick=()=>{{if(timer){{clearInterval(timer);timer=null;btn.innerHTML='&#9654;';return;}}
    btn.innerHTML='&#10073;&#10073;';
    timer=setInterval(()=>{{year=year>=Y1?Y0:year+1;slider.value=year;apply();}},600);}};
  apply();
}});
{_MODAL_JS}
</script></body></html>"""


# ---------------------------------------------------------------------------
# Rendering — the gene x species table
# ---------------------------------------------------------------------------


def _render_table_html(rows: list[dict], *, organism_count: int) -> str:
    rows_js = json.dumps(rows, separators=(",", ":"))
    class_opts = "".join(
        f'<option value="{escape(k)}">{escape(v["name"])}</option>'
        for k, v in org.DRUG_CLASSES.items()
    )
    call = "amr.workflows.BuildAmrAtlas"
    n_crit = sum(1 for r in rows if r["critical"])
    return f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Resistance gene &times; species</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
  body{{margin:0;font-family:system-ui,sans-serif;color:#222;background:#fafafa}}
  header{{padding:16px 20px 10px;background:#fff;border-bottom:1px solid #e5e5e5}}
  h1{{margin:0 0 4px;font-size:18px}} .sub{{color:#666;font-size:13px;max-width:760px}}
  .controls{{display:flex;flex-wrap:wrap;gap:10px;align-items:center;padding:10px 20px;background:#fff;
    border-bottom:1px solid #e5e5e5;position:sticky;top:0;z-index:2}}
  .controls input,.controls select{{padding:6px 8px;font-size:13px;border:1px solid #bbb;border-radius:5px}}
  .controls label{{font-size:13px;color:#444}}
  .wrap{{overflow-x:auto;padding:0 20px 60px}}
  table{{border-collapse:collapse;width:100%;background:#fff;margin-top:12px;font-size:13px}}
  th,td{{padding:6px 10px;text-align:left;border-bottom:1px solid #eee;white-space:nowrap}}
  th{{position:sticky;top:52px;background:#f3f3f3;cursor:pointer;user-select:none}}
  td.v{{text-align:right;font-variant-numeric:tabular-nums}}
  .pill{{display:inline-block;padding:1px 7px;border-radius:10px;font-size:11px;color:#fff}}
  .critrow{{background:#fff6f6}} .critbadge{{color:#b30000;font-weight:700}}
  .emerge{{color:#0a7d33;font-weight:600}}
  .infobtn{{margin-left:auto;padding:6px 10px;font-size:13px;cursor:pointer;background:#fff8e1;
    border:1px solid #f6c343;border-radius:5px;color:#5d4b00}}
  code{{background:#f0f0f0;padding:0 3px;border-radius:3px}}
  {_MODAL_CSS}
</style></head>
<body>
<header>
  <h1>Antibiotic-resistance genes by pathogen species</h1>
  <div class="sub">One row per resistance gene &times; species from NCBI Pathogen Detection
  ({organism_count} pathogen groups). <span class="critbadge">&#9888; {n_crit:,}</span> rows are
  critical mechanisms (carbapenemase, colistin <code>mcr</code>, ESBL, vancomycin). Co-occurrence
  is not proof of transfer — read the caveats.</div>
</header>
<div class="controls">
  <input id="q" placeholder="Search gene or species…" style="min-width:220px">
  <select id="clsf"><option value="">All drug classes</option>{class_opts}</select>
  <label><input type="checkbox" id="critonly"> critical only</label>
  <span id="count" style="color:#666;font-size:13px"></span>
  <button id="infobtn" class="infobtn">&#8505;&#65039; About this data</button>
</div>
<div class="wrap"><table id="t"><thead><tr>
  <th data-k="gene">Gene</th><th data-k="class_name">Drug class</th><th data-k="species">Species</th>
  <th data-k="first_year">First seen</th><th data-k="n">Isolates</th><th data-k="n_countries">Countries</th>
</tr></thead><tbody id="tb"></tbody></table></div>
{_modal("table")}
{_attribution(call)}
<script>
const ROWS={rows_js};
const tb=document.getElementById('tb'), q=document.getElementById('q'),
  clsf=document.getElementById('clsf'), critonly=document.getElementById('critonly'),
  count=document.getElementById('count');
let sortK='n', sortDir=-1;
const RECENT=(new Date().getFullYear?0:0); // no Date in build; "emerging" = first_year>=2015 client-side
function view(){{
  const term=q.value.trim().toLowerCase(), cf=clsf.value, co=critonly.checked;
  let r=ROWS.filter(x=>(!cf||x.class===cf)&&(!co||x.critical)
    &&(!term||x.gene.toLowerCase().includes(term)||x.species.toLowerCase().includes(term)));
  r.sort((a,b)=>{{const av=a[sortK],bv=b[sortK];
    if(av===null)return 1; if(bv===null)return -1;
    return (av<bv?-1:av>bv?1:0)*sortDir;}});
  count.textContent=r.length.toLocaleString()+' of '+ROWS.length.toLocaleString()+' rows';
  tb.innerHTML=r.slice(0,3000).map(x=>{{
    const fy=x.first_year===null?'&mdash;':x.first_year;
    const emerge=(x.first_year!==null&&x.first_year>=2015)?' <span class="emerge">new</span>':'';
    return `<tr class="${{x.critical?'critrow':''}}">`
      +`<td>${{x.critical?'<span class=critbadge>&#9888;</span> ':''}}${{x.gene}}</td>`
      +`<td>${{x.class_name}}</td><td><i>${{x.species}}</i></td>`
      +`<td class="v">${{fy}}${{emerge}}</td>`
      +`<td class="v">${{x.n.toLocaleString()}}</td><td class="v">${{x.n_countries}}</td></tr>`;
  }}).join('');
}}
document.querySelectorAll('th').forEach(th=>th.onclick=()=>{{
  const k=th.dataset.k; if(sortK===k)sortDir=-sortDir; else{{sortK=k;sortDir=(k==='gene'||k==='species'||k==='class_name')?1:-1;}}
  view();}});
q.oninput=view; clsf.onchange=view; critonly.onchange=view;
view();
{_MODAL_JS}
</script></body></html>"""
