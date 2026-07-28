"""Event facet handlers for the AMR domain — thin layers over ``_lib``."""

from __future__ import annotations

import os
from typing import Any

from .. import organisms as org
from .._lib import (
    build_gene_species_table,
    build_resistance_map,
    download_organism,
)

SRC = "amr.sources"
MAPS = "amr.maps"


def handle_list_organisms(params: dict[str, Any]) -> dict[str, Any]:
    """Return the curated pathogen group list (or the small quick demo set)."""
    step_log = params.get("_step_log")
    groups = org.QUICK_ORGANISMS if bool(params.get("quick")) else org.organism_groups()
    if step_log:
        step_log(f"ListOrganisms: {len(groups)} groups", level="success")
    return {"organisms": list(groups), "count": len(groups)}


def handle_download_organism(params: dict[str, Any]) -> dict[str, Any]:
    """Stream one NCBI Pathogen Detection organism group + cache its aggregate."""
    step_log = params.get("_step_log")
    group = params.get("group") or ""
    try:
        res = download_organism(
            group,
            force=bool(params.get("force")),
            max_rows=int(params.get("max_rows") or 0),
        )
        if step_log:
            step_log(
                f"DownloadOrganism[{group}]: {res.n_isolates} isolates, "
                f"{res.n_dated} dated, {res.n_geo_matched} geo-matched "
                f"({res.n_geo_unmatched} unmatched), {res.n_genes} genes",
                level="success",
            )
        return {
            "group": res.group,
            "isolate_count": res.n_isolates,
            "dated_count": res.n_dated,
            "geo_matched": res.n_geo_matched,
            "geo_unmatched": res.n_geo_unmatched,
            "gene_count": res.n_genes,
        }
    except Exception as exc:
        if step_log:
            step_log(f"DownloadOrganism[{group}]: {exc}", level="error")
        raise


def handle_build_resistance_map(params: dict[str, Any]) -> dict[str, Any]:
    """Fan-in: merge organism aggregates + render the year-slider world choropleth."""
    step_log = params.get("_step_log")
    try:
        res = build_resistance_map()
        if step_log:
            step_log(
                f"BuildResistanceMap: {res.country_count} countries, "
                f"{res.organism_count} organisms, {res.isolate_count} isolates, "
                f"{res.year_min}-{res.year_max} -> {res.html_path}",
                level="success",
            )
        return {
            "html_path": res.html_path,
            "geojson_path": res.geojson_path,
            "country_count": res.country_count,
            "organism_count": res.organism_count,
            "isolate_count": res.isolate_count,
            "year_min": res.year_min,
            "year_max": res.year_max,
        }
    except Exception as exc:
        if step_log:
            step_log(f"BuildResistanceMap: {exc}", level="error")
        raise


def handle_build_gene_species_table(params: dict[str, Any]) -> dict[str, Any]:
    """Fan-in: merge per-gene rollups + render the ranked gene x species table."""
    step_log = params.get("_step_log")
    try:
        res = build_gene_species_table()
        if step_log:
            step_log(
                f"BuildGeneSpeciesTable: {res.row_count} rows across "
                f"{res.organism_count} organisms -> {res.html_path}",
                level="success",
            )
        return {
            "html_path": res.html_path,
            "row_count": res.row_count,
            "organism_count": res.organism_count,
        }
    except Exception as exc:
        if step_log:
            step_log(f"BuildGeneSpeciesTable: {exc}", level="error")
        raise


_DISPATCH: dict[str, Any] = {
    f"{SRC}.ListOrganisms": handle_list_organisms,
    f"{SRC}.DownloadOrganism": handle_download_organism,
    f"{MAPS}.BuildResistanceMap": handle_build_resistance_map,
    f"{MAPS}.BuildGeneSpeciesTable": handle_build_gene_species_table,
}


def handle(payload: dict) -> dict:
    facet = payload["_facet_name"]
    handler = _DISPATCH.get(facet)
    if handler is None:
        raise ValueError(f"Unknown facet: {facet}")
    return handler(payload)


def register_handlers(runner) -> None:
    for facet_name in _DISPATCH:
        runner.register_handler(
            facet_name=facet_name,
            module_uri=f"file://{os.path.abspath(__file__)}",
            entrypoint="handle",
        )


def register_poller(poller) -> None:
    for facet_name, handler in _DISPATCH.items():
        poller.register(facet_name, handler)
