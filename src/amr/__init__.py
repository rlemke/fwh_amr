"""AMR domain — antibiotic-resistance gene spread from NCBI Pathogen Detection.

Builds two artifacts from NCBI Pathogen Detection's open per-isolate metadata
(``AMR_genotypes`` = AMRFinderPlus gene calls, ``collection_date``,
``geo_loc_name``, ``scientific_name``):

1. **BuildResistanceMap** — a world choropleth of resistance-gene *prevalence*
   (share of sequenced isolates carrying a gene of a chosen drug class) with a
   YEAR SLIDER and a drug-class dropdown (carbapenem / colistin / ESBL / …).
2. **BuildGeneSpeciesTable** — a ranked, explainable table of gene × species
   co-occurrences over time (first-seen year, isolate count, countries), so
   emerging gene/organism combinations surface.

The pipeline fans out one ``DownloadOrganism`` per pathogen group and fans in.
It ingests NCBI's *precomputed* AMRFinderPlus calls — no assemblies are
downloaded and AMRFinderPlus is not re-run. Discovered by the Facetwork runner
via the ``facetwork.domains`` entry point in pyproject.toml::

    [project.entry-points."facetwork.domains"]
    amr = "amr:domain"

Honesty boundaries (surfaced in every artifact's "About this data" modal):
co-occurrence is **not** horizontal-gene-transfer proof; sampling is heavily
surveillance-biased (US/EU over-represented); genotype presence is **not** a
phenotype; scope is NCBI's foodborne + healthcare pathogen panel, not all
bacteria.
"""

from __future__ import annotations

from pathlib import Path

from facetwork.domains import DomainPackage

from .handlers import register_all_registry_handlers

domain = DomainPackage(
    name="amr",
    ffl_dir=Path(__file__).parent / "ffl",
    register_handlers=register_all_registry_handlers,
)
