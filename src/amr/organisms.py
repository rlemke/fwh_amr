"""Curated pathogen catalog + AMR gene → drug-class classification.

The organism list is the fan-out unit set: each entry is one NCBI Pathogen
Detection *organism group* (the directory name under
``ftp.ncbi.nlm.nih.gov/pathogen/Results/<group>/``). We ship a curated panel of
foodborne + WHO-priority healthcare pathogens with good Pathogen Detection
coverage rather than all ~107 groups.

``classify_gene`` maps an AMRFinderPlus gene symbol (as it appears in the
``AMR_genotypes`` metadata column) to a drug class. Classification is by
established gene-nomenclature prefixes/families; it is deliberately coarse and
documented — allele-level ESBL-vs-narrow distinctions (e.g. blaTEM-1 vs an ESBL
TEM variant) are **not** resolvable from the symbol alone and we say so.
``CRITICAL_CLASSES`` marks the WHO "critical priority" mechanisms that the map
and table highlight.
"""

from __future__ import annotations

import re

# --- Curated pathogen panel (FTP group name, display label, foodborne?) -------
# Verified present under Results/ on 2026-07-27.
ORGANISMS: list[dict] = [
    {"group": "Salmonella", "label": "Salmonella", "kind": "foodborne"},
    {"group": "Escherichia_coli_Shigella", "label": "E. coli / Shigella", "kind": "foodborne"},
    {"group": "Campylobacter", "label": "Campylobacter", "kind": "foodborne"},
    {"group": "Listeria", "label": "Listeria", "kind": "foodborne"},
    {"group": "Vibrio_cholerae", "label": "Vibrio cholerae", "kind": "foodborne"},
    {"group": "Vibrio_parahaemolyticus", "label": "Vibrio parahaemolyticus", "kind": "foodborne"},
    {"group": "Klebsiella", "label": "Klebsiella", "kind": "healthcare"},
    {"group": "Acinetobacter", "label": "Acinetobacter", "kind": "healthcare"},
    {"group": "Pseudomonas_aeruginosa", "label": "Pseudomonas aeruginosa", "kind": "healthcare"},
    {"group": "Staphylococcus_aureus", "label": "Staphylococcus aureus", "kind": "healthcare"},
    {"group": "Enterococcus_faecium", "label": "Enterococcus faecium", "kind": "healthcare"},
    {"group": "Enterococcus_faecalis", "label": "Enterococcus faecalis", "kind": "healthcare"},
    {"group": "Neisseria_gonorrhoeae", "label": "Neisseria gonorrhoeae", "kind": "healthcare"},
    {"group": "Clostridioides_difficile", "label": "Clostridioides difficile", "kind": "healthcare"},
]

# Groups small enough for a quick example/CI run (metadata < ~30 MB).
QUICK_ORGANISMS = ["Vibrio_cholerae"]


def organism_groups() -> list[str]:
    return [o["group"] for o in ORGANISMS]


def organism_label(group: str) -> str:
    for o in ORGANISMS:
        if o["group"] == group:
            return o["label"]
    return group.replace("_", " ")


# --- Drug-class catalog -------------------------------------------------------
# key -> (display name, choropleth ramp end colour, critical?)
DRUG_CLASSES: dict[str, dict] = {
    "carbapenem":    {"name": "Carbapenemase", "color": "#67000d", "critical": True},
    "colistin":      {"name": "Colistin (mcr)", "color": "#3f007d", "critical": True},
    "esbl":          {"name": "ESBL (CTX-M)", "color": "#7a0177", "critical": True},
    "glycopeptide":  {"name": "Vancomycin (van)", "color": "#4d004b", "critical": True},
    "beta_lactam":   {"name": "Other β-lactam", "color": "#08519c", "critical": False},
    "aminoglycoside": {"name": "Aminoglycoside", "color": "#006d2c", "critical": False},
    "fluoroquinolone": {"name": "Fluoroquinolone", "color": "#993404", "critical": False},
    "macrolide":     {"name": "Macrolide", "color": "#54278f", "critical": False},
    "tetracycline":  {"name": "Tetracycline", "color": "#8c6d31", "critical": False},
    "phenicol":      {"name": "Phenicol", "color": "#525252", "critical": False},
    "sulfonamide":   {"name": "Sulfonamide", "color": "#a63603", "critical": False},
    "trimethoprim":  {"name": "Trimethoprim", "color": "#6a51a3", "critical": False},
    "fosfomycin":    {"name": "Fosfomycin", "color": "#238b45", "critical": False},
    "rifamycin":     {"name": "Rifamycin", "color": "#ce1256", "critical": False},
    "other":         {"name": "Other / unclassified", "color": "#999999", "critical": False},
}

CRITICAL_CLASSES = [k for k, v in DRUG_CLASSES.items() if v["critical"]]

# Carbapenemase β-lactamase families (critical). Matched against the bla allele.
_CARBAPENEMASE = re.compile(
    r"^bla(kpc|ndm|vim|imp|spm|sim|gim|dim|bic|aim|smb|tmb|fim)",
    re.I,
)
# OXA carbapenemases: the intrinsic/acquired carbapenemase OXA families
# (48-like in Enterobacterales; 23/24/40/58/143/235 in Acinetobacter; 51 intrinsic).
_OXA_CARBAPENEMASE = re.compile(
    r"^blaoxa-?(23|24|40|48|51|58|72|143|181|204|232|235|244|436)",
    re.I,
)


def classify_gene(symbol: str) -> str:
    """Map an AMRFinderPlus gene symbol to a drug-class key (see DRUG_CLASSES)."""
    s = symbol.strip().strip('"').strip()
    if not s:
        return "other"
    low = s.lower()

    # β-lactams: split the critical mechanisms out of the broad bla* bucket.
    if low.startswith("bla"):
        if _CARBAPENEMASE.match(low) or _OXA_CARBAPENEMASE.match(low):
            return "carbapenem"
        if low.startswith("blactx"):
            return "esbl"
        return "beta_lactam"

    if low.startswith("mcr-") or low.startswith("mcr.") or low.startswith("mcr_") or low == "mcr":
        return "colistin"
    if low.startswith("van"):
        return "glycopeptide"
    if low.startswith(("aac", "aad", "aph", "ant", "arm", "rmt", "str", "spc", "aph3", "aph6")):
        return "aminoglycoside"
    if low.startswith(("qnr", "oqx", "qep", "qnrvc")) or low == "aac(6')-ib-cr":
        return "fluoroquinolone"
    # QRDR point mutations reported by AMRFinderPlus as <gene>_<AAsub>, e.g.
    # gyrA_S83I, parC_S85L — the dominant fluoroquinolone-resistance mechanism.
    if low.startswith(("gyra_", "gyrb_", "parc_", "pare_")):
        return "fluoroquinolone"
    if low.startswith(("erm", "mph", "mef", "ere", "msr", "lnu", "lsa", "vga", "vat", "car")):
        return "macrolide"
    if low.startswith("tet"):
        return "tetracycline"
    if low.startswith(("cat", "cml", "flor", "fex", "cfr", "pp-flo", "estdl")):
        return "phenicol"
    if low.startswith("sul"):
        return "sulfonamide"
    if low.startswith("dfr"):
        return "trimethoprim"
    if low.startswith("fos"):
        return "fosfomycin"
    if low.startswith(("arr", "rph")):
        return "rifamycin"
    return "other"


def is_critical_gene(symbol: str) -> bool:
    return classify_gene(symbol) in CRITICAL_CLASSES


# --- geo_loc_name -> ISO3 (Natural Earth join) --------------------------------
# NCBI/INSDC country vocabulary -> ISO_A3. Covers the high-surveillance
# contributors and common names; unmatched countries are dropped from the map
# layer (still counted in global totals) and the drop count is logged/reported.
COUNTRY_ISO3: dict[str, str] = {
    "usa": "USA", "united states": "USA", "united states of america": "USA",
    "united kingdom": "GBR", "uk": "GBR", "england": "GBR", "scotland": "GBR",
    "wales": "GBR", "northern ireland": "GBR",
    "canada": "CAN", "mexico": "MEX", "brazil": "BRA", "argentina": "ARG",
    "chile": "CHL", "colombia": "COL", "peru": "PER", "ecuador": "ECU",
    "china": "CHN", "japan": "JPN", "south korea": "KOR", "korea": "KOR",
    "republic of korea": "KOR", "north korea": "PRK", "india": "IND",
    "pakistan": "PAK", "bangladesh": "BGD", "nepal": "NPL", "sri lanka": "LKA",
    "thailand": "THA", "viet nam": "VNM", "vietnam": "VNM", "cambodia": "KHM",
    "laos": "LAO", "myanmar": "MMR", "malaysia": "MYS", "singapore": "SGP",
    "indonesia": "IDN", "philippines": "PHL", "taiwan": "TWN", "hong kong": "HKG",
    "australia": "AUS", "new zealand": "NZL",
    "germany": "DEU", "france": "FRA", "italy": "ITA", "spain": "ESP",
    "portugal": "PRT", "netherlands": "NLD", "belgium": "BEL", "switzerland": "CHE",
    "austria": "AUT", "denmark": "DNK", "sweden": "SWE", "norway": "NOR",
    "finland": "FIN", "iceland": "ISL", "ireland": "IRL", "poland": "POL",
    "czech republic": "CZE", "czechia": "CZE", "slovakia": "SVK", "hungary": "HUN",
    "romania": "ROU", "bulgaria": "BGR", "greece": "GRC", "croatia": "HRV",
    "serbia": "SRB", "slovenia": "SVN", "russia": "RUS", "russian federation": "RUS",
    "ukraine": "UKR", "belarus": "BLR", "lithuania": "LTU", "latvia": "LVA",
    "estonia": "EST", "turkey": "TUR", "cyprus": "CYP", "luxembourg": "LUX",
    "israel": "ISR", "saudi arabia": "SAU", "united arab emirates": "ARE",
    "qatar": "QAT", "kuwait": "KWT", "oman": "OMN", "iran": "IRN", "iraq": "IRQ",
    "jordan": "JOR", "lebanon": "LBN", "egypt": "EGY", "morocco": "MAR",
    "tunisia": "TUN", "algeria": "DZA", "libya": "LBY",
    "south africa": "ZAF", "nigeria": "NGA", "kenya": "KEN", "tanzania": "TZA",
    "uganda": "UGA", "ghana": "GHA", "ethiopia": "ETH", "senegal": "SEN",
    "cameroon": "CMR", "cote d'ivoire": "CIV", "ivory coast": "CIV",
    "democratic republic of the congo": "COD", "zimbabwe": "ZWE", "zambia": "ZMB",
    "malawi": "MWI", "mozambique": "MOZ", "rwanda": "RWA", "botswana": "BWA",
    "mali": "MLI", "burkina faso": "BFA", "niger": "NER", "sudan": "SDN",
    "somalia": "SOM", "madagascar": "MDG", "gambia": "GMB",
    "haiti": "HTI", "dominican republic": "DOM", "cuba": "CUB", "jamaica": "JAM",
    "guatemala": "GTM", "honduras": "HND", "el salvador": "SLV", "nicaragua": "NIC",
    "costa rica": "CRI", "panama": "PAN", "bolivia": "BOL", "paraguay": "PRY",
    "uruguay": "URY", "venezuela": "VEN",
    "kazakhstan": "KAZ", "uzbekistan": "UZB", "georgia": "GEO", "armenia": "ARM",
    "azerbaijan": "AZE", "mongolia": "MNG", "afghanistan": "AFG", "yemen": "YEM",
}


def country_to_iso3(geo_loc_name: str) -> str | None:
    """Parse ``geo_loc_name`` ("Country" or "Country: region") to ISO_A3.

    Returns None when the value is empty/NULL or the country isn't in the map
    (caller counts these toward a global total but omits them from the map).
    """
    if not geo_loc_name:
        return None
    v = geo_loc_name.strip().strip('"').strip()
    if not v or v.upper() == "NULL" or v.lower() in ("not applicable", "missing", "unknown", "not collected"):
        return None
    country = v.split(":", 1)[0].strip().lower()
    return COUNTRY_ISO3.get(country)
