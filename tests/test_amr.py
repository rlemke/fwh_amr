"""Tests for the AMR domain — offline (no network; requests is mocked)."""

from __future__ import annotations

import json

import pytest

from amr import organisms as org


# --------------------------------------------------------------------------
# gene classification
# --------------------------------------------------------------------------


@pytest.mark.parametrize("sym,cls", [
    ("blaKPC-2", "carbapenem"),
    ("blaNDM-5", "carbapenem"),
    ("blaVIM-1", "carbapenem"),
    ("blaOXA-48", "carbapenem"),
    ("blaOXA-23", "carbapenem"),       # Acinetobacter carbapenemase
    ("blaOXA-61", "beta_lactam"),      # NOT a carbapenemase OXA
    ("blaCTX-M-15", "esbl"),
    ("blaTEM-1", "beta_lactam"),
    ("blaSHV-11", "beta_lactam"),
    ("mcr-1", "colistin"),
    ("mcr-9.1", "colistin"),
    ("vanA", "glycopeptide"),
    ("aac(6')-Ib", "aminoglycoside"),
    ("aph(3')-IIIa", "aminoglycoside"),
    ("armA", "aminoglycoside"),
    ("tet(O)", "tetracycline"),
    ("tetM", "tetracycline"),
    ("ermB", "macrolide"),
    ("mph(A)", "macrolide"),
    ("sul1", "sulfonamide"),
    ("dfrA17", "trimethoprim"),
    ("qnrS1", "fluoroquinolone"),
    ("gyrA_S83I", "fluoroquinolone"),
    ("parC_S85L", "fluoroquinolone"),
    ("fosA7", "fosfomycin"),
    ("floR", "phenicol"),
    ("catA1", "phenicol"),
    ("something_unknown", "other"),
])
def test_classify_gene(sym, cls):
    assert org.classify_gene(sym) == cls


def test_critical_flag():
    assert org.is_critical_gene("blaNDM-1")
    assert org.is_critical_gene("mcr-1")
    assert org.is_critical_gene("vanA")
    assert not org.is_critical_gene("tet(O)")
    assert set(org.CRITICAL_CLASSES) == {"carbapenem", "colistin", "esbl", "glycopeptide"}


# --------------------------------------------------------------------------
# country normalization
# --------------------------------------------------------------------------


@pytest.mark.parametrize("geo,iso", [
    ("USA: California", "USA"),
    ("USA", "USA"),
    ("United Kingdom", "GBR"),
    ("United Kingdom: England", "GBR"),
    ("Viet Nam", "VNM"),
    ("China: Beijing", "CHN"),
    ("NULL", None),
    ("", None),
    ("missing", None),
    ("Atlantis", None),
])
def test_country_to_iso3(geo, iso):
    assert org.country_to_iso3(geo) == iso


# --------------------------------------------------------------------------
# parse helpers
# --------------------------------------------------------------------------


def test_parse_year():
    from amr._lib import _parse_year
    assert _parse_year("2015-03-01") == 2015
    assert _parse_year("2015") == 2015
    assert _parse_year("NULL") is None
    assert _parse_year("") is None
    assert _parse_year("19") is None


def test_parse_genes():
    from amr._lib import _parse_genes
    assert _parse_genes('"blaOXA-61,tet(O)"') == ["blaOXA-61", "tet(O)"]
    assert _parse_genes('"blaKPC-2=COMPLETE,mcr-1=PARTIAL"') == ["blaKPC-2", "mcr-1"]
    assert _parse_genes("NULL") == []
    assert _parse_genes("") == []
    assert _parse_genes('"vanA"') == ["vanA"]


# --------------------------------------------------------------------------
# offline download (mocked requests) + fan-in builds
# --------------------------------------------------------------------------


_TSV = "\t".join([
    "collection_date", "geo_loc_name", "scientific_name", "AMR_genotypes",
])
_TSV_ROWS = [
    _TSV,
    "\t".join(["2016-05-01", "USA: California", "Vibrio cholerae", '"blaCTX-M-15,tet(O)"']),
    "\t".join(["2018", "India", "Vibrio cholerae", '"blaNDM-5,mcr-1,sul1"']),
    "\t".join(["NULL", "USA", "Vibrio cholerae", '"blaTEM-1"']),
    "\t".join(["2019-01-01", "USA: Texas", "Vibrio cholerae", "NULL"]),
    "\t".join(["2019-06-01", "Atlantis", "Vibrio cholerae", '"vanA"']),
]


class _FakeResp:
    def __init__(self, *, text="", lines=None):
        self.text = text
        self._lines = lines or []
        self.encoding = "utf-8"

    def raise_for_status(self):
        pass

    def iter_lines(self, decode_unicode=False):
        yield from self._lines


class _FakeRequests:
    def get(self, url, **kw):
        if url.endswith("/Metadata/"):
            return _FakeResp(text='<a href="PDG000000000.1.metadata.tsv">x</a>')
        if url.endswith(".metadata.tsv"):
            return _FakeResp(lines=_TSV_ROWS)
        raise AssertionError(f"unexpected URL {url}")


@pytest.fixture()
def cache_env(tmp_path, monkeypatch):
    monkeypatch.setenv("FW_AMR_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setenv("FW_AMR_OUTPUT_DIR", str(tmp_path / "out"))
    return tmp_path


def test_download_organism_offline(cache_env, monkeypatch):
    from amr import _lib
    monkeypatch.setattr(_lib, "requests", _FakeRequests())
    res = _lib.download_organism("Vibrio_cholerae", force=True)
    assert res.n_isolates == 5
    assert res.n_dated == 4                 # 4 rows have a parseable year
    assert res.n_geo_matched == 4           # USA x3, India x1 (Atlantis unmatched)
    assert res.n_geo_unmatched == 1         # Atlantis
    # genes across all rows: blaCTX-M-15, tet(O), blaNDM-5, mcr-1, sul1, blaTEM-1, vanA
    assert res.n_genes == 7

    from amr import storage as st
    agg_path = st.join(st.organisms_root(), "Vibrio_cholerae.json")
    with st.open_read(agg_path) as f:
        blob = json.load(f)
    # India 2018 isolate carries a carbapenemase + colistin
    assert blob["prevalence"]["IND"]["2018"][1]["carbapenem"] == 1
    assert blob["prevalence"]["IND"]["2018"][1]["colistin"] == 1
    assert blob["genes"]["mcr-1"]["class"] == "colistin"
    assert blob["genes"]["mcr-1"]["first_year"] == 2018


def _write_world(cache_env):
    from amr import storage as st
    world = {"type": "FeatureCollection", "features": [
        {"type": "Feature",
         "properties": {"NAME": "United States", "ISO_A3": "USA"},
         "geometry": {"type": "Polygon", "coordinates": [[[-100, 30], [-100, 40], [-90, 40], [-90, 30], [-100, 30]]]}},
        {"type": "Feature",
         "properties": {"NAME": "India", "ISO_A3": "IND"},
         "geometry": {"type": "Polygon", "coordinates": [[[70, 10], [70, 30], [88, 30], [88, 10], [70, 10]]]}},
    ]}
    with st.open_write(st.join(st.cache_root(), "world-countries.geojson"), "w") as f:
        json.dump(world, f)


def test_build_map_and_table(cache_env, monkeypatch):
    from amr import _lib
    monkeypatch.setattr(_lib, "requests", _FakeRequests())
    _lib.download_organism("Vibrio_cholerae", force=True)
    _write_world(cache_env)

    m = _lib.build_resistance_map()
    assert m.organism_count == 1
    assert m.country_count == 2            # USA + IND have in-window dated isolates
    with open(m.html_path, encoding="utf-8") as f:
        html = f.read()
    assert "About this data" in html
    assert "surveillance-biased" in html
    assert "carbapenem" in html

    t = _lib.build_gene_species_table()
    assert t.row_count == 7
    with open(t.html_path, encoding="utf-8") as f:
        thtml = f.read()
    assert "Vibrio cholerae" in thtml
    assert "mcr-1" in thtml
    assert "co-occurrence" in thtml.lower() or "not proof of transfer" in thtml.lower()


def test_build_requires_aggregates(cache_env):
    from amr import _lib
    with pytest.raises(RuntimeError, match="no organism aggregates"):
        _lib.build_resistance_map()
