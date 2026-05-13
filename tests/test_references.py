"""Unit tests for references.py and the bundled loci.json."""

from __future__ import annotations

import pytest

from linezolid_amr import references as ref_mod


def test_loci_loads():
    data = ref_mod.load_loci()
    assert "linezolid_positions_ecoli_23s" in data
    assert "organisms" in data
    assert "ecoli_reference" in data


def test_supported_organisms_exact():
    expected = {
        "Staphylococcus_aureus",
        "Enterococcus_faecalis",
        "Enterococcus_faecium",
        "Streptococcus_pneumoniae",
        "Mycobacterium_tuberculosis",
    }
    assert set(ref_mod.list_organisms()) == expected


def test_g2576_present():
    positions = ref_mod.get_linezolid_positions()
    g2576 = [p for p in positions if p.ecoli_position == 2576]
    assert len(g2576) == 1
    assert g2576[0].ref_base == "G"
    assert "T" in g2576[0].resistance_bases


def test_canonical_positions_complete():
    """All canonical E. coli linezolid positions must be present."""
    expected_positions = {2032, 2447, 2453, 2500, 2503, 2504, 2505, 2534, 2572, 2576, 2603}
    found = {p.ecoli_position for p in ref_mod.get_linezolid_positions()}
    assert expected_positions <= found, f"missing positions: {expected_positions - found}"


def test_organism_metadata_consistent():
    for organism in ref_mod.list_organisms():
        o = ref_mod.get_organism(organism)
        length = o.end - o.start + 1
        # Should be within a reasonable bacterial 23S length window
        assert 2800 <= length <= 3300, f"{organism}: unexpected length {length}"
        assert o.strand in ("+", "-")
        assert o.genome_accession.startswith("NC_")


def test_unknown_organism_raises():
    with pytest.raises(KeyError):
        ref_mod.get_organism("Escherichia_coli")


def test_cache_dir_respects_env(monkeypatch, tmp_path):
    monkeypatch.setenv("LINEZOLID_AMR_REFDIR", str(tmp_path))
    assert ref_mod.cache_dir() == tmp_path
