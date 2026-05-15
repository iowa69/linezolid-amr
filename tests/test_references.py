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
    }
    assert set(ref_mod.list_organisms()) == expected


def test_g2576_present():
    positions = ref_mod.get_linezolid_positions()
    g2576 = [p for p in positions if p.ecoli_position == 2576]
    assert len(g2576) == 1
    assert g2576[0].ref_base == "G"
    assert "T" in g2576[0].resistance_bases


def test_canonical_positions_complete():
    expected_positions = {2032, 2447, 2453, 2500, 2503, 2504, 2505, 2534, 2572, 2576, 2603}
    found = {p.ecoli_position for p in ref_mod.get_linezolid_positions()}
    assert expected_positions <= found, f"missing positions: {expected_positions - found}"


def test_every_position_has_pmid():
    data = ref_mod.load_loci()
    for item in data["linezolid_positions_ecoli_23s"]:
        pmids = item.get("pmid", [])
        assert pmids, f"position {item['ecoli_position']} has no PMID citation"
        for p in pmids:
            assert p.isdigit() and len(p) >= 7, f"invalid PMID {p}"


def test_mtb_removed():
    """M. tuberculosis was dropped in 0.1.2 — covered by domain-specific tools."""
    assert "Mycobacterium_tuberculosis" not in ref_mod.list_organisms()


def test_organism_metadata_consistent():
    for organism in ref_mod.list_organisms():
        o = ref_mod.get_organism(organism)
        length = o.end - o.start + 1
        assert 2800 <= length <= 3000, f"{organism}: unexpected length {length}"
        assert o.strand in ("+", "-")
        assert o.genome_accession.startswith("NC_")


def test_unknown_organism_raises():
    with pytest.raises(KeyError):
        ref_mod.get_organism("Escherichia_coli")


def test_cache_dir_respects_env(monkeypatch, tmp_path):
    monkeypatch.setenv("LINEZOLID_AMR_REFDIR", str(tmp_path))
    assert ref_mod.cache_dir() == tmp_path


def test_bundled_references_present(monkeypatch, tmp_path):
    """Every supported organism must have a bundled FASTA + BED + position map."""
    # Force the resolver to skip any user cache, falling back to bundled package data.
    monkeypatch.setenv("LINEZOLID_AMR_REFDIR", str(tmp_path / "no-such-dir"))
    for organism in ref_mod.list_organisms():
        fasta = ref_mod.organism_fasta_path(organism)
        bed = ref_mod.organism_bed_path(organism)
        pmap = ref_mod.organism_position_map_path(organism)
        assert fasta.exists(), f"missing bundled FASTA: {organism}"
        assert bed.exists(), f"missing bundled BED: {organism}"
        assert pmap.exists(), f"missing bundled position map: {organism}"
        # Sanity: fasta header + ≥2800 bp of sequence
        head = fasta.read_text().splitlines()
        assert head[0].startswith(">"), f"bad FASTA header: {fasta}"
        seq_len = sum(len(l) for l in head[1:] if not l.startswith(">"))
        assert seq_len > 2800, f"{organism} 23S length {seq_len} too short"


def test_ecoli_master_bundled(monkeypatch, tmp_path):
    monkeypatch.setenv("LINEZOLID_AMR_REFDIR", str(tmp_path / "no-such-dir"))
    assert ref_mod.ecoli_fasta_path().exists()
