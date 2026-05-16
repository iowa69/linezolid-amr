"""Unit tests for in-house MLST (linezolid_amr.internal_mlst)."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from linezolid_amr import internal_mlst as im
from linezolid_amr import mlst as legacy  # backwards-compat re-export


def test_pubmlst_organism_mapping_complete():
    expected = {
        "Staphylococcus_aureus",
        "Enterococcus_faecalis",
        "Enterococcus_faecium",
        "Streptococcus_pneumoniae",
    }
    assert set(im.PUBMLST_DBS.keys()) == expected


def test_each_scheme_has_a_seven_locus_layout():
    for org, info in im.PUBMLST_DBS.items():
        sdir = im.scheme_dir(info["name"])
        assert sdir.exists(), f"{org}: scheme dir missing"
        loci = (sdir / "loci.txt").read_text().splitlines()
        loci = [l.strip() for l in loci if l.strip()]
        assert len(loci) == 7, f"{org}: expected 7 loci, got {len(loci)}"


def test_each_scheme_has_alleles_and_profile():
    for org, info in im.PUBMLST_DBS.items():
        sdir = im.scheme_dir(info["name"])
        loci = (sdir / "loci.txt").read_text().splitlines()
        for locus in loci:
            if not locus.strip():
                continue
            fa = sdir / f"{locus}.fasta.gz"
            assert fa.exists(), f"{org}/{locus}: allele FASTA missing"
            # At least one allele header
            with gzip.open(fa, "rt") as fh:
                first = fh.readline()
                assert first.startswith(">"), f"{org}/{locus}: bad FASTA"
        # Profile table has ST + the 7 loci as columns
        header = (sdir / "profiles.tsv").read_text().splitlines()[0].split("\t")
        assert header[0] == "ST"
        for locus in loci:
            assert locus in header, f"{org}: locus {locus} missing from profiles"


def test_all_schemes_available_says_yes():
    assert im.all_schemes_available() is True


def test_mlst_result_organism_property():
    r = im.MlstResult(file="x.fa", scheme="efaecium", st="80",
                       alleles={"atpA": "1"}, raw="")
    assert r.organism == "Enterococcus_faecium"


def test_legacy_alias_reexports():
    """linezolid_amr.mlst is a back-compat alias for internal_mlst."""
    assert legacy.MlstResult is im.MlstResult
    assert legacy.run_mlst is im.run_internal_mlst
    assert legacy.MlstNoMatch is im.MlstNoMatch


def test_unsupported_organism_raises():
    with pytest.raises(im.MlstUnsupportedOrganism):
        im.run_internal_mlst(Path("/nonexistent.fa"), "Escherichia_coli")


def test_null_locus_matches_pubmlst_zero():
    """E. faecium ST1478 has pstS=0 (deleted). Our '-' (no hit) must match
    that profile entry so the ST gets assigned correctly."""
    sdir = im.scheme_dir("efaecium")
    loci, profiles = im._load_profiles(sdir)
    # ST1478 profile: atpA=9, ddl=1, gdh=1, purK=1, gyd=1, pstS=0, adk=1
    # Simulate a sample where pstS is a genuine deletion -> '-'
    locus_to_call = {"atpA": "9", "ddl": "1", "gdh": "1", "purK": "1",
                     "gyd": "1", "pstS": "-", "adk": "1"}
    tup = tuple(locus_to_call[l] for l in loci)
    st, _ = im._assign_st(tup, loci, profiles)
    assert st == "1478"


def test_partial_allele_still_unassigned():
    """A '~n' partial call must NOT trigger the null-equivalence path."""
    sdir = im.scheme_dir("efaecium")
    loci, profiles = im._load_profiles(sdir)
    locus_to_call = {"atpA": "9", "ddl": "1", "gdh": "1", "purK": "1",
                     "gyd": "1", "pstS": "~1", "adk": "1"}
    tup = tuple(locus_to_call[l] for l in loci)
    st, _ = im._assign_st(tup, loci, profiles)
    assert st == "-"
