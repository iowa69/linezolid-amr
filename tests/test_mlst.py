"""Unit tests for mlst.py (parser + organism mapping)."""

from __future__ import annotations

import pytest

from linezolid_amr import mlst as mlst_mod


def test_scheme_mapping_complete():
    expected_orgs = {
        "Staphylococcus_aureus",
        "Enterococcus_faecalis",
        "Enterococcus_faecium",
        "Streptococcus_pneumoniae",
    }
    assert set(mlst_mod.MLST_SCHEME_TO_ORGANISM.values()) == expected_orgs


def test_unsupported_schemes_known():
    assert "kpneumoniae" in mlst_mod.KNOWN_UNSUPPORTED_SCHEMES
    assert "ecoli" in mlst_mod.KNOWN_UNSUPPORTED_SCHEMES


def test_mlst_result_organism_property():
    r = mlst_mod.MlstResult(file="x.fa", scheme="efaecium", st="17",
                             alleles={"atpA": "1"}, raw="x.fa\tefaecium\t17\tatpA(1)")
    assert r.organism == "Enterococcus_faecium"


def test_mlst_result_unsupported_scheme():
    r = mlst_mod.MlstResult(file="x.fa", scheme="ecoli", st="-",
                             alleles={}, raw="")
    assert r.organism is None


def test_mlst_unsupported_exception_classes():
    assert issubclass(mlst_mod.MlstUnsupportedOrganism, ValueError)
    assert issubclass(mlst_mod.MlstNoMatch, ValueError)
