"""Unit tests for summary.py wide/long CSV builders."""

from __future__ import annotations

from pathlib import Path

from linezolid_amr.amrfinder import AmrHit
from linezolid_amr.rrna23s import PileupCall
from linezolid_amr import summary as summary_mod


def _h(symbol, klass, et="AMR", cov=100.0, ident=99.5):
    return AmrHit(
        contig="ctg", start=0, end=10, strand="+",
        gene_symbol=symbol, sequence_name=symbol, scope="core",
        element_type=et, element_subtype="AMR", class_=klass, subclass=klass,
        method="EXACTX", target_length=10, reference_length=10,
        coverage_pct=cov, identity_pct=ident, accession="NG_x",
    )


def _p(ecoli_pos, ref, alt_base, alt_af, depth, drug="linezolid"):
    return PileupCall(
        organism="Enterococcus_faecium",
        ref_contig="Enterococcus_faecium_23S",
        species_position=ecoli_pos + 14,
        ecoli_position=ecoli_pos,
        ref_base=ref, depth=depth,
        counts={"A": 0, "C": 0, "G": int(depth * (1 - alt_af)), "T": int(depth * alt_af),
                "N": 0, "DEL": 0},
        alt_alleles=[{"base": alt_base, "count": int(depth * alt_af), "af": alt_af, "resistance": True}],
        is_resistance=True, drug=drug,
    )


def test_wide_row_lzd_first_columns():
    amr = [_h("cfr(D)", "OXAZOLIDINONE"), _h("vanA", "GLYCOPEPTIDE"),
           _h("23S_G2576T", "OXAZOLIDINONE")]
    pile = [_p(2576, "G", "T", 0.66, 173)]
    row = summary_mod.build_wide_row("S1", "Enterococcus_faecium", "17", "efaecium", amr, pile, True)
    assert row["linezolid_call"] == "POS"
    assert "G2576T" in row["lzd_23S_mutations"]
    assert "cfr(D)" in row["lzd_genes"]
    # AMR_GLYCOPEPTIDE bucket present
    assert any(k.startswith("AMR_GLYCOPEPTIDE") for k in row)


def test_wide_neg_call():
    row = summary_mod.build_wide_row("S2", "Staphylococcus_aureus", "5", "saureus",
                                      [_h("blaZ", "BETA-LACTAM")], [], False)
    assert row["linezolid_call"] == "neg"
    assert row["lzd_23S_mutations"] == ""
    assert row["lzd_genes"] == ""


def test_long_rows_include_both():
    amr = [_h("cfr(D)", "OXAZOLIDINONE"), _h("vanA", "GLYCOPEPTIDE")]
    pile = [_p(2576, "G", "T", 0.66, 173)]
    rows = summary_mod.build_long_rows("S1", "Enterococcus_faecium", "17", "efaecium", amr, pile)
    kinds = {r["feature_kind"] for r in rows}
    assert "23S_LZD_mutation" in kinds
    # AmrHit is captured by its element_type, "AMR"
    assert "AMR" in kinds
    features = {r["feature"] for r in rows}
    assert "G2576T" in features
    assert "cfr(D)" in features
    assert "vanA" in features


def test_csv_writers_roundtrip(tmp_path: Path):
    amr = [_h("cfr(D)", "OXAZOLIDINONE")]
    pile = [_p(2576, "G", "T", 0.66, 173)]
    wide = summary_mod.build_wide_row("S1", "Enterococcus_faecium", "17", "efaecium", amr, pile, True)
    long = summary_mod.build_long_rows("S1", "Enterococcus_faecium", "17", "efaecium", amr, pile)
    summary_mod.write_wide_csv([wide], tmp_path / "wide.csv")
    summary_mod.write_long_csv(long, tmp_path / "long.csv")
    text_wide = (tmp_path / "wide.csv").read_text()
    text_long = (tmp_path / "long.csv").read_text()
    assert "linezolid_call" in text_wide.splitlines()[0]
    assert "feature_kind" in text_long.splitlines()[0]
    assert "POS" in text_wide
    assert "G2576T" in text_long
