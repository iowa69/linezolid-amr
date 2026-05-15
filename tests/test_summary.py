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


def test_wide_row_columns_match_new_format():
    amr = [_h("cfr(D)", "OXAZOLIDINONE"), _h("vanA", "GLYCOPEPTIDE"),
           _h("23S_G2576T", "OXAZOLIDINONE")]
    pile = [_p(2576, "G", "T", 0.66, 173)]
    row = summary_mod.build_wide_row("S1", "Enterococcus_faecium", "17", "efaecium",
                                       "atpA(9)|ddl(1)|gdh(1)|purK(1)|gyd(12)|pstS(1)|adk(1)",
                                       amr, pile, True)
    assert row["linezolid_call"] == "POS"
    assert row["ST"] == "17"
    assert "atpA(9)" in row["mlst_alleles"]
    assert int(row["lzd_n_23S_mutations"]) == 1
    assert float(row["lzd_max_23S_af"]) > 0.6
    # Linezolid genes pooled under LZD__
    assert "LZD__cfr(D)" in row
    assert "LZD__23S_G2576T" in row
    # 23S pileup mutation gets its own column with AF
    assert "LZD_23S_AF__G2576T" in row
    # Glycopeptide stays in its AMR class bucket
    assert any(k.startswith("AMR_GLYCOPEPTIDE__") for k in row)


def test_wide_neg_call():
    row = summary_mod.build_wide_row("S2", "Staphylococcus_aureus", "5", "saureus",
                                      "arcC(1)|aroE(1)|glpF(1)|gmk(1)|pta(1)|tpi(1)|yqiL(1)",
                                      [_h("blaZ", "BETA-LACTAM")], [], False)
    assert row["linezolid_call"] == "neg"
    assert row["lzd_n_23S_mutations"] == "0"
    assert row["lzd_max_23S_af"] == ""
    # No LZD columns when none detected
    assert not any(k.startswith("LZD__") for k in row)
    assert not any(k.startswith("LZD_23S_AF__") for k in row)
    # blaZ ends up in AMR_BETA-LACTAM
    assert any(k.startswith("AMR_BETA-LACTAM__") for k in row)


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
    wide = summary_mod.build_wide_row("S1", "Enterococcus_faecium", "17", "efaecium",
                                       "atpA(9)|ddl(1)|gdh(1)|purK(1)|gyd(12)|pstS(1)|adk(1)",
                                       amr, pile, True)
    long = summary_mod.build_long_rows("S1", "Enterococcus_faecium", "17", "efaecium", amr, pile)
    summary_mod.write_wide_csv([wide], tmp_path / "wide.csv")
    summary_mod.write_long_csv(long, tmp_path / "long.csv")
    text_wide = (tmp_path / "wide.csv").read_text()
    text_long = (tmp_path / "long.csv").read_text()
    # Wide CSV: identity prefix + LZD column + 23S AF column all present
    header_wide = text_wide.splitlines()[0]
    assert "linezolid_call" in header_wide
    assert "mlst_alleles" in header_wide
    assert "LZD__cfr(D)" in header_wide
    assert "LZD_23S_AF__G2576T" in header_wide
    # Cells
    assert "POS" in text_wide
    # Long CSV unchanged
    assert "feature_kind" in text_long.splitlines()[0]
    assert "G2576T" in text_long
