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


def _p(ecoli_pos, ref, alt_base, alt_af, depth, drug="linezolid", passes=True):
    return PileupCall(
        organism="Enterococcus_faecium",
        ref_contig="Enterococcus_faecium_23S",
        species_position=ecoli_pos + 14,
        ecoli_position=ecoli_pos,
        ref_base=ref, depth=depth,
        counts={"A": 0, "C": 0, "G": int(depth * (1 - alt_af)), "T": int(depth * alt_af),
                "N": 0, "DEL": 0},
        alt_alleles=[{"base": alt_base, "count": int(depth * alt_af),
                       "af": alt_af, "resistance": True, "passes_threshold": passes}],
        is_resistance=passes, drug=drug,
    )


def test_wide_row_bare_gene_names():
    amr = [_h("cfr(D)", "OXAZOLIDINONE"), _h("vanA", "GLYCOPEPTIDE"),
           _h("23S_G2576T", "OXAZOLIDINONE")]
    pile = [_p(2576, "G", "T", 0.66, 173, passes=True)]
    row = summary_mod.build_wide_row(
        "S1", "Enterococcus_faecium", "17", "efaecium",
        "atpA(9)|ddl(1)|gdh(1)|purK(1)|gyd(12)|pstS(1)|adk(1)",
        amr, pile, True,
    )
    assert row["linezolid_call"] == "POS"
    assert row["ST"] == "17"
    # Bare gene names — no LZD__ / AMR_ / VIRULENCE__ prefixes
    assert "cfr(D)" in row
    assert "23S_G2576T" in row
    assert "vanA" in row
    assert not any(k.startswith("LZD__") or k.startswith("AMR_") or k.startswith("VIRULENCE__") for k in row)
    # 23S pileup mutation present (also bare)
    assert "G2576T" in row
    assert float(row["G2576T"]) > 0.6


def test_wide_subthreshold_pileup_dropped():
    """Resistance alleles below --min-af must NOT appear as columns in the wide CSV."""
    pile = [_p(2576, "G", "T", 0.02, 200, passes=False)]  # 2% AF, below 15%
    row = summary_mod.build_wide_row(
        "S2", "Staphylococcus_aureus", "5", "saureus",
        "arcC(1)|aroE(1)|glpF(1)|gmk(1)|pta(1)|tpi(1)|yqiL(1)",
        [_h("blaZ", "BETA-LACTAM")], pile, False,
    )
    assert row["linezolid_call"] == "neg"
    assert "G2576T" not in row
    assert row["lzd_max_23S_af"] == ""
    # AMR gene present (bare)
    assert "blaZ" in row


def test_long_csv_keeps_sub_threshold():
    """Long CSV must keep sub-threshold AFs for transparency."""
    pile = [_p(2576, "G", "T", 0.02, 200, passes=False)]
    rows = summary_mod.build_long_rows("S2", "Staphylococcus_aureus", "5", "saureus",
                                        [], pile)
    assert any(r["feature"] == "G2576T" and r["passes_threshold"] == "NO" for r in rows)


def test_csv_writers_roundtrip(tmp_path: Path):
    amr = [_h("cfr(D)", "OXAZOLIDINONE")]
    pile = [_p(2576, "G", "T", 0.66, 173, passes=True)]
    wide = summary_mod.build_wide_row("S1", "Enterococcus_faecium", "17", "efaecium",
                                       "atpA(9)|ddl(1)|gdh(1)|purK(1)|gyd(12)|pstS(1)|adk(1)",
                                       amr, pile, True)
    long = summary_mod.build_long_rows("S1", "Enterococcus_faecium", "17", "efaecium", amr, pile)
    summary_mod.write_wide_csv([wide], tmp_path / "wide.csv")
    summary_mod.write_long_csv(long, tmp_path / "long.csv")
    header_wide = (tmp_path / "wide.csv").read_text().splitlines()[0]
    assert "linezolid_call" in header_wide
    assert "mlst_alleles" in header_wide
    assert "cfr(D)" in header_wide
    assert "G2576T" in header_wide
    assert "LZD__" not in header_wide  # no prefix
    text_long = (tmp_path / "long.csv").read_text()
    assert "feature_kind" in text_long.splitlines()[0]
    assert "G2576T" in text_long
