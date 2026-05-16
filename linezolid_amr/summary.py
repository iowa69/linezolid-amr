"""Stats-friendly CSV summaries.

Wide CSV: one row per sample, one column per detected gene/mutation. Column
names are the bare gene symbol (e.g. ``blaZ``, ``cfr(D)``, ``23S_G2576T``).
Linezolid 23S read-level mutations are reported only when AF ≥ ``--min-af``
(default 0.15); the AF cell holds the raw fraction. The long CSV keeps the
full per-feature view including sub-threshold AFs for transparency.
"""

from __future__ import annotations

import csv
from pathlib import Path

from linezolid_amr.amrfinder import AmrHit
from linezolid_amr.rrna23s import PileupCall


# -------------------------- helpers --------------------------

def _amr_column(h: AmrHit) -> str:
    """Wide-CSV column = the bare gene symbol. Class/element type are not
    encoded into the header (the AMRFinderPlus TSV carries that detail)."""
    return h.gene_symbol or h.sequence_name or "unknown"


def _passing_resistance_alts(p: PileupCall) -> list[tuple[str, float]]:
    """Only resistance-base alts that cleared the threshold get into the wide CSV."""
    return [
        (f"{p.ref_base}{p.ecoli_position}{a['base']}", a["af"])
        for a in p.alt_alleles
        if a["resistance"] and a.get("passes_threshold")
    ]


def _all_resistance_alts(p: PileupCall) -> list[tuple[str, float]]:
    """Every observed resistance-base alt (long CSV uses this — transparency)."""
    return [
        (f"{p.ref_base}{p.ecoli_position}{a['base']}", a["af"])
        for a in p.alt_alleles
        if a["resistance"]
    ]


def _amr_value(h: AmrHit) -> str:
    return f"{h.identity_pct:.1f}"


# Linezolid-relevant column ordering: cfr/optrA/poxtA/23S/L3/L4/L22 hits come
# first among the gene columns, then the actual LZD pileup mutations, then
# everything else alphabetical.
def _column_rank(col: str, lzd_amr_cols: set[str], lzd_pileup_cols: set[str]) -> tuple[int, str]:
    if col in lzd_amr_cols:
        return (0, col.lower())
    if col in lzd_pileup_cols:
        return (1, col.lower())
    return (2, col.lower())


# Fixed leading columns — always present in this order.
_LEADING = (
    "sample", "organism", "mlst_scheme", "ST", "mlst_alleles",
    "linezolid_call", "lzd_n_23S_mutations", "lzd_max_23S_af",
)


# -------------------------- builders --------------------------

def build_wide_row(
    sample: str,
    organism: str,
    st: str | None,
    mlst_scheme: str | None,
    mlst_alleles: str | None,
    amr_hits: list[AmrHit],
    pileup_calls: list[PileupCall],
    linezolid_call: bool,
) -> dict[str, str]:
    """Build the wide-CSV row for a single sample.

    Wide-CSV policy: only 23S resistance alleles AT OR ABOVE the threshold
    appear as columns. Sub-threshold AFs are still visible in the long CSV
    and the per-sample pileup TSV.
    """
    n_pos_positions = sum(1 for p in pileup_calls if p.is_resistance)
    passing_afs = [af for p in pileup_calls for _, af in _passing_resistance_alts(p)]
    max_af = max(passing_afs, default=0.0)
    row: dict[str, str] = {
        "sample": sample,
        "organism": organism or "",
        "mlst_scheme": mlst_scheme or "",
        "ST": st or "",
        "mlst_alleles": mlst_alleles or "",
        "linezolid_call": "POS" if linezolid_call else "neg",
        "lzd_n_23S_mutations": str(n_pos_positions),
        "lzd_max_23S_af": f"{max_af:.4f}" if max_af else "",
    }
    # AMRFinderPlus → one column per gene (bare gene name)
    for h in amr_hits:
        col = _amr_column(h)
        new_val = _amr_value(h)
        prev = row.get(col)
        if prev is None or float(new_val) > float(prev):
            row[col] = new_val
    # 23S pileup → only above-threshold mutations
    for p in pileup_calls:
        for mut, af in _passing_resistance_alts(p):
            col = mut
            prev_af = float(row[col]) if col in row else -1.0
            if af > prev_af:
                row[col] = f"{af:.4f}"
    return row


def build_long_rows(
    sample: str,
    organism: str,
    st: str | None,
    mlst_scheme: str | None,
    amr_hits: list[AmrHit],
    pileup_calls: list[PileupCall],
) -> list[dict[str, str]]:
    """Long-format rows: one per detected feature (sub-threshold AFs included)."""
    base = {"sample": sample, "organism": organism, "mlst_scheme": mlst_scheme or "", "ST": st or ""}
    out: list[dict[str, str]] = []
    for p in pileup_calls:
        for a in p.alt_alleles:
            if not a["resistance"]:
                continue
            out.append({
                **base,
                "feature_kind": "23S_LZD_mutation",
                "feature": f"{p.ref_base}{p.ecoli_position}{a['base']}",
                "class": "OXAZOLIDINONE",
                "subclass": "LINEZOLID",
                "evidence": f"E.coli pos {p.ecoli_position}; species pos {p.species_position}",
                "depth": str(p.depth),
                "alt_count": str(a["count"]),
                "alt_af": f"{a['af']:.4f}",
                "passes_threshold": "YES" if a.get("passes_threshold") else "NO",
                "coverage_pct": "",
                "identity_pct": "",
                "contig": p.ref_contig,
            })
    for h in amr_hits:
        out.append({
            **base,
            "feature_kind": h.element_type or "AMR",
            "feature": h.gene_symbol or h.sequence_name,
            "class": h.class_,
            "subclass": h.subclass,
            "evidence": h.method,
            "depth": "",
            "alt_count": "",
            "alt_af": "",
            "coverage_pct": f"{h.coverage_pct:.2f}",
            "identity_pct": f"{h.identity_pct:.2f}",
            "contig": h.contig,
        })
    return out


# -------------------------- writers --------------------------

def write_wide_csv(rows: list[dict[str, str]], path: Path) -> None:
    """Write rows as a CSV. Leading identity columns first, then LZD-relevant
    AMR hits, then 23S pileup mutations, then everything else alphabetical."""
    if not rows:
        path.write_text(",".join(_LEADING) + "\n")
        return
    # Identify LZD-relevant AMR columns vs LZD pileup mutation columns vs others
    lzd_amr_cols: set[str] = set()
    lzd_pileup_cols: set[str] = set()
    # LZD AMR hits are those whose AmrHit was flagged linezolid-relevant; we can't
    # introspect AmrHit from rows alone, so apply a lightweight rule on column name.
    LZD_GENE_TOKENS = ("cfr", "optr", "poxt", "23s_", "rplc", "rpld", "rplv")
    for r in rows:
        for k in r:
            if k in _LEADING:
                continue
            k_low = k.lower()
            # 23S pileup mutation pattern: <BASE><digits><BASE>
            if len(k) >= 4 and k[0] in "ACGT" and k[-1] in "ACGT" and k[1:-1].isdigit():
                lzd_pileup_cols.add(k)
            elif any(tok in k_low for tok in LZD_GENE_TOKENS):
                lzd_amr_cols.add(k)
    extras = sorted(
        {k for r in rows for k in r.keys() if k not in _LEADING},
        key=lambda c: _column_rank(c, lzd_amr_cols, lzd_pileup_cols),
    )
    fieldnames = list(_LEADING) + extras
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def write_long_csv(rows: list[dict[str, str]], path: Path) -> None:
    fieldnames = [
        "sample", "organism", "mlst_scheme", "ST",
        "feature_kind", "feature", "class", "subclass", "evidence",
        "depth", "alt_count", "alt_af",
        "coverage_pct", "identity_pct", "contig",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
