"""Stats-friendly CSV summaries.

We emit two CSV files per sample (and aggregated for ``folder`` mode):

* ``<sample>.summary_wide.csv`` — **one row per sample**, where every detected
  gene / mutation gets **its own column**. Column names are prefixed by group:

  ============================================  ============================
  Column                                        Meaning
  ============================================  ============================
  ``sample, organism, mlst_scheme, ST,``        Identity & MLST
  ``mlst_alleles, linezolid_call,``
  ``lzd_n_23S_mutations, lzd_max_23S_af``
  ``LZD__<gene_or_mutation>``                   Linezolid-relevant hit;
                                                value = identity (%)
  ``LZD_23S_AF__<ecoli_pos>``                   23S pileup AF at LZD site;
                                                value = AF in [0,1]
  ``AMR_<CLASS>__<gene>``                       AMR hit; value = identity (%)
  ``VIRULENCE__<gene>``                         Virulence hit; value = id (%)
  ``STRESS_<CLASS>__<gene>``                    Stress hit; value = id (%)
  ============================================  ============================

  Folder mode (``ALL_samples.summary_wide.csv``) takes the union of every
  detected column across the cohort so users get a clean sample × gene matrix
  ready for pandas / R / Excel pivot tables.

* ``<sample>.summary_long.csv`` — long format, one row per detected feature.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

from linezolid_amr.amrfinder import AmrHit, LINEZOLID_GENE_KEYWORDS, LINEZOLID_PROTEIN_KEYWORDS
from linezolid_amr.rrna23s import PileupCall


# -------------------------- column grouping helpers --------------------------

def _is_linezolid_relevant(h: AmrHit) -> bool:
    return h.is_linezolid_relevant


def _amr_column(h: AmrHit) -> str:
    """Return the ``<group>__<gene>`` column name for an AMRFinderPlus hit."""
    gene = h.gene_symbol or h.sequence_name or "unknown"
    cls = (h.class_ or "OTHER").upper().replace(" ", "_").replace("/", "_")
    et = (h.element_type or "AMR").upper()
    if _is_linezolid_relevant(h):
        # Pool every linezolid-relevant hit (cfr, optrA, poxtA, 23S_*, L3/L4/L22, …)
        # into a single LZD group so the linezolid block in the header is compact.
        return f"LZD__{gene}"
    if et == "VIRULENCE":
        return f"VIRULENCE__{gene}"
    if et == "STRESS":
        return f"STRESS_{cls}__{gene}"
    return f"AMR_{cls}__{gene}"


def _pileup_resistance_alts(p: PileupCall) -> list[tuple[str, float]]:
    """All (mutation_name, AF) pairs for resistance-known alts at this position,
    regardless of threshold (always-report-the-proportion policy).
    """
    out = []
    for a in p.alt_alleles:
        if not a["resistance"]:
            continue
        mut = f"{p.ref_base}{p.ecoli_position}{a['base']}"
        out.append((mut, a["af"]))
    return out


def _amr_value(h: AmrHit) -> str:
    """Cell value for an AMRFinderPlus hit — identity % formatted to 1 decimal."""
    return f"{h.identity_pct:.1f}"


# ---------------------------- main builders ----------------------------

# Fixed leading columns — always present in this order.
_LEADING = (
    "sample", "organism", "mlst_scheme", "ST", "mlst_alleles",
    "linezolid_call", "lzd_n_23S_mutations", "lzd_max_23S_af",
)


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

    Linezolid policy: ``linezolid_call=POS`` iff any resistance allele is at
    or above the min-AF threshold; *but* every observed resistance allele
    (even sub-threshold) gets its own ``LZD_23S_AF__<mutation>`` column with
    the raw AF so the reader can judge borderline cases.
    """
    n_pos_positions = sum(1 for p in pileup_calls if p.is_resistance)
    # Max AF over ALL observed resistance alleles (above-threshold or not),
    # because the user explicitly asked to "always show the proportion".
    all_res_afs = [af for p in pileup_calls for _, af in _pileup_resistance_alts(p)]
    max_af = max(all_res_afs, default=0.0)
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
    # AMRFinderPlus → one column per gene
    for h in amr_hits:
        col = _amr_column(h)
        prev = row.get(col)
        new_val = _amr_value(h)
        if prev is None or float(new_val) > float(prev):
            row[col] = new_val
    # 23S pileup → one column per observed resistance allele (always reported)
    for p in pileup_calls:
        for mut, af in _pileup_resistance_alts(p):
            col = f"LZD_23S_AF__{mut}"
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
    """Long-format rows: one per detected feature."""
    base = {"sample": sample, "organism": organism, "mlst_scheme": mlst_scheme or "", "ST": st or ""}
    out: list[dict[str, str]] = []
    # Emit every resistance-base alt seen at any canonical LZD position, even if
    # below the AF threshold (always-show-the-proportion). Threshold passage is
    # captured by the `passes_threshold` column so users can filter trivially.
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


# ---------------------------- writers ----------------------------

def _group_rank(col: str) -> tuple[int, str, str]:
    """Sort key that keeps groups together and orderly.

    Order: LZD → LZD_23S_AF → AMR_* (alphabetised by class) → VIRULENCE → STRESS_*.
    Within each group, genes are sorted alphabetically.
    """
    if "__" not in col:
        return (-1, col, "")
    prefix, gene = col.split("__", 1)
    bucket_order = [
        "LZD",
        "LZD_23S_AF",
    ]
    if prefix in bucket_order:
        return (bucket_order.index(prefix), prefix, gene.lower())
    if prefix.startswith("AMR_"):
        return (10, prefix, gene.lower())
    if prefix == "VIRULENCE":
        return (20, prefix, gene.lower())
    if prefix.startswith("STRESS_"):
        return (30, prefix, gene.lower())
    return (40, prefix, gene.lower())


def write_wide_csv(rows: list[dict[str, str]], path: Path) -> None:
    """Write rows as a CSV with union-of-keys columns. Identity prefix first."""
    if not rows:
        path.write_text(",".join(_LEADING) + "\n")
        return
    extras = sorted({k for r in rows for k in r.keys() if k not in _LEADING}, key=_group_rank)
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
