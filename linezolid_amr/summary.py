"""Kleborate-style wide CSV summary + long-format per-feature CSV.

Wide CSV: one row per sample, columns grouped:
    sample | organism | ST | linezolid_call | lzd_23S_mutations | lzd_genes
          | AMR_<class1> | AMR_<class2> | ... | virulence | stress | metal

Long CSV: one row per detected feature (gene/mutation), regardless of class.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

from linezolid_amr.amrfinder import AmrHit
from linezolid_amr.rrna23s import PileupCall


# Element types we promote into their own top-level "buckets".
# Anything labeled as AMR with a Class collapses into AMR_<Class>.
ELEMENT_BUCKETS = ("AMR", "VIRULENCE", "STRESS")


def _fmt_amr_hit(h: AmrHit) -> str:
    """A compact one-liner for a gene/mutation: name(cov:ident)."""
    return f"{h.gene_symbol}({h.coverage_pct:.0f}/{h.identity_pct:.1f})"


def _fmt_lzd_position(p: PileupCall) -> str:
    """Compact LZD position formatter for the wide CSV: e.g. G2576T:0.66(173)."""
    alts = [a for a in p.alt_alleles if a["resistance"]]
    if not alts:
        return ""
    top = max(alts, key=lambda a: a["af"])
    return f"{p.ref_base}{p.ecoli_position}{top['base']}:{top['af']:.4f}({p.depth})"


def _bucket_class(amr: AmrHit) -> str:
    """Bucket key used as a wide-CSV column header."""
    et = (amr.element_type or "").upper()
    cls = (amr.class_ or "OTHER").upper().replace("/", "_").replace(" ", "_")
    if et == "AMR":
        return f"AMR_{cls}"
    if et == "VIRULENCE":
        return "VIRULENCE"
    if et == "STRESS":
        return f"STRESS_{cls}" if amr.class_ else "STRESS"
    return f"{et}_{cls}" if et else cls


def build_wide_row(
    sample: str,
    organism: str,
    st: str | None,
    mlst_scheme: str | None,
    amr_hits: list[AmrHit],
    pileup_calls: list[PileupCall],
    linezolid_call: bool,
) -> dict[str, str]:
    """Return a single dict representing the sample's wide-CSV row."""
    # 1) Linezolid columns — always first after identification
    lzd_23s_mutations = [_fmt_lzd_position(p) for p in pileup_calls if p.is_resistance]
    lzd_genes = sorted({h.gene_symbol for h in amr_hits if h.is_linezolid_relevant and h.gene_symbol})

    row: dict[str, str] = {
        "sample": sample,
        "organism": organism or "",
        "mlst_scheme": mlst_scheme or "",
        "ST": st or "",
        "linezolid_call": "POS" if linezolid_call else "neg",
        "lzd_23S_mutations": ";".join(lzd_23s_mutations),
        "lzd_genes": ";".join(lzd_genes),
    }

    # 2) Bucket all other AMR / virulence / stress hits by class
    buckets: dict[str, list[str]] = {}
    for h in amr_hits:
        key = _bucket_class(h)
        buckets.setdefault(key, []).append(_fmt_amr_hit(h))
    for key, vals in sorted(buckets.items()):
        row[key] = ";".join(sorted(set(vals)))
    return row


def build_long_rows(
    sample: str,
    organism: str,
    st: str | None,
    mlst_scheme: str | None,
    amr_hits: list[AmrHit],
    pileup_calls: list[PileupCall],
) -> list[dict[str, str]]:
    """Long-format rows: one per detected feature (gene/mutation)."""
    base = {"sample": sample, "organism": organism, "mlst_scheme": mlst_scheme or "", "ST": st or ""}
    out: list[dict[str, str]] = []

    # 23S resistance positions (one row each, even if multiple alt alleles)
    for p in pileup_calls:
        if not p.is_resistance:
            continue
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
                "coverage_pct": "",
                "identity_pct": "",
                "contig": p.ref_contig,
            })

    # AMRFinderPlus hits
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


# ---------------- writers ---------------- #

def write_wide_csv(rows: list[dict[str, str]], path: Path) -> None:
    """Write rows as a CSV with union-of-keys columns. LZD columns first."""
    if not rows:
        path.write_text("sample,organism,mlst_scheme,ST,linezolid_call,lzd_23S_mutations,lzd_genes\n")
        return
    # Pinned-first columns
    fixed_first = [
        "sample", "organism", "mlst_scheme", "ST",
        "linezolid_call", "lzd_23S_mutations", "lzd_genes",
    ]
    extras = sorted({k for r in rows for k in r.keys() if k not in fixed_first})
    fieldnames = fixed_first + extras
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
