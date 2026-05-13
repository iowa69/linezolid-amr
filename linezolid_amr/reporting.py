"""Combine AMRFinderPlus hits and 23S pileup calls into a unified report."""

from __future__ import annotations

import json
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from linezolid_amr import __version__
from linezolid_amr.amrfinder import AmrHit
from linezolid_amr.rrna23s import PileupCall


def _amr_to_dict(h: AmrHit) -> dict:
    return {
        "contig": h.contig,
        "start": h.start,
        "end": h.end,
        "strand": h.strand,
        "gene_symbol": h.gene_symbol,
        "sequence_name": h.sequence_name,
        "scope": h.scope,
        "element_type": h.element_type,
        "element_subtype": h.element_subtype,
        "class": h.class_,
        "subclass": h.subclass,
        "method": h.method,
        "coverage_pct": h.coverage_pct,
        "identity_pct": h.identity_pct,
        "accession": h.accession,
        "linezolid_relevant": h.is_linezolid_relevant,
    }


def build_report(
    sample: str,
    organism: str,
    amr_hits: list[AmrHit] | None,
    pileup_calls: list[PileupCall] | None,
    vcf_path: Path | None,
    parameters: dict,
) -> dict:
    amr_hits = amr_hits or []
    pileup_calls = pileup_calls or []

    lzd_amr_hits = [h for h in amr_hits if h.is_linezolid_relevant]
    lzd_23s_hits = [c for c in pileup_calls if c.is_resistance]

    summary = {
        "linezolid_resistance_call": bool(lzd_amr_hits or lzd_23s_hits),
        "n_total_amr_hits": len(amr_hits),
        "n_linezolid_amr_hits": len(lzd_amr_hits),
        "n_23s_positions_evaluated": len(pileup_calls),
        "n_23s_positions_with_resistance_allele": len(lzd_23s_hits),
        "linezolid_amr_genes_detected": sorted({h.gene_symbol for h in lzd_amr_hits if h.gene_symbol}),
        "linezolid_23s_positions": [
            {
                "ecoli_position": c.ecoli_position,
                "species_position": c.species_position,
                "ref_base": c.ref_base,
                "resistance_alleles_present": [a["base"] for a in c.alt_alleles if a["resistance"]],
                "max_resistance_af": max(
                    (a["af"] for a in c.alt_alleles if a["resistance"]), default=0.0
                ),
                "depth": c.depth,
            }
            for c in lzd_23s_hits
        ],
    }

    return {
        "tool": "linezolid-amr",
        "version": __version__,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "sample": sample,
        "organism": organism,
        "parameters": parameters,
        "summary": summary,
        "amr_hits": [_amr_to_dict(h) for h in amr_hits],
        "rrna23s_pileup": [
            {
                **asdict(c),
                "counts": dict(c.counts),
            }
            for c in pileup_calls
        ],
        "rrna23s_vcf": str(vcf_path) if vcf_path else None,
    }


def write_json(report: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        json.dump(report, fh, indent=2, default=str)


def write_text_summary(report: dict, path: Path) -> None:
    s = report["summary"]
    lines = [
        f"# linezolid-amr report",
        f"sample: {report['sample']}",
        f"organism: {report['organism']}",
        f"generated: {report['generated_at']}",
        "",
        f"linezolid_resistance_call: {s['linezolid_resistance_call']}",
        f"total AMR hits: {s['n_total_amr_hits']}",
        f"linezolid AMR genes: {', '.join(s['linezolid_amr_genes_detected']) or '(none)'}",
        f"23S positions evaluated: {s['n_23s_positions_evaluated']}",
        f"23S positions with resistance allele: {s['n_23s_positions_with_resistance_allele']}",
        "",
        "## 23S resistance positions (E. coli numbering)",
    ]
    if not s["linezolid_23s_positions"]:
        lines.append("(none above thresholds)")
    else:
        lines.append("ecoli_pos\tref\talleles\tmax_af\tdepth")
        for p in s["linezolid_23s_positions"]:
            lines.append(
                f"{p['ecoli_position']}\t{p['ref_base']}\t"
                f"{','.join(p['resistance_alleles_present'])}\t"
                f"{p['max_resistance_af']:.4f}\t{p['depth']}"
            )
    path.write_text("\n".join(lines) + "\n")
