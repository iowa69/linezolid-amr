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


def _coverage_summary(pileup_calls: list[PileupCall], parameters: dict) -> dict:
    """Depth across the screened positions, and the resulting limit of detection.

    A negative call is only as trustworthy as the depth behind it. At 30x, the
    3-read minimum means nothing below ~10% could have been called no matter
    how real it was; reporting that number stops a shallow negative from being
    read as a confident one.
    """
    if not pileup_calls:
        return {}
    depths = sorted(c.base_depth or c.depth for c in pileup_calls)
    n = len(depths)
    median = depths[n // 2] if n % 2 else (depths[n // 2 - 1] + depths[n // 2]) / 2
    min_depth = depths[0]

    min_af = float(parameters.get("min_af", 0.15))
    min_alt_reads = int(parameters.get("min_alt_reads", 3))
    depth_floor = int(parameters.get("min_depth", 20))

    # The lowest fraction that could have produced a call at the shallowest
    # position: whichever of the frequency and read-count gates binds harder.
    lod = max(min_af, (min_alt_reads / min_depth) if min_depth else 1.0)
    return {
        "min_depth_observed": min_depth,
        "median_depth": median,
        "max_depth_observed": depths[-1],
        "positions_below_min_depth": sum(1 for d in depths if d < depth_floor),
        "limit_of_detection_af": round(min(1.0, lod), 4),
        "adequate_depth": min_depth >= depth_floor,
    }


def _mechanisms(lzd_amr_hits: list[AmrHit], positions: list[dict]) -> list[str]:
    """Name the resistance mechanisms behind a positive call.

    Acquired genes and target-site mutations have different clinical
    implications — one is transferable, the other is not — so a bare
    POSITIVE/negative is not enough to act on.
    """
    out: list[str] = []
    genes = sorted({h.gene_symbol for h in lzd_amr_hits if h.gene_symbol})
    if genes:
        out.append("acquired_gene:" + ",".join(genes))
    if positions:
        muts = ",".join(p["mutation"] for p in positions if p["mutation"])
        fixed = [p for p in positions if p["zygosity"] == "fixed"]
        hetero = [p for p in positions if p["zygosity"] == "heteroresistant"]
        if fixed:
            out.append("23S_target_mutation_fixed:" + muts)
        if hetero:
            out.append("23S_target_mutation_heteroresistant:" + muts)
    return out


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


# An allele at or above this fraction is present in essentially every rrn
# operon, so the isolate is uniformly resistant rather than heteroresistant.
FIXED_AF_THRESHOLD = 0.90


def classify_zygosity(af: float) -> str:
    """Describe what an allele fraction implies about the operon population."""
    if af >= FIXED_AF_THRESHOLD:
        return "fixed"
    return "heteroresistant"


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
    # Positions where a resistance allele was seen but did not survive filtering.
    # Reporting these is what makes a negative call interpretable. Report them
    # even where a second allele at the same position passed — "G2576T called,
    # G2576A rejected for strand bias" is more informative than hiding half of
    # what was observed.
    rejected = [c for c in pileup_calls if c.filtered_resistance]

    def _passing_alts(c: PileupCall) -> list[dict]:
        return [a for a in c.alt_alleles if a["resistance"] and a.get("passes_filters", True)]

    positions = []
    for c in lzd_23s_hits:
        alts = _passing_alts(c)
        max_af = max((a["af"] for a in alts), default=0.0)
        best = max(alts, key=lambda a: a["af"], default=None)
        positions.append({
            "ecoli_position": c.ecoli_position,
            "species_position": c.species_position,
            "ref_base": c.ref_base,
            "resistance_alleles_present": [a["base"] for a in alts],
            "max_resistance_af": max_af,
            "af_ci_low": best["af_ci_low"] if best else 0.0,
            "af_ci_high": best["af_ci_high"] if best else 0.0,
            "zygosity": classify_zygosity(max_af),
            "est_mutated_operons": best["est_mutated_operons"] if best else 0,
            "rrn_copy_number": best["rrn_copy_number"] if best else 0,
            "strand_bias_p": best["strand_bias_p"] if best else 1.0,
            "depth": c.depth,
            "base_depth": c.base_depth,
            "mutation": f"{c.ref_base}{c.ecoli_position}{best['base']}" if best else "",
        })

    heteroresistant = [p for p in positions if p["zygosity"] == "heteroresistant"]
    coverage = _coverage_summary(pileup_calls, parameters)
    summary = {
        "linezolid_resistance_call": bool(lzd_amr_hits or lzd_23s_hits),
        "resistance_mechanism": _mechanisms(lzd_amr_hits, positions),
        "heteroresistance_detected": bool(heteroresistant),
        "n_total_amr_hits": len(amr_hits),
        "n_linezolid_amr_hits": len(lzd_amr_hits),
        "n_23s_positions_evaluated": len(pileup_calls),
        "n_23s_positions_with_resistance_allele": len(lzd_23s_hits),
        "n_23s_positions_rejected_by_filters": len(rejected),
        "coverage_23s": coverage,
        "linezolid_amr_genes_detected": sorted({h.gene_symbol for h in lzd_amr_hits if h.gene_symbol}),
        "linezolid_23s_positions": positions,
        "rejected_23s_positions": [
            {
                "ecoli_position": c.ecoli_position,
                "ref_base": c.ref_base,
                "depth": c.depth,
                "rejections": c.filtered_resistance,
            }
            for c in rejected
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
    call = "POSITIVE" if s["linezolid_resistance_call"] else "negative"
    lines = [
        "# linezolid-amr report",
        f"sample: {report['sample']}",
        f"organism: {report['organism']}",
        f"tool version: {report['version']}",
        f"generated: {report['generated_at']}",
        "",
        f"linezolid_resistance_call: {call}",
    ]
    if s.get("resistance_mechanism"):
        lines.append("mechanism: " + "; ".join(s["resistance_mechanism"]))
    if s.get("heteroresistance_detected"):
        lines.append(
            "HETERORESISTANCE: a resistance allele is present in only part of the "
            "rrn operon population (see allele fractions below)."
        )
    lines += [
        f"total AMR hits: {s['n_total_amr_hits']}",
        f"linezolid AMR genes: {', '.join(s['linezolid_amr_genes_detected']) or '(none)'}",
        f"23S positions evaluated: {s['n_23s_positions_evaluated']}",
        f"23S positions with resistance allele: {s['n_23s_positions_with_resistance_allele']}",
    ]
    cov = s.get("coverage_23s") or {}
    if cov:
        lines.append(
            f"23S depth: median {cov['median_depth']:.0f}x, "
            f"minimum {cov['min_depth_observed']}x "
            f"(limit of detection ~AF {cov['limit_of_detection_af']:.2f})"
        )
        if not cov["adequate_depth"]:
            lines.append(
                f"CAUTION: {cov['positions_below_min_depth']} position(s) fell below the "
                f"minimum depth. A negative result at those positions is not informative."
            )
    lines += [
        "",
        "## 23S resistance positions (E. coli numbering)",
    ]
    if not s["linezolid_23s_positions"]:
        lines.append("(none above thresholds)")
    else:
        lines.append(
            "mutation\tref\talleles\tAF\t95% CI\toperons\tzygosity\tdepth"
        )
        for p in s["linezolid_23s_positions"]:
            operons = (
                f"{p['est_mutated_operons']}/{p['rrn_copy_number']}"
                if p.get("rrn_copy_number") else "-"
            )
            lines.append(
                f"{p.get('mutation') or p['ecoli_position']}\t{p['ref_base']}\t"
                f"{','.join(p['resistance_alleles_present'])}\t"
                f"{p['max_resistance_af']:.4f}\t"
                f"{p.get('af_ci_low', 0):.3f}-{p.get('af_ci_high', 0):.3f}\t"
                f"{operons}\t{p.get('zygosity', '')}\t{p['depth']}"
            )

    rejected = s.get("rejected_23s_positions") or []
    if rejected:
        lines += [
            "",
            "## Resistance alleles observed but rejected by filters",
            "These did not meet the evidence bar for a positive call. Inspect the",
            "*.23S_lzd_evidence.tsv file before dismissing them.",
            "ecoli_pos\tref\tdepth\treason",
        ]
        for r in rejected:
            lines.append(
                f"{r['ecoli_position']}\t{r['ref_base']}\t{r['depth']}\t"
                f"{','.join(r['rejections'])}"
            )
    path.write_text("\n".join(lines) + "\n")
