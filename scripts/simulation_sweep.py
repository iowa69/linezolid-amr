#!/usr/bin/env python3
"""Measure detection performance across allele fractions, organisms and depths.

Simulates 23S reads carrying G2576T at a known fraction, runs the real
mapping + pileup path, and reports how often the tool calls resistance. This
quantifies the sensitivity/specificity trade-off the default thresholds imply,
and is the harness to re-run when a threshold changes.

Usage:
    python scripts/simulation_sweep.py [--replicates 5] [--outdir DIR]

Requires minimap2 and samtools on PATH.
"""

from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from linezolid_amr import references as ref_mod  # noqa: E402
from linezolid_amr import rrna23s as rrna  # noqa: E402
from tests.simulate import read_fasta_single, simulate_reads  # noqa: E402

# Species coordinate of E. coli G2576, from the bundled BED files.
G2576_SPECIES_POS = {
    "Staphylococcus_aureus": 2603,
    "Enterococcus_faecalis": 2589,
    "Enterococcus_faecium": 2590,
    "Streptococcus_pneumoniae": 2578,
}

AF_GRID = [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.33, 0.50, 0.75, 1.0]
DEPTH_GRID = [30, 100, 300]


def run_one(organism: str, seq: str, fasta: Path, af: float, depth: int,
            seed: int, work: Path) -> tuple[bool, float, int]:
    r1 = work / "r1.fastq.gz"
    r2 = work / "r2.fastq.gz"
    pos = G2576_SPECIES_POS[organism]
    simulate_reads(
        seq, r1, r2,
        variants={pos: ("T", af)} if af > 0 else {},
        depth=depth, seed=seed,
    )
    bam = rrna.map_reads(fasta, r1, r2, work, threads=4, sample="s")
    calls = rrna.pileup_at_positions(bam, fasta, organism)
    call = next(c for c in calls if c.ecoli_position == 2576)
    observed = next((a["af"] for a in call.alt_alleles if a["base"] == "T"), 0.0)
    # Any positive anywhere counts as a sample-level call.
    any_pos = any(c.is_resistance for c in calls)
    return any_pos, observed, call.depth


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--organisms", nargs="*", default=list(G2576_SPECIES_POS))
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()

    for tool in ("minimap2", "samtools"):
        if shutil.which(tool) is None:
            print(f"error: {tool} not on PATH", file=sys.stderr)
            return 2

    rows = []
    tmp_root = args.outdir or Path(tempfile.mkdtemp(prefix="lzd-sweep-"))
    tmp_root.mkdir(parents=True, exist_ok=True)

    for organism in args.organisms:
        fasta_src = ref_mod.organism_fasta_path(organism)
        _name, seq = read_fasta_single(fasta_src)
        staged = ref_mod.stage_reference(organism, tmp_root / organism)
        copies = ref_mod.rrn_copy_number(organism)
        for depth in DEPTH_GRID:
            for af in AF_GRID:
                pos_calls = 0
                afs = []
                for rep in range(args.replicates):
                    work = tmp_root / organism / f"d{depth}_af{af}_r{rep}"
                    work.mkdir(parents=True, exist_ok=True)
                    called, observed, obs_depth = run_one(
                        organism, seq, staged, af, depth,
                        seed=hash((organism, depth, af, rep)) % 10_000, work=work,
                    )
                    pos_calls += called
                    afs.append(observed)
                rows.append({
                    "organism": organism, "rrn": copies, "depth": depth,
                    "true_af": af, "calls": pos_calls, "n": args.replicates,
                    "mean_obs_af": sum(afs) / len(afs),
                })
                print(
                    f"{organism:26s} rrn={copies} depth={depth:4d} "
                    f"true_af={af:.2f} called={pos_calls}/{args.replicates} "
                    f"mean_obs_af={sum(afs)/len(afs):.4f}"
                )

    print()
    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)
    neg = [r for r in rows if r["true_af"] == 0.0]
    fp = sum(r["calls"] for r in neg)
    print(f"Specificity: {sum(r['n'] for r in neg) - fp}/{sum(r['n'] for r in neg)} "
          f"true negatives correctly called ({fp} false positives)")

    for label, lo, hi in [
        ("sub-threshold (AF < 0.15)", 0.001, 0.1499),
        ("at/above threshold (AF >= 0.15)", 0.15, 1.01),
    ]:
        subset = [r for r in rows if lo <= r["true_af"] <= hi]
        if not subset:
            continue
        called = sum(r["calls"] for r in subset)
        total = sum(r["n"] for r in subset)
        print(f"{label}: {called}/{total} called ({100*called/total:.1f}%)")

    print()
    print("Detection by true allele fraction (all organisms, all depths):")
    for af in AF_GRID:
        subset = [r for r in rows if r["true_af"] == af]
        called = sum(r["calls"] for r in subset)
        total = sum(r["n"] for r in subset)
        bar = "#" * int(30 * called / total) if total else ""
        print(f"  AF {af:4.2f}: {called:3d}/{total:3d}  {bar}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
