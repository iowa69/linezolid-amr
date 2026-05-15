"""Compare our in-house MLST against Seemann's mlst on a folder of assemblies.

Usage:
    python scripts/validate_mlst_vs_seemann.py \
        --folder /home/iowa/Desktop/MORI/test \
        --seemann-mlst /home/iowa/miniconda3/envs/mlst/bin/mlst \
        --out validation.tsv

Reports allele-by-allele agreement. Exit code = number of mismatches.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from linezolid_amr import internal_mlst as im


def run_seemann(binary: Path, fasta: Path, env_name: str | None = None) -> tuple[str, str, dict[str, str]] | None:
    """Run Seemann's mlst. Returns (scheme, ST, alleles dict) or None on no match.

    If ``env_name`` is set, runs via ``conda run -n <env_name>`` so the
    correct PERL5LIB is loaded (mlst needs its own env's Perl modules).
    """
    if env_name:
        cmd = ["conda", "run", "-n", env_name, "mlst", "--quiet", "--nopath", str(fasta)]
    else:
        cmd = [str(binary), "--quiet", "--nopath", str(fasta)]
    out = subprocess.run(cmd, capture_output=True, text=True, check=False)
    line = out.stdout.strip()
    if not line:
        return None
    parts = line.split("\t")
    if len(parts) < 4:
        return None
    _file, scheme, st, *allele_fields = parts
    alleles = {}
    for af in allele_fields:
        if "(" in af and af.endswith(")"):
            g, v = af.split("(", 1)
            alleles[g] = v.rstrip(")")
    return scheme, st, alleles


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folder", required=True, type=Path)
    ap.add_argument("--seemann-mlst", required=True, type=Path,
                    help="Path to Seemann mlst binary (used if --seemann-env not given).")
    ap.add_argument("--seemann-env", default=None,
                    help="Conda env name where Seemann mlst lives (preferred — avoids Perl @INC mismatch).")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--threads", default=8, type=int)
    args = ap.parse_args()

    if not args.seemann_env and not args.seemann_mlst.exists():
        sys.exit(f"seemann mlst not found: {args.seemann_mlst}")

    rows = []
    n_match_st = 0
    n_match_alleles = 0
    n_total = 0
    mismatches = 0

    for fa in sorted(args.folder.glob("*.fasta")):
        n_total += 1
        sample = fa.stem
        # Seemann
        seem = run_seemann(args.seemann_mlst, fa, env_name=args.seemann_env)
        if seem is None:
            seem_scheme, seem_st, seem_alleles = "-", "-", {}
        else:
            seem_scheme, seem_st, seem_alleles = seem
        # Ours — infer organism first
        try:
            organism, ours = im.infer_organism(fa, threads=args.threads)
            ours_scheme = ours.scheme
            ours_st = ours.st
            ours_alleles = ours.alleles
        except Exception as e:  # noqa: BLE001
            ours_scheme, ours_st, ours_alleles = "-", "-", {}

        st_match = (seem_st == ours_st)
        alleles_match = (seem_alleles == ours_alleles)
        if st_match:
            n_match_st += 1
        if alleles_match:
            n_match_alleles += 1
        if not (st_match and alleles_match):
            mismatches += 1

        rows.append({
            "sample": sample,
            "seem_scheme": seem_scheme, "seem_ST": seem_st,
            "ours_scheme": ours_scheme, "ours_ST": ours_st,
            "st_match": "YES" if st_match else "NO",
            "alleles_match": "YES" if alleles_match else "NO",
            "seem_alleles": "|".join(f"{k}({v})" for k, v in seem_alleles.items()),
            "ours_alleles": "|".join(f"{k}({v})" for k, v in ours_alleles.items()),
        })
        print(f"  {sample}: seem={seem_scheme}/{seem_st}  ours={ours_scheme}/{ours_st}  "
              f"ST={'OK' if st_match else 'DIFF'}  alleles={'OK' if alleles_match else 'DIFF'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        fields = list(rows[0].keys()) if rows else []
        fh.write("\t".join(fields) + "\n")
        for r in rows:
            fh.write("\t".join(r[k] for k in fields) + "\n")

    print()
    print(f"Total samples: {n_total}")
    print(f"ST exact match:      {n_match_st} / {n_total} ({n_match_st/n_total*100:.1f}%)")
    print(f"Alleles exact match: {n_match_alleles} / {n_total} ({n_match_alleles/n_total*100:.1f}%)")
    print(f"Mismatches:          {mismatches}")
    print(f"Output: {args.out}")
    return mismatches


if __name__ == "__main__":
    sys.exit(main())
