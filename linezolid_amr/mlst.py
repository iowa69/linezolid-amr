"""MLST (Seemann's `mlst`) runner and organism inference.

`mlst` auto-detects the scheme. We map its scheme name to the AMRFinderPlus
`--organism` string used by AMRFinderPlus and by our 23S analysis.
"""

from __future__ import annotations

import csv
import io
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

# Map mlst scheme -> AMRFinderPlus organism (only the ones we support for 23S).
# Note `mlst` ships many schemes; for genera with several species (e.g. neisseria)
# you might need additional logic. We support the four 23S-LZD species.
MLST_SCHEME_TO_ORGANISM = {
    "saureus":       "Staphylococcus_aureus",
    "efaecalis":     "Enterococcus_faecalis",
    "efaecium":      "Enterococcus_faecium",
    "spneumoniae":   "Streptococcus_pneumoniae",
}

# Helpful mapping for the error message — which schemes mlst would *call* but we
# don't have a 23S reference for. Lets us tell users "you have an X — out of scope".
KNOWN_UNSUPPORTED_SCHEMES = {
    "ecoli":         "Escherichia coli",
    "klebsiella":    "Klebsiella spp.",
    "kpneumoniae":   "Klebsiella pneumoniae",
    "pseudomonas":   "Pseudomonas spp.",
    "paeruginosa":   "Pseudomonas aeruginosa",
    "abaumannii":    "Acinetobacter baumannii",
    "campylobacter": "Campylobacter spp.",
    "spyogenes":     "Streptococcus pyogenes",
    "sagalactiae":   "Streptococcus agalactiae",
}


@dataclass
class MlstResult:
    file: str
    scheme: str
    st: str
    alleles: dict[str, str]
    raw: str

    @property
    def organism(self) -> str | None:
        return MLST_SCHEME_TO_ORGANISM.get(self.scheme)


class MlstUnsupportedOrganism(ValueError):
    """MLST detected an organism we don't have an LZD 23S reference for."""


class MlstNoMatch(ValueError):
    """MLST could not determine any scheme."""


def mlst_available() -> bool:
    return shutil.which("mlst") is not None


def run_mlst(assembly: Path, threads: int = 4) -> MlstResult:
    """Run `mlst` on the assembly and parse its TSV output."""
    if not mlst_available():
        raise RuntimeError(
            "mlst not found on PATH. Install via `conda install -c bioconda mlst`."
        )
    cmd = [
        "mlst",
        "--quiet",
        "--threads", str(threads),
        "--nopath",
        str(assembly),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"mlst failed (exit {proc.returncode}): {proc.stderr.strip()}")
    line = proc.stdout.strip()
    if not line:
        raise MlstNoMatch("mlst returned no output")
    fields = line.split("\t")
    # mlst output: <file> <scheme> <ST> <allele1> <allele2> ...
    file_, scheme, st, *allele_fields = fields
    alleles: dict[str, str] = {}
    for af in allele_fields:
        if "(" in af and af.endswith(")"):
            gene, val = af.split("(", 1)
            alleles[gene] = val.rstrip(")")
        else:
            alleles[af] = ""
    return MlstResult(file=file_, scheme=scheme, st=st, alleles=alleles, raw=line)


def infer_organism(assembly: Path, threads: int = 4) -> tuple[str, MlstResult]:
    """Run MLST and return (AMRFinderPlus organism, MlstResult).

    Raises MlstNoMatch if mlst can't pick a scheme, or MlstUnsupportedOrganism
    if the scheme is not in our supported LZD 23S species list.
    """
    res = run_mlst(assembly, threads=threads)
    if res.scheme in ("-", "", "?"):
        raise MlstNoMatch(
            f"mlst could not determine a scheme for {assembly.name}. "
            f"Provide --organism manually."
        )
    organism = MLST_SCHEME_TO_ORGANISM.get(res.scheme)
    if not organism:
        friendly = KNOWN_UNSUPPORTED_SCHEMES.get(res.scheme, res.scheme)
        raise MlstUnsupportedOrganism(
            f"mlst detected scheme '{res.scheme}' ({friendly}) which is outside the "
            f"linezolid 23S catalogue. Supported organisms: "
            f"{', '.join(sorted(MLST_SCHEME_TO_ORGANISM.values()))}. "
            f"Override with --organism if you want to run AMRFinderPlus + 23S analysis anyway."
        )
    return organism, res
