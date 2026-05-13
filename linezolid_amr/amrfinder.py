"""AMRFinderPlus runner and result parser."""

from __future__ import annotations

import csv
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

# Substrings (lowercase) that flag a result as linezolid-relevant.
# Includes oxazolidinone resistance genes and ribosomal-protein mutation targets.
LINEZOLID_GENE_KEYWORDS = (
    "linezolid",
    "oxazolidinone",
    "cfr",
    "optra",
    "poxta",
)
LINEZOLID_PROTEIN_KEYWORDS = (
    "rplc",  # L3
    "rpld",  # L4
    "rplv",  # L22
    "23s",
)


@dataclass
class AmrHit:
    contig: str
    start: int
    end: int
    strand: str
    gene_symbol: str
    sequence_name: str
    scope: str
    element_type: str
    element_subtype: str
    class_: str
    subclass: str
    method: str
    target_length: int
    reference_length: int
    coverage_pct: float
    identity_pct: float
    accession: str
    raw: dict = field(default_factory=dict)

    @property
    def is_linezolid_relevant(self) -> bool:
        text = " ".join(
            [self.gene_symbol, self.sequence_name, self.class_, self.subclass]
        ).lower()
        if any(k in text for k in LINEZOLID_GENE_KEYWORDS):
            return True
        if any(k in self.gene_symbol.lower() for k in LINEZOLID_PROTEIN_KEYWORDS):
            return True
        return False


def amrfinder_available() -> bool:
    return shutil.which("amrfinder") is not None


def run_amrfinder(
    assembly: Path,
    organism: str,
    outdir: Path,
    threads: int = 4,
    extra_args: Iterable[str] | None = None,
) -> Path:
    """Run AMRFinderPlus on an assembly. Returns the path to the TSV output."""
    if not amrfinder_available():
        raise RuntimeError(
            "amrfinder not found on PATH. Install via 'conda install -c bioconda ncbi-amrfinderplus' "
            "or 'mamba install -c bioconda ncbi-amrfinderplus'."
        )
    outdir.mkdir(parents=True, exist_ok=True)
    tsv = outdir / "amrfinder.tsv"
    log = outdir / "amrfinder.log"

    cmd = [
        "amrfinder",
        "--nucleotide",
        str(assembly),
        "--organism",
        organism,
        "--output",
        str(tsv),
        "--threads",
        str(threads),
        "--name",
        assembly.stem,
    ]
    if extra_args:
        cmd.extend(extra_args)

    with log.open("w") as logfh:
        logfh.write("$ " + " ".join(cmd) + "\n")
        logfh.flush()
        proc = subprocess.run(cmd, stdout=logfh, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(
            f"amrfinder failed (exit {proc.returncode}). See {log}"
        )
    return tsv


# AMRFinderPlus TSV header is stable across recent versions; we tolerate column
# additions by reading via DictReader.
_FIELD_ALIASES = {
    "Contig id": ("Contig id", "Contig"),
    "Start": ("Start",),
    "Stop": ("Stop",),
    "Strand": ("Strand",),
    "Gene symbol": ("Gene symbol", "Element symbol"),
    "Sequence name": ("Sequence name", "Element name"),
    "Scope": ("Scope",),
    "Element type": ("Element type",),
    "Element subtype": ("Element subtype",),
    "Class": ("Class",),
    "Subclass": ("Subclass",),
    "Method": ("Method",),
    "Target length": ("Target length",),
    "Reference sequence length": ("Reference sequence length", "Reference length"),
    "% Coverage of reference sequence": (
        "% Coverage of reference sequence",
        "% Coverage of reference",
    ),
    "% Identity to reference sequence": (
        "% Identity to reference sequence",
        "% Identity to reference",
    ),
    "Accession of closest sequence": (
        "Accession of closest sequence",
        "Closest reference accession",
    ),
}


def _pick(row: dict, key: str, default: str = "") -> str:
    for alias in _FIELD_ALIASES.get(key, (key,)):
        if alias in row and row[alias] != "":
            return row[alias]
    return default


def parse_amrfinder_tsv(tsv_path: Path) -> list[AmrHit]:
    hits: list[AmrHit] = []
    if not tsv_path.exists():
        return hits
    with tsv_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                start = int(_pick(row, "Start") or 0)
                end = int(_pick(row, "Stop") or 0)
                target_len = int(_pick(row, "Target length") or 0)
                ref_len = int(_pick(row, "Reference sequence length") or 0)
                cov = float(_pick(row, "% Coverage of reference sequence") or 0)
                ident = float(_pick(row, "% Identity to reference sequence") or 0)
            except ValueError:
                start = end = target_len = ref_len = 0
                cov = ident = 0.0
            hits.append(
                AmrHit(
                    contig=_pick(row, "Contig id"),
                    start=start,
                    end=end,
                    strand=_pick(row, "Strand"),
                    gene_symbol=_pick(row, "Gene symbol"),
                    sequence_name=_pick(row, "Sequence name"),
                    scope=_pick(row, "Scope"),
                    element_type=_pick(row, "Element type"),
                    element_subtype=_pick(row, "Element subtype"),
                    class_=_pick(row, "Class"),
                    subclass=_pick(row, "Subclass"),
                    method=_pick(row, "Method"),
                    target_length=target_len,
                    reference_length=ref_len,
                    coverage_pct=cov,
                    identity_pct=ident,
                    accession=_pick(row, "Accession of closest sequence"),
                    raw=dict(row),
                )
            )
    return hits


def linezolid_relevant_hits(hits: Iterable[AmrHit]) -> list[AmrHit]:
    return [h for h in hits if h.is_linezolid_relevant]
