"""Reference selection and loci metadata for 23S rRNA linezolid resistance analysis.

References ship inside the package under data/references/bundled/ so the tool
works offline out of the box. Users can override individual files by setting
LINEZOLID_AMR_REFDIR to a directory containing same-named files; that
directory wins over bundled data when files exist there. The legacy
fetch-references command still works and writes into the override directory.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

LOCI_RESOURCE = "loci.json"
PACKAGE_REFDIR = "linezolid_amr.data.references"
PACKAGE_BUNDLED_REFDIR = "linezolid_amr.data.references.bundled"


@dataclass(frozen=True)
class LinezolidPosition:
    ecoli_position: int
    ref_base: str
    resistance_bases: tuple[str, ...]
    drug: str
    note: str | None = None


@dataclass(frozen=True)
class OrganismRef:
    organism: str
    amrfinder_organism: str
    description: str
    genome_accession: str
    start: int
    end: int
    strand: str
    expected_length_bp: int


@dataclass(frozen=True)
class EcoliRef:
    label: str
    description: str
    genome_accession: str
    start: int
    end: int
    strand: str
    expected_length_bp: int


def load_loci() -> dict:
    """Load the bundled loci.json metadata."""
    with resources.files(PACKAGE_REFDIR).joinpath(LOCI_RESOURCE).open("r") as fh:
        return json.load(fh)


def list_organisms() -> list[str]:
    return sorted(load_loci()["organisms"].keys())


def get_linezolid_positions() -> list[LinezolidPosition]:
    data = load_loci()
    out = []
    for item in data["linezolid_positions_ecoli_23s"]:
        out.append(
            LinezolidPosition(
                ecoli_position=item["ecoli_position"],
                ref_base=item["ref_base"],
                resistance_bases=tuple(item["resistance_bases"]),
                drug=item["drug"],
                note=item.get("note"),
            )
        )
    return out


def get_ecoli_reference() -> EcoliRef:
    e = load_loci()["ecoli_reference"]
    return EcoliRef(
        label=e["label"],
        description=e["description"],
        genome_accession=e["genome_accession"],
        start=e["start"],
        end=e["end"],
        strand=e["strand"],
        expected_length_bp=e["expected_length_bp"],
    )


def get_organism(organism: str) -> OrganismRef:
    data = load_loci()
    if organism not in data["organisms"]:
        raise KeyError(
            f"Organism '{organism}' is not supported for 23S analysis. "
            f"Supported: {', '.join(list_organisms())}"
        )
    o = data["organisms"][organism]
    return OrganismRef(
        organism=organism,
        amrfinder_organism=o["amrfinder_organism"],
        description=o["description"],
        genome_accession=o["genome_accession"],
        start=o["start"],
        end=o["end"],
        strand=o["strand"],
        expected_length_bp=o["expected_length_bp"],
    )


def override_dir() -> Path | None:
    """Optional override directory for users who want to swap in custom references.

    Resolution order:
      1) $LINEZOLID_AMR_REFDIR if set
      2) $XDG_DATA_HOME/linezolid-amr/references if it exists
      3) ~/.local/share/linezolid-amr/references if it exists

    Returns None if no override directory is set or present.
    """
    env = os.environ.get("LINEZOLID_AMR_REFDIR")
    if env:
        return Path(env).expanduser().resolve()
    xdg = os.environ.get("XDG_DATA_HOME")
    candidate = (Path(xdg).expanduser() if xdg else Path.home() / ".local" / "share")
    candidate = (candidate / "linezolid-amr" / "references").resolve()
    if candidate.exists():
        return candidate
    return None


# Back-compat alias (older code/tests used cache_dir()). Always returns a path —
# either the override location or the default cache path (whether it exists or not).
def cache_dir() -> Path:
    od = override_dir()
    if od is not None:
        return od
    xdg = os.environ.get("XDG_DATA_HOME")
    base = Path(xdg).expanduser() if xdg else Path.home() / ".local" / "share"
    return (base / "linezolid-amr" / "references").resolve()


def _bundled_path(name: str) -> Path:
    """Return a filesystem Path to a bundled reference file inside the package."""
    return Path(str(resources.files(PACKAGE_BUNDLED_REFDIR).joinpath(name)))


def _resolve(name: str) -> Path:
    """Override directory wins if it has the file; otherwise package data."""
    od = override_dir()
    if od is not None:
        candidate = od / name
        if candidate.exists():
            return candidate
    return _bundled_path(name)


def organism_fasta_path(organism: str) -> Path:
    return _resolve(f"{organism}_23S.fasta")


def organism_bed_path(organism: str) -> Path:
    return _resolve(f"{organism}_23S_lzd_positions.bed")


def organism_position_map_path(organism: str) -> Path:
    return _resolve(f"{organism}_23S_position_map.tsv")


def ecoli_fasta_path() -> Path:
    return _resolve("ecoli_K12_23S_rrlB.fasta")


def ensure_references_available(organism: str) -> tuple[Path, Path]:
    """Return (fasta, bed). Bundled references are always present."""
    fasta = organism_fasta_path(organism)
    bed = organism_bed_path(organism)
    if not fasta.exists() or not bed.exists():
        raise FileNotFoundError(
            f"References for '{organism}' not found. Bundled set may be missing — "
            f"reinstall linezolid-amr or run 'linezolid-amr fetch-references' "
            f"to populate an override directory."
        )
    return fasta, bed
