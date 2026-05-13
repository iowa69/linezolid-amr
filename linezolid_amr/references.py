"""Reference selection and loci metadata for 23S rRNA linezolid resistance analysis."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Iterable

LOCI_RESOURCE = "loci.json"
PACKAGE_REFDIR = "linezolid_amr.data.references"


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


def cache_dir() -> Path:
    """Resolve the local reference cache directory.

    Resolution order:
      1) $LINEZOLID_AMR_REFDIR if set
      2) $XDG_DATA_HOME/linezolid-amr/references
      3) ~/.local/share/linezolid-amr/references
    """
    env = os.environ.get("LINEZOLID_AMR_REFDIR")
    if env:
        return Path(env).expanduser().resolve()
    xdg = os.environ.get("XDG_DATA_HOME")
    base = Path(xdg).expanduser() if xdg else Path.home() / ".local" / "share"
    return (base / "linezolid-amr" / "references").resolve()


def organism_fasta_path(organism: str) -> Path:
    return cache_dir() / f"{organism}_23S.fasta"


def organism_bed_path(organism: str) -> Path:
    """BED of canonical LZD positions in species-specific coordinates (computed by fetch-references)."""
    return cache_dir() / f"{organism}_23S_lzd_positions.bed"


def organism_position_map_path(organism: str) -> Path:
    """TSV mapping E. coli positions to species-specific positions in the bundled reference."""
    return cache_dir() / f"{organism}_23S_position_map.tsv"


def ecoli_fasta_path() -> Path:
    return cache_dir() / "ecoli_K12_23S_rrlB.fasta"


def ensure_references_available(organism: str) -> tuple[Path, Path]:
    """Return (fasta, bed). Raises if references not yet fetched."""
    fasta = organism_fasta_path(organism)
    bed = organism_bed_path(organism)
    if not fasta.exists() or not bed.exists():
        raise FileNotFoundError(
            f"References for '{organism}' not found in cache ({cache_dir()}). "
            f"Run 'linezolid-amr fetch-references' first."
        )
    return fasta, bed
