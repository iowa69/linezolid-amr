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
import shutil
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

LOCI_RESOURCE = "loci.json"
PACKAGE_REFDIR = "linezolid_amr.data.references"
PACKAGE_BUNDLED_REFDIR = "linezolid_amr.data.references.bundled"


# Evidence tiers, strongest first. A position drives a positive resistance call
# only if its tier is at or above the run's threshold; everything is reported
# either way. This keeps lab-only substitutions visible without letting them
# manufacture clinical positives.
TIER_ESTABLISHED = "established"
TIER_ASSOCIATED = "associated"
TIER_EXPERIMENTAL = "experimental"
EVIDENCE_TIERS = (TIER_ESTABLISHED, TIER_ASSOCIATED, TIER_EXPERIMENTAL)
DEFAULT_EVIDENCE_TIER = TIER_ASSOCIATED


def tier_rank(tier: str) -> int:
    """Lower rank means stronger evidence."""
    try:
        return EVIDENCE_TIERS.index(tier)
    except ValueError:
        return len(EVIDENCE_TIERS)


def tier_is_callable(tier: str, threshold: str = DEFAULT_EVIDENCE_TIER) -> bool:
    """True when *tier* is strong enough to drive a positive call at *threshold*."""
    return tier_rank(tier) <= tier_rank(threshold)


@dataclass(frozen=True)
class LinezolidPosition:
    ecoli_position: int
    ref_base: str
    resistance_bases: tuple[str, ...]
    drug: str
    note: str | None = None
    # Literature designation, which may spell the wild-type letter using a
    # different organism's base (e.g. "C2534U" where E. coli carries A).
    published_as: str = ""
    evidence_tier: str = TIER_ESTABLISHED
    organisms: tuple[str, ...] = ()
    pmid: tuple[str, ...] = ()
    numbering_caveat: str | None = None


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
    # Number of rrn operons in a typical genome of this species. Sets the grid
    # of allele fractions a heteroresistant isolate can produce (k/n operons).
    rrn_copy_number: int = 0


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


def get_linezolid_positions(
    min_tier: str | None = None,
) -> list[LinezolidPosition]:
    """Return the curated linezolid positions, strongest evidence first.

    ``min_tier`` filters to positions at or above that evidence tier; the
    default returns every curated position.
    """
    data = load_loci()
    out = []
    for item in data["linezolid_positions_ecoli_23s"]:
        tier = item.get("evidence_tier", TIER_ESTABLISHED)
        if min_tier is not None and not tier_is_callable(tier, min_tier):
            continue
        out.append(
            LinezolidPosition(
                # ecoli_ref_base is authoritative; ref_base is the pre-v3 name.
                ecoli_position=item["ecoli_position"],
                ref_base=item.get("ecoli_ref_base", item["ref_base"]),
                resistance_bases=tuple(item["resistance_bases"]),
                drug=item["drug"],
                note=item.get("note"),
                published_as=item.get("published_as", ""),
                evidence_tier=tier,
                organisms=tuple(item.get("organisms", ())),
                pmid=tuple(item.get("pmid", ())),
                numbering_caveat=item.get("numbering_caveat"),
            )
        )
    out.sort(key=lambda p: (tier_rank(p.evidence_tier), p.ecoli_position))
    return out


def position_tiers() -> dict[int, str]:
    """Map each curated E. coli position to its evidence tier."""
    return {p.ecoli_position: p.evidence_tier for p in get_linezolid_positions()}


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
        rrn_copy_number=int(o.get("rrn_copy_number", 0)),
    )


def rrn_copy_number(organism: str) -> int:
    """rrn operon count for *organism*, or 0 when unknown."""
    try:
        return get_organism(organism).rrn_copy_number
    except KeyError:
        return 0


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


# Bundled-only accessors. The normal lookups above let a user override
# directory win, which is right at runtime but wrong for tooling that must
# read or write the files actually shipped in the package — a developer with
# a stale override would otherwise regenerate their own copy and leave the
# release data untouched.

def bundled_organism_fasta_path(organism: str) -> Path:
    return _bundled_path(f"{organism}_23S.fasta")


def bundled_organism_bed_path(organism: str) -> Path:
    return _bundled_path(f"{organism}_23S_lzd_positions.bed")


def bundled_organism_position_map_path(organism: str) -> Path:
    return _bundled_path(f"{organism}_23S_position_map.tsv")


def bundled_ecoli_fasta_path() -> Path:
    return _bundled_path("ecoli_K12_23S_rrlB.fasta")


def ecoli_fasta_path() -> Path:
    return _resolve("ecoli_K12_23S_rrlB.fasta")


def stage_reference(organism: str, workdir: Path) -> Path:
    """Copy the organism 23S FASTA into *workdir* and return the copy.

    Indexing a FASTA writes a sibling ``.fai``, and both pysam and bcftools need
    one. The bundled references live inside the installed package, which is
    read-only in a conda/system install and should not accumulate build
    artifacts even when it is writable. Staging a copy per run keeps the
    package pristine and makes concurrent runs independent.
    """
    src = organism_fasta_path(organism)
    if not src.exists():
        raise FileNotFoundError(
            f"Bundled 23S reference for '{organism}' not found at {src}."
        )
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    dest = workdir / src.name
    if not dest.exists() or dest.stat().st_size != src.stat().st_size:
        shutil.copyfile(src, dest)
        # A stale index next to a freshly copied FASTA would be silently wrong.
        fai = Path(str(dest) + ".fai")
        if fai.exists():
            fai.unlink()
    return dest


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
