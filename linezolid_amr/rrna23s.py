"""23S rRNA read mapping, pileup at canonical LZD positions, and full-locus VCF."""

from __future__ import annotations

import csv
import shutil
import subprocess
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

import pysam

from linezolid_amr.references import (
    OrganismRef,
    get_linezolid_positions,
    get_organism,
    organism_bed_path,
    organism_fasta_path,
    organism_position_map_path,
)


# Allele frequency below which we consider the position effectively wild type.
# Linezolid heteroresistance literature commonly uses 1-3% as the lower bound
# because (a) Illumina error around 0.5-1% per position, and (b) most rRNA
# operon copies must remain wild type for fitness.
DEFAULT_MIN_AF = 0.01
DEFAULT_MIN_DEPTH = 20


def _need(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(
            f"'{tool}' not found on PATH. Install via 'conda install -c bioconda {tool}'."
        )


@dataclass
class PileupCall:
    organism: str
    ref_contig: str
    species_position: int
    ecoli_position: int
    ref_base: str
    depth: int
    counts: dict[str, int]
    alt_alleles: list[dict]
    is_resistance: bool
    drug: str
    note: str | None = None

    def to_row(self) -> dict:
        d = asdict(self)
        d["counts"] = ";".join(f"{b}={c}" for b, c in self.counts.items() if c)
        d["alt_alleles"] = ";".join(
            f"{a['base']}:{a['count']}:{a['af']:.4f}{'*' if a['resistance'] else ''}"
            for a in self.alt_alleles
        )
        return d


def map_reads(
    reference_fasta: Path,
    r1: Path,
    r2: Path | None,
    outdir: Path,
    threads: int = 4,
    sample: str = "sample",
) -> Path:
    """minimap2 sr preset -> samtools sort -> indexed BAM."""
    for tool in ("minimap2", "samtools"):
        _need(tool)
    outdir.mkdir(parents=True, exist_ok=True)
    bam = outdir / f"{sample}.23S.bam"

    mm_cmd = [
        "minimap2",
        "-ax",
        "sr",
        "-t",
        str(threads),
        "-R",
        f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA",
        str(reference_fasta),
        str(r1),
    ]
    if r2:
        mm_cmd.append(str(r2))

    sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", str(bam), "-"]
    log = outdir / f"{sample}.23S.map.log"

    with log.open("w") as logfh:
        logfh.write("$ " + " ".join(mm_cmd) + " | " + " ".join(sort_cmd) + "\n")
        logfh.flush()
        mm = subprocess.Popen(mm_cmd, stdout=subprocess.PIPE, stderr=logfh)
        sort = subprocess.Popen(sort_cmd, stdin=mm.stdout, stdout=logfh, stderr=logfh)
        mm.stdout.close()
        sort.communicate()
        mm.wait()
    if mm.returncode != 0 or sort.returncode != 0:
        raise RuntimeError(f"mapping failed (minimap2={mm.returncode}, samtools={sort.returncode}); see {log}")
    pysam.index(str(bam))
    return bam


def ensure_fasta_index(fasta: Path) -> None:
    if not Path(str(fasta) + ".fai").exists():
        pysam.faidx(str(fasta))


def _read_position_map(map_path: Path) -> dict[int, int]:
    """species_position (1-based) -> ecoli_position (1-based)."""
    out: dict[int, int] = {}
    with map_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sp = row["species_position"]
            if sp == "." or sp == "":
                continue
            out[int(sp)] = int(row["ecoli_position"])
    return out


def pileup_at_positions(
    bam: Path,
    reference_fasta: Path,
    organism: str,
    min_af: float = DEFAULT_MIN_AF,
    min_depth: int = DEFAULT_MIN_DEPTH,
) -> list[PileupCall]:
    """Pileup at canonical LZD positions; report allele frequencies."""
    bed = organism_bed_path(organism)
    pos_map = organism_position_map_path(organism)
    if not bed.exists() or not pos_map.exists():
        raise FileNotFoundError(
            f"Missing BED or position map for {organism}. Run 'linezolid-amr fetch-references'."
        )
    sp_to_ecoli = _read_position_map(pos_map)

    # Parse BED -> list of (species_pos_1b, ref_base, resistance_bases, drug)
    targets: list[tuple[str, int, str, tuple[str, ...], str]] = []
    with bed.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            contig = parts[0]
            species_pos = int(parts[2])  # 1-based end == 1-based position
            ref_base = parts[4]
            resistance = tuple(parts[5].split(","))
            drug = parts[6] if len(parts) > 6 else "linezolid"
            targets.append((contig, species_pos, ref_base, resistance, drug))

    ensure_fasta_index(reference_fasta)
    bamfh = pysam.AlignmentFile(str(bam), "rb")
    fa = pysam.FastaFile(str(reference_fasta))

    contig_name = fa.references[0]
    calls: list[PileupCall] = []

    for _bed_contig, species_pos_1b, ref_base, resistance_bases, drug in targets:
        pos0 = species_pos_1b - 1
        counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "DEL": 0}
        depth = 0
        # truncate=True restricts to single position; min_base_quality=13 is samtools default
        for col in bamfh.pileup(
            contig=contig_name,
            start=pos0,
            end=pos0 + 1,
            truncate=True,
            min_base_quality=13,
            stepper="samtools",
            ignore_overlaps=True,
        ):
            for pread in col.pileups:
                if pread.is_del or pread.is_refskip:
                    counts["DEL"] += 1
                    depth += 1
                    continue
                qpos = pread.query_position
                if qpos is None:
                    continue
                base = pread.alignment.query_sequence[qpos].upper()
                if base in counts:
                    counts[base] += 1
                else:
                    counts["N"] += 1
                depth += 1

        ecoli_pos = sp_to_ecoli.get(species_pos_1b, -1)
        alt_alleles = []
        total_non_n = sum(counts[b] for b in "ACGT")
        for base in "ACGT":
            if base == ref_base.upper():
                continue
            c = counts[base]
            if c == 0:
                continue
            af = c / total_non_n if total_non_n else 0.0
            if af < min_af and depth >= min_depth:
                continue
            alt_alleles.append(
                {
                    "base": base,
                    "count": c,
                    "af": af,
                    "resistance": base in [r.upper() for r in resistance_bases],
                }
            )

        is_resistance = any(a["resistance"] and depth >= min_depth for a in alt_alleles)
        calls.append(
            PileupCall(
                organism=organism,
                ref_contig=contig_name,
                species_position=species_pos_1b,
                ecoli_position=ecoli_pos,
                ref_base=ref_base,
                depth=depth,
                counts=counts,
                alt_alleles=alt_alleles,
                is_resistance=is_resistance,
                drug=drug,
            )
        )

    bamfh.close()
    fa.close()
    return calls


def call_full_vcf(
    bam: Path,
    reference_fasta: Path,
    outdir: Path,
    sample: str = "sample",
    threads: int = 4,
) -> Path:
    """Run bcftools mpileup | call -mv to produce a full 23S VCF (any variant)."""
    for tool in ("bcftools",):
        _need(tool)
    ensure_fasta_index(reference_fasta)
    outdir.mkdir(parents=True, exist_ok=True)
    vcf = outdir / f"{sample}.23S.vcf.gz"
    log = outdir / f"{sample}.23S.vcf.log"

    cmd = (
        f"bcftools mpileup -f {reference_fasta} -a AD,DP -d 100000 --threads {threads} {bam} "
        f"| bcftools call -mv -Oz --ploidy 1 -o {vcf} && bcftools index {vcf}"
    )
    with log.open("w") as logfh:
        logfh.write("$ " + cmd + "\n")
        logfh.flush()
        proc = subprocess.run(cmd, shell=True, stdout=logfh, stderr=logfh)
    if proc.returncode != 0:
        raise RuntimeError(f"bcftools failed (exit {proc.returncode}); see {log}")
    return vcf


def write_pileup_tsv(calls: Iterable[PileupCall], path: Path) -> None:
    rows = [c.to_row() for c in calls]
    if not rows:
        path.write_text(
            "organism\tref_contig\tspecies_position\tecoli_position\tref_base\tdepth\t"
            "counts\talt_alleles\tis_resistance\tdrug\tnote\n"
        )
        return
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
