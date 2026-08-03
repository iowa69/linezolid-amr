"""23S rRNA read mapping, pileup at canonical LZD positions, and full-locus VCF.

The pileup is the heart of the tool. rRNA is encoded by several near-identical
rrn operons, so reads from every copy collapse onto the single reference copy we
map against. A resistance mutation present in only some operons therefore shows
up as an intermediate allele fraction rather than a clean homozygous call —
that intermediate fraction is exactly the heteroresistance signal we are after.

Because the signal of interest lives at low allele fractions, it shares that
space with sequencing and alignment artifacts. Each candidate allele is
therefore annotated with the evidence needed to tell the two apart: strand
distribution, base qualities, a binomial confidence interval on the fraction,
and the number of rrn operons the fraction implies.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Iterable

import pysam

from linezolid_amr import references as references_mod
from linezolid_amr import stats
from linezolid_amr.references import (
    OrganismRef,
    get_linezolid_positions,
    get_organism,
    organism_bed_path,
    organism_fasta_path,
    organism_position_map_path,
    rrn_copy_number,
    stage_reference,
)


# Minimum allele frequency required to call a positive linezolid resistance.
# Set to 15 % so that for the worst-case organism (S. aureus, 5 rrn operons)
# a mutation must be present in roughly one full operon copy (1/5 = 20 %)
# before we flag the sample. Lower observations are still reported in the
# output (with their AF) but do NOT trigger a POS call.
DEFAULT_MIN_AF = 0.15
DEFAULT_MIN_DEPTH = 20

# Read-level quality gates. min_base_quality is enforced inside htslib, which
# drops failing bases from the pileup column entirely.
DEFAULT_MIN_BASEQ = 13
DEFAULT_MIN_MAPQ = 20

# An allele needs a few independent reads before its fraction means anything;
# 3 is the count a 15 % fraction produces at the 20x minimum depth.
DEFAULT_MIN_ALT_READS = 3

# Strand-bias rejection. Both conditions must hold, because at high depth a
# Fisher test alone flags biologically ordinary imbalance: the p-value must be
# significant AND essentially all alt reads must sit on one strand.
DEFAULT_STRAND_BIAS_P = 1e-3
DEFAULT_STRAND_BIAS_MINOR_FRAC = 0.05

# Matches the -d passed to bcftools mpileup on the VCF path, so the two
# outputs describe the same reads.
DEFAULT_MAX_DEPTH = 100_000


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
    # Reads whose base passed the quality gates and is a real A/C/G/T — this is
    # the denominator every allele fraction is computed against. `depth` is the
    # broader count including deletions, kept for continuity with earlier output.
    base_depth: int = 0
    ref_fwd: int = 0
    ref_rev: int = 0
    evidence_tier: str = "established"
    # Populated when a resistance allele is present but rejected by a filter,
    # so a negative call can always be explained.
    filtered_resistance: list[str] = field(default_factory=list)

    def to_row(self) -> dict:
        d = asdict(self)
        d["counts"] = ";".join(f"{b}={c}" for b, c in self.counts.items() if c)
        d["alt_alleles"] = ";".join(
            f"{a['base']}:{a['count']}:{a['af']:.4f}{'*' if a['resistance'] else ''}"
            for a in self.alt_alleles
        )
        d["filtered_resistance"] = ";".join(self.filtered_resistance)
        return d


@dataclass(frozen=True)
class PileupFilters:
    """Tunable gates applied to every candidate allele."""

    min_af: float = DEFAULT_MIN_AF
    min_depth: int = DEFAULT_MIN_DEPTH
    min_baseq: int = DEFAULT_MIN_BASEQ
    min_mapq: int = DEFAULT_MIN_MAPQ
    min_alt_reads: int = DEFAULT_MIN_ALT_READS
    strand_bias_p: float = DEFAULT_STRAND_BIAS_P
    strand_bias_minor_frac: float = DEFAULT_STRAND_BIAS_MINOR_FRAC
    apply_strand_filter: bool = True
    # Ceiling on pileup depth. rRNA loci exceed pysam's 8000 default.
    max_depth: int = DEFAULT_MAX_DEPTH
    # Weakest evidence tier still allowed to produce a positive call.
    evidence_tier: str = references_mod.DEFAULT_EVIDENCE_TIER


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
        # Close our copy so minimap2 sees SIGPIPE if samtools dies early.
        mm.stdout.close()
        sort.communicate()
        mm.wait()
    if mm.returncode != 0 or sort.returncode != 0:
        raise RuntimeError(
            f"mapping failed (minimap2={mm.returncode}, samtools={sort.returncode}); see {log}"
        )
    pysam.index(str(bam))
    return bam


def ensure_fasta_index(fasta: Path) -> None:
    """Index a FASTA in place.

    Callers should hand this a *staged* copy (see references.stage_reference);
    the bundled references live inside the installed package, which is commonly
    read-only and must not accumulate .fai files.
    """
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


@dataclass(frozen=True)
class BedTarget:
    contig: str
    species_position: int  # 1-based
    ref_base: str
    resistance_bases: tuple[str, ...]
    drug: str
    evidence_tier: str = references_mod.TIER_ESTABLISHED


def _parse_bed(bed: Path) -> list[BedTarget]:
    """Parse the LZD target BED.

    The evidence-tier column became mandatory in v0.2. A BED without it is
    refused rather than assumed: guessing "established" would let archaeal and
    engineered substitutions drive clinical positives and silently turn
    ``--evidence-tier established`` into a no-op, while guessing the weakest
    tier would suppress every call. Neither wrong answer is acceptable when the
    fix — regenerating the file — is one command.
    """
    targets: list[BedTarget] = []
    with bed.open() as fh:
        for lineno, line in enumerate(fh, 1):
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                raise ValueError(
                    f"{bed}:{lineno} has {len(parts)} columns; the evidence-tier "
                    f"column has been required since v0.2. This file was written by "
                    f"an older version — regenerate it with 'linezolid-amr "
                    f"fetch-references', or delete the override directory to fall "
                    f"back on the references bundled in the package."
                )
            tier = parts[7].strip()
            if tier not in references_mod.EVIDENCE_TIERS:
                raise ValueError(
                    f"{bed}:{lineno} declares unknown evidence tier {tier!r}; "
                    f"expected one of {', '.join(references_mod.EVIDENCE_TIERS)}."
                )
            targets.append(
                BedTarget(
                    contig=parts[0],
                    species_position=int(parts[2]),  # 1-based end == 1-based position
                    ref_base=parts[4],
                    resistance_bases=tuple(parts[5].split(",")),
                    drug=parts[6] or "linezolid",
                    evidence_tier=tier,
                )
            )
    return targets


def _collect_observations(
    bamfh: pysam.AlignmentFile,
    contig: str,
    pos0: int,
    filters: PileupFilters,
    fastafile: pysam.FastaFile | None = None,
) -> tuple[dict[str, int], dict[str, dict[str, int]], int]:
    """Tally the bases seen at a single reference position.

    Returns ``(counts, strand_counts, raw_depth)`` where ``strand_counts`` maps
    each base to ``{"fwd": n, "rev": n}``.

    htslib applies ``min_base_quality`` while building the column, so bases
    failing it never reach us. Mapping quality is filtered here.
    """
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "DEL": 0}
    strand: dict[str, dict[str, int]] = {b: {"fwd": 0, "rev": 0} for b in "ACGT"}
    raw_depth = 0

    for col in bamfh.pileup(
        contig=contig,
        start=pos0,
        end=pos0 + 1,
        truncate=True,
        min_base_quality=filters.min_baseq,
        min_mapping_quality=filters.min_mapq,
        stepper="samtools",
        ignore_overlaps=True,
        # rRNA is multi-copy and every operon's reads land here, so this locus
        # reaches depths that pysam's 8000 default silently truncates — which
        # would under-report depth, widen the confidence interval, and disagree
        # with the companion VCF (bcftools runs at -d 100000).
        max_depth=filters.max_depth,
        # Supplying the reference enables BAQ, which down-weights bases in
        # poorly-anchored alignments near indels. bcftools applies it on the VCF
        # path by default; without it here the two shipped outputs could
        # disagree about the same position.
        fastafile=fastafile,
    ):
        for pread in col.pileups:
            aln = pread.alignment
            if aln.mapping_quality < filters.min_mapq:
                continue
            if pread.is_del or pread.is_refskip:
                counts["DEL"] += 1
                raw_depth += 1
                continue
            qpos = pread.query_position
            if qpos is None:
                # Base failed the htslib quality gate.
                continue
            base = aln.query_sequence[qpos].upper()
            raw_depth += 1
            if base in "ACGT":
                counts[base] += 1
                strand[base]["rev" if aln.is_reverse else "fwd"] += 1
            else:
                counts["N"] += 1
    return counts, strand, raw_depth


def _evaluate_alt(
    base: str,
    counts: dict[str, int],
    strand: dict[str, dict[str, int]],
    ref_base: str,
    base_depth: int,
    raw_depth: int,
    is_known_resistance: bool,
    copy_number: int,
    filters: PileupFilters,
) -> dict:
    """Build the annotated record for one candidate alternate allele."""
    count = counts[base]
    af = count / base_depth if base_depth else 0.0
    alt_fwd = strand[base]["fwd"]
    alt_rev = strand[base]["rev"]
    ref_fwd = strand.get(ref_base, {}).get("fwd", 0)
    ref_rev = strand.get(ref_base, {}).get("rev", 0)

    sb_p, minor_frac = stats.strand_bias(ref_fwd, ref_rev, alt_fwd, alt_rev)
    ci_low, ci_high = stats.wilson_interval(count, base_depth)
    operons, expected_af = stats.nearest_operon_count(af, copy_number)

    # Strand bias is only meaningful where the position actually has coverage on
    # both strands; single-end libraries and one-sided tiling legitimately give
    # one-strand data for every allele, reference included.
    #
    # Judge that on the whole column, not on the reference allele alone. A
    # near-complete artifact leaves no reference read on the artifact's strand,
    # so a reference-only test would see ref_fwd == 0, conclude the data is
    # single-stranded and skip the filter — making it weaker precisely as the
    # artifact grows more extreme. Summing both alleles keeps genuinely
    # single-end data exempt (one strand is empty overall) while holding a
    # two-stranded library to the test.
    two_sided_coverage = (ref_fwd + alt_fwd) > 0 and (ref_rev + alt_rev) > 0

    reasons: list[str] = []
    if raw_depth < filters.min_depth:
        reasons.append("low_depth")
    if af < filters.min_af:
        reasons.append("low_af")
    if count < filters.min_alt_reads:
        reasons.append("few_alt_reads")
    if (
        filters.apply_strand_filter
        and two_sided_coverage
        and sb_p < filters.strand_bias_p
        and minor_frac < filters.strand_bias_minor_frac
    ):
        reasons.append("strand_bias")

    passes_threshold = af >= filters.min_af and raw_depth >= filters.min_depth
    return {
        "base": base,
        "count": count,
        "af": af,
        "resistance": is_known_resistance,
        "passes_threshold": passes_threshold,
        "passes_filters": not reasons,
        "filters": reasons,
        "fwd": alt_fwd,
        "rev": alt_rev,
        "strand_bias_p": sb_p,
        "strand_minor_frac": minor_frac,
        "af_ci_low": ci_low,
        "af_ci_high": ci_high,
        "est_mutated_operons": operons,
        "expected_af_for_operons": expected_af,
        "rrn_copy_number": copy_number,
    }


def pileup_at_positions(
    bam: Path,
    reference_fasta: Path,
    organism: str,
    min_af: float = DEFAULT_MIN_AF,
    min_depth: int = DEFAULT_MIN_DEPTH,
    filters: PileupFilters | None = None,
) -> list[PileupCall]:
    """Pileup at canonical LZD positions; report annotated allele frequencies.

    ``min_af``/``min_depth`` are honoured for backwards compatibility; pass a
    :class:`PileupFilters` to control the full gate set.
    """
    if filters is None:
        filters = PileupFilters(min_af=min_af, min_depth=min_depth)

    bed = organism_bed_path(organism)
    pos_map = organism_position_map_path(organism)
    if not bed.exists() or not pos_map.exists():
        raise FileNotFoundError(
            f"Missing BED or position map for {organism}. Run 'linezolid-amr fetch-references'."
        )
    sp_to_ecoli = _read_position_map(pos_map)
    targets = _parse_bed(bed)
    copy_number = rrn_copy_number(organism)

    ensure_fasta_index(reference_fasta)
    bamfh = pysam.AlignmentFile(str(bam), "rb")
    fa = pysam.FastaFile(str(reference_fasta))

    available = set(fa.references)
    calls: list[PileupCall] = []

    for target in targets:
        # Honour the contig named in the BED. Falling back to the first record
        # only makes sense for the single-record references we ship.
        if target.contig in available:
            contig_name = target.contig
        elif len(fa.references) == 1:
            contig_name = fa.references[0]
        else:
            raise ValueError(
                f"BED contig '{target.contig}' not found in {reference_fasta} "
                f"(has: {', '.join(sorted(available))})."
            )

        species_pos_1b = target.species_position
        ref_base = target.ref_base
        pos0 = species_pos_1b - 1
        counts, strand, raw_depth = _collect_observations(
            bamfh, contig_name, pos0, filters, fastafile=fa
        )

        base_depth = sum(counts[b] for b in "ACGT")
        ecoli_pos = sp_to_ecoli.get(species_pos_1b, -1)
        resistance_upper = {r.upper() for r in target.resistance_bases}
        # RNA literature writes uracil; the sequencing alphabet has thymine.
        if "U" in resistance_upper:
            resistance_upper.add("T")

        # Evidence tier gates whether this position may drive a positive call.
        tier_callable = references_mod.tier_is_callable(
            target.evidence_tier, filters.evidence_tier
        )

        alt_alleles: list[dict] = []
        for base in "ACGT":
            if base == ref_base.upper():
                continue
            if counts[base] == 0:
                continue
            is_known_resistance = base in resistance_upper
            rec = _evaluate_alt(
                base, counts, strand, ref_base.upper(), base_depth, raw_depth,
                is_known_resistance, copy_number, filters,
            )
            rec["evidence_tier"] = target.evidence_tier
            if is_known_resistance and not tier_callable:
                # Reported with its fraction, but not evidence enough to call.
                rec["filters"] = rec["filters"] + ["evidence_tier"]
                rec["passes_filters"] = False
            # Unrelated low-frequency noise is dropped from the output; known
            # resistance alleles are always kept so the reader sees the fraction
            # even when it is sub-threshold.
            if not is_known_resistance and rec["af"] < filters.min_af:
                continue
            alt_alleles.append(rec)

        is_resistance = any(a["resistance"] and a["passes_filters"] for a in alt_alleles)
        # Only surface rejections that are actually informative. A resistance
        # base sitting at 0.1% of an 800x pileup is ordinary sequencing noise;
        # announcing it as a "rejected resistance allele" would cry wolf on
        # essentially every sample. What deserves attention is an allele that
        # cleared the frequency and depth bar and was then rejected on quality
        # or evidence grounds.
        _ABUNDANCE_GATES = {"low_af", "low_depth", "few_alt_reads"}
        filtered_resistance = [
            f"{a['base']}({','.join(a['filters'])})"
            for a in alt_alleles
            if a["resistance"]
            and not a["passes_filters"]
            and not (_ABUNDANCE_GATES & set(a["filters"]))
        ]

        calls.append(
            PileupCall(
                organism=organism,
                ref_contig=contig_name,
                species_position=species_pos_1b,
                ecoli_position=ecoli_pos,
                ref_base=ref_base,
                depth=raw_depth,
                counts=counts,
                alt_alleles=alt_alleles,
                is_resistance=is_resistance,
                drug=target.drug,
                base_depth=base_depth,
                ref_fwd=strand.get(ref_base.upper(), {}).get("fwd", 0),
                ref_rev=strand.get(ref_base.upper(), {}).get("rev", 0),
                filtered_resistance=filtered_resistance,
                evidence_tier=target.evidence_tier,
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
    _need("bcftools")
    ensure_fasta_index(reference_fasta)
    outdir.mkdir(parents=True, exist_ok=True)
    vcf = outdir / f"{sample}.23S.vcf.gz"
    log = outdir / f"{sample}.23S.vcf.log"

    mpileup = [
        "bcftools", "mpileup",
        "-f", str(reference_fasta),
        "-a", "AD,DP",
        "-d", "100000",
        "--threads", str(threads),
        str(bam),
    ]
    call = [
        "bcftools", "call",
        "-mv", "-Oz", "--ploidy", "1",
        "-o", str(vcf),
    ]

    with log.open("w") as logfh:
        logfh.write("$ " + " ".join(mpileup) + " | " + " ".join(call) + "\n")
        logfh.flush()
        # No shell: paths containing spaces or shell metacharacters are safe.
        mp = subprocess.Popen(mpileup, stdout=subprocess.PIPE, stderr=logfh)
        cl = subprocess.Popen(call, stdin=mp.stdout, stdout=logfh, stderr=logfh)
        mp.stdout.close()
        cl.communicate()
        mp.wait()
        if mp.returncode != 0 or cl.returncode != 0:
            raise RuntimeError(
                f"bcftools failed (mpileup={mp.returncode}, call={cl.returncode}); see {log}"
            )
        idx = subprocess.run(["bcftools", "index", "-f", str(vcf)], stdout=logfh, stderr=logfh)
    if idx.returncode != 0:
        raise RuntimeError(f"bcftools index failed (exit {idx.returncode}); see {log}")
    return vcf


_PILEUP_COLUMNS = [
    "organism", "ref_contig", "species_position", "ecoli_position", "ref_base",
    "depth", "base_depth", "counts", "alt_alleles", "is_resistance", "drug",
    "evidence_tier", "ref_fwd", "ref_rev", "filtered_resistance", "note",
]


def write_pileup_tsv(calls: Iterable[PileupCall], path: Path) -> None:
    rows = [c.to_row() for c in calls]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=_PILEUP_COLUMNS, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


_EVIDENCE_COLUMNS = [
    "organism", "ecoli_position", "species_position", "ref_base", "alt_base",
    "is_resistance_allele", "evidence_tier", "depth", "base_depth", "alt_count", "af",
    "af_ci_low", "af_ci_high", "alt_fwd", "alt_rev", "ref_fwd", "ref_rev",
    "strand_bias_p", "strand_minor_frac", "est_mutated_operons",
    "rrn_copy_number", "expected_af_for_operons", "passes_filters", "filters",
]


def write_evidence_tsv(calls: Iterable[PileupCall], path: Path) -> None:
    """One row per observed alternate allele, with the full statistical evidence.

    This is the file to inspect when deciding whether a borderline fraction is
    real heteroresistance or an artifact.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=_EVIDENCE_COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for c in calls:
            for a in c.alt_alleles:
                w.writerow({
                    "organism": c.organism,
                    "ecoli_position": c.ecoli_position,
                    "species_position": c.species_position,
                    "ref_base": c.ref_base,
                    "alt_base": a["base"],
                    "is_resistance_allele": "YES" if a["resistance"] else "NO",
                    "evidence_tier": a.get("evidence_tier", c.evidence_tier),
                    "depth": c.depth,
                    "base_depth": c.base_depth,
                    "alt_count": a["count"],
                    "af": f"{a['af']:.4f}",
                    "af_ci_low": f"{a['af_ci_low']:.4f}",
                    "af_ci_high": f"{a['af_ci_high']:.4f}",
                    "alt_fwd": a["fwd"],
                    "alt_rev": a["rev"],
                    "ref_fwd": c.ref_fwd,
                    "ref_rev": c.ref_rev,
                    "strand_bias_p": f"{a['strand_bias_p']:.3g}",
                    "strand_minor_frac": f"{a['strand_minor_frac']:.4f}",
                    "est_mutated_operons": a["est_mutated_operons"],
                    "rrn_copy_number": a["rrn_copy_number"],
                    "expected_af_for_operons": f"{a['expected_af_for_operons']:.4f}",
                    "passes_filters": "YES" if a["passes_filters"] else "NO",
                    "filters": ";".join(a["filters"]) or "PASS",
                })
