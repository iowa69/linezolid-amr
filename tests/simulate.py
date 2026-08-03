"""Deterministic paired-end read simulator for 23S heteroresistance testing.

Generates FASTQ pairs from a reference 23S sequence in which a chosen subset of
reads carry a substitution at a given position. This lets the test suite assert
that the pileup recovers a *known* alt-allele frequency, which is the property
the whole heteroresistance claim rests on.

The simulator is intentionally simple (uniform fragment placement, independent
per-base errors) — it is a correctness harness, not an attempt to reproduce
Illumina error structure.
"""

from __future__ import annotations

import gzip
import random
from pathlib import Path

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def read_fasta_single(path: Path) -> tuple[str, str]:
    """Read a single-record FASTA; return (name, sequence)."""
    name = ""
    chunks: list[str] = []
    with Path(path).open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    break
                name = line[1:].split()[0]
            else:
                chunks.append(line.upper())
    return name, "".join(chunks)


def _phred(q: int) -> str:
    return chr(min(q, 93) + 33)


def simulate_reads(
    reference: str,
    out_r1: Path,
    out_r2: Path,
    *,
    variants: dict[int, tuple[str, float]] | None = None,
    depth: int = 200,
    read_len: int = 150,
    frag_len: int = 350,
    frag_sd: int = 40,
    error_rate: float = 0.001,
    base_quality: int = 37,
    error_quality: int = 37,
    seed: int = 1234,
    sample: str = "sim",
) -> tuple[Path, Path]:
    """Simulate paired-end reads over *reference*.

    ``variants`` maps a **1-based** reference position to ``(alt_base, fraction)``.
    Each generated fragment independently carries the alt base with probability
    ``fraction``, so the observed alt-allele frequency converges on ``fraction``.

    ``error_quality`` is the quality score assigned to bases corrupted by the
    error model; set it low to emulate the low-confidence calls that a
    base-quality filter is supposed to discard.
    """
    variants = variants or {}
    rng = random.Random(seed)
    ref_len = len(reference)
    n_fragments = max(1, int(depth * ref_len / (2 * read_len)))

    # Decide per-fragment genotype up front so the alt fraction is exact-ish and
    # independent of which fragments happen to span the site.
    alt_calls: dict[int, list[bool]] = {
        pos: [rng.random() < frac for _ in range(n_fragments)]
        for pos, (_b, frac) in variants.items()
    }

    r1_records: list[str] = []
    r2_records: list[str] = []

    for i in range(n_fragments):
        flen = max(read_len, int(rng.gauss(frag_len, frag_sd)))
        flen = min(flen, ref_len)
        start = rng.randint(0, ref_len - flen)  # 0-based
        frag = list(reference[start : start + flen])

        # Apply variant genotypes to this fragment
        for pos, (alt_base, _frac) in variants.items():
            idx = (pos - 1) - start
            if 0 <= idx < flen and alt_calls[pos][i]:
                frag[idx] = alt_base

        quals = [base_quality] * flen
        # Independent sequencing errors
        if error_rate > 0:
            for j in range(flen):
                if rng.random() < error_rate:
                    frag[j] = rng.choice([b for b in "ACGT" if b != frag[j]])
                    quals[j] = error_quality

        frag_seq = "".join(frag)
        r1_seq = frag_seq[:read_len]
        r1_qual = "".join(_phred(q) for q in quals[:read_len])
        r2_seq = revcomp(frag_seq[-read_len:])
        r2_qual = "".join(_phred(q) for q in quals[-read_len:][::-1])

        name = f"{sample}_frag{i}"
        r1_records.append(f"@{name}/1\n{r1_seq}\n+\n{r1_qual}\n")
        r2_records.append(f"@{name}/2\n{r2_seq}\n+\n{r2_qual}\n")

    out_r1 = Path(out_r1)
    out_r2 = Path(out_r2)
    out_r1.parent.mkdir(parents=True, exist_ok=True)
    _write_fastq(out_r1, r1_records)
    _write_fastq(out_r2, r2_records)
    return out_r1, out_r2


def _write_fastq(path: Path, records: list[str]) -> None:
    data = "".join(records)
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(data)
    else:
        path.write_text(data)


def simulate_strand_biased_reads(
    reference: str,
    out_r1: Path,
    out_r2: Path,
    *,
    position: int,
    alt_base: str,
    fraction: float,
    depth: int = 200,
    read_len: int = 150,
    seed: int = 99,
    sample: str = "sb",
) -> tuple[Path, Path]:
    """Simulate an artifact: the alt allele appears on R1 (forward) reads only.

    Real heteroresistance is carried by the template and therefore appears on
    both strands. A one-strand-only alt is the classic signature of a sequencing
    or alignment artifact, and is what the strand-bias filter must reject.
    """
    rng = random.Random(seed)
    ref_len = len(reference)
    n_fragments = max(1, int(depth * ref_len / (2 * read_len)))
    r1_records: list[str] = []
    r2_records: list[str] = []

    for i in range(n_fragments):
        flen = min(350, ref_len)
        start = rng.randint(0, ref_len - flen)
        frag = list(reference[start : start + flen])
        idx = (position - 1) - start

        r1_span = range(0, min(read_len, flen))
        carries = rng.random() < fraction
        # Mutate only within the R1 footprint so R2 never sees the alt allele.
        if carries and 0 <= idx < flen and idx in r1_span:
            frag[idx] = alt_base

        frag_seq = "".join(frag)
        qual = "".join(_phred(37) for _ in range(read_len))
        r1_records.append(f"@{sample}_f{i}/1\n{frag_seq[:read_len]}\n+\n{qual}\n")
        r2_records.append(f"@{sample}_f{i}/2\n{revcomp(frag_seq[-read_len:])}\n+\n{qual}\n")

    _write_fastq(Path(out_r1), r1_records)
    _write_fastq(Path(out_r2), r2_records)
    return Path(out_r1), Path(out_r2)
