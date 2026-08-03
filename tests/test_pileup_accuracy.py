"""End-to-end accuracy tests for the 23S heteroresistance pileup.

These run the *real* pipeline — minimap2 mapping followed by the pysam pileup —
against simulated reads whose alt-allele fraction is known exactly. They are the
only tests that can substantiate the tool's central claim: that a resistance
allele present in a minority of rrn operons is recovered at the right frequency.

Skipped automatically when minimap2/samtools are absent.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from linezolid_amr import references as ref_mod
from linezolid_amr import rrna23s as rrna
from tests.simulate import read_fasta_single, simulate_reads, simulate_strand_biased_reads

pytestmark = pytest.mark.integration

HAVE_TOOLS = shutil.which("minimap2") is not None and shutil.which("samtools") is not None
requires_tools = pytest.mark.skipif(
    not HAVE_TOOLS, reason="minimap2/samtools not on PATH"
)

ORGANISM = "Enterococcus_faecium"
# Species coordinate of E. coli G2576 for E. faecium, from the bundled BED.
G2576_SPECIES_POS = 2590


@pytest.fixture(scope="module")
def reference() -> tuple[Path, str, str]:
    fasta = ref_mod.organism_fasta_path(ORGANISM)
    name, seq = read_fasta_single(fasta)
    return fasta, name, seq


def _call_at(calls, ecoli_pos: int):
    for c in calls:
        if c.ecoli_position == ecoli_pos:
            return c
    raise AssertionError(f"no call at E. coli position {ecoli_pos}")


def _alt_af(call, base: str) -> float:
    for a in call.alt_alleles:
        if a["base"] == base:
            return a["af"]
    return 0.0


def _run(tmp_path: Path, reference, *, variants, depth=200, seed=1, **kw):
    fasta, _name, seq = reference
    r1 = tmp_path / "sim_R1.fastq.gz"
    r2 = tmp_path / "sim_R2.fastq.gz"
    simulate_reads(seq, r1, r2, variants=variants, depth=depth, seed=seed, **kw)
    bam = rrna.map_reads(fasta, r1, r2, tmp_path, threads=2, sample="sim")
    return rrna.pileup_at_positions(bam, fasta, ORGANISM, min_af=0.15, min_depth=20)


# --------------------------------------------------------------------------
# Core claim: a minority resistance allele is recovered at the right frequency
# --------------------------------------------------------------------------

@requires_tools
@pytest.mark.parametrize("truth_af", [0.0, 0.10, 0.20, 0.35, 0.50, 0.80, 1.0])
def test_recovered_af_matches_truth(tmp_path, reference, truth_af):
    """Observed AF at G2576T must track the simulated fraction within tolerance."""
    calls = _run(
        tmp_path,
        reference,
        variants={G2576_SPECIES_POS: ("T", truth_af)} if truth_af else {},
        depth=300,
        seed=int(truth_af * 100) + 7,
    )
    call = _call_at(calls, 2576)
    observed = _alt_af(call, "T")
    assert abs(observed - truth_af) < 0.06, (
        f"simulated AF={truth_af}, recovered AF={observed:.4f} "
        f"(depth={call.depth}, counts={call.counts})"
    )


@requires_tools
def test_heteroresistance_is_called_positive(tmp_path, reference):
    """A 25% G2576T minority population must yield a POSITIVE resistance call."""
    calls = _run(tmp_path, reference, variants={G2576_SPECIES_POS: ("T", 0.25)}, depth=300)
    call = _call_at(calls, 2576)
    assert call.is_resistance, f"25% G2576T not called: {call.to_row()}"


@requires_tools
def test_wildtype_is_called_negative(tmp_path, reference):
    """A clean wild-type sample must produce no positive call at any position."""
    calls = _run(tmp_path, reference, variants={}, depth=300)
    positives = [c for c in calls if c.is_resistance]
    assert not positives, f"false positives on wild-type input: {[c.to_row() for c in positives]}"


@requires_tools
def test_all_canonical_positions_are_evaluated(tmp_path, reference):
    """Every position in loci.json must appear in the pileup output."""
    calls = _run(tmp_path, reference, variants={}, depth=100)
    expected = {p.ecoli_position for p in ref_mod.get_linezolid_positions()}
    assert {c.ecoli_position for c in calls} == expected


@requires_tools
def test_depth_is_reported_and_nonzero(tmp_path, reference):
    calls = _run(tmp_path, reference, variants={}, depth=150)
    for c in calls:
        assert c.depth > 0, f"zero depth at E. coli {c.ecoli_position}"


# --------------------------------------------------------------------------
# Noise rejection
# --------------------------------------------------------------------------

@requires_tools
def test_low_quality_bases_do_not_inflate_af(tmp_path, reference):
    """Errors emitted with low base quality must be filtered out of the AF.

    A 5% error rate at Q2 would otherwise push several positions over the 15%
    threshold purely from sequencing noise.
    """
    fasta, _name, seq = reference
    r1 = tmp_path / "lq_R1.fastq.gz"
    r2 = tmp_path / "lq_R2.fastq.gz"
    simulate_reads(
        seq, r1, r2, variants={}, depth=300, seed=5,
        error_rate=0.05, error_quality=2,
    )
    bam = rrna.map_reads(fasta, r1, r2, tmp_path, threads=2, sample="lq")
    calls = rrna.pileup_at_positions(bam, fasta, ORGANISM, min_af=0.15, min_depth=20)
    positives = [c for c in calls if c.is_resistance]
    assert not positives, (
        "low-quality sequencing noise produced a resistance call: "
        f"{[c.to_row() for c in positives]}"
    )


@requires_tools
def test_subthreshold_allele_is_reported_but_not_called(tmp_path, reference):
    """Below --min-af, the allele must still be visible but must not be POS.

    This is the transparency contract: clinicians see the frequency and decide.
    """
    calls = _run(tmp_path, reference, variants={G2576_SPECIES_POS: ("T", 0.05)}, depth=400)
    call = _call_at(calls, 2576)
    assert not call.is_resistance
    assert _alt_af(call, "T") > 0.0, "sub-threshold resistance allele was hidden entirely"


# --------------------------------------------------------------------------
# Artifact rejection — the precision half of the problem
# --------------------------------------------------------------------------

@requires_tools
def test_strand_biased_artifact_is_not_called(tmp_path, reference):
    """An alt allele seen on one strand only must not produce a POSITIVE call.

    Before the strand filter existed this exact input was reported as a
    confident 20% G2576T heteroresistant population.
    """
    fasta, _name, seq = reference
    r1 = tmp_path / "sb_R1.fastq.gz"
    r2 = tmp_path / "sb_R2.fastq.gz"
    simulate_strand_biased_reads(
        seq, r1, r2, position=G2576_SPECIES_POS, alt_base="T", fraction=0.45, depth=300
    )
    bam = rrna.map_reads(fasta, r1, r2, tmp_path, threads=2, sample="sb")
    calls = rrna.pileup_at_positions(bam, fasta, ORGANISM)
    call = _call_at(calls, 2576)

    assert _alt_af(call, "T") > 0.15, "precondition: the artifact should clear the AF threshold"
    assert not call.is_resistance, "strand-biased artifact was called as resistance"
    assert call.filtered_resistance, "rejection happened but was not explained in the output"
    assert any("strand_bias" in f for f in call.filtered_resistance)


@requires_tools
def test_strand_filter_can_be_disabled(tmp_path, reference):
    """--no-strand-filter restores the previous, more permissive behaviour."""
    fasta, _name, seq = reference
    r1 = tmp_path / "sb2_R1.fastq.gz"
    r2 = tmp_path / "sb2_R2.fastq.gz"
    simulate_strand_biased_reads(
        seq, r1, r2, position=G2576_SPECIES_POS, alt_base="T", fraction=0.45, depth=300
    )
    bam = rrna.map_reads(fasta, r1, r2, tmp_path, threads=2, sample="sb2")
    calls = rrna.pileup_at_positions(
        bam, fasta, ORGANISM, filters=rrna.PileupFilters(apply_strand_filter=False)
    )
    assert _call_at(calls, 2576).is_resistance


@requires_tools
def test_genuine_heteroresistance_survives_the_strand_filter(tmp_path, reference):
    """The filter must not cost sensitivity: real minority alleles still pass.

    Guards the tool's headline property — 100% sensitivity including
    sub-consensus populations — against over-filtering.
    """
    for truth_af in (0.18, 0.25, 0.42):
        calls = _run(
            tmp_path, reference,
            variants={G2576_SPECIES_POS: ("T", truth_af)},
            depth=300, seed=int(truth_af * 1000),
        )
        call = _call_at(calls, 2576)
        assert call.is_resistance, (
            f"genuine {truth_af:.0%} heteroresistance rejected: "
            f"filters={call.filtered_resistance}, alt={call.alt_alleles}"
        )


@requires_tools
def test_single_end_reads_do_not_trip_the_strand_filter(tmp_path, reference):
    """With single-end data every read is forward; that is not strand bias.

    Without the two-sided-coverage guard, single-end input would reject every
    allele it found.
    """
    fasta, _name, seq = reference
    r1 = tmp_path / "se_R1.fastq.gz"
    r2 = tmp_path / "se_R2.fastq.gz"
    simulate_reads(seq, r1, r2, variants={G2576_SPECIES_POS: ("T", 0.40)}, depth=300, seed=11)
    # Map R1 only — the reverse mate is simply not supplied.
    bam = rrna.map_reads(fasta, r1, None, tmp_path, threads=2, sample="se")
    calls = rrna.pileup_at_positions(bam, fasta, ORGANISM)
    call = _call_at(calls, 2576)
    assert call.is_resistance, (
        f"single-end input misclassified as strand bias: {call.filtered_resistance}"
    )


# --------------------------------------------------------------------------
# Evidence annotation
# --------------------------------------------------------------------------

@requires_tools
def test_allele_fraction_confidence_interval_covers_truth(tmp_path, reference):
    truth_af = 0.30
    calls = _run(tmp_path, reference, variants={G2576_SPECIES_POS: ("T", truth_af)}, depth=400)
    call = _call_at(calls, 2576)
    alt = next(a for a in call.alt_alleles if a["base"] == "T")
    assert alt["af_ci_low"] <= truth_af <= alt["af_ci_high"], (
        f"95% CI ({alt['af_ci_low']:.3f}, {alt['af_ci_high']:.3f}) excludes truth {truth_af}"
    )


@requires_tools
def test_operon_estimate_is_sensible(tmp_path, reference):
    """E. faecium carries 6 rrn operons, so a 50% allele implies about 3 of them."""
    calls = _run(tmp_path, reference, variants={G2576_SPECIES_POS: ("T", 0.50)}, depth=400)
    alt = next(a for a in _call_at(calls, 2576).alt_alleles if a["base"] == "T")
    assert alt["rrn_copy_number"] == 6
    assert alt["est_mutated_operons"] == 3


@requires_tools
def test_evidence_tsv_is_written_and_parsable(tmp_path, reference):
    import csv

    calls = _run(tmp_path, reference, variants={G2576_SPECIES_POS: ("T", 0.30)}, depth=300)
    out = tmp_path / "evidence.tsv"
    rrna.write_evidence_tsv(calls, out)
    with out.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    g2576 = [r for r in rows if r["ecoli_position"] == "2576" and r["alt_base"] == "T"]
    assert g2576, "G2576T missing from the evidence table"
    row = g2576[0]
    assert row["is_resistance_allele"] == "YES"
    assert row["filters"] == "PASS"
    assert int(row["alt_fwd"]) > 0 and int(row["alt_rev"]) > 0


@requires_tools
def test_bundled_reference_is_not_modified_by_a_run(tmp_path, reference):
    """A run must not write .fai files into the installed package directory."""
    from linezolid_amr import references as rm

    bundled = rm.organism_fasta_path(ORGANISM)
    fai = Path(str(bundled) + ".fai")
    existed_before = fai.exists()

    staged = rm.stage_reference(ORGANISM, tmp_path / "staged")
    r1 = tmp_path / "st_R1.fastq.gz"
    r2 = tmp_path / "st_R2.fastq.gz"
    _fasta, _name, seq = reference
    simulate_reads(seq, r1, r2, variants={}, depth=60, seed=2)
    bam = rrna.map_reads(staged, r1, r2, tmp_path, threads=2, sample="st")
    rrna.pileup_at_positions(bam, staged, ORGANISM)

    assert Path(str(staged) + ".fai").exists(), "staged copy should carry the index"
    if not existed_before:
        assert not fai.exists(), f"run created an index inside the package: {fai}"
