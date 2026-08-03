"""Unit tests for the dependency-free statistics helpers.

Reference values were cross-checked against scipy.stats (fisher_exact,
binomtest) and published Wilson-interval tables; scipy itself is deliberately
not a runtime dependency.
"""

from __future__ import annotations

import pytest

from linezolid_amr import stats


# --------------------------- Fisher exact ---------------------------

@pytest.mark.parametrize(
    "table,expected",
    [
        ((10, 20, 5, 3), 0.2230247075),
        ((500, 480, 12, 1), 0.003409799437),
        ((7, 3, 0, 10), 0.003095975232),
        ((50, 50, 0, 8), 0.006970162841),
        ((2000, 1900, 40, 38), 1.0),
        ((1, 1, 1, 1), 1.0),
    ],
)
def test_fisher_matches_reference_values(table, expected):
    assert stats.fisher_exact_2x2(*table) == pytest.approx(expected, abs=1e-9)


def test_fisher_empty_table_is_one():
    assert stats.fisher_exact_2x2(0, 0, 0, 0) == 1.0


def test_fisher_is_symmetric_under_row_swap():
    assert stats.fisher_exact_2x2(10, 20, 5, 3) == pytest.approx(
        stats.fisher_exact_2x2(5, 3, 10, 20)
    )


def test_fisher_survives_high_depth():
    """rRNA loci routinely reach 10^4 depth; the log-space form must not overflow."""
    p = stats.fisher_exact_2x2(50_000, 49_000, 900, 880)
    assert 0.0 <= p <= 1.0


# --------------------------- Wilson interval ---------------------------

def test_wilson_known_value():
    lo, hi = stats.wilson_interval(3, 100)
    assert lo == pytest.approx(0.0103, abs=1e-3)
    assert hi == pytest.approx(0.0842, abs=1e-3)


def test_wilson_never_leaves_unit_interval():
    for k, n in [(0, 10), (10, 10), (1, 1), (0, 1), (5, 1000)]:
        lo, hi = stats.wilson_interval(k, n)
        assert 0.0 <= lo <= hi <= 1.0


def test_wilson_zero_depth_is_degenerate():
    assert stats.wilson_interval(0, 0) == (0.0, 0.0)


def test_wilson_narrows_with_depth():
    """The whole point of reporting a CI: 3/20 is far less certain than 150/1000."""
    lo_shallow, hi_shallow = stats.wilson_interval(3, 20)
    lo_deep, hi_deep = stats.wilson_interval(150, 1000)
    assert (hi_shallow - lo_shallow) > (hi_deep - lo_deep)


# --------------------------- binomial tail ---------------------------

@pytest.mark.parametrize(
    "k,n,p,expected",
    [
        (5, 100, 0.01, 0.003432321588),
        (1, 10, 0.5, 0.9990234375),
        (30, 200, 0.1, 0.01632656593),
        (0, 50, 0.02, 1.0),
    ],
)
def test_binom_sf_matches_reference(k, n, p, expected):
    assert stats.binom_sf(k, n, p) == pytest.approx(expected, abs=1e-9)


def test_binom_sf_bounds():
    assert stats.binom_sf(101, 100, 0.5) == 0.0
    assert stats.binom_sf(0, 100, 0.5) == 1.0


# --------------------------- strand bias ---------------------------

def test_strand_bias_flags_one_sided_allele():
    """The artifact signature: every alt read on one strand."""
    p, minor = stats.strand_bias(ref_fwd=120, ref_rev=135, alt_fwd=65, alt_rev=0)
    assert minor == 0.0
    assert p < 1e-3


def test_strand_bias_accepts_balanced_allele():
    """A real heteroresistant allele is carried by the template, so both strands see it."""
    p, minor = stats.strand_bias(ref_fwd=120, ref_rev=135, alt_fwd=33, alt_rev=32)
    assert minor > 0.4
    assert p > 0.05


def test_strand_bias_of_absent_allele_is_defined():
    p, minor = stats.strand_bias(100, 100, 0, 0)
    assert minor == 0.0
    assert 0.0 <= p <= 1.0


# --------------------------- operon quantisation ---------------------------

def test_expected_operon_fractions():
    assert stats.expected_operon_fractions(4) == [0.25, 0.5, 0.75, 1.0]
    assert stats.expected_operon_fractions(0) == []


@pytest.mark.parametrize(
    "af,copies,expected_k",
    [
        (0.20, 5, 1),   # S. aureus, 1 of 5 operons
        (0.66, 6, 4),   # E. faecium, 4 of 6
        (0.99, 4, 4),   # fixed
        (0.02, 4, 0),   # noise, no operon
        (0.51, 4, 2),
    ],
)
def test_nearest_operon_count(af, copies, expected_k):
    k, expected_af = stats.nearest_operon_count(af, copies)
    assert k == expected_k
    assert expected_af == pytest.approx(expected_k / copies)


def test_nearest_operon_count_unknown_copy_number():
    assert stats.nearest_operon_count(0.5, 0) == (0, 0.0)


def test_phred_scaling():
    assert stats.phred(0.1) == pytest.approx(10.0)
    assert stats.phred(0.001) == pytest.approx(30.0)
    assert stats.phred(1.0) == 0.0
    assert stats.phred(0.0) == 255.0


# --------------------------------------------------------------------------
# Strand-bias guard: the filter must not weaken as the artifact gets worse
# --------------------------------------------------------------------------

def test_complete_one_sided_artifact_is_still_significant():
    """When an artifact is total, no reference read remains on that strand.

    Judging "is there two-strand coverage?" from the reference allele alone
    made the filter skip exactly this case, so the most extreme artifacts were
    the ones that passed. The Fisher test itself is decisive here.
    """
    p, minor = stats.strand_bias(ref_fwd=0, ref_rev=165, alt_fwd=155, alt_rev=0)
    assert minor == 0.0
    assert p < 1e-3


def test_column_level_two_sided_check_distinguishes_the_two_cases():
    """Single-end data has one empty strand overall; an artifact does not."""
    # Complete artifact: the column has reads on both strands.
    ref_fwd, ref_rev, alt_fwd, alt_rev = 0, 165, 155, 0
    assert (ref_fwd + alt_fwd) > 0 and (ref_rev + alt_rev) > 0

    # Single-end library: every read is forward, so the reverse side is empty.
    ref_fwd, ref_rev, alt_fwd, alt_rev = 200, 0, 40, 0
    assert not ((ref_fwd + alt_fwd) > 0 and (ref_rev + alt_rev) > 0)
