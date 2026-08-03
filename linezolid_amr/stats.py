"""Small, dependency-free statistics used by the 23S variant filters.

Everything here is implemented in terms of :mod:`math` so the package keeps a
light dependency footprint (no scipy/numpy) — important for a bioconda recipe
that already pulls in pysam and the AMRFinderPlus stack.

All routines work in log space so they stay accurate at the read depths
(10^2–10^5) that rRNA loci reach in a typical WGS run.
"""

from __future__ import annotations

from math import exp, lgamma, log, sqrt

# 97.5th percentile of the standard normal — two-sided 95% coverage.
Z_95 = 1.959963984540054


def _lchoose(n: int, k: int) -> float:
    """log(n choose k)."""
    if k < 0 or k > n or n < 0:
        return float("-inf")
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def _log_hypergeom_pmf(a: int, b: int, c: int, d: int) -> float:
    """log P of the 2x2 table [[a, b], [c, d]] under fixed margins."""
    n = a + b + c + d
    if n == 0:
        return 0.0
    return _lchoose(a + b, a) + _lchoose(c + d, c) - _lchoose(n, a + c)


def fisher_exact_2x2(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact test p-value for the table [[a, b], [c, d]].

    Uses the conventional "sum of tables no more likely than the observed one"
    definition, matching R's ``fisher.test`` and scipy's ``fisher_exact``.
    """
    n = a + b + c + d
    if n == 0:
        return 1.0
    row1 = a + b
    row2 = c + d
    col1 = a + c
    lo = max(0, col1 - row2)
    hi = min(row1, col1)
    logp_obs = _log_hypergeom_pmf(a, b, c, d)
    # Relative tolerance guards against float noise making an equally-likely
    # table look strictly more likely than the observed one.
    tol = 1e-7
    total = 0.0
    for x in range(lo, hi + 1):
        lp = _log_hypergeom_pmf(x, row1 - x, col1 - x, row2 - (col1 - x))
        if lp <= logp_obs + tol:
            total += exp(lp)
    return min(1.0, total)


def wilson_interval(k: int, n: int, z: float = Z_95) -> tuple[float, float]:
    """Wilson score confidence interval for a binomial proportion.

    Preferred over the normal approximation because allele fractions of
    interest here sit near the boundaries (a 2/200 observation must not yield a
    negative lower bound).
    """
    if n <= 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p). Upper-tail probability."""
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    # Sum the smaller tail for numerical stability, then complement.
    if k > n * p:
        total = 0.0
        for i in range(k, n + 1):
            total += exp(_lchoose(n, i) + i * log(p) + (n - i) * log1p_neg(p))
        return min(1.0, total)
    total = 0.0
    for i in range(0, k):
        total += exp(_lchoose(n, i) + i * log(p) + (n - i) * log1p_neg(p))
    return max(0.0, 1.0 - total)


def log1p_neg(p: float) -> float:
    """log(1 - p), accurate for small p."""
    from math import log1p

    return log1p(-p)


def phred(p: float, cap: float = 255.0) -> float:
    """Convert a probability to a Phred-scaled quality, capped."""
    if p <= 0.0:
        return cap
    if p >= 1.0:
        return 0.0
    return min(cap, -10.0 * log(p, 10))


def strand_bias(
    ref_fwd: int, ref_rev: int, alt_fwd: int, alt_rev: int
) -> tuple[float, float]:
    """Return ``(p_value, minor_strand_fraction)`` for an alt allele.

    ``p_value`` is a two-sided Fisher exact test of the strand distribution of
    the alt allele against that of the reference allele. ``minor_strand_fraction``
    is the share of alt-supporting reads on whichever strand carries fewer of
    them — 0.0 means the allele is seen on one strand only.
    """
    p = fisher_exact_2x2(ref_fwd, ref_rev, alt_fwd, alt_rev)
    alt_total = alt_fwd + alt_rev
    minor = (min(alt_fwd, alt_rev) / alt_total) if alt_total else 0.0
    return p, minor


def expected_operon_fractions(copy_number: int) -> list[float]:
    """Allele fractions expected when *k* of *copy_number* rrn operons are mutated.

    rRNA genes are present in multiple near-identical operons, so a
    heteroresistant isolate produces allele fractions clustered near k/n rather
    than at arbitrary values. Returns fractions for k = 1..copy_number.
    """
    if copy_number <= 0:
        return []
    return [k / copy_number for k in range(1, copy_number + 1)]


def nearest_operon_count(af: float, copy_number: int) -> tuple[int, float]:
    """Estimate how many rrn operons carry the allele.

    Returns ``(k, expected_af)`` for the k in 0..copy_number whose expected
    fraction k/copy_number is closest to the observed ``af``.
    """
    if copy_number <= 0:
        return (0, 0.0)
    best_k = min(range(copy_number + 1), key=lambda k: abs(k / copy_number - af))
    return (best_k, best_k / copy_number)
