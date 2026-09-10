"""Z2 linear algebra for the Bell difference sampling learner, and the paper's formulas.

A Bell difference sample of an n-qubit stabilizer state is a uniformly random
element of its unsigned stabilizer group, written as a 2n-bit vector (a | b) with
Z-part a and X-part b. Learning the group is collecting n linearly independent
samples; Lemma 5 of the paper gives the exact probability of that happening after
k draws, and that single formula is where the whole cloning lower bound comes from.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np


# --- decoding Q# results -----------------------------------------------------------------


def bell_difference(results, n: int) -> np.ndarray:
    """XOR of the two Bell labels returned by `BellDifferenceSample`: a 2n vector (a | b)."""
    bits = np.array([1 if str(r) in ("One", "1") else 0 for r in results], dtype=np.uint8)
    return bits[: 2 * n] ^ bits[2 * n :]


def to_paulis(vec: np.ndarray, n: int) -> list:
    """(a | b) -> per-qubit Pauli codes 0=I, 1=X, 2=Y, 3=Z for Q#'s Measure."""
    a, b = vec[:n], vec[n:]
    return [int({(0, 0): 0, (0, 1): 1, (1, 1): 2, (1, 0): 3}[(int(a[i]), int(b[i]))]) for i in range(n)]


def symplectic_product(u: np.ndarray, v: np.ndarray, n: int) -> int:
    return int((u[:n] @ v[n:] + u[n:] @ v[:n]) % 2)


# --- incremental basis over Z2 ------------------------------------------------------------


class Span:
    """Row-reduced basis of the samples seen so far; rank grows by one per independent sample."""

    def __init__(self, width: int):
        self.width = width
        self.rows: list[np.ndarray] = []      # reduced rows, each with a distinct pivot
        self.pivots: list[int] = []

    def reduce(self, v: np.ndarray) -> np.ndarray:
        v = v.astype(np.uint8).copy()
        for row, p in zip(self.rows, self.pivots):
            if v[p]:
                v ^= row
        return v

    def add(self, v: np.ndarray) -> bool:
        """Adds v; returns True if it was independent of the current span."""
        v = self.reduce(v)
        if not v.any():
            return False
        p = int(np.nonzero(v)[0][0])
        for k, row in enumerate(self.rows):
            if row[p]:
                self.rows[k] = row ^ v
        self.rows.append(v)
        self.pivots.append(p)
        return True

    @property
    def rank(self) -> int:
        return len(self.rows)

    def basis(self) -> np.ndarray:
        return np.array(self.rows, dtype=np.uint8)


# --- the paper's formulas ------------------------------------------------------------------


def rank_distribution(k: int, n: int) -> np.ndarray:
    """Lemma 5: P[exactly alpha of k uniform draws from Z2^n are independent], alpha = 0..n.

    Computed as a Markov chain on the rank: a fresh uniform draw is independent of
    a rank-r subspace with probability 1 - 2^(r-n).
    """
    p = np.zeros(n + 1)
    p[0] = 1.0
    for _ in range(k):
        q = np.zeros(n + 1)
        for r in range(n + 1):
            if p[r] == 0:
                continue
            stay = 2.0 ** (r - n)
            q[r] += p[r] * stay
            if r < n:
                q[r + 1] += p[r] * (1 - stay)
        p = q
    return p


def p_learned(k: int, n: int) -> float:
    """Probability that k samples span the whole group (Lemma 5 with alpha = n)."""
    return float(rank_distribution(k, n)[n])


def p_span_limit() -> float:
    """lim_n prod_{i>=1} (1 - 2^-i) = 0.2888..., the constant behind the paper's 0.14."""
    return float(np.prod([1 - 2.0 ** (-i) for i in range(1, 60)]))


def expected_samples(n: int) -> float:
    """E[number of uniform draws until they span Z2^n] = sum_r 1/(1 - 2^(r-n))."""
    return float(sum(1.0 / (1 - 2.0 ** (r - n)) for r in range(n)))


@lru_cache(maxsize=None)
def log_lagrangian_completions(n: int, r: int) -> float:
    """log of the number of Lagrangian subspaces of Z2^(2n) containing a fixed isotropic
    r-dim one, prod_{i=1}^{n-r} (2^i + 1). These are the stabilizer groups consistent
    with r independent samples. Kept in log space: the count exceeds float range early."""
    return float(sum(np.log1p(2.0**i) + i * np.log(2.0) - np.log(2.0**i) for i in range(1, n - r + 1)))


def inv_completions(n: int, r: int) -> float:
    return float(np.exp(-log_lagrangian_completions(n, r)))


def consistent_learner_success(k: int, n: int) -> float:
    """P[a uniformly random consistent hypothesis is the true group | k samples] = E[1/N(rank)].

    This is exactly the acceptance probability of the paper's distinguisher
    (Algorithm 1 / 3) on true samples.
    """
    p = rank_distribution(k, n)
    return float(sum(p[r] * inv_completions(n, r) for r in range(n + 1)))


def amplification_advantage(t: int, n: int) -> float:
    """Advantage of the paper's distinguisher against the best Bell-sampling amplifier
    (t-1 true samples plus one replayed from their span): success(t) - success(t-1).
    Theorem 6 / Corollary 7 lower-bound the optimal error at t = n by this quantity
    and show it is >= 0.14 for every n."""
    return consistent_learner_success(t, n) - consistent_learner_success(t - 1, n)


def advantage_from_ranks(ranks_t: np.ndarray, ranks_tm1: np.ndarray, n: int) -> float:
    """Same advantage, from measured rank samples instead of Lemma 5."""
    f = lambda rs: float(np.mean([inv_completions(n, int(r)) for r in rs]))
    return f(ranks_t) - f(ranks_tm1)
