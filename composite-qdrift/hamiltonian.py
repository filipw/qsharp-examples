"""Skewed-norm local Hamiltonian, plus the exact classical reference it is scored against.

Amplitude ordering matches Q#'s state dump: qubit 0 is the most significant bit,
so every many-qubit operator is built by kron-ing qubit 0 first.
"""

from __future__ import annotations

import numpy as np

I2 = np.eye(2, dtype=complex)
PAULI = [
    I2,
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex),
]
PAULI_NAME = "IXYZ"


def pauli_string_matrix(term: tuple[int, ...]) -> np.ndarray:
    out = np.array([[1.0 + 0j]])
    for p in term:
        out = np.kron(out, PAULI[p])
    return out


def dense_hamiltonian(terms, coeffs, n: int) -> np.ndarray:
    h = np.zeros((2**n, 2**n), dtype=complex)
    for term, c in zip(terms, coeffs):
        h += c * pauli_string_matrix(term)
    return h


def build_hamiltonian(n: int, n_small: int, lambda_small: float, seed: int):
    """Transverse-field Ising backbone (norms ~1) plus a dense bath of weak 2-local couplings.

    The point of the skew is that the bath contributes many *terms* but little
    1-norm: a product formula pays for the count, qDRIFT pays for the norm.

    The bath coefficients are rescaled so that sum |h_j| over the small terms is
    exactly `lambda_small`. That matters because the two costs move in opposite
    directions with the number of bath terms -- Trotterizing the tail costs O(n_small)
    gates per step, while qDRIFT-ing it costs O(lambda_small^2 / eps). Holding the
    1-norm fixed while the count grows is precisely the skew a composite simulator
    is supposed to exploit.
    """
    rng = np.random.default_rng(seed)
    terms: list[tuple[int, ...]] = []
    coeffs: list[float] = []

    def add(term, c):
        terms.append(tuple(term))
        coeffs.append(float(c))

    # Large terms: nearest-neighbour ZZ couplings and transverse X fields.
    for i in range(n - 1):
        term = [0] * n
        term[i] = term[i + 1] = 3
        add(term, 1.0 + 0.3 * rng.standard_normal())
    for i in range(n):
        term = [0] * n
        term[i] = 1
        add(term, 0.8 + 0.3 * rng.standard_normal())

    n_big = len(terms)

    # Small terms: random weight-2 Paulis on random pairs, no duplicates.
    # There are only C(n, 2) * 9 distinct weight-2 strings, so asking for more
    # would spin forever in the rejection loop below.
    available = 9 * n * (n - 1) // 2
    if n_small > available - n_big:
        raise ValueError(
            f"n_small={n_small} exceeds the {available - n_big} distinct weight-2 "
            f"Pauli strings left on {n} qubits"
        )

    seen = {t for t in terms}
    while len(terms) - n_big < n_small:
        i, j = rng.choice(n, size=2, replace=False)
        pi, pj = rng.integers(1, 4, size=2)
        term = [0] * n
        term[i], term[j] = int(pi), int(pj)
        key = tuple(term)
        if key in seen:
            continue
        seen.add(key)
        add(term, rng.standard_normal())

    big_terms, big_coeffs = terms[:n_big], np.array(coeffs[:n_big])
    small_terms, small_coeffs = terms[n_big:], np.array(coeffs[n_big:])
    small_coeffs *= lambda_small / np.abs(small_coeffs).sum()
    return big_terms, big_coeffs, small_terms, small_coeffs


def all_pauli_strings(n: int, max_weight: int) -> list[tuple[int, ...]]:
    """Every non-identity Pauli string on n qubits with weight <= max_weight."""
    out = []
    for code in range(1, 4**n):
        term = tuple((code >> (2 * (n - 1 - k))) & 3 for k in range(n))
        if 0 < sum(1 for p in term if p) <= max_weight:
            out.append(term)
    return out


def build_powerlaw(n: int, n_terms: int, alpha: float, seed: int, max_weight: int = 3):
    """H = sum_j a_j P_j with |a_j| proportional to j^(-alpha), sum_j |a_j| = 1.

    This is the coefficient profile the paper singles out ("power-law
    interactions"): the tail mass lambda_K = sum_{j>K} |a_j| decays only
    polynomially, so the optimal Trotter/qDRIFT split K*(eps) keeps moving as eps
    shrinks and minimising the paper's objective gives a gate count scaling as
    eps^(-1/(2 alpha - 1)) (our calculation; the paper does not write the exponent out).

    Terms are random Pauli strings of weight <= max_weight so that every term
    costs O(1) two-qubit gates, as the paper's gate model assumes. The list comes
    back sorted by |a_j| descending, which is the order the split K refers to.
    """
    rng = np.random.default_rng(seed)
    pool = all_pauli_strings(n, max_weight)
    if n_terms > len(pool):
        raise ValueError(f"only {len(pool)} Pauli strings of weight <= {max_weight} on {n} qubits")
    picks = rng.choice(len(pool), size=n_terms, replace=False)
    terms = [pool[i] for i in picks]
    mags = np.arange(1, n_terms + 1, dtype=float) ** (-alpha)
    mags /= mags.sum()
    signs = rng.choice([-1.0, 1.0], size=n_terms)
    return terms, mags * signs


def term_label(term: tuple[int, ...]) -> str:
    return "".join(PAULI_NAME[p] for p in term)


def pauli_weight(term: tuple[int, ...]) -> int:
    return sum(1 for p in term if p != 0)


# --- classical mirror of Q#'s PrepareInput -------------------------------------


def _single_qubit_op(n: int, i: int, m: np.ndarray) -> np.ndarray:
    out = np.array([[1.0 + 0j]])
    for k in range(n):
        out = np.kron(out, m if k == i else I2)
    return out


def _cnot(n: int, control: int, target: int) -> np.ndarray:
    dim = 2**n
    out = np.zeros((dim, dim), dtype=complex)
    for idx in range(dim):
        bits = [(idx >> (n - 1 - k)) & 1 for k in range(n)]
        if bits[control] == 1:
            bits[target] ^= 1
        j = sum(b << (n - 1 - k) for k, b in enumerate(bits))
        out[j, idx] = 1.0
    return out


def input_angles(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.uniform(0.2, np.pi - 0.2, size=2 * n)


def input_state(n: int, angles: np.ndarray) -> np.ndarray:
    psi = np.zeros(2**n, dtype=complex)
    psi[0] = 1.0
    for i in range(n):
        th = angles[i]
        ry = np.array([[np.cos(th / 2), -np.sin(th / 2)], [np.sin(th / 2), np.cos(th / 2)]], dtype=complex)
        psi = _single_qubit_op(n, i, ry) @ psi
    for i in range(n - 1):
        psi = _cnot(n, i, i + 1) @ psi
    for i in range(n):
        th = angles[n + i]
        rz = np.diag([np.exp(-1j * th / 2), np.exp(1j * th / 2)])
        psi = _single_qubit_op(n, i, rz) @ psi
    return psi


def exact_evolved_state(h: np.ndarray, time: float, psi0: np.ndarray) -> np.ndarray:
    w, v = np.linalg.eigh(h)
    return v @ (np.exp(-1j * w * time) * (v.conj().T @ psi0))
