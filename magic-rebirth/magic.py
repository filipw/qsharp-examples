"""Density-matrix side of the demo.

Partial traces of the dumped Q# states, the real GHZ-X manifold read-off with the paper's
closed-form robustness of magic (Theorem 4), its dual witness, negativity, concurrence,
the single-qubit octahedron, a brute-force stabilizer-polytope LP that certifies the
closed form on small registers, and the paper's threshold formulas.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq, linprog

MANIFOLD_TOL = 1e-8


# ----------------------------------------------------------------------------- states

def reduced_states(psi: np.ndarray, n_first: int, n_rest: int) -> tuple[np.ndarray, np.ndarray]:
    """Reduced states of the first n_first and the remaining n_rest qubits of a pure state.
    The QDK's dump orders amplitudes with qubit 0 as the most significant bit."""
    M = psi.reshape(2**n_first, 2**n_rest)
    return M @ M.conj().T, M.T @ M.conj()


def cat_state(alpha: float, n: int) -> np.ndarray:
    psi = np.zeros(2**n, dtype=complex)
    psi[0] = alpha
    psi[-1] = np.sqrt(1 - alpha**2)
    return psi


def amplitude_damp(rho: np.ndarray, gamma: float, n: int) -> np.ndarray:
    """E_gamma^(x)n applied to an n-qubit density matrix via the Kraus operators of Eq. (40)."""
    E0 = np.diag([1.0, np.sqrt(1 - gamma)])
    E1 = np.array([[0.0, np.sqrt(gamma)], [0.0, 0.0]])
    out = rho
    for i in range(n):
        d_l, d_r = 2**i, 2**(n - i - 1)
        acc = np.zeros_like(out)
        for K in (E0, E1):
            Kf = np.kron(np.kron(np.eye(d_l), K), np.eye(d_r))
            acc += Kf @ out @ Kf.conj().T
        out = acc
    return out


def dephase(rho: np.ndarray, p: float, n: int) -> np.ndarray:
    """D_p^(x)n, D_p(rho) = (1 - p) rho + p Z rho Z on every qubit."""
    out = rho
    for i in range(n):
        Z = np.kron(np.kron(np.eye(2**i), np.diag([1.0, -1.0])), np.eye(2**(n - i - 1)))
        out = (1 - p) * out + p * Z @ out @ Z
    return out


# ------------------------------------------------------------------- GHZ-X manifold

def ghzx(rho: np.ndarray) -> tuple[float, float, float, float]:
    """(p0, p1, c, residual): the endpoint populations, the (0^n, 1^n) coherence, and the
    largest matrix element that a real GHZ-X state must not have. residual ~ 0 certifies
    that the state is in the manifold, where Theorem 4 applies."""
    off = rho - np.diag(np.diag(rho))
    off[0, -1] = 0
    off[-1, 0] = 0
    resid = max(float(np.abs(off).max()), abs(float(rho[0, -1].imag)), float(np.abs(np.diag(rho).imag).max()))
    return float(rho[0, 0].real), float(rho[-1, -1].real), float(rho[0, -1].real), resid


def robustness_ghzx(rho: np.ndarray) -> float:
    """Signed-decomposition robustness of magic on the real GHZ-X manifold, Eq. (53):
    R = 1 + 2 max{0, |c| - p0, |c| - p1}; R = 1 exactly on the stabilizer polytope."""
    p0, p1, c, resid = ghzx(rho)
    if resid > MANIFOLD_TOL:
        raise ValueError(f"state is not in the real GHZ-X manifold (residual {resid:.1e})")
    return 1 + 2 * max(0.0, abs(c) - p0, abs(c) - p1)


def robustness_certificate(rho: np.ndarray) -> tuple[float, float, bool]:
    """Certifies the robustness without the closed form, for any n.

    Lower bound: the dual witness W_{s,j} of Eq. (60) evaluated on rho (its feasibility,
    |Tr(W sigma)| <= 1 for every stabilizer state, is the three-line argument of Eqs. 61-63;
    verify.py brute-forces it for n <= 4). Upper bound: the explicit signed decomposition of
    Eqs. (59) and (70) into computational-basis states and the two GHZ states, checked here
    by reconstructing rho from it with non-negative convex weights. Returns
    (lower, upper, decomposition_valid); lower == upper pins R down."""
    d = rho.shape[0]
    p0, p1, c, _ = ghzx(rho)
    s = 1 if c >= 0 else -1
    m = min(p0, p1)
    t = max(0.0, abs(c) - m)
    lower = 1 + 2 * abs(c) - 2 * m

    def ghz(sign):
        v = np.zeros(d)
        v[0], v[-1] = 1 / np.sqrt(2), sign / np.sqrt(2)
        return np.outer(v, v)

    omega = (rho + t * ghz(-s)) / (1 + t)
    cw = float(omega[0, -1].real)
    weights = np.diag(omega).real.copy()
    weights[0] -= abs(cw)
    weights[-1] -= abs(cw)
    recon = 2 * abs(cw) * ghz(1 if cw >= 0 else -1) + np.diag(weights)
    valid = bool(weights.min() >= -1e-12 and np.abs(recon - omega).max() < 1e-12
                 and np.abs((1 + t) * omega - t * ghz(-s) - rho).max() < 1e-12)
    upper = (1 + t) * (2 * abs(cw) + weights.sum()) + t
    return lower, float(upper), valid


def witness_value(p0: float, p1: float, c: float) -> float:
    """Best of the two dual witnesses W_{s,j} of Eq. (60): Tr(W rho) = 1 + 2|c| - 2 min(p0, p1).
    Every stabilizer state has |Tr(W sigma)| <= 1, so a value above 1 certifies magic."""
    return 1 + 2 * abs(c) - 2 * min(p0, p1)


def decoded_qubit(rho: np.ndarray) -> tuple[np.ndarray, float]:
    """Parity-syndrome extraction, exactly: project onto span{|0^n>, |1^n>}, decode to one qubit.
    Returns (rho_tilde, success probability), Eq. (19)."""
    p0, p1, c, _ = ghzx(rho)
    ps = p0 + p1
    return np.array([[p0, c], [c, p1]]) / ps, ps


# ---------------------------------------------------------------- entanglement

def negativity(rho: np.ndarray, n: int, m: int) -> float:
    """Bipartite negativity (||rho^{T_A}||_1 - 1)/2 across the first m qubits | the rest."""
    d1, d2 = 2**m, 2**(n - m)
    pt = rho.reshape(d1, d2, d1, d2).transpose(2, 1, 0, 3).reshape(d1 * d2, d1 * d2)
    ev = np.linalg.eigvalsh((pt + pt.conj().T) / 2)
    return float((np.abs(ev).sum() - 1) / 2)


def concurrence(rho: np.ndarray) -> float:
    """Wootters concurrence of a two-qubit state."""
    sy = np.array([[0, -1j], [1j, 0]])
    yy = np.kron(sy, sy)
    ev = np.linalg.eigvals(rho @ yy @ rho.conj() @ yy)
    lam = np.sort(np.sqrt(np.abs(ev)))[::-1]
    return float(max(0.0, lam[0] - lam[1] - lam[2] - lam[3]))


# ------------------------------------------------------------------ one qubit

def bloch(rho: np.ndarray) -> tuple[float, float, float]:
    return float(2 * rho[0, 1].real), float(-2 * rho[0, 1].imag), float((rho[0, 0] - rho[1, 1]).real)


def robustness_qubit(rho: np.ndarray) -> float:
    """R = 1 + max{0, |x| + |y| + |z| - 1}: the stabilizer octahedron is |x| + |y| + |z| <= 1."""
    x, y, z = bloch(rho)
    return 1 + max(0.0, abs(x) + abs(y) + abs(z) - 1)


# ------------------------------------------------ brute-force stabilizer polytope

def stabilizer_states(n: int) -> np.ndarray:
    """Every pure n-qubit stabilizer state, as the orbit of |0^n> under H, S and CNOT.
    Rows are state vectors with the global phase fixed. Counts: 6, 60, 1080, 36720."""
    d = 2**n
    H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    S = np.diag([1, 1j])

    def one(g, i):
        return np.kron(np.kron(np.eye(2**i), g), np.eye(2**(n - i - 1)))

    def cnot(c, t):
        U = np.zeros((d, d))
        for x in range(d):
            bits = [(x >> (n - 1 - k)) & 1 for k in range(n)]
            if bits[c]:
                bits[t] ^= 1
            U[sum(b << (n - 1 - k) for k, b in enumerate(bits)), x] = 1
        return U

    gates = [one(H, i) for i in range(n)] + [one(S, i) for i in range(n)]
    gates += [cnot(c, t) for c in range(n) for t in range(n) if c != t]

    def key(v):
        idx = np.flatnonzero(np.abs(v) > 1e-9)[0]
        w = v / (v[idx] / abs(v[idx]))
        return tuple(np.round(w.real, 6)) + tuple(np.round(w.imag, 6))

    start = np.zeros(d, dtype=complex)
    start[0] = 1
    seen = {key(start): start}
    frontier = [start]
    while frontier:
        F = np.array(frontier)
        new = []
        for G in gates:
            for v in F @ G.T:
                k = key(v)
                if k not in seen:
                    seen[k] = v
                    new.append(v)
        frontier = new
    return np.array(list(seen.values()))


def constant_weight_support(state: np.ndarray, n: int) -> bool:
    """Proposition 12's magic-insulator test: every computational string in the support has
    the same Hamming weight."""
    supp = np.flatnonzero(np.abs(state) > 1e-9)
    weights = {bin(int(x)).count("1") for x in supp}
    return len(weights) == 1


def _features(M: np.ndarray, iu) -> np.ndarray:
    return np.concatenate([np.diag(M).real, M[iu].real, M[iu].imag])


def robustness_lp(rho: np.ndarray, states: np.ndarray) -> float:
    """Signed-decomposition robustness by linear programming over all pure stabilizer states:
    min sum |q_j| subject to sum_j q_j |phi_j><phi_j| = rho (Eq. 47)."""
    N, d = states.shape
    iu = np.triu_indices(d, 1)
    A = np.empty((d * d, N))
    for j, s in enumerate(states):
        A[:, j] = _features(np.outer(s, s.conj()), iu)
    b = _features(rho, iu)
    res = linprog(np.ones(2 * N), A_eq=np.hstack([A, -A]), b_eq=b, bounds=(0, None), method="highs")
    if res.status != 0:
        raise RuntimeError(f"LP failed: {res.message}")
    return float(res.fun)


# --------------------------------------------------------------------- theory

def cat_ghzx(alpha: float, gamma: float, n: int) -> tuple[float, float, float]:
    """(P_0, P_n, c) of the amplitude-damped cat, Eq. (1)."""
    beta = np.sqrt(1 - alpha**2)
    return alpha**2 + beta**2 * gamma**n, beta**2 * (1 - gamma) ** n, alpha * beta * (1 - gamma) ** (n / 2)


def thresholds(alpha: float, n: int, eta: float = 1.0) -> tuple[float, float, float]:
    """(gamma_-, gamma_+, gamma_e) for the cat alpha|0^n> + beta|1^n> under amplitude damping,
    optionally composed with a phase flip of strength p, eta = (1 - 2p)^2 (Appendix D):
    gamma_+ = 1 - eta r^{2/n}, gamma_e = eta r^{2/n}, gamma_- the root of P_0 = c.
    gamma_- is NaN when there is no stabilizer window (r >= eta^{n/2})."""
    beta = np.sqrt(1 - alpha**2)
    r = alpha / beta
    if r >= 1:                                   # no stabilizer window at all (Theorem 1)
        return float("nan"), float("nan"), float("nan")
    g_plus = 1 - eta * r ** (2 / n)
    g_e = eta * r ** (2 / n)

    def f(g):
        p0, _, c = cat_ghzx(alpha, g, n)
        return p0 - eta ** (n / 2) * c

    g_minus = brentq(f, 0.0, g_plus) if f(0.0) < 0 < f(g_plus) else float("nan")
    return g_minus, g_plus, g_e


def regime_boundaries(n: int) -> tuple[float, float]:
    """alpha_1^(n) < alpha_2^(n) of Corollary 1 (Eqs. 99-103): below alpha_1 entanglement dies
    before magic does (Regime I), between them magic dies first but is reborn only after
    entanglement is gone (II), above alpha_2 rebirth happens while still entangled (III)."""
    r1 = (1 + 2 ** (2 / n)) ** (-n / 2)
    r2 = 2 ** (-n / 2)
    return r1 / np.sqrt(1 + r1**2), r2 / np.sqrt(1 + r2**2)


def yield_curve(alpha: float, s: float) -> float:
    """Expected robustness yield per register on the reborn branch, Eq. (22), s = (1-gamma)^{n/2}.
    Maximal, alpha^2 / 2, at s = r / 2."""
    beta = np.sqrt(1 - alpha**2)
    return 2 * s * (alpha * beta - beta**2 * s)


def dephasing_for_exponent(gamma: float, a: float) -> float:
    """Phase-flip strength that turns amplitude damping (coherence factor q^{1/2}, q = 1 - gamma)
    into the concurrent-dephasing semigroup of Theorem 2 with coherence factor q^a,
    a = 1/2 + Gamma_phi / kappa: (1 - 2p) q^{1/2} = q^a."""
    q = 1 - gamma
    return (1 - q ** (a - 0.5)) / 2
