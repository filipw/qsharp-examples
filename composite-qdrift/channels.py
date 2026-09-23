"""Exact classical evaluation of the simulation channels the Q# code implements.

The randomized simulators (qDRIFT, composite) are channels, not unitaries. Their
error is a property of the *average* output, which a handful of sampled circuits
estimates only noisily. On 6 qubits the channels can instead be applied exactly:
a qDRIFT segment is the m-fold power of one linear map on 64x64 density matrices,

    E1(rho) = cos^2(tau) rho - i cos(tau) sin(tau) [B/lambda_B, rho]
              + sin^2(tau) sum_j p_j P_j rho P_j,

and the deterministic factors are unitary conjugations. Everything here mirrors
the gate order of the Q# operations exactly; `verify.py` checks that agreement.

Amplitude ordering follows Q# (qubit 0 = most significant bit) throughout.
"""

from __future__ import annotations

import numpy as np

import hamiltonian as ham

DIM_CACHE: dict = {}


# --- Pauli algebra -------------------------------------------------------------------


def pauli_matrices(terms) -> np.ndarray:
    """Stack of the full-width Pauli-string matrices, shape (L, d, d)."""
    return np.array([ham.pauli_string_matrix(t) for t in terms])


def term_exponential(pauli: np.ndarray, theta: float) -> np.ndarray:
    """exp(-i theta P) for a Pauli string P (P^2 = I)."""
    d = pauli.shape[0]
    return np.cos(theta) * np.eye(d) - 1j * np.sin(theta) * pauli


def suzuki_u(order: int) -> float:
    return 1.0 / (4.0 - 4.0 ** (1.0 / (order - 1)))


def stages(order: int) -> int:
    """Sweeps over the term list in one product-formula step: 1, 2, 2*5^(k-1)."""
    if order == 1:
        return 1
    return 2 * 5 ** (order // 2 - 1)


# --- deterministic product formulas -----------------------------------------------------


class ProductFormula:
    """Mirror of Q#'s `SuzukiStep`/`TrotterEvolve` on a fixed term list, with caching."""

    def __init__(self, paulis: np.ndarray, coeffs: np.ndarray):
        self.paulis = paulis
        self.coeffs = np.asarray(coeffs, dtype=float)
        self.d = paulis.shape[1]
        self._cache: dict[tuple[int, float], np.ndarray] = {}

    def step(self, order: int, time: float) -> np.ndarray:
        """Unitary of one order-`order` step over duration `time` (Q# gate order)."""
        key = (order, float(time))
        if key in self._cache:
            return self._cache[key]
        eye = np.eye(self.d)
        if order == 1:
            u = eye
            for p, c in zip(self.paulis, self.coeffs):
                u = term_exponential(p, c * time) @ u
        elif order == 2:
            u = eye
            half = time / 2.0
            for p, c in zip(self.paulis, self.coeffs):
                u = term_exponential(p, c * half) @ u
            for p, c in zip(self.paulis[::-1], self.coeffs[::-1]):
                u = term_exponential(p, c * half) @ u
        else:
            s = suzuki_u(order)
            u = eye
            for frac in (s, s, 1.0 - 4.0 * s, s, s):
                u = self.step(order - 2, frac * time) @ u
        self._cache[key] = u
        return u

    def evolve(self, order: int, time: float, steps: int) -> np.ndarray:
        return np.linalg.matrix_power(self.step(order, time / steps), steps)


# --- qDRIFT as an exact channel -----------------------------------------------------------


class QDriftChannel:
    """The single-sample qDRIFT map for a term list, applied exactly.

    E1(tau) = cos^2 I - i cos sin (B/lambda (x) I - I (x) (B/lambda)^T) + sin^2 M2 on
    row-major vec(rho), where M2 = sum_j p_j P_j (x) P_j^T does not depend on tau and
    is built once. A dense d^2 x d^2 matvec is ~4x faster than the batched 64x64
    conjugation sum at L = 400, and repeated squaring of the same matrix handles
    the 10^5-10^6 samples pure qDRIFT needs.
    """

    def __init__(self, paulis: np.ndarray, coeffs: np.ndarray):
        coeffs = np.asarray(coeffs, dtype=float)
        self.paulis = paulis
        self.d = paulis.shape[1]
        self.lam = float(np.abs(coeffs).sum())
        self.probs = np.abs(coeffs) / self.lam
        self.bhat = np.tensordot(coeffs / self.lam, paulis, axes=(0, 0))   # B / lambda_B
        self._m2 = None
        self._comm = None

    def _parts(self):
        if self._m2 is None:
            d = self.d
            eye = np.eye(d)
            m2 = np.zeros((d * d, d * d), dtype=complex)
            for p, w in zip(self.paulis, self.probs):
                m2 += w * np.kron(p, p.T)
            self._m2 = m2
            self._comm = np.kron(self.bhat, eye) - np.kron(eye, self.bhat.T)
        return self._m2, self._comm

    def dense(self, tau: float) -> np.ndarray:
        """E1(tau) as a d^2 x d^2 matrix on row-major vec(rho)."""
        m2, comm = self._parts()
        c, s = np.cos(tau), np.sin(tau)
        return c * c * np.eye(self.d * self.d) - 1j * c * s * comm + s * s * m2

    def apply(self, rho: np.ndarray, tau: float, samples: int, method: str = "dense") -> np.ndarray:
        """E1(tau)^samples applied to rho."""
        if samples == 0:
            return rho
        if method == "direct":
            c, s = np.cos(tau), np.sin(tau)
            for _ in range(samples):
                twirl = np.tensordot(self.probs, np.matmul(np.matmul(self.paulis, rho), self.paulis), axes=(0, 0))
                rho = c * c * rho - 1j * c * s * (self.bhat @ rho - rho @ self.bhat) + s * s * twirl
            return rho
        e = self.dense(tau)
        v = rho.reshape(-1)
        for _ in range(samples):
            v = e @ v
        return v.reshape(self.d, self.d)

    def apply_power(self, rho: np.ndarray, tau: float, samples: int) -> np.ndarray:
        """Same map as `apply`, via repeated squaring; for very large sample counts."""
        e = np.linalg.matrix_power(self.dense(tau), samples)
        return (e @ rho.reshape(-1)).reshape(self.d, self.d)


# --- the composite channel -----------------------------------------------------------------


class Composite:
    """Hagan-Wiebe composite channel: mirror of Q#'s `CompositeEvolve`."""

    def __init__(self, a_paulis, a_coeffs, b_paulis, b_coeffs):
        self.a = ProductFormula(a_paulis, a_coeffs)
        self.b = QDriftChannel(b_paulis, b_coeffs)

    def _step(self, rho, outer: int, inner: int, time: float, m: int) -> np.ndarray:
        if outer == 2:
            ua = self.a.step(inner, time / 2.0)
            rho = ua @ rho @ ua.conj().T
            if m > 0:
                rho = self.b.apply(rho, self.b.lam * time / m, m)
            return ua @ rho @ ua.conj().T
        s = suzuki_u(outer)
        for frac in (s, s, 1.0 - 4.0 * s, s, s):
            rho = self._step(rho, outer - 2, inner, frac * time, m)
        return rho

    def evolve(self, rho0: np.ndarray, order: int, time: float, steps: int, m: int) -> np.ndarray:
        rho = rho0
        for _ in range(steps):
            rho = self._step(rho, order, order, time / steps, m)
        return rho


# --- error metrics -------------------------------------------------------------------------


def trace_distance(rho: np.ndarray, psi: np.ndarray) -> float:
    """1/2 || rho - |psi><psi| ||_1, the paper's error metric on this input."""
    diff = rho - np.outer(psi, psi.conj())
    return 0.5 * float(np.abs(np.linalg.eigvalsh((diff + diff.conj().T) / 2)).sum())


def infidelity(rho: np.ndarray, psi: np.ndarray) -> float:
    return float(1.0 - np.real(psi.conj() @ rho @ psi))


# --- gate accounting -------------------------------------------------------------------------


def two_qubit_gates_per_exp(weights) -> np.ndarray:
    """CNOTs Q#'s `Exp` spends on a weight-w Pauli: Rzz plus a SpreadZ ladder, 2(w-1)."""
    w = np.asarray(weights)
    return np.where(w >= 2, 2 * (w - 1), 0)


def trotter_cost(order: int, steps: int, weights) -> tuple[int, int]:
    """(rotations, two-qubit gates) for `steps` steps of an order-`order` formula."""
    sweeps = stages(order) * steps
    return sweeps * len(weights), int(sweeps * two_qubit_gates_per_exp(weights).sum())


def qdrift_cost(samples: int, weights, probs) -> tuple[int, float]:
    """(rotations, expected two-qubit gates) for `samples` qDRIFT samples."""
    return samples, float(samples * (probs * two_qubit_gates_per_exp(weights)).sum())


def composite_cost(order: int, steps: int, m: int, a_weights, b_weights, b_probs) -> tuple[int, float]:
    """Per the nested construction: 5^(k-1) base steps per outer step, each with two
    inner A-factors of `stages(order)` sweeps and one B-segment of m samples."""
    base_steps = 5 ** (order // 2 - 1) * steps
    a_rot, a_2q = trotter_cost(order, 2 * base_steps, a_weights)
    b_rot, b_2q = qdrift_cost(base_steps * m, b_weights, b_probs)
    return a_rot + b_rot, a_2q + b_2q
