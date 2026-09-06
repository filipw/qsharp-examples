"""Classical stabilizer tableau (Aaronson-Gottesman) used as ground truth.

Each row is a signed Pauli (-1)^r * prod_j P(x_j, z_j) with P(1,1) = Y. Only the n
stabilizer rows are tracked, since that is all the demo needs to check what the
learner recovers.
"""

from __future__ import annotations

import numpy as np


class Tableau:
    def __init__(self, n: int):
        self.n = n
        self.x = np.zeros((n, n), dtype=np.uint8)
        self.z = np.eye(n, dtype=np.uint8)
        self.r = np.zeros(n, dtype=np.uint8)

    def h(self, q: int) -> None:
        self.r ^= self.x[:, q] & self.z[:, q]
        self.x[:, q], self.z[:, q] = self.z[:, q].copy(), self.x[:, q].copy()

    def s(self, q: int) -> None:
        self.r ^= self.x[:, q] & self.z[:, q]
        self.z[:, q] ^= self.x[:, q]

    def cnot(self, c: int, t: int) -> None:
        self.r ^= self.x[:, c] & self.z[:, t] & (self.x[:, t] ^ self.z[:, c] ^ 1)
        self.x[:, t] ^= self.x[:, c]
        self.z[:, c] ^= self.z[:, t]

    def apply(self, gates, a, b) -> "Tableau":
        for g, i, j in zip(gates, a, b):
            if g == 0:
                self.h(i)
            elif g == 1:
                self.s(i)
            else:
                self.cnot(i, j)
        return self

    def symplectic(self) -> np.ndarray:
        """Generators as rows [z | x] over Z2, matching the Bell-label convention (a, b)."""
        return np.concatenate([self.z, self.x], axis=1)

    def sign_of(self, z: np.ndarray, x: np.ndarray) -> int | None:
        """Sign (0 or 1) of the group element with Z-part z and X-part x, or None if the
        Pauli is not in the stabilizer group. Multiplies rows with full phase tracking."""
        coeffs = solve_z2(self.symplectic(), np.concatenate([z, x]))
        if coeffs is None:
            return None
        px = np.zeros(self.n, dtype=np.uint8)
        pz = np.zeros(self.n, dtype=np.uint8)
        phase = 0                                   # exponent of i, mod 4
        for k in np.nonzero(coeffs)[0]:
            phase = (phase + 2 * int(self.r[k]) + _g_sum(px, pz, self.x[k], self.z[k])) % 4
            px ^= self.x[k]
            pz ^= self.z[k]
        assert phase % 2 == 0, "stabilizer group element with imaginary phase"
        return (phase // 2) % 2


def _g_sum(x1, z1, x2, z2) -> int:
    """Sum over qubits of the AG phase function g for the product P(x1,z1) * P(x2,z2)."""
    total = 0
    for a, b, c, d in zip(x1, z1, x2, z2):
        if a == 0 and b == 0:
            continue
        if a == 1 and b == 1:
            total += int(d) - int(c)
        elif a == 1:
            total += int(d) * (2 * int(c) - 1)
        else:
            total += int(c) * (1 - 2 * int(d))
    return total % 4


def solve_z2(rows: np.ndarray, target: np.ndarray):
    """Coefficients c with c @ rows = target over Z2, or None if target is not in the span."""
    m, w = rows.shape
    aug = np.concatenate([rows.astype(np.uint8), np.eye(m, dtype=np.uint8)], axis=1)
    piv = []
    row = 0
    for col in range(w):
        pr = next((k for k in range(row, m) if aug[k, col]), None)
        if pr is None:
            continue
        aug[[row, pr]] = aug[[pr, row]]
        for k in range(m):
            if k != row and aug[k, col]:
                aug[k] ^= aug[row]
        piv.append(col)
        row += 1
    t = target.astype(np.uint8).copy()
    coeffs = np.zeros(m, dtype=np.uint8)
    for k, col in enumerate(piv):
        if t[col]:
            t ^= aug[k, :w]
            coeffs ^= aug[k, w:]
    return coeffs if not t.any() else None
