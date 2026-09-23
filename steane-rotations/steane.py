"""The classical side of the demo: the [[7,1,3]] Steane code in numpy, the rotation-and-correct
protocol as exact linear algebra on 128-dimensional operators, the read-off of the logical
channel per syndrome (angle, dephasing, probability, leakage), the decoder for transversal
readout strings, tomography from expectation values, and the closed forms of arXiv:2608.20676.

Conventions follow the QDK dump: qubit 0 is the most significant bit of the basis index.
Rz(theta) = exp(-i theta Z / 2), as in the paper and in Q#.
"""

from __future__ import annotations

import numpy as np

N = 7
DIM = 2**N
H = np.array([[0, 0, 0, 1, 1, 1, 1],
              [0, 1, 1, 0, 0, 1, 1],
              [1, 0, 1, 0, 1, 0, 1]], dtype=np.uint8)
ROWS = [tuple(np.flatnonzero(r)) for r in H]

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)


# ------------------------------------------------------------------ bit strings

def bits(index: int) -> np.ndarray:
    """Basis index -> 7 bits, qubit 0 first."""
    return np.array([(index >> (N - 1 - j)) & 1 for j in range(N)], dtype=np.uint8)


def index(b) -> int:
    return int(sum(int(x) << (N - 1 - j) for j, x in enumerate(b)))


ALL_BITS = np.array([bits(i) for i in range(DIM)])              # (128, 7)
WEIGHTS = ALL_BITS.sum(axis=1).astype(int)


def flagged(s: int) -> int:
    """Qubit the decoder flips for syndrome s, or -1."""
    return s - 1


# ------------------------------------------------------------------- the code

def _span(rows) -> list[np.ndarray]:
    out = []
    for m in range(2 ** len(rows)):
        v = np.zeros(N, dtype=np.uint8)
        for k, r in enumerate(rows):
            if (m >> k) & 1:
                v ^= r
        out.append(v)
    return out


CODEWORDS0 = _span([H[i] for i in range(3)])                     # 8 strings of weight 0 or 4
CODEWORDS1 = [c ^ 1 for c in CODEWORDS0]                         # weight 7 or 3


def _uniform(words) -> np.ndarray:
    v = np.zeros(DIM, dtype=complex)
    for w in words:
        v[index(w)] = 1
    return v / np.sqrt(len(words))


ZERO_L = _uniform(CODEWORDS0)
ONE_L = _uniform(CODEWORDS1)
PLUS_L = (ZERO_L + ONE_L) / np.sqrt(2)
PLUS_I_L = (ZERO_L + 1j * ONE_L) / np.sqrt(2)
LOGICAL_INPUTS = {0: ZERO_L, 1: ONE_L, 2: PLUS_L, 3: PLUS_I_L}
INPUT_BLOCH = {0: (0, 0, 1), 1: (0, 0, -1), 2: (1, 0, 0), 3: (0, 1, 0)}
V = np.column_stack([ZERO_L, ONE_L])                             # (128, 2) isometry onto the code space


def on_qubits(single: np.ndarray, qubits) -> np.ndarray:
    """Tensor product with `single` on the listed qubits and identity elsewhere."""
    out = np.array([[1.0 + 0j]])
    for j in range(N):
        out = np.kron(out, single if j in qubits else I2)
    return out


X_STABILIZERS = [on_qubits(X, r) for r in ROWS]
Z_STABILIZERS = [on_qubits(Z, r) for r in ROWS]
X_L = on_qubits(X, range(N))
Z_L = on_qubits(Z, range(N))


def syndrome_projector(s: int) -> np.ndarray:
    """Projector onto the joint X-stabilizer eigenspace with syndrome s = s1 s2 s3 (binary)."""
    P = np.eye(DIM, dtype=complex)
    for r, Xr in enumerate(X_STABILIZERS):
        sign = -1 if (s >> (2 - r)) & 1 else 1
        P = P @ (np.eye(DIM) + sign * Xr) / 2
    return P


PROJECTORS = [syndrome_projector(s) for s in range(8)]


def correction(s: int) -> np.ndarray:
    return np.eye(DIM, dtype=complex) if s == 0 else on_qubits(Z, [flagged(s)])


CORRECTIONS = [correction(s) for s in range(8)]


# --------------------------------------------------------------- the channel

def rz_transversal(theta: float) -> np.ndarray:
    """Rz(theta)^{x7} is diagonal: exp(-i theta (7 - 2|x|) / 2) on the string x."""
    return np.diag(np.exp(-1j * theta * (N - 2 * WEIGHTS) / 2))


def dephase(rho: np.ndarray, p: float) -> np.ndarray:
    """D_p on every qubit, D_p(rho) = (1 - p) rho + p Z rho Z, Kraus form."""
    out = rho
    for j in range(N):
        Zj = on_qubits(Z, [j])
        out = (1 - p) * out + p * Zj @ out @ Zj
    return out


def physical_channel(rho: np.ndarray, theta: float, p: float) -> np.ndarray:
    """N_{p,theta}^{x7} of Eq. 2: dephase, then rotate (the two commute)."""
    R = rz_transversal(theta)
    return R @ dephase(rho, p) @ R.conj().T


def dephasing_mask(theta: float, p: float) -> np.ndarray:
    """Eq. 11 as an elementwise factor: lambda^{|x xor y|} e^{i theta (|x| - |y|)}."""
    lam = 1 - 2 * p
    ham = (ALL_BITS[:, None, :] != ALL_BITS[None, :, :]).sum(axis=2)
    return lam**ham * np.exp(1j * theta * (WEIGHTS[:, None] - WEIGHTS[None, :]))


def branch(op: np.ndarray, s: int) -> np.ndarray:
    """C_s Pi_s op Pi_s C_s: the subnormalized branch of syndrome s, back in the code space."""
    return CORRECTIONS[s] @ PROJECTORS[s] @ op @ PROJECTORS[s] @ CORRECTIONS[s]


def logical_block(op: np.ndarray) -> tuple[np.ndarray, float]:
    """(V^dagger op V, leakage): the 2x2 block of an operator in the {|0_L>, |1_L>} basis and
    the weight it has outside the code space."""
    block = V.conj().T @ op @ V
    return block, float((np.trace(op) - np.trace(block)).real)


# ------------------------------------------------------ the logical channel

class LogicalChannel:
    """A (subnormalized) single-qubit channel as a 4x4 matrix on vec(rho) in the
    {|0><0|, |0><1|, |1><0|, |1><1|} basis, plus the leakage per input."""

    def __init__(self, matrix: np.ndarray, leakage: float = 0.0):
        self.matrix = np.asarray(matrix, dtype=complex)
        self.leakage = leakage

    @staticmethod
    def from_pairs(inputs: list[np.ndarray], outputs: list[np.ndarray], leakage: float = 0.0) -> "LogicalChannel":
        """Linear inversion from input/output pairs of 2x2 operators (4 independent inputs)."""
        A = np.column_stack([r.reshape(-1) for r in inputs])
        B = np.column_stack([r.reshape(-1) for r in outputs])
        return LogicalChannel(B @ np.linalg.inv(A), leakage)

    def apply(self, rho: np.ndarray) -> np.ndarray:
        return (self.matrix @ rho.reshape(-1)).reshape(2, 2)

    @property
    def probability(self) -> float:
        return float(np.trace(self.apply(np.diag([1.0, 0.0]))).real)

    @property
    def eta(self) -> complex:
        """<0| E(|0><1|) |1>, the coherence factor of Eq. 7."""
        return complex(self.apply(np.array([[0, 1], [0, 0]], dtype=complex))[0, 1])

    @property
    def angle(self) -> float:
        """Logical rotation angle, Eq. 9, modulo 2 pi (0 by convention when p_s = 0)."""
        return 0.0 if self.probability < 1e-12 else float(wrap(-np.angle(self.eta)))

    @property
    def dephasing(self) -> float:
        """Logical dephasing rate, Eq. 10."""
        ps = self.probability
        return 0.0 if ps < 1e-12 else float((1 - abs(self.eta) / ps) / 2)

    def covariance_residual(self) -> float:
        """How far the channel is from 'populations stay, coherences get one complex factor':
        the largest entry of E(|a><b|) outside its expected place, plus |p(|0>) - p(|1>)|."""
        e00 = self.apply(np.diag([1.0, 0.0]))
        e11 = self.apply(np.diag([0.0, 1.0]))
        e01 = self.apply(np.array([[0, 1], [0, 0]], dtype=complex))
        e10 = self.apply(np.array([[0, 0], [1, 0]], dtype=complex))
        p0, p1 = e00[0, 0].real, e11[1, 1].real
        stray = max(abs(e00[0, 1]), abs(e00[1, 0]), abs(e00[1, 1]), abs(e11[0, 1]), abs(e11[1, 0]), abs(e11[0, 0]),
                    abs(e01[0, 0]), abs(e01[1, 0]), abs(e01[1, 1]), abs(e10[0, 0]), abs(e10[0, 1]), abs(e10[1, 1]),
                    abs(e01[0, 1] - np.conj(e10[1, 0])))
        return float(max(stray, abs(p0 - p1)))

    def unitarity_residual(self) -> float:
        """0 when the normalized channel is a unitary Z rotation: |eta| = p_s."""
        return float(abs(abs(self.eta) - self.probability))

    def bloch_map(self) -> tuple[np.ndarray, np.ndarray]:
        """(M, t) of the normalized channel on Bloch vectors, r -> M r + t."""
        ps = self.probability
        paulis = [X, Y, Z]
        M = np.zeros((3, 3))
        t = np.zeros(3)
        out_id = self.apply(np.eye(2, dtype=complex)) / ps
        for a, Pa in enumerate(paulis):
            t[a] = np.trace(Pa @ out_id).real / 2
            for b, Pb in enumerate(paulis):
                M[a, b] = np.trace(Pa @ self.apply(Pb)).real / 2 / ps
        return M, t


def exact_channel(theta: float, p: float, s: int, second: tuple[float, float, int] | None = None) -> LogicalChannel:
    """The logical channel of syndrome s after one round (theta, p), computed by applying the
    protocol to the four operators |a_L><b_L|. With `second` = (theta2, p2, s2), a second round
    follows the first: E_{s2}(theta2, p2) after E_s(theta, p), both branches exact."""
    basis = [np.outer(V[:, a], V[:, b].conj()) for a in range(2) for b in range(2)]
    outs, leak = [], 0.0
    for op in basis:
        out = branch(physical_channel(op, theta, p), s)
        if second is not None:
            th2, p2, s2 = second
            out = branch(physical_channel(out, th2, p2), s2)
        blk, lk = logical_block(out)
        outs.append(blk)
        leak = max(leak, abs(lk))
    ins = [np.outer(np.eye(2)[a], np.eye(2)[b]) for a in range(2) for b in range(2)]
    return LogicalChannel.from_pairs(ins, outs, leak)


def channel_from_states(outputs: dict[int, np.ndarray]) -> LogicalChannel:
    """Process tomography from the four logical inputs: `outputs[i]` is the (subnormalized)
    128x128 branch state for input i of LOGICAL_INPUTS; the channel is the linear inversion."""
    ins, outs, leak = [], [], 0.0
    for i in range(4):
        v = np.array([[1, 0], [0, 1], [1, 1], [1, 1j]][i], dtype=complex) / (1 if i < 2 else np.sqrt(2))
        ins.append(np.outer(v, v.conj()))
        blk, lk = logical_block(outputs[i])
        outs.append(blk)
        leak = max(leak, abs(lk))
    return LogicalChannel.from_pairs(ins, outs, leak)


def reduced_data_state(psi: np.ndarray, n_env: int) -> np.ndarray:
    """Trace the last n_env qubits out of a pure state of 7 + n_env qubits."""
    M = psi.reshape(DIM, 2**n_env)
    return M @ M.conj().T


def project_data(rho: np.ndarray) -> dict[int, np.ndarray]:
    """All eight corrected branches C_s Pi_s rho Pi_s C_s of a data state."""
    return {s: branch(rho, s) for s in range(8)}


# ------------------------------------------------------------ shot decoding

def decode(strings: np.ndarray) -> np.ndarray:
    """Transversal readout strings (shots x 7 bits) -> logical eigenvalues +-1: the nearest
    Hamming code word (flip the qubit flagged by H x), then the parity of its weight."""
    x = np.asarray(strings, dtype=np.uint8).copy()
    syn = (x @ H.T) % 2
    j = 4 * syn[:, 0] + 2 * syn[:, 1] + syn[:, 2]
    rows = np.flatnonzero(j)
    x[rows, j[rows] - 1] ^= 1
    return 1 - 2 * (x.sum(axis=1) % 2).astype(int)


def bloch_from_expectations(ex: dict[int, dict[str, float]]) -> tuple[np.ndarray, np.ndarray]:
    """Affine Bloch map (M, t) from the mean output Pauli of every input state,
    ex[input] = {"x": <X>, "y": <Y>, "z": <Z>}, by solving r_out = M r_in + t exactly
    (four inputs, twelve unknowns)."""
    A = np.array([list(INPUT_BLOCH[i]) + [1.0] for i in range(4)])      # (4, 4)
    B = np.array([[ex[i]["x"], ex[i]["y"], ex[i]["z"]] for i in range(4)])
    sol = np.linalg.solve(A, B)                                          # rows: M^T then t
    return sol[:3].T, sol[3]


def fit_rotation(M: np.ndarray) -> tuple[float, float]:
    """The paper's least-squares split E = E' o Rz(phi): phi maximises Tr(M R(-phi)), and the
    residual E' has average gate infidelity (3 - Tr M R(-phi)) / 6. Returns (phi, infidelity)."""
    phi = float(np.arctan2(M[1, 0] - M[0, 1], M[0, 0] + M[1, 1]))
    c, s = np.cos(phi), np.sin(phi)
    Rm = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
    return phi, float((3 - np.trace(M @ Rm)) / 6)


# ------------------------------------------------------------------- theory

def eta_paper(theta: float, p: float, trivial: bool) -> complex:
    """Eqs. 14 and 15."""
    lam = 1 - 2 * p
    a = 3 + np.exp(-4j * theta)
    b = 7 + np.exp(-8j * theta)
    if trivial:
        return np.exp(1j * theta) / 64 * (14 * lam**3 * a + lam**7 * b)
    return np.exp(1j * theta) / 64 * (2 * lam**3 * a - lam**7 * b)


def prob_paper(theta: float, p: float, trivial: bool) -> float:
    """Eqs. 16 and 17."""
    lam = 1 - 2 * p
    pt = 1 / 8 + 7 / 32 * lam**4 * (3 + np.cos(4 * theta))
    return pt if trivial else (1 - pt) / 7


def angle_paper(theta: float, p: float, trivial: bool) -> float:
    return float(-np.angle(eta_paper(theta, p, trivial)))


def dephasing_paper(theta: float, p: float, trivial: bool) -> float:
    ps = prob_paper(theta, p, trivial)
    return 0.0 if ps == 0 else float((1 - abs(eta_paper(theta, p, trivial)) / ps) / 2)


def ideal_angle(theta: float, trivial: bool) -> float:
    """Eqs. 22 and 23 (the nontrivial angle 3 theta wrapped to (-pi, pi])."""
    if trivial:
        return float(-np.angle(np.exp(1j * theta) * (7 + np.exp(-4j * theta)) ** 2))
    return float(np.angle(np.exp(3j * theta)))


def ideal_prob(theta: float, trivial: bool) -> float:
    """Eqs. 20 and 21."""
    return (25 + 7 * np.cos(4 * theta)) / 32 if trivial else (1 - np.cos(4 * theta)) / 32


def trivial_angle_closed(theta: float) -> float:
    """tan(phi_t / 2) = -(7 c^4 s^3 + s^7) / (c^7 + 7 c^3 s^4), c = cos(theta/2), s = sin(theta/2):
    the trivial-syndrome angle written out in the half-angle, an independent form of Eq. 22."""
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return float(2 * np.arctan2(-(7 * c**4 * s**3 + s**7), c**7 + 7 * c**3 * s**4))


def wrap(a):
    """Angle modulo 2 pi, as numpy's arg convention (-pi, pi]."""
    return np.angle(np.exp(1j * np.asarray(a, dtype=float)))


def angle_diff(a: float, b: float) -> float:
    return float(abs(wrap(a - b)))
