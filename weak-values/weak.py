"""Numpy mirror of the Q# experiment and the closed forms it is compared with.

Units: positions in grid steps, momenta in radians per grid step. The pointer is a Gaussian of
position standard deviation sigma (momentum standard deviation 1 / (2 sigma)) and the kick is
+-dp, so the one number that says how weak the measurement is is w = dp * sigma, the kick
measured against the pointer's coherence momentum 1 / sigma (Sec. II of arXiv:2609.37822).
"""

from __future__ import annotations

import numpy as np

# The spin: |0> = |+z>, |1> = |-z>.
PLUS_X = np.array([1, 1]) / np.sqrt(2)
MINUS_X = np.array([1, -1]) / np.sqrt(2)
SIGMA_Z = np.diag([1.0, -1.0])


def ry(theta: float) -> np.ndarray:
    """exp(-i theta Y / 2), the Q# Ry."""
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([[c, -s], [s, c]])


def final_state(delta: float) -> np.ndarray:
    """The post-selected spin state |f(delta)> = Ry(delta)|-x>."""
    return ry(delta) @ MINUS_X


def orthogonal_state(delta: float) -> np.ndarray:
    """|n(delta)> = Ry(delta)|+x>, the outcome the post-selection throws away."""
    return ry(delta) @ PLUS_X


def delta_for(weak: float) -> float:
    """The post-selection angle at which the weak value of sigma_z is `weak`: cot(delta/2) = weak."""
    return 2 * np.arctan(1 / weak)


# ------------------------------------------------------------------ closed forms

def weak_value(delta: float) -> float:
    """<f|sigma_z|i> / <f|i> for |i> = |+x>, |f> = Ry(delta)|-x>: cos(delta/2) / sin(delta/2)."""
    return 1 / np.tan(delta / 2)


def kappa(w: float) -> float:
    """Overlap of the two kicked pointer states, <psi e^{+i dp z} | psi e^{-i dp z}> = exp(-2 w^2)
    for the Gaussian pointer (the characteristic function of N(0, sigma^2) at 2 dp). For another
    real pointer the formulas below hold with kappa = int |psi|^2 cos(2 dp z) dz."""
    return float(np.exp(-2 * w * w))


def amplification(delta, w: float):
    """Mean momentum of the post-selected pointer divided by the kick dp, for any strength:
        A(delta) = sin(delta) / (1 - kappa cos(delta)).
    It tends to the weak value cot(delta/2) as w -> 0, equals +-1 at delta = +-pi/2 for every w
    (the post-selection then picks |+-z> and the pointer just reads the eigenvalue), and passes
    through 0 at delta = 0. Derivation in the README."""
    return np.sin(delta) / (1 - kappa(w) * np.cos(delta))


def amplification_orthogonal(delta, w: float):
    """The same for the discarded outcome |n(delta)>: -sin(delta) / (1 + kappa cos(delta))."""
    return -np.sin(delta) / (1 + kappa(w) * np.cos(delta))


def postselection_probability(delta, w: float):
    """P(f) = (1 - kappa cos(delta)) / 2; sin^2(delta/2) = |<f|i>|^2 in the weak limit."""
    return (1 - kappa(w) * np.cos(delta)) / 2


def peak(w: float) -> tuple[float, float]:
    """(delta*, A_max): dA/d delta = 0 at cos(delta*) = kappa, where A = 1 / sqrt(1 - kappa^2).
    For small w that is delta* ~ 2w and A_max ~ 1 / (2w), so the largest shift A_max * dp tends
    to 1 / (2 sigma) however weak the kick: set by the pointer's position spread, and equal to
    its momentum spread only because this Gaussian has minimum uncertainty."""
    k = kappa(w)
    return float(np.arccos(k)), float(1 / np.sqrt(1 - k * k))


# ------------------------------------------------------------------ the grid model

def grid(n: int) -> np.ndarray:
    """Positions z_j = j - (N - 1)/2 of the n-qubit pointer register."""
    size = 2**n
    return np.arange(size) - (size - 1) / 2


def momenta(n: int) -> np.ndarray:
    """Momentum of read-out index k: 2 pi k / N with k >= N/2 taken as k - N."""
    size = 2**n
    k = np.arange(size)
    return 2 * np.pi * np.where(k >= size // 2, k - size, k) / size


def pointer(n: int, sigma: float) -> np.ndarray:
    """The Gaussian pointer as PreparePointer builds it."""
    z = grid(n)
    a = np.exp(-z * z / (4 * sigma * sigma))
    return a / np.linalg.norm(a)


def kicked(n: int, sigma: float, dp: float) -> np.ndarray:
    """(2, N) array: the pointer amplitude attached to |+z> (kicked by +dp) and |-z> (-dp), with
    the spin in |+x> before the kick."""
    z = grid(n)
    psi = pointer(n, sigma)
    return np.array([psi * np.exp(1j * dp * z), psi * np.exp(-1j * dp * z)]) / np.sqrt(2)


def branch(n: int, sigma: float, dp: float, spin: np.ndarray) -> np.ndarray:
    """The unnormalised pointer state left by projecting the spin onto `spin`, in momentum
    representation (index k as the Q# read-out leaves it)."""
    pos = spin.conj() @ kicked(n, sigma, dp)
    return np.fft.fft(pos) / np.sqrt(2**n)


def experiment(n: int, sigma: float, dp: float, delta: float) -> dict:
    """Probability and mean momentum (in units of dp) of both outcomes, and the momentum
    distribution of each, from the grid model."""
    p = momenta(n)
    out = {}
    for name, spin in (("f", final_state(delta)), ("n", orthogonal_state(delta))):
        amp = branch(n, sigma, dp, spin)
        prob = np.abs(amp) ** 2
        total = float(prob.sum())
        out[name] = dict(prob=total, shift=float((p * prob).sum() / total / dp), dist=prob / total)
    return out
