"""Thin bridge between the Q# project and numpy: dumped states become reduced density
matrices, shots become bit arrays."""

from __future__ import annotations

import os

import numpy as np
from qdk import code, qsharp

import magic

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

# `qsharp.run(..., seed=s)` seeds shot i with s + i, so two runs whose seed ranges overlap
# share random numbers (seeds s and s + 1 give the same outcomes shifted by one shot). Runs
# that must be independent are therefore spaced by more than any shot count used here.
SEED_STRIDE = 1_000_000


def init() -> None:
    qsharp.init(project_root=PROJECT_ROOT)


# ------------------------------------------------------------------ exact states

def dumped_state(callable_, *args) -> np.ndarray:
    result = qsharp.run(callable_, 1, *args, save_events=True)
    return np.array(result[0]["dumps"][0].as_dense_state(), dtype=complex)


def damped_cat(n: int, alpha: float, gamma: float) -> tuple[np.ndarray, np.ndarray]:
    """(rho_S, rho_E): the damped cat and the state its environment is left in."""
    psi = dumped_state(code.Main.SimDampedCat, n, alpha, gamma)
    return magic.reduced_states(psi, n, n)


def damped_dephased_cat(n: int, alpha: float, gamma: float, p: float) -> np.ndarray:
    psi = dumped_state(code.Main.SimDampedDephasedCat, n, alpha, gamma, p)
    return magic.reduced_states(psi, n, 2 * n)[0]


def damped_circuit(n: int, circuit, gamma: float) -> np.ndarray:
    gates, a, b = circuit
    psi = dumped_state(code.Main.SimDampedCircuit, n, gates, a, b, gamma)
    return magic.reduced_states(psi, n, n)[0]


def cat(n: int, alpha: float) -> np.ndarray:
    return dumped_state(code.Main.SimCat, n, alpha)


# Gate codes (0 = H, 1 = S, 2 = CNOT, 3 = X) for the stabilizer inputs of Sec. VII.
CIRCUITS = {
    "Phi+": (2, ([0, 2], [0, 0], [0, 1])),                    # (|00> + |11>)/sqrt2
    "Psi+": (2, ([0, 2, 3], [0, 0, 1], [0, 1, 0])),           # (|01> + |10>)/sqrt2
    "GHZ3": (3, ([0, 2, 2], [0, 0, 0], [0, 1, 2])),           # (|000> + |111>)/sqrt2
    "Phi+ x 0": (3, ([0, 2], [0, 0], [0, 1])),                # (|00> + |11>)|0>/sqrt2
    "00+": (3, ([0], [2], [0])),                              # |00>|+>
}


# ------------------------------------------------------------------------ shots

def _bits(results) -> np.ndarray:
    return np.array([[1 if str(r) == "One" else 0 for r in shot] for shot in results], dtype=np.uint8)


def _seed(seed: int | None, k: int) -> int | None:
    return None if seed is None else seed + k * SEED_STRIDE


def measure_cat(n: int, alpha: float, gamma: float, basis: int, shots: int, seed: int | None = None) -> np.ndarray:
    """`shots` rows of n bits: Z-basis outcomes (basis 0) or X-basis outcomes (basis 1)."""
    return _bits(qsharp.run(code.Main.MeasureCat, shots, n, alpha, gamma, basis, seed=seed))


def witness_from_shots(n: int, alpha: float, gamma: float, shots: int, seed: int | None = None) -> dict:
    """The dual witness from two measurement settings: P_0 and P_n from the Z strings, 2c from
    the parity of the X strings. Standard errors are binomial; the witness error adds the
    error of whichever endpoint population is the smaller one. The damped cat's coherence is
    non-negative, so c is not folded through an absolute value here."""
    z = measure_cat(n, alpha, gamma, 0, shots, _seed(seed, 0))
    x = measure_cat(n, alpha, gamma, 1, shots, _seed(seed, 1))
    p0 = float(np.mean(~z.any(axis=1)))
    pn = float(np.mean(z.all(axis=1)))
    xpar = 1 - 2 * (x.sum(axis=1).astype(int) % 2)
    c = float(xpar.mean()) / 2
    se_p0 = float(np.sqrt(p0 * (1 - p0) / shots))
    se_pn = float(np.sqrt(pn * (1 - pn) / shots))
    se_c = float(xpar.std(ddof=1) / np.sqrt(shots)) / 2
    w = 1 + 2 * c - 2 * min(p0, pn)
    se_w = float(np.sqrt(4 * se_c**2 + 4 * (se_p0**2 if p0 < pn else se_pn**2)))
    return dict(p0=p0, pn=pn, c=c, witness=w, se=se_w, se_p0=se_p0, se_pn=se_pn, se_c=se_c)


def extract(n: int, alpha: float, gamma: float, basis: int, shots: int, seed: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """(syndromes, readout): n - 1 parity bits and the decoded qubit's outcome per shot."""
    bits = _bits(qsharp.run(code.Main.ExtractMagic, shots, n, alpha, gamma, basis, seed=seed))
    return bits[:, :-1], bits[:, -1]


def decoded_bloch_from_shots(n: int, alpha: float, gamma: float, shots: int, seed: int | None = None) -> dict:
    """Post-selected Bloch vector of the extracted qubit and the success probability, with
    standard errors, from three independent runs (Z, X, Y readout)."""
    out = {}
    accepted = 0
    total = 0
    for k, (basis, name) in enumerate(((0, "z"), (1, "x"), (2, "y"))):
        syn, r = extract(n, alpha, gamma, basis, shots, _seed(seed, k))
        keep = ~syn.any(axis=1)
        accepted += int(keep.sum())
        total += shots
        vals = 1 - 2 * r[keep].astype(float)
        out[name] = float(vals.mean())
        out["se_" + name] = float(vals.std(ddof=1) / np.sqrt(len(vals)))
    out["p_success"] = accepted / total
    out["se_p_success"] = float(np.sqrt(out["p_success"] * (1 - out["p_success"]) / total))
    return out
