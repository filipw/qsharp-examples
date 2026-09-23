"""Thin bridge between the Q# project and numpy: dumped states become data-register density
matrices, shots become syndromes and decoded logical eigenvalues."""

from __future__ import annotations

import os

import numpy as np
from qdk import code, qsharp

import steane

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

# `qsharp.run(..., seed=s)` seeds shot i with s + i, so two runs whose seed ranges overlap
# share random numbers (seeds s and s + 1 give the same outcomes shifted by one shot). Runs
# that must be independent are therefore spaced by more than any shot count used here.
SEED_STRIDE = 1_000_000


def init() -> None:
    qsharp.init(project_root=PROJECT_ROOT)


def seed_for(seed: int | None, k: int) -> int | None:
    return None if seed is None else seed + k * SEED_STRIDE


# ------------------------------------------------------------------ exact states

def _dense(dump) -> np.ndarray:
    return np.array(dump.as_dense_state(), dtype=complex)


def encoded(input_state: int) -> np.ndarray:
    """The 7-qubit state vector of logical input 0..3."""
    r = qsharp.run(code.Main.SimEncoded, 1, input_state, save_events=True)
    return _dense(r[0]["dumps"][0])


def rotated(input_state: int, theta: float, p: float) -> np.ndarray:
    """Data density matrix (128 x 128) after the transversal rotation and dephasing, before any
    syndrome measurement: the environment of the dilation is traced out."""
    r = qsharp.run(code.Main.SimRotated, 1, input_state, theta, p, save_events=True)
    return steane.reduced_data_state(_dense(r[0]["dumps"][0]), 7)


def corrected(input_state: int, theta: float, p: float, shots: int, seed: int | None = None) -> list[tuple[int, np.ndarray]]:
    """`shots` runs of the full round with the Steane block and the correction applied in Q#:
    (measured syndrome, data density matrix after the correction) per shot."""
    r = qsharp.run(code.Main.SimOneRound, shots, input_state, theta, p, save_events=True, seed=seed)
    return [(_syndrome(x["result"]), steane.reduced_data_state(_dense(x["dumps"][0]), 7)) for x in r]


def two_rounds(theta: float, p: float, shots: int, seed: int | None = None) -> list[tuple[int, np.ndarray]]:
    """(round-1 syndrome, data density matrix after the second rotation) per shot, from the
    21-qubit dumps; the state is determined by the syndrome, so repeats are checks."""
    r = qsharp.run(code.Main.SimTwoRounds, shots, theta, p, save_events=True, seed=seed)
    return [(_syndrome(x["result"]), steane.reduced_data_state(_dense(x["dumps"][0]), 14)) for x in r]


# ------------------------------------------------------------------------ shots

def _bits(results) -> np.ndarray:
    return np.array([[1 if str(r) == "One" else 0 for r in shot] for shot in results], dtype=np.uint8)


def _syndrome(res) -> int:
    b = [1 if str(r) == "One" else 0 for r in res]
    return 4 * b[0] + 2 * b[1] + b[2]


def _syndromes(bits: np.ndarray) -> np.ndarray:
    return (4 * bits[:, 0] + 2 * bits[:, 1] + bits[:, 2]).astype(int)


def one_round(input_state: int, theta: float, p: float, basis: int, shots: int, seed: int | None = None, variant: str = "block") -> tuple[np.ndarray, np.ndarray]:
    """(syndrome 0..7, logical eigenvalue +-1) per shot; basis 0 = Z, 1 = X, 2 = Y. The variant
    picks the Steane block (default), joint measurements ("direct") or DrawRandomBool flips
    ("drawn"), the last two for checks."""
    op = {"block": code.Main.RunOneRound, "direct": code.Main.RunOneRoundDirect, "drawn": code.Main.RunOneRoundDrawn}[variant]
    b = _bits(qsharp.run(op, shots, input_state, theta, p, basis, seed=seed))
    return _syndromes(b[:, :3]), steane.decode(b[:, 3:])


def draw_and_measure(shots: int, seed: int | None = None) -> np.ndarray:
    """Rows of [coin 1, coin 2, measurement 1, measurement 2] from DrawAndMeasure."""
    return _bits(qsharp.run(code.Main.DrawAndMeasure, shots, seed=seed))


def two_round_shots(theta: float, p: float, basis: int, shots: int, seed: int | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(s1, s2, logical eigenvalue) per shot for the +theta / -theta protocol on |+_L>."""
    b = _bits(qsharp.run(code.Main.RunTwoRounds, shots, theta, p, basis, seed=seed))
    return _syndromes(b[:, :3]), _syndromes(b[:, 3:6]), steane.decode(b[:, 6:])


# ------------------------------------------------------------------- estimators

def class_stats(syndromes: np.ndarray, values: np.ndarray, trivial: bool) -> dict:
    """Mean +-1 outcome and its standard error within a syndrome class, plus the class rate."""
    keep = (syndromes == 0) if trivial else (syndromes != 0)
    n = int(keep.sum())
    rate = n / len(syndromes)
    if n == 0:
        return dict(n=0, rate=0.0, se_rate=0.0, mean=float("nan"), se=float("nan"))
    v = values[keep].astype(float)
    return dict(n=n, rate=rate, se_rate=float(np.sqrt(rate * (1 - rate) / len(syndromes))),
                mean=float(v.mean()), se=float(v.std(ddof=1) / np.sqrt(n)) if n > 1 else float("nan"))


def ramsey(theta: float, p: float, shots: int, seed: int | None = None) -> dict:
    """One X-readout run on |+_L>: P(X_L = +1 | trivial) and P(X_L = +1 | nontrivial) with
    their standard errors, and the two class probabilities."""
    syn, val = one_round(2, theta, p, 1, shots, seed)
    out = {}
    for name, triv in (("t", True), ("n", False)):
        c = class_stats(syn, val, triv)
        out[name] = dict(p_plus=(1 + c["mean"]) / 2, se=c["se"] / 2, rate=c["rate"], se_rate=c["se_rate"], n=c["n"])
    return out


def tomography(theta: float, p: float, shots: int, seed: int | None = None) -> dict:
    """Twelve runs (4 logical inputs x 3 readout bases), each with its own seed block, and the
    per-class output expectations: ex[cls][input] = {x, y, z, se_x, se_y, se_z}."""
    ex = {"t": {}, "n": {}}
    rates = {"t": [], "n": []}
    k = 0
    for inp in range(4):
        for cls in ex:
            ex[cls][inp] = {}
        for basis, name in ((1, "x"), (2, "y"), (0, "z")):
            syn, val = one_round(inp, theta, p, basis, shots, seed_for(seed, k))
            k += 1
            for cls, triv in (("t", True), ("n", False)):
                c = class_stats(syn, val, triv)
                ex[cls][inp][name] = c["mean"]
                ex[cls][inp]["se_" + name] = c["se"]
                rates[cls].append(c["rate"])
    return dict(ex=ex, rate={c: float(np.mean(v)) for c, v in rates.items()}, runs=k)


def two_round_bloch(theta: float, p: float, shots: int, seed: int | None = None) -> dict:
    """<X_L>, <Y_L> of the output of the +theta / -theta protocol on |+_L> per syndrome pair
    tt, tn, nt, nn, from two runs (X and Y readout), with the angle, the contraction and the
    dephasing of the pair's channel and delta-method standard errors."""
    ax, bx, vx = two_round_shots(theta, p, 1, shots, seed_for(seed, 0))
    ay, by, vy = two_round_shots(theta, p, 2, shots, seed_for(seed, 1))
    out = {}
    for name, (f1, f2) in {"tt": (True, True), "tn": (True, False), "nt": (False, True), "nn": (False, False)}.items():
        kx = ((ax == 0) == f1) & ((bx == 0) == f2)
        ky = ((ay == 0) == f1) & ((by == 0) == f2)
        nx, ny = int(kx.sum()), int(ky.sum())
        x = vx[kx].astype(float)
        y = vy[ky].astype(float)
        mx, my = (x.mean() if nx else float("nan")), (y.mean() if ny else float("nan"))
        sex = x.std(ddof=1) / np.sqrt(nx) if nx > 1 else float("nan")
        sey = y.std(ddof=1) / np.sqrt(ny) if ny > 1 else float("nan")
        c = float(np.hypot(mx, my))
        se_phi = float(np.sqrt((my * sex) ** 2 + (mx * sey) ** 2) / c**2) if c > 0 else float("nan")
        se_c = float(np.sqrt((mx * sex) ** 2 + (my * sey) ** 2) / c) if c > 0 else float("nan")
        rate = (nx + ny) / (2 * shots)
        out[name] = dict(x=float(mx), y=float(my), se_x=float(sex), se_y=float(sey), phi=float(np.arctan2(my, mx)), se_phi=se_phi,
                         contraction=c, se_contraction=se_c, q=(1 - c) / 2, se_q=se_c / 2,
                         rate=rate, se_rate=float(np.sqrt(rate * (1 - rate) / (2 * shots))), n=nx + ny)
    return out
