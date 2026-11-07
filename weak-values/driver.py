"""Thin bridge between the Q# project and numpy: dumped states become joint (spin outcome,
momentum) probabilities, shots become arrays of outcomes and momenta."""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from qdk import code, qsharp

import weak

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

# `qsharp.run(..., seed=s)` seeds shot i with s + i, so two runs whose seed ranges overlap
# share random numbers. Runs that must be independent are spaced by more than any shot count.
SEED_STRIDE = 10_000_000
CHUNK = 5_000


def init() -> None:
    qsharp.init(project_root=PROJECT_ROOT)


def seed_for(seed: int | None, k: int) -> int | None:
    return None if seed is None else seed + k * SEED_STRIDE


def bit_reverse(n: int) -> np.ndarray:
    """rev[i] = i with its n bits reversed."""
    i = np.arange(2**n)
    out = np.zeros_like(i)
    for b in range(n):
        out |= ((i >> b) & 1) << (n - 1 - b)
    return out


# ------------------------------------------------------------------ exact states

def _dense(dump) -> np.ndarray:
    return np.array(dump.as_dense_state(), dtype=complex)


def pointer_state(n: int, sigma: float) -> np.ndarray:
    """The prepared pointer as a function of j (the dump is big-endian, like the register)."""
    r = qsharp.run(code.Main.SimPointer, 1, n, sigma, save_events=True)
    return _dense(r[0]["dumps"][0])


def kicked_state(n: int, sigma: float, dp: float) -> np.ndarray:
    """(2, N): the pointer amplitude attached to spin |0> = |+z> and |1> = |-z> after the kick.
    The spin is the most significant qubit of the dump."""
    r = qsharp.run(code.Main.SimKicked, 1, n, sigma, dp, save_events=True)
    return _dense(r[0]["dumps"][0]).reshape(2, 2**n)


def joint(n: int, sigma: float, dp: float, delta: float) -> np.ndarray:
    """(2, N) amplitudes of (spin outcome, read-out index k) just before the measurements:
    row 1 is the post-selected outcome |f>, row 0 the orthogonal |n>. The pointer qubits of the
    dump are big-endian while the read-out k is little-endian, so the columns are bit-reversed
    back into k order."""
    r = qsharp.run(code.Main.SimExperiment, 1, n, sigma, dp, delta, save_events=True)
    amp = _dense(r[0]["dumps"][0]).reshape(2, 2**n)
    return amp[:, bit_reverse(n)]


def exact(n: int, sigma: float, dp: float, delta: float) -> dict:
    """Outcome probabilities and mean momenta (in units of dp) read off the Q# state."""
    amp = joint(n, sigma, dp, delta)
    p = weak.momenta(n)
    out = {}
    for name, row in (("f", 1), ("n", 0)):
        prob = np.abs(amp[row]) ** 2
        total = float(prob.sum())
        if total > 0:
            shift = float((p * prob).sum() / total / dp) if dp != 0 else float("nan")
            dist = prob / total
        else:   # an outcome that cannot occur (|f> at delta = 0 with no kick)
            shift, dist = float("nan"), np.full_like(prob, np.nan)
        out[name] = dict(prob=total, shift=shift, dist=dist)
    return out


# ------------------------------------------------------------------------ shots

def _run_chunk(args: tuple) -> tuple[np.ndarray, np.ndarray]:
    n, sigma, dp, delta, count, seed = args
    res = qsharp.run(code.Main.RunExperiment, count, n, sigma, dp, delta, seed=seed)
    hit = np.array([str(s) == "One" for s, _ in res], dtype=bool)
    k = np.array([k for _, k in res], dtype=int)
    return hit, k


def shots(n: int, sigma: float, dp: float, delta: float, count: int, seed: int | None = None, pool=None) -> tuple[np.ndarray, np.ndarray]:
    """(post-selected?, read-out index k) per shot; weak.momenta(n)[k] is the momentum. With a process pool the run is cut into chunks of
    CHUNK shots, chunk c seeded seed + c * CHUNK; since qsharp.run seeds shot i with seed + i,
    that is exactly the stream of one serial run with the same seed (verify.py checks it)."""
    if pool is None:
        hit, k = _run_chunk((n, sigma, dp, delta, count, seed))
    else:
        jobs = [(n, sigma, dp, delta, min(CHUNK, count - start), None if seed is None else seed + start)
                for start in range(0, count, CHUNK)]
        parts = list(pool.map(_run_chunk, jobs))
        hit = np.concatenate([h for h, _ in parts])
        k = np.concatenate([x for _, x in parts])
    return hit, k


def make_pool(workers: int | None = None) -> ProcessPoolExecutor:
    """Worker processes with the Q# project loaded, one per logical core up to 8 by default."""
    return ProcessPoolExecutor(max_workers=workers or min(8, os.cpu_count() or 1), initializer=init)


def plane_wave(n: int, k0: int, count: int = 3) -> list[int]:
    return list(qsharp.run(code.Main.RunPlaneWave, count, n, k0))


def mean_shift(p: np.ndarray, dp: float) -> tuple[float, float, int]:
    """Mean momentum in units of dp, its standard error and the number of shots."""
    m = len(p)
    if m == 0:
        return float("nan"), float("nan"), 0
    se = float(p.std(ddof=1) / np.sqrt(m) / dp) if m > 1 else float("nan")
    return float(p.mean() / dp), se, m
