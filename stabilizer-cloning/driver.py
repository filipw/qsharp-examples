"""Q# <-> Python bridge: compile the Bell-difference circuit to QIR, run it on the
QDK stabilizer simulator, decode outcomes."""

from __future__ import annotations

import os

import numpy as np
from qdk import qsharp, code
from qdk.simulation import run_qir

import learner

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
_PAULI = [qsharp.Pauli.I, qsharp.Pauli.X, qsharp.Pauli.Y, qsharp.Pauli.Z]


def init() -> None:
    qsharp.init(project_root=PROJECT_ROOT, target_profile=qsharp.TargetProfile.Base)


def random_clifford_circuit(n: int, rng: np.random.Generator, depth: int | None = None):
    """A random sequence of H, S and CNOT gates; depth defaults to 3n. Not a uniformly
    random Clifford, but every statistic measured here depends only on the stabilizer
    group being n-dimensional, which any Clifford circuit guarantees."""
    depth = 3 * n if depth is None else depth
    gates, a, b = [], [], []
    for _ in range(depth):
        g = int(rng.integers(3))
        i = int(rng.integers(n))
        j = int(rng.integers(n))
        while j == i:
            j = int(rng.integers(n))
        gates.append(g)
        a.append(i)
        b.append(j)
    return gates, a, b


def bell_difference_samples(n: int, circuit, shots: int, seed: int) -> np.ndarray:
    """`shots` independent Bell difference samples of the state, as rows (a | b)."""
    gates, a, b = circuit
    qir = qsharp.compile(code.Main.BellDifferenceSample, n, gates, a, b)
    results = run_qir(qir, shots, None, seed, "clifford")
    return np.array([learner.bell_difference(r, n) for r in results], dtype=np.uint8)


def measure_generators(n: int, circuit, generators: np.ndarray, seed: int) -> np.ndarray:
    """Signs (0 = +1, 1 = -1) of the given group elements on one fresh copy."""
    gates, a, b = circuit
    paulis = [[_PAULI[c] for c in learner.to_paulis(g, n)] for g in generators]
    qir = qsharp.compile(code.Main.MeasureGenerators, n, gates, a, b, paulis)
    result = run_qir(qir, 1, None, seed, "clifford")[0]
    return np.array([1 if str(r) in ("One", "1") else 0 for r in result], dtype=np.uint8)
