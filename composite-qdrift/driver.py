"""Thin bridge between the Q# project and the numpy reference."""

from __future__ import annotations

import os

import numpy as np
import qsharp

import hamiltonian as ham

_PAULI = [qsharp.Pauli.I, qsharp.Pauli.X, qsharp.Pauli.Y, qsharp.Pauli.Z]
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))


def init(seed: int = 1234) -> None:
    qsharp.init(project_root=PROJECT_ROOT)
    qsharp.set_classical_seed(seed)
    qsharp.set_quantum_seed(seed)


def qs_paulis(terms) -> list[list]:
    return [[_PAULI[p] for p in t] for t in terms]


def run_states(callable_, shots: int, *args, seed_base: int | None = None) -> np.ndarray:
    """Runs a Sim* entry point and returns the dumped state vector of every shot.

    `qsharp.set_classical_seed` is re-applied at the start of *each* shot, so a
    multi-shot run of a randomized simulator would return the same realization
    `shots` times over. Randomized entry points therefore pass `seed_base` and get
    one independent, individually reproducible realization per seed.
    """
    if seed_base is None:
        results = qsharp.run(callable_, shots, *args, save_events=True)
        return np.array([r["dumps"][0].as_dense_state() for r in results], dtype=complex)

    states = []
    for m in range(shots):
        qsharp.set_classical_seed(seed_base + m)
        result = qsharp.run(callable_, 1, *args, save_events=True)
        states.append(result[0]["dumps"][0].as_dense_state())
    return np.array(states, dtype=complex)


def fidelity_stats(states: np.ndarray, target: np.ndarray) -> tuple[float, float]:
    """Mean infidelity 1 - <target|rho|target> and its standard error over shots."""
    overlaps = np.abs(states.conj() @ target) ** 2
    infid = 1.0 - overlaps
    mean = float(infid.mean())
    stderr = float(infid.std(ddof=1) / np.sqrt(len(infid))) if len(infid) > 1 else 0.0
    return mean, stderr
