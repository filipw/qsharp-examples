"""Fault-tolerant cost of the three strategies at matched accuracy.

The sweep measures algorithm-level cost (two-qubit gates). This pushes the cheapest
configuration of each strategy that reaches a target trace distance through the
full resource estimator, so the comparison lands in physical qubits, T states and
runtime. The estimator's error budget is set to the target, holding rotation
synthesis to the same standard as the simulation error.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import qsharp
from qsharp.estimator import EstimatorParams, QubitParams, QECScheme

import driver
import experiment as exp
import hamiltonian as ham
import predict as P

TARGETS = [1e-2, 1e-3, 1e-4]


def cheapest(rows, family, target):
    hit = [r for r in rows if r["family"] == family and r["trace_distance"] <= target]
    return min(hit, key=lambda r: r["twoq"]) if hit else None


def main(alpha: float) -> None:
    rows = P.load_sweep(alpha)
    driver.init(seed=1)
    terms, coeffs = ham.build_powerlaw(exp.N_QUBITS, exp.N_TERMS, alpha, exp.HAM_SEED, exp.MAX_WEIGHT)
    angles = list(ham.input_angles(exp.N_QUBITS, exp.STATE_SEED))
    q_all, c_all = driver.qs_paulis(terms), list(coeffs)

    def call(row):
        n, t = exp.N_QUBITS, exp.TIME
        if row["family"] == "trotter":
            return qsharp.code.Main.EstTrotter, (n, angles, q_all, c_all, t, row["order"], row["steps"])
        if row["family"] == "qdrift":
            return qsharp.code.Main.EstQDrift, (n, angles, q_all, c_all, t, row["samples"])
        K = row["K"]
        return qsharp.code.Main.EstComposite, (n, angles, driver.qs_paulis(terms[:K]), list(coeffs[:K]),
                                               driver.qs_paulis(terms[K:]), list(coeffs[K:]), t,
                                               row["order"], row["steps"], row["m"])

    lines = [f"alpha = {alpha}, L = {exp.N_TERMS}, t = {exp.TIME}; qubit model {QubitParams.GATE_NS_E4}, surface code"]
    for target in TARGETS:
        params = EstimatorParams()
        params.error_budget = target
        params.qubit_params.name = QubitParams.GATE_NS_E4
        params.qec_scheme.name = QECScheme.SURFACE_CODE
        lines.append(f"\ntarget trace distance {target:g}")
        lines.append(f"  {'strategy':<10} {'configuration':<26} {'2q gates':>9} {'D':>10} {'phys qubits':>12} {'T states':>11} {'runtime':>12}")
        for family in ("trotter", "qdrift", "composite"):
            row = cheapest(rows, family, target)
            if row is None:
                lines.append(f"  {family:<10} {'not reached in the swept grid':<26}")
                continue
            res = qsharp.estimate(*call(row)[:1], params, *call(row)[1])
            pc = res["physicalCounts"]
            cfg = (f"K={row['K']} o={row['order']} r={row['steps']}" + (f" m={row['m']}" if family == "composite" else "")
                   if family != "qdrift" else f"N={row['samples']}")
            lines.append(f"  {family:<10} {cfg:<26} {row['twoq']:>9.0f} {row['trace_distance']:>10.2e} "
                         f"{pc['physicalQubits']:>12,} {pc['breakdown']['numTstates']:>11,} {pc['runtime'] / 1e6:>9.2f} ms")
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(P.RESULTS, f"estimates_alpha{alpha}.txt"), "w") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 2.0)
