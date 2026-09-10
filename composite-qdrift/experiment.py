"""Measure the paper's cost objective on power-law Hamiltonians.

arXiv:2607.19852 (Zlokapa, Allen, Harrow): for H = sum_j a_j h_j with a_1 >= a_2 >= ...
and sum_j a_j = 1, the Hagan-Wiebe composite channel (K largest terms Trotterized, the
rest qDRIFT-sampled, K chosen optimally) reaches trace-distance eps with

    G = O( min_{0<=K<=L} [ K t + t^2 lambda_K^2 / eps ] )        lambda_K = sum_{j>K} a_j

two-qubit gates for every such Hamiltonian, and for every coefficient profile there exist
Hamiltonians on which no circuit can do better (their Theorem 1), so the objective is the
optimal complexity of the class. K = L is plain Trotter, K = 0 plain qDRIFT. For power-law
coefficients a_j ~ j^(-alpha) the optimum slides with eps and minimising the objective
gives a gate count scaling as eps^(-1/(2 alpha - 1)). A single instance can only test the
upper-bound side: that the measured cost follows the objective and its minimisation.

This sweep measures, for several alpha, every (K, order, steps, samples) on a grid:
the exact trace distance of the channel on a fixed entangled input (no Monte Carlo)
and its two-qubit-gate cost. Q# provides the circuits and the deterministic states;
`channels.py` evaluates the randomized channels exactly and is cross-checked against
the Q# circuits by `verify.py`.
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time

import numpy as np
import qsharp

import channels as ch
import driver
import hamiltonian as ham

N_QUBITS = 6
N_TERMS = 400
MAX_WEIGHT = 3
ALPHAS = [1.5, 2.0, 3.0]
TIME = 2.0                     # sum_j |a_j| = 1, so this is the paper's t
STATE_SEED = 3
HAM_SEED = 11

TROTTER_STEPS = {2: [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128], 4: [1, 2, 3, 4, 6, 8, 12, 16]}
QDRIFT_SAMPLES = [int(round(x)) for x in np.logspace(2, 6, 9)]
SPLITS = [1, 2, 4, 8, 16, 32, 64, 128, 256]
COMPOSITE_STEPS = {2: [1, 2, 4, 8, 16], 4: [1, 2, 4]}
SAMPLES_PER_SEGMENT = lambda K: sorted({max(1, K // 4), K, 4 * K})

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
FIELDS = ["alpha", "family", "K", "order", "steps", "m", "samples", "rotations", "twoq",
          "trace_distance", "infidelity", "seconds"]


def run_alpha(alpha: float, writer, log) -> None:
    terms, coeffs = ham.build_powerlaw(N_QUBITS, N_TERMS, alpha, HAM_SEED, MAX_WEIGHT)
    weights = np.array([ham.pauli_weight(t) for t in terms])
    paulis = ch.pauli_matrices(terms)
    L = len(terms)
    angles = list(ham.input_angles(N_QUBITS, STATE_SEED))
    psi0 = ham.input_state(N_QUBITS, np.array(angles))
    rho0 = np.outer(psi0, psi0.conj())
    psi_exact = ham.exact_evolved_state(ham.dense_hamiltonian(terms, coeffs, N_QUBITS), TIME, psi0)
    np.save(os.path.join(RESULTS, f"coeffs_alpha{alpha}.npy"), np.abs(coeffs))
    np.save(os.path.join(RESULTS, f"weights_alpha{alpha}.npy"), weights)

    def emit(family, K, order, steps, m, samples, rot, twoq, rho_or_psi, t0):
        if rho_or_psi.ndim == 1:
            psi = rho_or_psi
            d = float(np.sqrt(max(0.0, 1.0 - abs(np.vdot(psi_exact, psi)) ** 2)))
            inf = float(1.0 - abs(np.vdot(psi_exact, psi)) ** 2)
        else:
            d = ch.trace_distance(rho_or_psi, psi_exact)
            inf = ch.infidelity(rho_or_psi, psi_exact)
        row = dict(alpha=alpha, family=family, K=K, order=order, steps=steps, m=m, samples=samples,
                   rotations=rot, twoq=twoq, trace_distance=d, infidelity=inf, seconds=time.time() - t0)
        writer.writerow(row)
        log.flush()
        print(f"  a={alpha} {family:<9} K={K:>3} o={order} r={steps:>3} m={m:>5}  2q={twoq:>9.0f}  "
              f"D={d:.3e}  ({row['seconds']:.1f}s)", file=log, flush=True)

    # K = L: deterministic product formula on every term, state straight from Q#.
    q_all, c_all = driver.qs_paulis(terms), list(coeffs)
    for order, steps_list in TROTTER_STEPS.items():
        for steps in steps_list:
            t0 = time.time()
            psi = driver.run_states(qsharp.code.Main.SimTrotter, 1, N_QUBITS, angles, q_all, c_all, TIME, order, steps)[0]
            rot, twoq = ch.trotter_cost(order, steps, weights)
            emit("trotter", L, order, steps, 0, 0, rot, twoq, psi, t0)

    # K = 0: pure qDRIFT, exact channel by repeated squaring.
    qd_all = ch.QDriftChannel(paulis, coeffs)
    for n in QDRIFT_SAMPLES:
        t0 = time.time()
        rho = qd_all.apply_power(rho0, qd_all.lam * TIME / n, n)
        rot, twoq = ch.qdrift_cost(n, weights, qd_all.probs)
        emit("qdrift", 0, 1, 1, n, n, rot, twoq, rho, t0)
    del qd_all

    # 0 < K < L: the composite channel, exact.
    for K in SPLITS:
        if K >= L:
            continue
        comp = ch.Composite(paulis[:K], coeffs[:K], paulis[K:], coeffs[K:])
        for order, steps_list in COMPOSITE_STEPS.items():
            for steps in steps_list:
                for m in SAMPLES_PER_SEGMENT(K):
                    t0 = time.time()
                    rho = comp.evolve(rho0, order, TIME, steps, m)
                    rot, twoq = ch.composite_cost(order, steps, m, weights[:K], weights[K:], comp.b.probs)
                    emit("composite", K, order, steps, m, 5 ** (order // 2 - 1) * steps * m, rot, twoq, rho, t0)
        del comp


def main(alphas=ALPHAS) -> None:
    os.makedirs(RESULTS, exist_ok=True)
    driver.init(seed=1)
    meta = dict(n_qubits=N_QUBITS, n_terms=N_TERMS, max_weight=MAX_WEIGHT, alphas=alphas, time=TIME,
                splits=SPLITS, trotter_steps=TROTTER_STEPS, composite_steps=COMPOSITE_STEPS,
                qdrift_samples=QDRIFT_SAMPLES)
    with open(os.path.join(RESULTS, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    for alpha in alphas:
        path = os.path.join(RESULTS, f"sweep_alpha{alpha}.csv")
        with open(path, "w", newline="") as fh, open(os.path.join(RESULTS, f"sweep_alpha{alpha}.log"), "w") as log:
            writer = csv.DictWriter(fh, fieldnames=FIELDS)
            writer.writeheader()
            t0 = time.time()
            run_alpha(alpha, writer, log)
            print(f"alpha={alpha} done in {time.time() - t0:.0f}s -> {path}", flush=True)


if __name__ == "__main__":
    main([float(a) for a in sys.argv[1:]] if len(sys.argv) > 1 else ALPHAS)
