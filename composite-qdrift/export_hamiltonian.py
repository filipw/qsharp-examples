"""Exports the sweep's Hamiltonians, input state and estimate table to JSON.

Everything `experiment.py` evolves comes out of numpy's `default_rng`: the 400 Pauli
strings and their signs (`build_powerlaw`, seed HAM_SEED) and the input-state angles
(`input_angles`, seed STATE_SEED). A port that cannot run numpy's generator (the
Copenhagen app recomputes Trotter states and the exact reference in TypeScript)
needs those draws as data. This writes them, using the project's own
`hamiltonian.py` and constants, and checks the result against the arrays the sweep
saved, so the export is the Hamiltonian the recorded numbers were measured on.

    results/hamiltonian.json        terms, signs, weights, |a_j| per alpha, input angles
    results/estimates_alpha2.0.json results/estimates_alpha2.0.txt as rows

Run with `python3 export_hamiltonian.py` (after `experiment.py` and `estimate_table.py`).
"""

from __future__ import annotations

import json
import os
import re

import numpy as np

import experiment as exp
import hamiltonian as ham

RESULTS = exp.RESULTS


def export_hamiltonian() -> str:
    reference = None
    magnitudes = {}
    for alpha in exp.ALPHAS:
        terms, coeffs = ham.build_powerlaw(exp.N_QUBITS, exp.N_TERMS, alpha, exp.HAM_SEED, exp.MAX_WEIGHT)
        labels = [ham.term_label(t) for t in terms]
        signs = [int(s) for s in np.sign(coeffs)]
        weights = [ham.pauli_weight(t) for t in terms]
        # Only the magnitudes depend on alpha: the strings and signs are the same draws.
        if reference is None:
            reference = (labels, signs, weights)
        elif (labels, signs, weights) != reference:
            raise SystemExit(f"alpha={alpha}: terms or signs differ from alpha={exp.ALPHAS[0]}")
        saved = np.load(os.path.join(RESULTS, f"coeffs_alpha{alpha}.npy"))
        if not np.array_equal(saved, np.abs(coeffs)):
            raise SystemExit(f"alpha={alpha}: |a_j| differ from results/coeffs_alpha{alpha}.npy")
        saved_w = np.load(os.path.join(RESULTS, f"weights_alpha{alpha}.npy"))
        if not np.array_equal(saved_w, weights):
            raise SystemExit(f"alpha={alpha}: Pauli weights differ from results/weights_alpha{alpha}.npy "
                             "(numpy's generator drew different terms than the recorded sweep)")
        magnitudes[str(alpha)] = [float(x) for x in np.abs(coeffs)]

    labels, signs, weights = reference
    angles = [float(x) for x in ham.input_angles(exp.N_QUBITS, exp.STATE_SEED)]
    out = dict(
        n_qubits=exp.N_QUBITS,
        n_terms=exp.N_TERMS,
        max_weight=exp.MAX_WEIGHT,
        time=exp.TIME,
        ham_seed=exp.HAM_SEED,
        state_seed=exp.STATE_SEED,
        alphas=exp.ALPHAS,
        note="a_j = signs[j] * magnitudes[alpha][j] on terms[j]; qubit 0 is the leftmost letter. "
             "Terms are sorted by |a_j| descending, the order the split K refers to.",
        terms=labels,
        signs=signs,
        weights=weights,
        magnitudes=magnitudes,
        angles=angles,
    )
    path = os.path.join(RESULTS, "hamiltonian.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    return path


HEADER = re.compile(r"^alpha = (?P<alpha>[\d.]+), L = (?P<L>\d+), t = (?P<t>[\d.]+); qubit model (?P<qubit>\S+), (?P<qec>.+)$")
TARGET = re.compile(r"^target trace distance (?P<target>\S+)$")
ROW = re.compile(
    r"^\s+(?P<family>trotter|qdrift|composite)\s+(?P<config>.+?)\s{2,}(?P<twoq>\d+)\s+(?P<d>\S+)"
    r"\s+(?P<qubits>[\d,]+)\s+(?P<tstates>[\d,]+)\s+(?P<runtime>[\d.]+) ms$"
)


def export_estimates(alpha: float = 2.0) -> str | None:
    src = os.path.join(RESULTS, f"estimates_alpha{alpha}.txt")
    if not os.path.exists(src):
        return None
    head, rows, target = None, [], None
    for line in open(src):
        line = line.rstrip("\n")
        if head is None and (m := HEADER.match(line)):
            head = m.groupdict()
        elif m := TARGET.match(line):
            target = float(m["target"])
        elif m := ROW.match(line):
            cfg = {k: int(v) for k, v in re.findall(r"(\w+)=(\d+)", m["config"])}
            order, steps = cfg.get("o", 1), cfg.get("r", 1)
            # The sweep CSV's columns: qDRIFT has m = samples = N; a composite's
            # samples are its total qDRIFT draws, 5^(order/2 - 1) * steps * m.
            per = cfg.get("m", cfg.get("N", 0))
            samples = cfg.get("N", 0) if m["family"] == "qdrift" else (
                5 ** (order // 2 - 1) * steps * per if m["family"] == "composite" else 0)
            rows.append(dict(
                target=target,
                family=m["family"],
                K=cfg.get("K", 0),
                order=order,
                steps=steps,
                m=per,
                samples=samples,
                twoq=int(m["twoq"]),
                trace_distance=float(m["d"]),
                physical_qubits=int(m["qubits"].replace(",", "")),
                t_states=int(m["tstates"].replace(",", "")),
                runtime_ms=float(m["runtime"]),
            ))
    if head is None or not rows:
        raise SystemExit(f"{src}: no table found")
    out = dict(alpha=float(head["alpha"]), n_terms=int(head["L"]), time=float(head["t"]),
               qubit_params=head["qubit"], qec_scheme=head["qec"].replace(" ", "_"),
               error_budget="equal to the target", rows=rows)
    path = os.path.join(RESULTS, f"estimates_alpha{alpha}.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    return path


if __name__ == "__main__":
    print("wrote", export_hamiltonian())
    path = export_estimates()
    if path:
        print("wrote", path)
