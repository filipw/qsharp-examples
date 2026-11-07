"""The sweeps.

A. exact:    the amplification A = <p>_f / dp of the post-selected pointer, read off the Q# state
             (no shots), over the post-selection angle delta for five measurement strengths
             w = dp * sigma, and over w at the angle where the weak value of S_z is 7. Also the
             exact momentum distributions behind the plots.
B. fig2:     the paper's Fig. 2 experiment from shots: delta from -pi/2 to pi/2 at w = 0.06, the strength
             fitted to the paper's green curve.
C. seven:    one long run at the "7" angle with w = 0.01, the result in the title.
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time

import numpy as np

import driver
import weak

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
SEED = 2026
N_QUBITS = 6            # pointer register: 64 grid points
SIGMA = 2.0             # pointer width in grid steps; the grid is 32 sigma long
SZ_TARGET = 7.0         # weak value of S_z = sigma_z / 2 at the post-selection angle
DELTA_SEVEN = weak.delta_for(2 * SZ_TARGET)
W_CURVES = [0.005, 0.02, 0.06, 0.2, 1.0]
_GEOM = np.geomspace(1e-4, np.pi / 2, 150)
DELTA_EXACT = np.concatenate([-_GEOM[::-1], [0.0], _GEOM])
W_SWEEP = np.geomspace(1e-3, 2.0, 61)
W_DIST = 0.01
FIG2_W = 0.06          # fitted to the Gaussian (green) curve of the paper's Fig. 2
FIG2_DELTA = [-np.pi / 2, -0.8, -0.4, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.4, 0.8, np.pi / 2]
FIG2_SHOTS = 40_000
SEVEN_W = 0.01
SEVEN_SHOTS = 1_000_000


def dp_of(w: float) -> float:
    return w / SIGMA


def write_csv(name: str, rows: list[dict]) -> None:
    with open(os.path.join(RESULTS, name), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)


# ------------------------------------------------------------------------ A. exact

def sweep_exact(pool=None) -> None:
    t = time.time()
    rows = []
    for w in W_CURVES:
        for d in DELTA_EXACT:
            e = driver.exact(N_QUBITS, SIGMA, dp_of(w), float(d))
            rows.append(dict(w=w, delta=float(d), a_f=e["f"]["shift"], a_n=e["n"]["shift"], p_f=e["f"]["prob"],
                             a_closed=float(weak.amplification(d, w)), a_n_closed=float(weak.amplification_orthogonal(d, w)),
                             p_closed=float(weak.postselection_probability(d, w))))
        print(f"  exact w={w}: {len(DELTA_EXACT)} angles", flush=True)
    write_csv("exact.csv", rows)

    rows = []
    for w in W_SWEEP:
        e = driver.exact(N_QUBITS, SIGMA, dp_of(float(w)), DELTA_SEVEN)
        d_star, a_max = weak.peak(float(w))
        e_star = driver.exact(N_QUBITS, SIGMA, dp_of(float(w)), d_star)
        rows.append(dict(w=float(w), a_seven=e["f"]["shift"], a_seven_closed=float(weak.amplification(DELTA_SEVEN, w)),
                         a_max=e_star["f"]["shift"], a_max_closed=a_max))
    write_csv("strength.csv", rows)
    print(f"  strength: {len(W_SWEEP)} values of w at delta = {DELTA_SEVEN:.5f} and at the peak", flush=True)

    p = weak.momenta(N_QUBITS)
    d_star, _ = weak.peak(W_DIST)
    seven = driver.exact(N_QUBITS, SIGMA, dp_of(W_DIST), DELTA_SEVEN)["f"]["dist"]
    at_peak = driver.exact(N_QUBITS, SIGMA, dp_of(W_DIST), d_star)["f"]["dist"]
    rest = driver.exact(N_QUBITS, SIGMA, 0.0, np.pi / 2)["f"]["dist"]
    write_csv("distributions.csv", [dict(w=W_DIST, p=float(p[i]), rest=float(rest[i]), seven=float(seven[i]), peak=float(at_peak[i]))
                                    for i in np.argsort(p)])
    print(f"A. exact done in {time.time() - t:.0f}s")


# ---------------------------------------------------------------------- B. Fig. 2

def sweep_fig2(pool=None) -> None:
    t = time.time()
    rows = []
    for i, d in enumerate(FIG2_DELTA):
        hit, k = driver.shots(N_QUBITS, SIGMA, dp_of(FIG2_W), float(d), FIG2_SHOTS, driver.seed_for(SEED, i), pool)
        p = weak.momenta(N_QUBITS)[k]
        a_f, se_f, n_f = driver.mean_shift(p[hit], dp_of(FIG2_W))
        a_n, se_n, _ = driver.mean_shift(p[~hit], dp_of(FIG2_W))
        rate = n_f / FIG2_SHOTS
        rows.append(dict(delta=float(d), w=FIG2_W, shots=FIG2_SHOTS, n_f=n_f, rate=rate, se_rate=float(np.sqrt(rate * (1 - rate) / FIG2_SHOTS)),
                         a_f=a_f, se_f=se_f, a_n=a_n, se_n=se_n,
                         exact_rate=float(weak.postselection_probability(d, FIG2_W)), exact_a_f=float(weak.amplification(d, FIG2_W)),
                         exact_a_n=float(weak.amplification_orthogonal(d, FIG2_W)),
                         weak_value=float(weak.weak_value(d)) if d != 0 else float("nan")))
        print(f"  fig2 delta={d:+.3f}: {n_f} of {FIG2_SHOTS} post-selected, S_z = {a_f / 2:.2f} +- {se_f / 2:.2f} (exact {weak.amplification(d, FIG2_W) / 2:.2f})", flush=True)
    write_csv("fig2.csv", rows)
    print(f"B. fig2 done in {time.time() - t:.0f}s")


# ---------------------------------------------------------------------- C. the 7

def sweep_seven(pool=None) -> None:
    t = time.time()
    dp = dp_of(SEVEN_W)
    hit, k = driver.shots(N_QUBITS, SIGMA, dp, DELTA_SEVEN, SEVEN_SHOTS, driver.seed_for(SEED, 100), pool)
    p = weak.momenta(N_QUBITS)[k]
    a_f, se_f, n_f = driver.mean_shift(p[hit], dp)
    a_n, se_n, _ = driver.mean_shift(p[~hit], dp)
    rate = n_f / SEVEN_SHOTS
    row = dict(delta=DELTA_SEVEN, w=SEVEN_W, shots=SEVEN_SHOTS, n_f=n_f, rate=rate, se_rate=float(np.sqrt(rate * (1 - rate) / SEVEN_SHOTS)),
               a_f=a_f, se_f=se_f, sz=a_f / 2, se_sz=se_f / 2, a_n=a_n, se_n=se_n,
               exact_rate=float(weak.postselection_probability(DELTA_SEVEN, SEVEN_W)),
               exact_a_f=float(weak.amplification(DELTA_SEVEN, SEVEN_W)), exact_a_n=float(weak.amplification_orthogonal(DELTA_SEVEN, SEVEN_W)),
               weak_value=float(weak.weak_value(DELTA_SEVEN)))
    write_csv("seven.csv", [row])
    momenta = weak.momenta(N_QUBITS)
    counts = np.bincount(k[hit], minlength=len(momenta))
    write_csv("seven_counts.csv", [dict(p=float(momenta[i]), f=int(counts[i])) for i in np.argsort(momenta)])
    print(f"  seven: {n_f} of {SEVEN_SHOTS} post-selected, A = {a_f:.2f} +- {se_f:.2f}, S_z = {a_f / 2:.2f} +- {se_f / 2:.2f}"
          f" (exact {weak.amplification(DELTA_SEVEN, SEVEN_W) / 2:.3f}, weak value {weak.weak_value(DELTA_SEVEN) / 2:.3f})", flush=True)
    print(f"C. seven done in {time.time() - t:.0f}s")


STAGES = dict(exact=sweep_exact, fig2=sweep_fig2, seven=sweep_seven)


def main() -> None:
    """`python experiment.py` runs every stage; `python experiment.py fig2` only that one.
    meta.json keeps the settings and the seconds of the last run of every stage."""
    os.makedirs(RESULTS, exist_ok=True)
    names = sys.argv[1:] or list(STAGES)
    unknown = [n for n in names if n not in STAGES]
    if unknown:
        sys.exit(f"unknown stage(s) {unknown}; choose from {list(STAGES)}")
    driver.init()
    path = os.path.join(RESULTS, "meta.json")
    meta = {}
    if os.path.exists(path):
        with open(path) as fh:
            meta = json.load(fh)
    meta.update(seed=SEED, n_qubits=N_QUBITS, sigma=SIGMA, sz_target=SZ_TARGET, delta_seven=DELTA_SEVEN, w_curves=W_CURVES,
                n_delta_exact=len(DELTA_EXACT), w_sweep=[float(W_SWEEP[0]), float(W_SWEEP[-1]), len(W_SWEEP)], w_dist=W_DIST,
                fig2_w=FIG2_W, fig2_delta=FIG2_DELTA, fig2_shots=FIG2_SHOTS, seven_w=SEVEN_W, seven_shots=SEVEN_SHOTS,
                chunk=driver.CHUNK)
    seconds = dict(meta.get("seconds", {}))
    t0 = time.time()
    with driver.make_pool() as pool:
        for name in names:
            t = time.time()
            STAGES[name](pool)
            seconds[name] = round(time.time() - t, 1)
    meta["seconds"] = seconds
    with open(path, "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
