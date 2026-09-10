"""Measure the sample complexity of learning, and hence cloning, stabilizer states.

For each n: random stabilizer states, Bell difference samples in sequence, the rank
of the collected group after every sample, the first k at which it reaches n, the
signs read from one extra copy checked against a classical tableau, and the total
copies the learner consumed (4 per sample + 1 for signs). The rank records are also
what the paper's cloning lower bound is made of (see plot.py).
"""

from __future__ import annotations

import csv
import json
import os
import time

import numpy as np

import driver
import learner
from tableau import Tableau

N_VALUES = [4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128]
STATES = 200            # random stabilizer states per n
EXTRA = 8               # samples recorded beyond n
SEED = 2026
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


def main() -> None:
    os.makedirs(RESULTS, exist_ok=True)
    driver.init()
    rng = np.random.default_rng(SEED)
    rows = []
    t0 = time.time()
    for n in N_VALUES:
        ranks = np.zeros((STATES, n + EXTRA), dtype=np.int32)
        tn = time.time()
        for s in range(STATES):
            circuit = driver.random_clifford_circuit(n, rng)
            tab = Tableau(n).apply(*circuit)
            samples = driver.bell_difference_samples(n, circuit, n + EXTRA, int(rng.integers(2**31)))
            span = learner.Span(2 * n)
            k_star = 0
            for k, v in enumerate(samples, start=1):
                span.add(v)
                ranks[s, k - 1] = span.rank
                if span.rank == n and k_star == 0:
                    k_star = k
            sign_ok = False
            negatives = -1
            if k_star:
                gens = span.basis()
                signs = driver.measure_generators(n, circuit, gens, int(rng.integers(2**31)))
                truth = [tab.sign_of(g[:n], g[n:]) for g in gens]
                sign_ok = all(t is not None and t == int(b) for t, b in zip(truth, signs))
                negatives = int(signs.sum())
            rows.append(dict(n=n, state=s, k_star=k_star, learned=int(bool(k_star) and sign_ok),
                             copies=(4 * k_star + 1) if k_star else -1, negatives=negatives))
        np.save(os.path.join(RESULTS, f"ranks_n{n}.npy"), ranks)
        done = [r for r in rows if r["n"] == n]
        learned = [r for r in done if r["learned"]]
        print(f"n={n:4d}: learned {len(learned)}/{STATES} within {n + EXTRA} samples, "
              f"mean copies {np.mean([r['copies'] for r in learned]):.1f} "
              f"(4n+7.4 = {4 * n + 7.4:.1f}), {time.time() - tn:.1f}s", flush=True)

    with open(os.path.join(RESULTS, "learning.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(RESULTS, "meta.json"), "w") as fh:
        json.dump(dict(n_values=N_VALUES, states=STATES, extra=EXTRA, seed=SEED,
                       seconds=time.time() - t0), fh, indent=2)
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
