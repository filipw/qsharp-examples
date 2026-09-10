"""Checks that the Q# Bell difference sampler and the learner do what the analysis assumes.

Run inside the demo's venv: `source .venv/bin/activate && python verify.py`.
"""

from __future__ import annotations

import sys

import numpy as np

import driver
import learner
from tableau import Tableau

failures: list[str] = []


def check(name: str, ok: bool, detail: str) -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        failures.append(name)


def main() -> int:
    driver.init()
    rng = np.random.default_rng(5)

    print("\n1. decoding convention: |0...0> has a Z-only stabilizer group")
    n = 5
    samples = driver.bell_difference_samples(n, ([], [], []), 200, 1)
    check("X-part is zero for |0^n>", not samples[:, n:].any(), f"max X-bit = {samples[:, n:].max()}")
    check("Z-part varies (nontrivial group elements)", samples[:, :n].any(), f"{len(np.unique(samples, axis=0))} distinct outcomes")

    print("\n2. every Bell difference sample lies in the unsigned stabilizer group")
    for n in (3, 8, 24):
        circuit = driver.random_clifford_circuit(n, rng)
        gens = Tableau(n).apply(*circuit).symplectic()
        samples = driver.bell_difference_samples(n, circuit, 300, 11 + n)
        commute = all(learner.symplectic_product(s, g, n) == 0 for s in samples for g in gens)
        check(f"n={n}", commute, "300/300 samples commute with all generators" if commute else "some sample anticommutes")

    print("\n3. samples are uniform over the group (n=4, 16 elements, 4000 shots)")
    n = 4
    circuit = driver.random_clifford_circuit(n, rng)
    samples = driver.bell_difference_samples(n, circuit, 4000, 3)
    _, counts = np.unique(samples, axis=0, return_counts=True)
    expected = 4000 / 16
    chi2 = float(((counts - expected) ** 2 / expected).sum()) if len(counts) == 16 else np.inf
    check("16 distinct outcomes, chi-square (15 dof) < 30", len(counts) == 16 and chi2 < 30,
          f"{len(counts)} outcomes, chi2 = {chi2:.1f}")

    print("\n4. rank after k samples follows Lemma 5")
    n = 8
    trials = 400
    ranks = np.zeros((trials, n + 4), dtype=int)
    for t in range(trials):
        circuit = driver.random_clifford_circuit(n, rng)
        samples = driver.bell_difference_samples(n, circuit, n + 4, 100 + t)
        span = learner.Span(2 * n)
        for k, s in enumerate(samples):
            span.add(s)
            ranks[t, k] = span.rank
    worst = 0.0
    for k in (n - 1, n, n + 1, n + 3):
        emp = float(np.mean(ranks[:, k - 1] == n))
        th = learner.p_learned(k, n)
        se = np.sqrt(th * (1 - th) / trials)
        worst = max(worst, abs(emp - th) / max(se, 1e-9))
        print(f"      k={k:2d}: P[rank=n] measured {emp:.3f}, Lemma 5 {th:.3f}")
    check("within 3 sigma at every k", worst < 3.0, f"worst deviation {worst:.1f} sigma")

    print("\n5. signs read from one copy match the tableau, and the learned state is the true one")
    ok_all = True
    for n in (4, 10, 20):
        circuit = driver.random_clifford_circuit(n, rng)
        tab = Tableau(n).apply(*circuit)
        span = learner.Span(2 * n)
        samples = driver.bell_difference_samples(n, circuit, 4 * n, 500 + n)
        for s in samples:
            span.add(s)
            if span.rank == n:
                break
        gens = span.basis()
        signs = driver.measure_generators(n, circuit, gens, 900 + n)
        truth = [tab.sign_of(g[:n], g[n:]) for g in gens]
        match = span.rank == n and all(t is not None and t == int(s) for t, s in zip(truth, signs))
        ok_all &= match
        print(f"      n={n:2d}: rank {span.rank}, signs {'match' if match else 'MISMATCH'} ({int(signs.sum())} negative)")
    check("learned signed generators == true stabilizer group", ok_all, "all three sizes")

    print("\n6. signs agree on three independent paths (stabilizer sim, full-state sim, tableau)")
    from qdk import qsharp, code
    from qdk.simulation import run_qir
    pauli = [qsharp.Pauli.I, qsharp.Pauli.X, qsharp.Pauli.Y, qsharp.Pauli.Z]
    mism = total = negs = 0
    for s in range(20):
        n = int(rng.integers(3, 11))
        circuit = driver.random_clifford_circuit(n, rng)
        tab = Tableau(n).apply(*circuit)
        gmat = tab.symplectic()
        gens = []
        for _ in range(n):                      # random group elements, so negatives are common
            c = rng.integers(2, size=n)
            if not c.any():
                c[0] = 1
            gens.append((c @ gmat) % 2)
        gens = np.array(gens, dtype=np.uint8)
        gates, a, b = circuit
        qir = qsharp.compile(code.Main.MeasureGenerators, n, gates, a, b,
                             [[pauli[k] for k in learner.to_paulis(g, n)] for g in gens])
        cl = [1 if str(r) == "One" else 0 for r in run_qir(qir, 1, None, 5, "clifford")[0]]
        fs = [1 if str(r) == "One" else 0 for r in run_qir(qir, 1, None, 5, "cpu")[0]]
        tb = [tab.sign_of(g[:n], g[n:]) for g in gens]
        for x, y, z in zip(cl, fs, tb):
            total += 1
            negs += z
            mism += int(x != y or y != z)
    check("no mismatches", mism == 0, f"{total} elements, {negs} negative, {mism} mismatches")

    print("\n7. a graph state's canonical generators are all positive (why learned bases look positive)")
    n = 8
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if rng.random() < 0.5]
    gates, a, b = [0] * n, list(range(n)), [(q + 1) % n for q in range(n)]
    for (i, j) in edges:                        # CZ = H CNOT H on the target
        gates += [0, 2, 0]
        a += [j, i, j]
        b += [i, j, i]
    tab = Tableau(n).apply(gates, a, b)
    canon = []
    for i in range(n):
        v = np.zeros(2 * n, dtype=np.uint8)
        v[n + i] = 1
        for (u, w) in edges:
            if u == i:
                v[w] = 1
            if w == i:
                v[u] = 1
        canon.append(v)
    canon = np.array(canon, dtype=np.uint8)
    t_signs = [tab.sign_of(g[:n], g[n:]) for g in canon]
    q_signs = list(map(int, driver.measure_generators(n, (gates, a, b), canon, 3)))
    check("X_i Z_N(i) all +1 on both paths", not any(t_signs) and not any(q_signs), f"tableau {t_signs}, Q# {q_signs}")

    print("\n8. the hypothesis count N(n, r) = prod(2^i + 1), brute force at n = 3")
    import itertools
    m = 3
    vecs = [np.array(v, dtype=np.uint8) for v in itertools.product([0, 1], repeat=2 * m) if any(v)]
    sp = lambda u, v: int((u[:m] @ v[m:] + u[m:] @ v[:m]) % 2)
    def rref(M):
        M = M.copy() % 2
        r = 0
        for c in range(2 * m):
            piv = next((i for i in range(r, M.shape[0]) if M[i, c]), None)
            if piv is None:
                continue
            M[[r, piv]] = M[[piv, r]]
            for i in range(M.shape[0]):
                if i != r and M[i, c]:
                    M[i] ^= M[r]
            r += 1
        return tuple(map(tuple, M[:r]))
    lag = set()
    for u in vecs:
        for v in vecs:
            if sp(u, v):
                continue
            for w in vecs:
                if sp(u, w) or sp(v, w):
                    continue
                R = rref(np.array([u, v, w]))
                if len(R) == 3:
                    lag.add(R)
    fixed = np.array(next(iter(lag)), dtype=np.uint8)
    counts = []
    for r in range(4):
        counts.append(sum(1 for R in lag if len(rref(np.vstack([np.array(R, dtype=np.uint8), fixed[:r]]))) == 3))
    expect = [round(np.exp(learner.log_lagrangian_completions(m, r))) for r in range(4)]
    check("completions for r = 0..3", counts == expect, f"brute force {counts}, formula {expect}")

    print("\n9. the constant behind the bound")
    p = learner.p_span_limit()
    check("prod(1 - 2^-i) = 0.2888 and p/2 >= 0.14", abs(p - 0.28879) < 1e-4 and p / 2 >= 0.14, f"p = {p:.5f}, p/2 = {p/2:.4f}")

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED: {', '.join(failures)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
