"""The sweeps.

A. Trajectories: magic (robustness), entanglement (negativity), the environment's copies of
   both, and the decoded extraction qubit, along gamma, from exact Q# states.
B. Thresholds gamma_-, gamma_+, gamma_e for n = 2..8 by bisection on Q# states, plus the same
   thresholds read off the environment register.
C. The dual witness from shots (two Pauli settings) along gamma.
D. The extracted qubit from post-selected shots.
E. Stabilizer inputs: Bell splitting and three-qubit magic-generators, robustness by LP.
F. Channel conditions: amplitude damping with concurrent dephasing, rebirth iff T_2 > T_1.
"""

from __future__ import annotations

import csv
import json
import os
import time

import numpy as np

import driver
import magic

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
SEED = 2026
GRID = np.linspace(0.0, 1.0, 401)
CURVES = [(2, 0.2), (2, 0.4), (2, 0.5), (3, 0.4), (4, 0.2), (6, 0.2)]
N_VALUES = list(range(2, 9))
ALPHAS = [0.2, 0.3, 0.4]
WITNESS = [(2, 0.4), (4, 0.2)]
WITNESS_POINTS = np.linspace(0.02, 0.98, 25)
WITNESS_SHOTS = 20000
EXTRACT = [(2, 0.2), (4, 0.2), (6, 0.2)]
EXTRACT_SHOTS = 10000
STAB_GRID = np.linspace(0.0, 1.0, 41)
EXPONENTS = [0.5, 0.75, 1.0, 1.25]      # a = 1/2 + Gamma_phi / kappa, T_2 / T_1 = 1 / a
CHANNEL = (3, 0.4)
KT = np.linspace(0.0, 5.0, 201)


def point_seed(point: int) -> int:
    """A base seed per shot-sweep point; each point uses up to four runs spaced by SEED_STRIDE."""
    return SEED + point * 4 * driver.SEED_STRIDE


def pt_min_eig(rho: np.ndarray, n: int) -> float:
    d1, d2 = 2, 2 ** (n - 1)
    pt = rho.reshape(d1, d2, d1, d2).transpose(2, 1, 0, 3).reshape(2**n, 2**n)
    return float(np.linalg.eigvalsh((pt + pt.conj().T) / 2)[0])


def bisect(pred, lo: float, hi: float, iters: int = 40) -> float:
    plo = pred(lo)
    if pred(hi) == plo:
        return float("nan")
    for _ in range(iters):
        mid = (lo + hi) / 2
        if pred(mid) == plo:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def write_csv(name: str, rows: list[dict]) -> None:
    with open(os.path.join(RESULTS, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def evaluate(n: int, alpha: float, gamma: float) -> dict:
    rs, re = driver.damped_cat(n, alpha, gamma)
    p0, pn, c, _ = magic.ghzx(rs)
    rt, ps = magic.decoded_qubit(rs)
    x, _, z = magic.bloch(rt)
    return dict(n=n, alpha=alpha, gamma=gamma, p0=p0, pn=pn, c=c,
                R=magic.robustness_ghzx(rs) - 1, neg=magic.negativity(rs, n, 1),
                R_env=magic.robustness_ghzx(re) - 1, neg_env=magic.negativity(re, n, 1),
                purity=float(np.trace(rs @ rs).real), purity_env=float(np.trace(re @ re).real),
                x=x, z=z, p_success=ps, witness=magic.witness_value(p0, pn, c) - 1)


def sweep_curves() -> None:
    rows = []
    t = time.time()
    for n, alpha in CURVES:
        rows += [evaluate(n, alpha, float(g)) for g in GRID]
        print(f"  curves n={n} alpha={alpha}: {len(GRID)} states", flush=True)
    write_csv("curves.csv", rows)
    print(f"A. curves done in {time.time() - t:.0f}s")


def measured_thresholds(n: int, alpha: float, env: bool = False) -> tuple[float, float, float]:
    """(gamma_-, gamma_+, gamma_e) of the system, or of its environment register."""
    def state(g):
        rs, re = driver.damped_cat(n, alpha, g)
        return re if env else rs

    def magic_lo(g):
        p0, pn, c, _ = magic.ghzx(state(g))
        return c > p0

    def magic_hi(g):
        p0, pn, c, _ = magic.ghzx(state(g))
        return c > pn

    def entangled(g):
        return pt_min_eig(state(g), n) < -1e-13

    # The rebirth facet c = P_n is crossed where the |1^n> amplitude is beta (1-gamma)^{n/2} = beta r,
    # so the bracket ends where that amplitude is still 0.01: the simulator drops amplitudes
    # far below that, which would make c and P_n both read as zero. (The environment's
    # copy of the same facet sits at small gamma, so its bracket is mirrored.)
    eps = 1e-6
    edge = 1 - 0.01 ** (2 / n)
    g_e = bisect(entangled, eps, 1 - eps)
    if alpha >= 1 / np.sqrt(2):                             # r >= 1: no stabilizer window (Theorem 1)
        return float("nan"), float("nan"), g_e
    if env:
        g_plus = bisect(magic_hi, 1 - edge, 1 - eps)      # environment reborn: expect gamma_e
        g_minus = bisect(magic_lo, g_plus, 1 - eps)       # environment magic dies going down: expect 1 - gamma_-
    else:
        g_plus = bisect(magic_hi, eps, edge)
        g_minus = bisect(magic_lo, eps, g_plus)
    return g_minus, g_plus, g_e


def sweep_thresholds() -> None:
    rows = []
    t = time.time()
    for alpha in ALPHAS:
        for n in N_VALUES:
            gm, gp, ge = measured_thresholds(n, alpha)
            em, ep, ee = measured_thresholds(n, alpha, env=True)
            tm, tp, te = magic.thresholds(alpha, n)
            a1, a2 = magic.regime_boundaries(n)
            regime = "I" if alpha < a1 else ("II" if alpha < a2 else "III")
            rows.append(dict(n=n, alpha=alpha, regime=regime, g_minus=gm, g_plus=gp, g_e=ge, g_sum=ge + gp,
                             th_minus=tm, th_plus=tp, th_e=te, env_minus=em, env_plus=ep, env_e=ee))
            print(f"  n={n} alpha={alpha} ({regime}): gamma_- {gm:.5f} ({tm:.5f})  gamma_+ {gp:.5f} ({tp:.5f})  "
                  f"gamma_e {ge:.5f} ({te:.5f})  sum {ge + gp:.6f}  env: reborn at {ep:.5f}, entangled from {ee:.5f}", flush=True)
    write_csv("thresholds.csv", rows)
    print(f"B. thresholds done in {time.time() - t:.0f}s")


def sweep_witness() -> None:
    rows = []
    t = time.time()
    for k, (n, alpha) in enumerate(WITNESS):
        for j, gamma in enumerate(WITNESS_POINTS):
            w = driver.witness_from_shots(n, alpha, float(gamma), WITNESS_SHOTS, point_seed(k * len(WITNESS_POINTS) + j))
            exact = magic.witness_value(*magic.cat_ghzx(alpha, float(gamma), n))
            rows.append(dict(n=n, alpha=alpha, gamma=float(gamma), witness=w["witness"] - 1, se=w["se"], exact=exact - 1,
                             p0=w["p0"], pn=w["pn"], c=w["c"], certified=int(w["witness"] - 1 > 2 * w["se"])))
        print(f"  witness n={n} alpha={alpha}: {len(WITNESS_POINTS)} points x 2 settings x {WITNESS_SHOTS} shots", flush=True)
    write_csv("witness.csv", rows)
    print(f"C. witness done in {time.time() - t:.0f}s")


def sweep_extraction() -> None:
    rows = []
    t = time.time()
    for k, (n, alpha) in enumerate(EXTRACT):
        gm, gp, _ = magic.thresholds(alpha, n)
        gammas = sorted(set(np.round(np.concatenate([np.linspace(0.05, gm, 3), np.linspace(gp, 0.97, 6)]), 4)))
        for j, gamma in enumerate(gammas):
            d = driver.decoded_bloch_from_shots(n, alpha, float(gamma), EXTRACT_SHOTS, point_seed(100 + 20 * k + j))
            rs, _ = driver.damped_cat(n, alpha, float(gamma))
            rt, ps = magic.decoded_qubit(rs)
            x, _, z = magic.bloch(rt)
            rows.append(dict(n=n, alpha=alpha, gamma=float(gamma), x=d["x"], z=d["z"], y=d["y"], se_x=d["se_x"], se_z=d["se_z"], se_y=d["se_y"],
                             p_success=d["p_success"], se_p=d["se_p_success"], exact_x=x, exact_z=z, exact_p=ps,
                             quality=abs(d["x"]) + abs(d["z"]), exact_quality=abs(x) + abs(z),
                             yield_=ps * (magic.robustness_qubit(rt) - 1), R_joint=magic.robustness_ghzx(rs) - 1))
        print(f"  extraction n={n} alpha={alpha}: {len(gammas)} points x 3 bases x {EXTRACT_SHOTS} shots", flush=True)
    write_csv("extraction.csv", rows)
    print(f"D. extraction done in {time.time() - t:.0f}s")


def sweep_stabilizer_inputs() -> None:
    rows = []
    t = time.time()
    states = {2: magic.stabilizer_states(2), 3: magic.stabilizer_states(3)}
    closed = {
        "Phi+": lambda g: g * (1 - g),
        "Psi+": lambda g: 0.0,
        "GHZ3": lambda g: (1 - g) ** 1.5 - (1 - g) ** 3,
        "Phi+ x 0": lambda g: g * (1 - g),
        "00+": lambda g: np.sqrt(1 - g) - (1 - g),
    }
    for name, (n, circ) in driver.CIRCUITS.items():
        for gamma in STAB_GRID:
            rho = driver.damped_circuit(n, circ, float(gamma))
            row = dict(state=name, n=n, gamma=float(gamma), R=magic.robustness_lp(rho, states[n]) - 1, closed=closed[name](float(gamma)))
            row["concurrence"] = magic.concurrence(rho) if n == 2 else float("nan")
            rows.append(row)
        print(f"  {name}: {len(STAB_GRID)} LPs over {len(states[n])} stabilizer states", flush=True)
    write_csv("stabilizer_inputs.csv", rows)
    print(f"E. stabilizer inputs done in {time.time() - t:.0f}s")


def sweep_channel() -> None:
    n, alpha = CHANNEL
    beta = np.sqrt(1 - alpha**2)
    r = alpha / beta
    rows, thr = [], []
    t = time.time()
    for a in EXPONENTS:
        def state(kt):
            gamma = 1 - np.exp(-kt)
            return driver.damped_dephased_cat(n, alpha, float(gamma), magic.dephasing_for_exponent(float(gamma), a))

        for kt in KT:
            rho = state(float(kt))
            p0, pn, c, _ = magic.ghzx(rho)
            rows.append(dict(a=a, kt=float(kt), gamma=1 - np.exp(-kt), R=magic.robustness_ghzx(rho) - 1, neg=magic.negativity(rho, n, 1), margin=c - pn))

        def reborn(kt):
            p0, pn, c, _ = magic.ghzx(state(kt))
            return c > pn

        kt_plus = bisect(reborn, 1e-6, KT[-1])
        q_plus_th = r ** (1 / ((1 - a) * n)) if a < 1 else float("nan")
        thr.append(dict(a=a, T2_over_T1=1 / a, kt_plus=kt_plus, q_plus=np.exp(-kt_plus) if not np.isnan(kt_plus) else float("nan"), q_plus_theory=q_plus_th,
                        kt_plus_theory=-np.log(q_plus_th) if a < 1 else float("nan"),
                        max_margin_after_death=max(rw["margin"] for rw in rows if rw["a"] == a and rw["kt"] > 0.5)))
        print(f"  a={a} (T2/T1 = {1 / a:.2f}): rebirth at kappa t = {kt_plus:.4f} (theory {thr[-1]['kt_plus_theory']:.4f})", flush=True)
    write_csv("channel.csv", rows)
    write_csv("channel_thresholds.csv", thr)
    print(f"F. channel conditions done in {time.time() - t:.0f}s")


STAGES = dict(curves=sweep_curves, thresholds=sweep_thresholds, witness=sweep_witness, extraction=sweep_extraction,
              stabilizer=sweep_stabilizer_inputs, channel=sweep_channel)


def main() -> None:
    """`python experiment.py` runs every sweep; `python experiment.py thresholds channel` only those."""
    import sys

    os.makedirs(RESULTS, exist_ok=True)
    driver.init()
    t0 = time.time()
    for name in (sys.argv[1:] or list(STAGES)):
        STAGES[name]()
    with open(os.path.join(RESULTS, "meta.json"), "w") as fh:
        json.dump(dict(seed=SEED, curves=CURVES, n_values=N_VALUES, alphas=ALPHAS, witness=WITNESS, witness_shots=WITNESS_SHOTS,
                       extract=EXTRACT, extract_shots=EXTRACT_SHOTS, exponents=EXPONENTS, channel=CHANNEL,
                       seconds=time.time() - t0), fh, indent=2)
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
