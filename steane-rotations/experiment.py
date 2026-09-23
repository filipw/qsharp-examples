"""The sweeps.

A. channel:    the exact logical channel of every syndrome, by process tomography on the four
               logical inputs dumped from Q#, over theta and the physical dephasing p.
B. ramsey:     logical Ramsey fringes from shots, P(X_L = +1 | syndrome class), and the class
               probabilities, over theta at p = 0 and p = 0.02.
C. tomography: shot-based process tomography of one round at a few theta (p = 0.02): the
               fitted logical angle, the residual infidelity, the logical dephasing.
D. tworounds:  the +theta / -theta protocol from shots, per syndrome pair: angle, dephasing
               and pair probability against the exact composition.
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time

import numpy as np

import driver
import steane as st

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
SEED = 2026
THETA_EXACT = np.linspace(0.0, np.pi / 2, 65)
P_EXACT = [0.0, 0.01, 0.02, 0.05]
THETA_RAMSEY = np.linspace(0.0, np.pi / 2, 25)
P_RAMSEY = [0.0, 0.02]
RAMSEY_SHOTS = 10000
THETA_TOMO = [0.05 * np.pi, 0.10 * np.pi, 0.15 * np.pi, 0.20 * np.pi, 0.25 * np.pi]
P_TOMO = 0.02
TOMO_SHOTS = 5000
BOOTSTRAP = 400
THETA_TWO = np.linspace(0.0, np.pi / 4, 11)
P_TWO = [0.0, 0.02]
TWO_SHOTS = 15000
CLASSES = {"t": 0, "n": 1}
PAIRS = {"tt": (0, 0), "tn": (0, 1), "nt": (1, 0), "nn": (1, 1)}


def point_seed(point: int) -> int:
    """A base seed per shot-sweep point; a point uses at most 16 runs spaced by SEED_STRIDE."""
    return SEED + point * 16 * driver.SEED_STRIDE


def write_csv(name: str, rows: list[dict]) -> None:
    with open(os.path.join(RESULTS, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


# ------------------------------------------------------------------ A. exact channel

def exact_from_qsharp(theta: float, p: float) -> dict[int, st.LogicalChannel]:
    """Process tomography of every syndrome branch from the four dumped logical inputs."""
    rotated = {i: driver.rotated(i, theta, p) for i in range(4)}
    branches = {i: st.project_data(rotated[i]) for i in range(4)}
    return {s: st.channel_from_states({i: branches[i][s] for i in range(4)}) for s in range(8)}


def sweep_channel() -> None:
    rows = []
    t = time.time()
    for p in P_EXACT:
        for theta in THETA_EXACT:
            chs = exact_from_qsharp(float(theta), p)
            nsym = max(float(np.abs(chs[s].matrix - chs[1].matrix).max()) for s in range(2, 8))
            for name, s in CLASSES.items():
                ch = chs[s]
                rows.append(dict(theta=float(theta), p=p, cls=name, prob=ch.probability, phi=ch.angle, q=ch.dephasing,
                                 leak=ch.leakage, covariance=ch.covariance_residual(), unitarity=ch.unitarity_residual(),
                                 nontrivial_spread=nsym, prob_paper=st.prob_paper(theta, p, s == 0),
                                 phi_paper=st.angle_paper(theta, p, s == 0), q_paper=st.dephasing_paper(theta, p, s == 0),
                                 phi_ideal=st.ideal_angle(theta, s == 0), prob_ideal=st.ideal_prob(theta, s == 0)))
        print(f"  channel p={p}: {len(THETA_EXACT)} angles x 4 inputs dumped", flush=True)
    write_csv("channel.csv", rows)
    print(f"A. channel done in {time.time() - t:.0f}s")


# ------------------------------------------------------------------------ B. Ramsey

def sweep_ramsey() -> None:
    rows = []
    t = time.time()
    k = 0
    for p in P_RAMSEY:
        for theta in THETA_RAMSEY:
            r = driver.ramsey(float(theta), p, RAMSEY_SHOTS, point_seed(k))
            k += 1
            for name, s in CLASSES.items():
                ch = st.exact_channel(float(theta), p, s)
                exact = (1 + (1 - 2 * ch.dephasing) * np.cos(ch.angle)) / 2 if ch.probability > 1e-9 else float("nan")
                rows.append(dict(theta=float(theta), p=p, cls=name, p_plus=r[name]["p_plus"], se=r[name]["se"], n=r[name]["n"],
                                 rate=r[name]["rate"], se_rate=r[name]["se_rate"], exact=exact,
                                 exact_rate=ch.probability if s == 0 else 1 - st.exact_channel(float(theta), p, 0).probability))
        print(f"  ramsey p={p}: {len(THETA_RAMSEY)} angles x {RAMSEY_SHOTS} shots", flush=True)
    write_csv("ramsey.csv", rows)
    print(f"B. ramsey done in {time.time() - t:.0f}s")


# -------------------------------------------------------------------- C. tomography

def fit_from_expectations(ex: dict) -> tuple[float, float, float]:
    """(fitted angle, residual infidelity, dephasing of the residual) from the Bloch map."""
    M, _ = st.bloch_from_expectations(ex)
    phi, r = st.fit_rotation(M)
    return phi, r, float((1 - np.hypot(M[0, 0] + M[1, 1], M[1, 0] - M[0, 1]) / 2) / 2)


def bootstrap_fit(ex: dict, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parametric bootstrap: resample every expectation from its standard error."""
    phis, rs, qs = [], [], []
    for _ in range(BOOTSTRAP):
        sample = {i: {a: float(np.clip(rng.normal(ex[i][a], ex[i]["se_" + a]), -1, 1)) for a in "xyz"} for i in range(4)}
        phi, r, q = fit_from_expectations(sample)
        phis.append(phi)
        rs.append(r)
        qs.append(q)
    return np.array(phis), np.array(rs), np.array(qs)


def sweep_tomography() -> None:
    rows = []
    rng = np.random.default_rng(SEED)
    t = time.time()
    for k, theta in enumerate(THETA_TOMO):
        tomo = driver.tomography(float(theta), P_TOMO, TOMO_SHOTS, point_seed(100 + k))
        for name, s in CLASSES.items():
            ex = tomo["ex"][name]
            phi, r, q = fit_from_expectations(ex)
            bphi, br, bq = bootstrap_fit(ex, rng)
            ch = st.exact_channel(float(theta), P_TOMO, s)
            Mx, _ = ch.bloch_map()
            phi_x, r_x = st.fit_rotation(Mx)
            rows.append(dict(theta=float(theta), p=P_TOMO, cls=name, phi=phi, se_phi=float(np.std(st.wrap(bphi - phi))), infidelity=r,
                             se_infidelity=float(np.std(br)), q=q, se_q=float(np.std(bq)), rate=tomo["rate"][name],
                             x_plus=ex[2]["x"], y_plus=ex[2]["y"], se_x_plus=ex[2]["se_x"], se_y_plus=ex[2]["se_y"],
                             exact_phi=ch.angle, exact_q=ch.dephasing, exact_infidelity=r_x, exact_prob=ch.probability, exact_rate=ch.probability * (1 if s == 0 else 7),
                             exact_x_plus=(1 - 2 * ch.dephasing) * np.cos(ch.angle), exact_y_plus=(1 - 2 * ch.dephasing) * np.sin(ch.angle),
                             phi_paper=st.angle_paper(theta, P_TOMO, s == 0), q_paper=st.dephasing_paper(theta, P_TOMO, s == 0)))
        print(f"  tomography theta={theta / np.pi:.2f}pi: {tomo['runs']} runs x {TOMO_SHOTS} shots", flush=True)
    write_csv("tomography.csv", rows)
    print(f"C. tomography done in {time.time() - t:.0f}s")


# -------------------------------------------------------------------- D. two rounds

def sweep_tworounds() -> None:
    rows = []
    t = time.time()
    k = 0
    for p in P_TWO:
        for theta in THETA_TWO:
            b = driver.two_round_bloch(float(theta), p, TWO_SHOTS, point_seed(200 + k))
            k += 1
            one = {s: st.exact_channel(float(theta), p, s) for s in (0, 1)}
            for name, (s1, s2) in PAIRS.items():
                ch = st.exact_channel(float(theta), p, s1, second=(-float(theta), p, s2))
                prob = ch.probability * (1 if s1 == 0 else 7) * (1 if s2 == 0 else 7)
                ideal = st.wrap(st.ideal_angle(theta, s1 == 0) - st.ideal_angle(theta, s2 == 0))
                rows.append(dict(theta=float(theta), p=p, pair=name, phi=b[name]["phi"], se_phi=b[name]["se_phi"], q=b[name]["q"], se_q=b[name]["se_q"],
                                 x=b[name]["x"], y=b[name]["y"], se_x=b[name]["se_x"], se_y=b[name]["se_y"],
                                 rate=b[name]["rate"], se_rate=b[name]["se_rate"], n=b[name]["n"],
                                 exact_phi=ch.angle, exact_q=ch.dephasing, exact_rate=prob,
                                 exact_x=(1 - 2 * ch.dephasing) * np.cos(ch.angle), exact_y=(1 - 2 * ch.dephasing) * np.sin(ch.angle),
                                 sum_of_angles=float(st.wrap(one[s1].angle - one[s2].angle)), ideal_phi=float(ideal),
                                 composed_q=(1 - (1 - 2 * one[s1].dephasing) * (1 - 2 * one[s2].dephasing)) / 2))
        print(f"  two rounds p={p}: {len(THETA_TWO)} angles x 2 readouts x {TWO_SHOTS} shots", flush=True)
    write_csv("tworounds.csv", rows)
    print(f"D. two rounds done in {time.time() - t:.0f}s")


STAGES = dict(channel=sweep_channel, ramsey=sweep_ramsey, tomography=sweep_tomography, tworounds=sweep_tworounds)


def main() -> None:
    """`python experiment.py` runs every sweep; `python experiment.py ramsey tworounds` only those.
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
    meta.update(seed=SEED, theta_exact=[float(x) for x in THETA_EXACT[[0, -1]]], n_theta_exact=len(THETA_EXACT), p_exact=P_EXACT,
                n_theta_ramsey=len(THETA_RAMSEY), p_ramsey=P_RAMSEY, ramsey_shots=RAMSEY_SHOTS,
                theta_tomo=[float(x) for x in THETA_TOMO], p_tomo=P_TOMO, tomo_shots=TOMO_SHOTS, bootstrap=BOOTSTRAP,
                n_theta_two=len(THETA_TWO), p_two=P_TWO, two_shots=TWO_SHOTS)
    seconds = dict(meta.get("seconds", {})) if isinstance(meta.get("seconds"), dict) else {}
    t0 = time.time()
    for name in names:
        t = time.time()
        STAGES[name]()
        seconds[name] = round(time.time() - t, 1)
    meta["seconds"] = seconds
    with open(path, "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
