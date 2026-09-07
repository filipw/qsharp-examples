"""Checks that the Q# circuits implement the channels the paper analyses, and that the
classical read-off (closed-form robustness, witness, negativity, extraction) is right.

Run with `python verify.py` (~1 min; the n = 4 stabilizer LP is most of it).
"""

from __future__ import annotations

import sys

import numpy as np

import driver
import magic

failures: list[str] = []


def check(name: str, ok: bool, detail: str) -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        failures.append(name)


def pt_min_eig(rho: np.ndarray, n: int) -> float:
    d1, d2 = 2, 2 ** (n - 1)
    pt = rho.reshape(d1, d2, d1, d2).transpose(2, 1, 0, 3).reshape(2**n, 2**n)
    return float(np.linalg.eigvalsh((pt + pt.conj().T) / 2)[0])


def bisect(pred, lo: float, hi: float, iters: int = 40) -> float:
    plo = pred(lo)
    assert pred(hi) != plo, "predicate does not change sign on the bracket"
    for _ in range(iters):
        mid = (lo + hi) / 2
        if pred(mid) == plo:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def main() -> int:
    driver.init()

    print("\n1. state preparation")
    for n, alpha in ((2, 0.4), (5, 0.2)):
        got = driver.cat(n, alpha)
        ov = abs(np.vdot(magic.cat_state(alpha, n), got))
        check(f"PrepareCat n={n}, alpha={alpha} == numpy", abs(ov - 1) < 1e-12, f"|<numpy|Q#>| = {ov:.14f}")

    print("\n2. the dilation implements amplitude damping, and its environment the complementary channel")
    worst_s = worst_e = worst_spec = 0.0
    for n, alpha, gamma in ((2, 0.4, 0.3), (3, 0.2, 0.7), (4, 0.6, 0.05), (6, 0.3, 0.5)):
        psi = magic.cat_state(alpha, n)
        rho0 = np.outer(psi, psi.conj())
        rs, re = driver.damped_cat(n, alpha, gamma)
        worst_s = max(worst_s, float(np.abs(rs - magic.amplitude_damp(rho0, gamma, n)).max()))
        worst_e = max(worst_e, float(np.abs(re - magic.amplitude_damp(rho0, 1 - gamma, n)).max()))
        ss, se = np.linalg.eigvalsh(rs), np.linalg.eigvalsh(re)
        worst_spec = max(worst_spec, float(np.abs(np.sort(ss)[::-1] - np.sort(se)[::-1]).max()))
    check("system == Kraus channel E_gamma (Eq. 40)", worst_s < 1e-12, f"max |diff| = {worst_s:.1e}")
    check("environment == E_{1-gamma} of the input (Eq. 87)", worst_e < 1e-12, f"max |diff| = {worst_e:.1e}")
    check("system and environment are isospectral (Eq. 88)", worst_spec < 1e-12, f"max |diff| = {worst_spec:.1e}")
    for name, (n, circ) in driver.CIRCUITS.items():
        rho_q = driver.damped_circuit(n, circ, 0.0)
        rho_d = driver.damped_circuit(n, circ, 0.35)
        diff = float(np.abs(rho_d - magic.amplitude_damp(rho_q, 0.35, n)).max())
        check(f"gate-code input {name}: damped Q# state == channel", diff < 1e-12, f"max |diff| = {diff:.1e}")
    p = 0.15
    rho_dd = driver.damped_dephased_cat(3, 0.4, 0.3, p)
    ref = magic.dephase(magic.amplitude_damp(np.outer(magic.cat_state(0.4, 3), magic.cat_state(0.4, 3)), 0.3, 3), p, 3)
    diff = float(np.abs(rho_dd - ref).max())
    p0, p1, c, _ = magic.ghzx(rho_dd)
    P0, Pn, C = magic.cat_ghzx(0.4, 0.3, 3)
    check("phase-flip dilation: coherence x (1-2p)^n, populations unchanged", diff < 1e-12 and abs(c - (1 - 2 * p) ** 3 * C) < 1e-12 and abs(p0 - P0) < 1e-12,
          f"max |diff| = {diff:.1e}, c/(1-2p)^n C = {c / ((1 - 2 * p) ** 3 * C):.12f}")

    print("\n3. the damped cat stays in the real GHZ-X manifold, with the paper's populations and coherence")
    worst = 0.0
    for n, alpha, gamma in ((2, 0.4, 0.3), (5, 0.2, 0.6), (8, 0.3, 0.25)):
        rs, _ = driver.damped_cat(n, alpha, gamma)
        p0, p1, c, resid = magic.ghzx(rs)
        P0, Pn, C = magic.cat_ghzx(alpha, gamma, n)
        worst = max(worst, resid, abs(p0 - P0), abs(p1 - Pn), abs(c - C))
    check("P_0, P_n, c of Eq. (1), no other coherences", worst < 1e-12, f"max deviation {worst:.1e}")

    print("\n4. closed-form robustness (Theorem 4) == brute-force stabilizer-polytope LP")
    states = {n: magic.stabilizer_states(n) for n in (2, 3, 4)}
    counts = {n: len(s) for n, s in states.items()}
    check("stabilizer-state counts 60, 1080, 36720", counts == {2: 60, 3: 1080, 4: 36720}, str(counts))
    ins = {n: int(sum(magic.constant_weight_support(s, n) for s in states[n])) for n in (2, 3, 4)}
    check("magic-insulators (constant-weight support) 8, 32, 220 (Corollary 7)", ins == {2: 8, 3: 32, 4: 220}, str(ins))
    worst = 0.0
    for n, alpha, gammas in ((2, 0.4, (0.1, 0.4, 0.8)), (3, 0.3, (0.05, 0.45, 0.9)), (4, 0.3, (0.2, 0.42, 0.7))):
        for gamma in gammas:
            rs, _ = driver.damped_cat(n, alpha, gamma)
            lp = magic.robustness_lp(rs, states[n])
            cf = magic.robustness_ghzx(rs)
            worst = max(worst, abs(lp - cf))
            print(f"      n={n} alpha={alpha} gamma={gamma}: LP {lp:.10f}  closed form {cf:.10f}")
    check("LP == 1 + 2 max{0, c - P_0, c - P_n} on both magic branches and in the window", worst < 1e-7, f"max |diff| = {worst:.1e}")
    worst_phi = worst_psi = 0.0
    for gamma in (0.1, 0.3, 0.5, 0.7, 0.9):
        worst_phi = max(worst_phi, abs(magic.robustness_lp(driver.damped_circuit(2, driver.CIRCUITS["Phi+"][1], gamma), states[2]) - (1 + gamma * (1 - gamma))))
        worst_psi = max(worst_psi, abs(magic.robustness_lp(driver.damped_circuit(2, driver.CIRCUITS["Psi+"][1], gamma), states[2]) - 1))
    check("Bell splitting (Eq. 30): LP gives R(Phi+) = 1 + gamma(1-gamma), R(Psi+) = 1", worst_phi < 1e-7 and worst_psi < 1e-7,
          f"max |diff| Phi+ {worst_phi:.1e}, Psi+ {worst_psi:.1e}")
    worst = 0.0
    for n in (2, 3, 4):
        d = 2**n
        for s in (1, -1):
            for j in (0, 1):
                W = np.eye(d, dtype=complex)
                W[0, -1] += s
                W[-1, 0] += s
                W[0 if j == 0 else -1, 0 if j == 0 else -1] -= 2
                worst = max(worst, float(np.abs(np.einsum("ki,ij,kj->k", states[n].conj(), W, states[n]).real).max()))
    check("the dual witness W_{s,j} is feasible: max |Tr(W sigma)| over every stabilizer state == 1", abs(worst - 1) < 1e-12, f"max = {worst:.12f}")
    worst_gap = 0.0
    all_valid = True
    for n in (5, 6, 7, 8):
        for gamma in (0.05, 0.3, 0.6, 0.9):
            rs, _ = driver.damped_cat(n, 0.3, gamma)
            lo, up, valid = magic.robustness_certificate(rs)
            all_valid &= valid
            worst_gap = max(worst_gap, abs(up - lo), abs(lo - magic.robustness_ghzx(rs)))
    check("n = 5..8: explicit decomposition (upper bound) meets the witness (lower bound), and both equal the closed form",
          all_valid and worst_gap < 1e-10, f"decompositions valid, max gap {worst_gap:.1e}")

    print("\n5. entanglement: every bipartition loses its negativity at the same gamma_e = r^{2/n}")

    def pt_min_eig_cut(rho, n, m):
        d1, d2 = 2**m, 2 ** (n - m)
        pt = rho.reshape(d1, d2, d1, d2).transpose(2, 1, 0, 3).reshape(2**n, 2**n)
        return float(np.linalg.eigvalsh((pt + pt.conj().T) / 2)[0])

    n, alpha = 4, 0.3
    deaths = [bisect(lambda g, m=m: pt_min_eig_cut(driver.damped_cat(n, alpha, g)[0], n, m) < -1e-13, 0.0, 0.999) for m in (1, 2, 3)]
    check("n=4: death threshold across 1|3, 2|2 and 3|1 cuts agree (Prop. 2, the 2x2 block determinant is cut-independent)",
          max(deaths) - min(deaths) < 1e-6, f"cuts give {deaths[0]:.7f}, {deaths[1]:.7f}, {deaths[2]:.7f}")
    worst = 0.0
    for n, alpha in ((2, 0.4), (3, 0.2), (5, 0.3)):
        ge = bisect(lambda g: pt_min_eig(driver.damped_cat(n, alpha, g)[0], n) < -1e-13, 0.0, 0.999)
        _, _, ge_th = magic.thresholds(alpha, n)
        worst = max(worst, abs(ge - ge_th))
        print(f"      n={n} alpha={alpha}: measured gamma_e {ge:.7f}, r^(2/n) = {ge_th:.7f}")
    check("measured gamma_e == r^(2/n)", worst < 1e-6, f"max |diff| = {worst:.1e}")

    print("\n6. shot-based witness from two settings agrees with the exact value")
    worst_sig = 0.0
    for n, alpha, gamma, seed in ((2, 0.4, 0.15, 11), (2, 0.4, 0.45, 12), (2, 0.4, 0.8, 13), (4, 0.2, 0.7, 14)):
        w = driver.witness_from_shots(n, alpha, gamma, 20000, seed)
        exact = magic.witness_value(*magic.cat_ghzx(alpha, gamma, n))
        sig = abs(w["witness"] - exact) / w["se"]
        worst_sig = max(worst_sig, sig)
        print(f"      n={n} alpha={alpha} gamma={gamma}: Tr(W rho) = {w['witness']:.4f} +- {w['se']:.4f}, exact {exact:.4f} ({sig:.1f} sigma)")
    check("within 3 sigma", worst_sig < 3.0, f"worst {worst_sig:.1f} sigma")
    z = driver.measure_cat(2, 0.4, 0.45, 0, 20000, 21)
    x = driver.measure_cat(2, 0.4, 0.45, 1, 20000, 22)
    zi = float(np.mean(1 - 2 * z[:, 0].astype(int)))
    zz = float(np.mean((1 - 2 * z[:, 0].astype(int)) * (1 - 2 * z[:, 1].astype(int))))
    xx = float(np.mean((1 - 2 * x[:, 0].astype(int)) * (1 - 2 * x[:, 1].astype(int))))
    lhs, rhs = 2 * abs(zi) + 2 * abs(xx), 1 + zz
    P0, Pn, C = magic.cat_ghzx(0.4, 0.45, 2)
    check("n=2 three-Pauli form (Eq. 12): 2|<ZI>| + 2|<XX>| <= 1 + <ZZ> inside the window",
          lhs <= rhs and C <= min(P0, Pn), f"lhs {lhs:.4f} vs rhs {rhs:.4f} at gamma=0.45 (in the window [0.324, 0.564])")

    print("\n7. parity-syndrome extraction (Theorem 3)")
    worst_exact = worst_lossless = 0.0
    for n, alpha in ((2, 0.2), (4, 0.2), (6, 0.3)):
        for gamma in np.linspace(0.05, 0.95, 19):
            rs, _ = driver.damped_cat(n, alpha, gamma)
            P0, Pn, C = magic.cat_ghzx(alpha, gamma, n)
            rt, ps = magic.decoded_qubit(rs)
            ref = np.array([[P0, C], [C, Pn]]) / (P0 + Pn)
            worst_exact = max(worst_exact, float(np.abs(rt - ref).max()), abs(ps - P0 - Pn))
            worst_lossless = max(worst_lossless, abs(ps * (magic.robustness_qubit(rt) - 1) - (magic.robustness_ghzx(rs) - 1)))
    check("decoded state == Eq. (19), success probability == P_0 + P_n", worst_exact < 1e-12, f"max |diff| = {worst_exact:.1e}")
    check("lossless identity (Eq. 21): (P_0 + P_n)[R(rho~) - 1] == R(rho_n) - 1", worst_lossless < 1e-12, f"max |diff| = {worst_lossless:.1e}")
    worst_sig = 0.0
    for n, alpha, gamma, seed in ((2, 0.2, 0.85, 31), (4, 0.2, 0.7, 34), (6, 0.3, 0.5, 37)):
        d = driver.decoded_bloch_from_shots(n, alpha, gamma, 10000, seed)
        rs, _ = driver.damped_cat(n, alpha, gamma)
        rt, ps = magic.decoded_qubit(rs)
        x, y, zc = magic.bloch(rt)
        sig = max(abs(d["x"] - x) / d["se_x"], abs(d["z"] - zc) / d["se_z"], abs(d["y"] - y) / d["se_y"], abs(d["p_success"] - ps) / d["se_p_success"])
        worst_sig = max(worst_sig, sig)
        print(f"      n={n} alpha={alpha} gamma={gamma}: shots x={d['x']:.3f} z={d['z']:.3f} p={d['p_success']:.3f}; exact x={x:.3f} z={zc:.3f} p={ps:.3f} ({sig:.1f} sigma)")
    check("post-selected Bloch vector and success probability from shots, within 3 sigma", worst_sig < 3.0, f"worst {worst_sig:.1f} sigma")

    print("\n8. dephased amplitude damping keeps the reflection (Appendix D)")
    n, alpha, p = 2, 0.3, 0.1127  # eta = (1 - 2p)^2 ~ 0.6
    eta = (1 - 2 * p) ** 2

    def st(g):
        return driver.damped_dephased_cat(n, alpha, g, p)

    g_plus = bisect(lambda g: (lambda q: q[2] > q[1])(magic.ghzx(st(g))), 0.0, 0.999)
    g_e = bisect(lambda g: pt_min_eig(st(g), n) < -1e-13, 0.0, 0.999)
    _, gp_th, ge_th = magic.thresholds(alpha, n, eta)
    check(f"eta={eta:.3f}: gamma_+ = 1 - eta r^(2/n), gamma_e = eta r^(2/n), sum = 1",
          abs(g_plus - gp_th) < 1e-6 and abs(g_e - ge_th) < 1e-6,
          f"measured gamma_+ {g_plus:.6f} (formula {gp_th:.6f}), gamma_e {g_e:.6f} (formula {ge_th:.6f}), sum {g_plus + g_e:.6f}")

    print("\n9. shots: reproducible from a seed, and the driver's runs are independent of each other")
    a = driver.measure_cat(3, 0.4, 0.3, 0, 500, 5)
    b = driver.measure_cat(3, 0.4, 0.3, 0, 500, 5)
    c2 = driver.measure_cat(3, 0.4, 0.3, 0, 500, 6)
    d2 = driver.measure_cat(3, 0.4, 0.3, 0, 500, 5 + driver.SEED_STRIDE)
    shifted = bool((a[1:] == c2[:-1]).all())
    check("same seed -> same outcomes; seeds s and s+1 are the same stream shifted by one shot (QDK seeds shot i with s+i)",
          bool((a == b).all()) and shifted and len(np.unique(a, axis=0)) > 1, f"shifted copy: {shifted}, {len(np.unique(a, axis=0))} distinct strings in 500 shots")
    check("seeds one SEED_STRIDE apart are independent", float((a[1:] == d2[:-1]).all(axis=1).mean()) < 0.9,
          f"shifted agreement {float((a[1:] == d2[:-1]).all(axis=1).mean()):.2f} (would be 1.00 for a shared stream)")
    sz, _ = driver.extract(4, 0.2, 0.7, 0, 2000, 77)
    sx, _ = driver.extract(4, 0.2, 0.7, 1, 2000, 77 + driver.SEED_STRIDE)
    agree = float((sz == sx).all(axis=1).mean())
    check("the three readout runs of the extraction do not share their syndromes", agree < 0.9, f"syndrome agreement between the Z and X runs {agree:.2f}")

    print("\n10. the paper's numbers")
    a1, a2 = magic.regime_boundaries(2)
    gm, gp, ge = magic.thresholds(0.4, 2)
    check("n=2 regime boundaries 1/sqrt10, 1/sqrt5; alpha=0.4 window [0.324, 0.564]",
          abs(a1 - 1 / np.sqrt(10)) < 1e-12 and abs(a2 - 1 / np.sqrt(5)) < 1e-12 and abs(gm - 0.324) < 5e-4 and abs(gp - 0.564) < 5e-4,
          f"alpha_1 = {a1:.4f}, alpha_2 = {a2:.4f}, window [{gm:.3f}, {gp:.3f}], gamma_e = {ge:.3f}")

    print("\n11. the pointwise mirror at n=2 (Eq. 13): magic at gamma == reflected concurrence at 1 - gamma")
    worst = 0.0
    for gamma in (0.6, 0.7, 0.8, 0.9):
        rs, _ = driver.damped_cat(2, 0.4, gamma)
        rm, _ = driver.damped_cat(2, 0.4, 1 - gamma)
        lhs = magic.robustness_ghzx(rs) - 1
        rhs = (1 - gamma) / gamma * magic.concurrence(rm)
        worst = max(worst, abs(lhs - rhs))
        print(f"      gamma={gamma}: R(gamma) - 1 = {lhs:.6f}, (1-gamma)/gamma * C(1-gamma) = {rhs:.6f}")
    check("R_A(gamma) - 1 == (1-gamma)/gamma C(1-gamma) on the reborn branch", worst < 1e-12, f"max |diff| = {worst:.1e}")

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED: {', '.join(failures)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
