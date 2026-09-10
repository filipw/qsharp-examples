"""Correctness checks: the Q# circuits against the exact classical channels.

Every number in the experiment rests on two agreements: (1) the numpy mirrors in
`channels.py` apply gates in exactly the order the Q# operations do, and (2) the
randomized Q# operations sample from exactly the channel `channels.py` evaluates.
Run with `python3 verify.py`.
"""

from __future__ import annotations

import sys

import numpy as np
import qsharp

import channels as ch
import driver
import hamiltonian as ham

N = 5
TIME = 1.5
TOL = 1e-10
failures: list[str] = []


def check(name: str, ok: bool, detail: str) -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        failures.append(name)


def mc_check(name, states, rho_exact, psi_exact):
    """Monte-Carlo mean infidelity from Q# vs the exact channel's, within 3 sigma."""
    mean, se = driver.fidelity_stats(states, psi_exact)
    exact = ch.infidelity(rho_exact, psi_exact)
    z = abs(mean - exact) / max(se, 1e-15)
    check(name, z < 3.0, f"Q# {mean:.4e} +- {se:.1e} vs exact {exact:.4e}  ({z:.1f} sigma)")


def main() -> int:
    driver.init(seed=99)
    angles = list(ham.input_angles(N, 5))
    psi0 = ham.input_state(N, np.array(angles))
    rho0 = np.outer(psi0, psi0.conj())
    terms, coeffs = ham.build_powerlaw(N, 60, 2.0, 7)
    P = ch.pauli_matrices(terms)
    q_all = driver.qs_paulis(terms)
    c_all = list(coeffs)
    K = 8
    q_a, q_b = driver.qs_paulis(terms[:K]), driver.qs_paulis(terms[K:])
    c_a, c_b = list(coeffs[:K]), list(coeffs[K:])
    exact = ham.exact_evolved_state(ham.dense_hamiltonian(terms, coeffs, N), TIME, psi0)

    print("\n1. state preparation")
    got = driver.run_states(qsharp.code.Main.SimInputState, 1, N, angles)[0]
    check("PrepareInput == numpy mirror", abs(abs(np.vdot(psi0, got)) - 1) < TOL,
          f"|<numpy|Q#>| = {abs(np.vdot(psi0, got)):.14f}")

    print("\n2. product formulas: Q# state == numpy mirror, all orders")
    pf = ch.ProductFormula(P, coeffs)
    for order, steps in ((1, 1), (2, 3), (4, 2)):
        st = driver.run_states(qsharp.code.Main.SimTrotter, 1, N, angles, q_all, c_all, TIME, order, steps)[0]
        ref = pf.evolve(order, TIME, steps) @ psi0
        check(f"order {order}, r={steps}", abs(abs(np.vdot(ref, st)) - 1) < TOL,
              f"|<numpy|Q#>| = {abs(np.vdot(ref, st)):.14f}")

    print("\n3. product formulas converge at their advertised order (exact reference)")
    for order, expected, steps in ((1, 1.0, [2, 4, 8, 16]), (2, 2.0, [2, 4, 8, 16]), (4, 4.0, [1, 2, 3, 4])):
        errs = [np.sqrt(max(1 - abs(np.vdot(exact, pf.evolve(order, TIME, r) @ psi0)) ** 2, 1e-16)) for r in steps]
        slope = -np.polyfit(np.log(steps), np.log(errs), 1)[0]
        check(f"order {order}", abs(slope - expected) < 0.25, f"state error ~ r^-{slope:.2f} (expected {expected:.0f})")

    print("\n4. qDRIFT channel: direct map == dense superoperator == Q# Monte Carlo")
    qd = ch.QDriftChannel(P, coeffs)
    rho_r = np.random.default_rng(3).standard_normal((2**N, 2**N)) + 0j
    rho_r = rho_r @ rho_r.conj().T
    rho_r /= np.trace(rho_r)
    tau = qd.lam * TIME / 50
    a, b, c = (qd.apply(rho_r, tau, 7, method="direct"), qd.apply(rho_r, tau, 7), qd.apply_power(rho_r, tau, 7))
    diff = max(np.abs(a - b).max(), np.abs(a - c).max())
    check("direct == dense == squared", diff < 1e-10, f"max |diff| = {diff:.1e}")
    samples = 300
    rho_q = qd.apply(rho0, qd.lam * TIME / samples, samples)
    states = driver.run_states(qsharp.code.Main.SimQDrift, 96, N, angles, q_all, c_all, TIME, samples, seed_base=4242)
    mc_check("Q# qDRIFT == exact channel", states, rho_q, exact)
    errs = [ch.trace_distance(qd.apply(rho0, qd.lam * TIME / n, n), exact) for n in (200, 800, 3200)]
    slope = -np.polyfit(np.log([200, 800, 3200]), np.log(errs), 1)[0]
    check("exact trace distance ~ 1/N", abs(slope - 1.0) < 0.15, f"D ~ N^-{slope:.2f}")

    print("\n5. composite channel: nested structure and sampling")
    comp = ch.Composite(P[:K], coeffs[:K], P[K:], coeffs[K:])
    for order, steps in ((2, 3), (4, 2)):
        # m = 0 empties every B-segment, so the composite is a deterministic nest of
        # A-factors; Q# and numpy must then agree exactly, which pins the outer/inner
        # recursion and its gate order.
        st = driver.run_states(qsharp.code.Main.SimComposite, 1, N, angles, q_a, c_a, q_b, c_b, TIME, order, steps, 0)[0]
        rho_det = comp.evolve(rho0, order, TIME, steps, 0)
        ov = np.real(st.conj() @ rho_det @ st)
        check(f"order {order} skeleton (m=0)", abs(ov - 1) < TOL, f"<Q#|rho_numpy|Q#> = {ov:.14f}")
    for order, steps, m in ((2, 2, 6), (4, 1, 4)):
        rho_c = comp.evolve(rho0, order, TIME, steps, m)
        states = driver.run_states(qsharp.code.Main.SimComposite, 96, N, angles, q_a, c_a, q_b, c_b, TIME, order, steps, m, seed_base=777)
        mc_check(f"order {order}, r={steps}, m={m}: Q# == exact channel", states, rho_c, exact)

    print("\n6. gate accounting == qsharp.logical_counts")
    prep = qsharp.logical_counts(qsharp.code.Main.EstInputState, N, angles)["rotationCount"]
    w_all = [ham.pauli_weight(t) for t in terms]
    for order, steps in ((2, 3), (4, 2)):
        got = qsharp.logical_counts(qsharp.code.Main.EstTrotter, N, angles, q_all, c_all, TIME, order, steps)["rotationCount"] - prep
        exp_ = ch.trotter_cost(order, steps, w_all)[0]
        check(f"Trotter order {order}, r={steps}", got == exp_, f"Q# {got} vs formula {exp_}")
    qb = ch.QDriftChannel(P[K:], coeffs[K:])
    for order, steps, m in ((2, 3, 5), (4, 2, 3)):
        got = qsharp.logical_counts(qsharp.code.Main.EstComposite, N, angles, q_a, c_a, q_b, c_b, TIME, order, steps, m)["rotationCount"] - prep
        exp_ = ch.composite_cost(order, steps, m, w_all[:K], w_all[K:], qb.probs)[0]
        check(f"composite order {order}, r={steps}, m={m}", got == exp_, f"Q# {got} vs formula {exp_}")

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED: {', '.join(failures)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
