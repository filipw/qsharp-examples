"""Checks that the Q# circuits implement the weak Stern-Gerlach experiment the paper analyses,
that the pointer shift read off the Q# state is what the closed form says, that the closed form
is what a continuum calculation gives (no grid, no Fourier transform), and that the shots agree
with the exact values.

Run with `python verify.py` (~1.5 min, mostly shots).
"""

from __future__ import annotations

import sys

import numpy as np
from scipy import integrate, stats

import driver
import weak

failures: list[str] = []
# Shot runs here are seeded SEED + k * SEED_STRIDE (qsharp.run seeds shot i with seed + i), half
# a stride away from experiment.py's seeds, so the two scripts never share random numbers.
SEED = 5_000_000
N, SIGMA = 6, 2.0


def nmax(*xs) -> float:
    """max that propagates NaN (Python's max(0.0, nan) is 0.0), so a NaN fails the check."""
    return float(np.max(np.asarray(xs, dtype=float)))


def check(name: str, ok: bool, detail: str) -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        failures.append(name)


def info(name: str, detail: str) -> None:
    print(f"  [INFO] {name}: {detail}")


def z_score(est: float, se: float, exact: float) -> float:
    return abs(est - exact) / se if se > 0 else (0.0 if abs(est - exact) < 1e-12 else float("inf"))


def continuum_shift(delta: float, dp: float, density, half_width: float) -> tuple[float, float]:
    """(P(f), <p>_f / dp) for a real pointer wave function with |psi(z)|^2 = density(z) on the
    real line, by quadrature: the post-selected pointer is psi(z) g(z) with
    g = (<f|+z> e^{i dp z} + <f|-z> e^{-i dp z}) / sqrt(2), and for real psi the momentum
    expectation needs only g': <p> = int |psi|^2 Im(g* g') / int |psi|^2 |g|^2. No grid and no
    Fourier transform, so this is independent of both the Q# read-out and the numpy mirror."""
    a, b = weak.final_state(delta).conj() / np.sqrt(2)

    def g(z):
        return a * np.exp(1j * dp * z) + b * np.exp(-1j * dp * z)

    def dg(z):
        return 1j * dp * (a * np.exp(1j * dp * z) - b * np.exp(-1j * dp * z))

    opts = dict(limit=400, epsabs=1e-14, epsrel=1e-12)
    norm = integrate.quad(lambda z: density(z) * abs(g(z)) ** 2, -half_width, half_width, **opts)[0]
    mom = integrate.quad(lambda z: density(z) * (np.conj(g(z)) * dg(z)).imag, -half_width, half_width, **opts)[0]
    total = integrate.quad(density, -half_width, half_width, **opts)[0]
    return norm / total, mom / norm / dp


def main() -> int:
    driver.init()
    size = 2**N

    print("\n1. the pointer")
    psi_q = driver.pointer_state(N, SIGMA)
    psi = weak.pointer(N, SIGMA)
    z = weak.grid(N)
    p = weak.momenta(N)
    rest = np.abs(np.fft.fft(psi) / np.sqrt(size)) ** 2
    check("PreparePointer is the numpy Gaussian", np.abs(psi_q - psi).max() < 1e-9, f"max |diff| = {np.abs(psi_q - psi).max():.1e}")
    pos_sd = float(np.sqrt((z * z * psi**2).sum()))
    mom_sd = float(np.sqrt((p * p * rest).sum()))
    check("position standard deviation sigma, momentum standard deviation 1 / (2 sigma), zero mean momentum",
          abs(pos_sd - SIGMA) < 1e-10 and abs(mom_sd - 1 / (2 * SIGMA)) < 1e-10 and abs((p * rest).sum()) < 1e-14,
          f"sd(z) = {pos_sd:.12f}, sd(p) = {mom_sd:.12f}, <p> = {(p * rest).sum():.1e}")
    check("the grid holds the packet: weight beyond 12 sigma in position and beyond 10 widths in momentum (the grid ends at 12.6) below 1e-15",
          float((psi[np.abs(z) > 12 * SIGMA] ** 2).sum()) < 1e-15 and float(rest[np.abs(p) > 10 / (2 * SIGMA)].sum()) < 1e-15,
          f"{(psi[np.abs(z) > 12 * SIGMA] ** 2).sum():.1e} and {rest[np.abs(p) > 10 / (2 * SIGMA)].sum():.1e}")

    print("\n2. the coupling exp(i dp sigma_z z)")
    worst, means = 0.0, []
    for dp in (0.0025, 0.1, -0.3, 1.0):
        kq = driver.kicked_state(N, SIGMA, dp)
        worst = nmax(worst, float(np.abs(kq - weak.kicked(N, SIGMA, dp)).max()))
        for row, sign in ((0, 1), (1, -1)):
            d = np.abs(np.fft.fft(kq[row])) ** 2
            means.append(abs((p * d).sum() / d.sum() - sign * dp))
    check("WeakKick leaves psi e^{+i dp z} on |+z> and psi e^{-i dp z} on |-z> (four kicks, one negative, one strong)", worst < 1e-9, f"max |diff| = {worst:.1e}")
    check("the two components carry mean momentum +dp and -dp", nmax(*means) < 1e-12, f"max |<p> -+ dp| = {nmax(*means):.1e}")

    print("\n3. the momentum read-out")
    reads = {k0: driver.plane_wave(N, k0) for k0 in (0, 1, 7, -5, 31, -32)}
    check("a plane wave of momentum 2 pi k0 / N is read as k0 mod N, every time", all(r == [k0 % size] * 3 for k0, r in reads.items()),
          ", ".join(f"{k0} -> {r[0]}" for k0, r in reads.items()))
    worst_d, worst_p = 0.0, 0.0
    for delta in (-1.2, 0.0, 0.3, np.pi / 2):
        e = driver.exact(N, SIGMA, 0.0, delta)
        for name in ("f", "n"):
            if e[name]["prob"] > 1e-12:   # at delta = 0 with no kick |f> never occurs
                worst_d = nmax(worst_d, float(np.abs(e[name]["dist"] - rest).max()))
        worst_p = nmax(worst_p, abs(e["f"]["prob"] - np.sin(delta / 2) ** 2))
    check("with no kick both outcomes leave the pointer at rest, and P(f) = |<f|+x>|^2 = sin^2(delta/2)", worst_d < 1e-10 and worst_p < 1e-12,
          f"max |dist - rest| = {worst_d:.1e}, max |P(f) - sin^2| = {worst_p:.1e}")

    print("\n4. the spin algebra")
    ds = np.linspace(-np.pi / 2, np.pi / 2, 41)
    ok = all(abs(weak.final_state(d) @ weak.PLUS_X - np.sin(d / 2)) < 1e-14
             and abs(weak.final_state(d) @ weak.SIGMA_Z @ weak.PLUS_X - np.cos(d / 2)) < 1e-14
             and abs(weak.final_state(d) @ weak.orthogonal_state(d)) < 1e-14 for d in ds)
    check("<f|i> = sin(delta/2), <f|sigma_z|i> = cos(delta/2), <n|f> = 0, so the weak value is cot(delta/2)", ok, "41 angles")
    d7 = weak.delta_for(14)
    check("delta = 2 atan(1/14) puts the weak value of sigma_z at 14, of S_z = sigma_z / 2 at 7", abs(weak.weak_value(d7) - 14) < 1e-12,
          f"delta = {d7:.6f} rad = {np.degrees(d7):.3f} deg, |<f|i>|^2 = {np.sin(d7 / 2) ** 2:.6f}")
    check("delta = +-pi/2 post-selects |+z> / |-z>, a strong Stern-Gerlach", np.allclose(weak.final_state(np.pi / 2), [1, 0]) and np.allclose(weak.final_state(-np.pi / 2), [0, -1]),
          f"|f(pi/2)> = {np.round(weak.final_state(np.pi / 2), 12)}, |f(-pi/2)> = {np.round(weak.final_state(-np.pi / 2), 12)}")

    print("\n5. Q# state against the numpy mirror and the closed form")
    rng = np.random.default_rng(7)
    pts = [(float(10 ** rng.uniform(-3, 0.3)), float(rng.uniform(-np.pi / 2, np.pi / 2))) for _ in range(40)]
    worst_amp = worst_a = worst_n = worst_p = 0.0
    for w, delta in pts:
        dp = w / SIGMA
        amp_q = driver.joint(N, SIGMA, dp, delta)
        amp_np = np.array([weak.branch(N, SIGMA, dp, weak.orthogonal_state(delta)), weak.branch(N, SIGMA, dp, weak.final_state(delta))])
        phase = np.vdot(amp_np.ravel(), amp_q.ravel())
        worst_amp = nmax(worst_amp, float(np.abs(amp_q - amp_np * phase / abs(phase)).max()), abs(abs(phase) - 1))
        e = driver.exact(N, SIGMA, dp, delta)
        worst_a = nmax(worst_a, abs(e["f"]["shift"] - weak.amplification(delta, w)) / max(1, abs(weak.amplification(delta, w))))
        worst_n = nmax(worst_n, abs(e["n"]["shift"] - weak.amplification_orthogonal(delta, w)))
        worst_p = nmax(worst_p, abs(e["f"]["prob"] - weak.postselection_probability(delta, w)))
    check("the joint (outcome, momentum) amplitudes Q# leaves equal the numpy mirror, up to one global phase", worst_amp < 1e-9, f"40 random (w, delta), max |diff| = {worst_amp:.1e}")
    check("<p>_f / dp = sin(delta) / (1 - kappa cos(delta)), kappa = e^{-2 w^2}", worst_a < 1e-10, f"max relative |diff| = {worst_a:.1e}")
    check("the discarded outcome reads -sin(delta) / (1 + kappa cos(delta))", worst_n < 1e-10, f"max |diff| = {worst_n:.1e}")
    check("P(f) = (1 - kappa cos(delta)) / 2", worst_p < 1e-12, f"max |diff| = {worst_p:.1e}")

    print("\n6. the closed form from a continuum calculation (no grid, no Fourier transform)")
    gauss = lambda zz: np.exp(-zz * zz / (2 * SIGMA**2))
    worst_c = 0.0
    for w, delta in [(0.005, d7), (0.01, d7), (0.05, 0.1), (0.2, -0.4), (1.0, 1.0), (2.0, 0.3)]:
        prob, a = continuum_shift(delta, w / SIGMA, gauss, 40 * SIGMA)
        worst_c = nmax(worst_c, abs(a - weak.amplification(delta, w)) / abs(weak.amplification(delta, w)), abs(prob - weak.postselection_probability(delta, w)))
    check("Gaussian pointer on the real line: P(f) and <p>_f / dp match the closed form", worst_c < 1e-8, f"6 points, max relative |diff| = {worst_c:.1e}")
    radius = 3.0
    tf = lambda zz: np.maximum(0.0, 1 - (zz / radius) ** 2)
    worst_tf, weak_tf = 0.0, 0.0
    for dp, delta in [(0.002, d7), (0.02, 0.2), (0.2, -0.5), (0.6, 1.2)]:
        u = 2 * dp * radius
        kap = 3 * (np.sin(u) - u * np.cos(u)) / u**3
        closed = np.sin(delta) / (1 - kap * np.cos(delta))
        prob, a = continuum_shift(delta, dp, tf, radius)
        worst_tf = nmax(worst_tf, abs(a - closed) / abs(closed), abs(prob - (1 - kap * np.cos(delta)) / 2))
    _, a_weak = continuum_shift(d7, 1e-5, tf, radius)
    check("1D Thomas-Fermi pointer (the paper's TF model): same formula with kappa its characteristic function 3(sin u - u cos u)/u^3, u = 2 dp R",
          worst_tf < 1e-8, f"4 points, max relative |diff| = {worst_tf:.1e}")
    check("and in the weak limit it reads the weak value too, whatever the pointer's shape", abs(a_weak / 14 - 1) < 1e-6, f"<p>_f / dp = {a_weak:.8f} at dp R = 3e-5")

    print("\n7. limits")
    ratio = {d: weak.amplification(d, 1e-4) / weak.weak_value(d) - 1 for d in (0.05, d7, 0.5, 1.2)}
    lead = {d: -2e-8 * np.cos(d) / (1 - np.cos(d)) for d in ratio}
    check("weak limit: A / cot(delta/2) - 1 = -2 w^2 cos(delta) / (1 - cos(delta)) + O(w^4), so 'weak' means w << sin(delta/2) = |<f|i>|",
          nmax(*(abs(ratio[d] / lead[d] - 1) for d in ratio)) < 1e-3, "w = 1e-4: " + ", ".join(f"delta = {d:.3f}: {ratio[d]:.3e} (leading term {lead[d]:.3e})" for d in ratio))
    e_edge = [driver.exact(N, SIGMA, w / SIGMA, s * np.pi / 2)["f"]["shift"] - s for w in (0.005, 0.2, 2.0) for s in (1, -1)]
    e_zero = [driver.exact(N, SIGMA, w / SIGMA, 0.0)["f"]["shift"] for w in (0.005, 0.2, 2.0)]
    check("A = +-1 at delta = +-pi/2 and A = 0 at delta = 0, for every strength (Q# state)", nmax(*map(abs, e_edge + e_zero)) < 1e-9,
          f"max |diff| = {nmax(*map(abs, e_edge + e_zero)):.1e}")
    strong = nmax(*(abs(driver.exact(N, SIGMA, 3.0 / SIGMA, float(d))["f"]["shift"]) for d in np.linspace(-np.pi / 2, np.pi / 2, 61)))
    check("strong measurement (w = 3): no amplification anywhere, |A| <= 1", strong <= 1 + 1e-9, f"max |A| over 61 angles = {strong:.6f}")
    worst_pk = 0.0
    for w in (0.005, 0.05, 0.3):
        dd = np.linspace(1e-6, np.pi / 2, 400001)
        a = weak.amplification(dd, w)
        d_star, a_max = weak.peak(w)
        worst_pk = nmax(worst_pk, abs(dd[np.argmax(a)] - d_star), abs(a.max() / a_max - 1))
    check("the peak is at cos(delta*) = kappa with A_max = 1 / sqrt(1 - kappa^2)", worst_pk < 1e-5, f"dense grid, max |diff| = {worst_pk:.1e}")
    sat = [weak.peak(w)[1] * w * 2 for w in (1e-2, 1e-3, 1e-4)]
    check("the largest shift A_max dp tends to 1 / (2 sigma), the coherence momentum (one momentum width for this Gaussian)", abs(sat[-1] - 1) < 1e-7, ", ".join(f"{s:.8f}" for s in sat))
    dd = np.linspace(-np.pi / 2, np.pi / 2, 20001)
    snr = nmax(*(float((np.abs(weak.amplification(dd, w)) * np.sqrt(weak.postselection_probability(dd, w))).max()) for w in (1e-3, 0.01, 0.05, 0.3, 2.0)))
    weak_snr = nmax(*(abs(abs(weak.amplification(d, 1e-5)) * np.sqrt(weak.postselection_probability(d, 1e-5)) - np.cos(d / 2)) for d in (0.05, d7, 0.5, 1.2)))
    check("shift x sqrt(post-selection rate) <= 1, the mean-shift form of the Fisher bound below (up to the change in the pointer's width)",
          snr <= 1 + 1e-12 and weak_snr < 1e-6, f"max A sqrt(P(f)) = {snr:.6f} over 5 strengths; it is cos(delta/2) in the weak limit (max |diff| {weak_snr:.1e})")
    p_dp = 1e-6
    lo = np.abs(np.fft.fft(weak.kicked(N, SIGMA, -p_dp)[0])) ** 2 * 2 / size
    hi = np.abs(np.fft.fft(weak.kicked(N, SIGMA, p_dp)[0])) ** 2 * 2 / size
    mid = np.abs(np.fft.fft(weak.kicked(N, SIGMA, 0.0)[0])) ** 2 * 2 / size
    keep = mid > 1e-12   # the tail bins carry nothing but rounding noise
    fisher = float((((hi - lo)[keep] / (2 * p_dp)) ** 2 / mid[keep]).sum())
    check("the bare kick (spin |+z>, read the momentum) has Fisher information 4 sigma^2 per shot about dp, the quantum limit 4 Var(sigma_z z)",
          abs(fisher / (4 * SIGMA**2) - 1) < 1e-6, f"{fisher:.8f} against 4 sigma^2 = {4 * SIGMA**2:g}, so no post-selection can beat it")
    info("at the peak the reading is half the weak value", ", ".join(f"w = {w:g}: A_max / cot(delta*/2) = {weak.peak(w)[1] / weak.weak_value(weak.peak(w)[0]):.4f}" for w in (0.05, 0.005)))

    print("\n8. shots")
    with driver.make_pool() as pool:
        h1, p1 = driver.shots(N, SIGMA, 0.05 / SIGMA, 0.3, 8_000, driver.seed_for(SEED, 0))
        h2, p2 = driver.shots(N, SIGMA, 0.05 / SIGMA, 0.3, 8_000, driver.seed_for(SEED, 0), pool)
        h3, p3 = driver.shots(N, SIGMA, 0.05 / SIGMA, 0.3, 8_000, driver.seed_for(SEED, 0), pool)
        check("parallel chunks reproduce the serial stream exactly, and a seeded run repeats", np.array_equal(h1, h2) and np.array_equal(p1, p2)
              and np.array_equal(h2, h3) and np.array_equal(p2, p3), "8000 shots in 2 chunks")
        worst_z, worst_se, k = 0.0, 0.0, 1
        for w, delta, count in ((0.05, 0.1, 120_000), (0.05, 0.4, 60_000), (0.2, -0.4, 40_000), (2.0, np.pi / 2, 20_000), (2.0, 0.6, 20_000)):
            dp = w / SIGMA
            hit, kk = driver.shots(N, SIGMA, dp, delta, count, driver.seed_for(SEED, k), pool)
            mom = weak.momenta(N)[kk]
            k += 1
            rate = hit.mean()
            zr = z_score(rate, np.sqrt(rate * (1 - rate) / count), weak.postselection_probability(delta, w))
            a_f, se_f, _ = driver.mean_shift(mom[hit], dp)
            a_n, se_n, _ = driver.mean_shift(mom[~hit], dp)
            za, zn = z_score(a_f, se_f, weak.amplification(delta, w)), z_score(a_n, se_n, weak.amplification_orthogonal(delta, w))
            exact = driver.exact(N, SIGMA, dp, delta)
            ex = exact["f"]["dist"]
            for name, se, n_shots in (("f", se_f, int(hit.sum())), ("n", se_n, int((~hit).sum()))):
                d_ex = exact[name]["dist"]
                spread = np.sqrt((weak.momenta(N) ** 2 * d_ex).sum() - (weak.momenta(N) * d_ex).sum() ** 2)
                worst_se = nmax(worst_se, abs(se / (spread / np.sqrt(n_shots) / dp) - 1))
            counts = np.bincount(kk[hit], minlength=size)
            expect = ex * hit.sum()
            keep = expect >= 5
            obs = np.append(counts[keep], counts[~keep].sum())
            exp = np.append(expect[keep], expect[~keep].sum())
            chi = stats.chisquare(obs, exp * obs.sum() / exp.sum())
            print(f"    w = {w}, delta = {delta:+.3f}: {hit.sum()} of {count} post-selected, rate {zr:.2f} sigma, "
                  f"A_f = {a_f:.3f} +- {se_f:.3f} ({za:.2f} sigma), A_n {zn:.2f} sigma, histogram chi^2 p = {chi.pvalue:.3f}")
            worst_z = nmax(worst_z, zr, za, zn)
            check(f"post-selected momentum histogram matches the exact distribution (w = {w}, delta = {delta:+.3f})", chi.pvalue > 1e-3, f"p = {chi.pvalue:.3f}")
    check("post-selection rates and both outcome means within 3.5 sigma of the closed form", worst_z < 3.5, f"worst {worst_z:.2f} sigma over 15 comparisons")
    check("the standard errors are the exact spread of each outcome's distribution over sqrt(shots), within 10%", worst_se < 0.1,
          f"worst relative difference {worst_se:.3f} over 10 means (so the error bars are neither too small nor too large)")

    print("\n9. the paper's Eqs. 5 and 6")
    worst5 = 0.0
    for w in (0.005, 0.06, 0.3, 1.0):
        dp = w / SIGMA
        dp_paper, sigma_paper = 2 * dp, np.sqrt(2) * SIGMA     # relative kick, psi ~ exp(-z^2 / 2 sigma^2)
        for delta in (-1.2, 0.1, d7, np.pi / 2):
            eq5 = dp_paper * np.cos(delta / 2) * np.sin(delta / 2) / (1 - np.exp(-dp_paper**2 * sigma_paper**2 / 4) * np.cos(delta))
            worst5 = nmax(worst5, abs(eq5 / dp - weak.amplification(delta, w)))
    check("Eq. 5 is the closed form here once its dp is the relative kick 2 dp and its sigma is sqrt(2) sigma", worst5 < 1e-12,
          f"16 points, max |diff| = {worst5:.1e}")
    worst6 = 0.0
    for w in (1e-4, 0.005, 0.06, 0.3, 1.0):
        dp = w / SIGMA
        dp_paper, sigma_paper = 2 * dp, np.sqrt(2) * SIGMA
        eq6 = dp_paper / (2 * np.sqrt(1 - np.exp(-dp_paper**2 * sigma_paper**2 / 2)))
        worst6 = nmax(worst6, abs(eq6 / (weak.peak(w)[1] * dp) - 1))
    check("Eq. 6, with e^(-dp^2 sigma^2 / 2) under the root, is the largest shift A_max dp at every strength, and 1/(2 sigma) for a small kick",
          worst6 < 1e-6, f"5 strengths, max relative |diff| = {worst6:.1e}")
    with np.errstate(invalid="ignore"):
        printed = dp_paper / (2 * np.sqrt(1 - np.exp(dp_paper**2 * sigma_paper**2 / 2)))
    info("Eq. 6 as printed, with e^(+dp^2 sigma^2 / 2)", f"the square root of a negative number for any kick: {printed}")

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED: {failures}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
