"""Figures and results/summary.txt."""

from __future__ import annotations

import csv
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import steane as st

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
INK, INK2, INK3, SURF = "#0b0b0b", "#52514e", "#8a8880", "#fcfcfb"
BLUE, ORANGE, AQUA, RED, PURPLE = "#2a78d6", "#eb6834", "#1baf7a", "#d63a3a", "#8e5bd6"
CLS = {"t": (BLUE, "trivial syndrome"), "n": (ORANGE, "nontrivial syndrome")}
PAIR = {"tt": (BLUE, "trivial, trivial"), "tn": (ORANGE, "trivial, nontrivial"), "nt": (AQUA, "nontrivial, trivial"), "nn": (RED, "nontrivial, nontrivial")}
PLUS_L = r"$|\overline{+}\rangle$"


def style(ax):
    ax.set_facecolor(SURF)
    ax.grid(True, color="#e6e5e0", lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#d8d7d1")
    ax.tick_params(colors=INK2, labelsize=8.5)


def load(name: str) -> list[dict]:
    with open(os.path.join(RESULTS, name)) as fh:
        rows = []
        for r in csv.DictReader(fh):
            out = {}
            for k, v in r.items():
                try:
                    out[k] = float(v)
                except ValueError:
                    out[k] = v
            rows.append(out)
        return rows


def col(rows, key, **where):
    return np.array([r[key] for r in rows if all(abs(r[k] - v) < 1e-9 if isinstance(v, float) else r[k] == v for k, v in where.items())])


def upper_bound(ax, x, q, se, c, marker, label):
    """Points on a log axis: a normal error bar where the value is resolved, and a downward
    triangle at the 2σ upper limit where it is consistent with zero. Both are drawn with
    errorbar so the legend keeps them in call order; `label` is used for the resolved points
    and, when none is resolved, for the bounds."""
    x, q, se = (np.asarray(v, dtype=float) for v in (x, q, se))
    res = q - 2 * se > 0
    hidden = label == "_nolegend_"
    if res.any():
        ax.errorbar(x[res], q[res], yerr=2 * se[res], fmt=marker, ms=5, mfc=SURF, mec=c, ecolor=c, elinewidth=1, capsize=2, label=label)
    if (~res).any():
        lab = "_nolegend_" if hidden else (label.replace("(±2σ)", "(2σ upper bound)") if not res.any() else "same, 2σ upper bound where consistent with 0")
        ax.errorbar(x[~res], (q + 2 * se)[~res], fmt="v", ms=5, mfc=SURF, mec=c, mew=1, label=lab)


def unwrap_pi(phi, prob):
    """Angle series in units of pi, unwrapped, with branches of vanishing probability masked."""
    phi = np.array(phi, dtype=float)
    phi[prob < 1e-9] = np.nan
    ok = ~np.isnan(phi)
    out = np.full_like(phi, np.nan)
    out[ok] = np.unwrap(phi[ok])
    return out / np.pi


def rewrap(a, to):
    """`a` shifted by a multiple of 2 pi onto the branch of `to`."""
    a, to = np.asarray(a, dtype=float), np.asarray(to, dtype=float)
    return a + 2 * np.pi * np.round((to - a) / (2 * np.pi))


def figure_rotations(channel, ramsey, meta, lines):
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.6), dpi=200)
    fig.patch.set_facecolor(SURF)

    # (A) the logical angle per syndrome against theta
    ax = axes[0]
    style(ax)
    th = col(channel, "theta", p=0.0, cls="t") / np.pi
    for name, (c, lab) in CLS.items():
        prob = col(channel, "prob", p=0.0, cls=name)
        ax.plot(th, unwrap_pi(col(channel, "phi", p=0.0, cls=name), prob), color=c, lw=2.4, label=f"{lab}, exact Q# channel, p = 0")
        ax.plot(th, unwrap_pi(col(channel, "phi_ideal", p=0.0, cls=name), prob), color=INK, lw=1.1, ls=(0, (5, 3)),
                label=r"Eqs. 22, 23: $-\arg[e^{i\theta}(7 + e^{-4i\theta})^2]$ and 3θ" if name == "t" else None)
        pn = col(channel, "prob", p=0.05, cls=name)
        ax.plot(th, unwrap_pi(col(channel, "phi", p=0.05, cls=name), pn), color=c, lw=1.0, ls=":", label=f"{lab}, p = 0.05")
    ax.plot([0.25, 0.25], [-0.25, 0.75], "o", color=INK, ms=5, mfc=SURF, zorder=5)
    ax.annotate(r"θ = π/4: φ_t = −π/4 ($T^\dagger$), probability 9/16", (0.25, -0.25), textcoords="offset points", xytext=(6, 4), fontsize=7.5, color=INK2)
    ax.annotate(r"θ = π/4: φ_n = 3π/4 = π − π/4 ($\bar{Z}T^\dagger$), 7/16", (0.25, 0.75), textcoords="offset points", xytext=(8, -4), fontsize=7.5, color=INK2)
    ax.plot([0.5], [-0.5], "s", color=INK, ms=5, mfc=SURF, zorder=5)
    ax.annotate(r"θ = π/2: φ_t = −π/2 ($S^{\otimes 7} = \bar{S}^\dagger$), trivial syndrome certain", (0.5, -0.5), textcoords="offset points", xytext=(-215, -13), fontsize=7.5, color=INK2)
    ax.axhline(0, color=INK3, lw=0.8)
    ax.set_xlim(0, 0.5)
    ax.set_ylim(-0.6, 1.6)
    ax.set_xlabel("physical rotation θ / π", color=INK2, fontsize=9.5)
    ax.set_ylabel("logical rotation φ_s / π", color=INK2, fontsize=9.5)
    ax.set_title("Logical angle φ_s against the physical angle θ, per syndrome", color=INK, fontsize=10.5, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.2, labelcolor=INK2, loc="upper left")
    for p in meta["p_exact"]:
        worst = dict(prob=0.0, phi=0.0, q=0.0, leak=0.0, cov=0.0, unit=0.0, spread=0.0)
        for r in channel:
            if abs(r["p"] - p) > 1e-9 or r["prob"] < 1e-9:
                continue
            worst["prob"] = max(worst["prob"], abs(r["prob"] - r["prob_paper"]))
            worst["phi"] = max(worst["phi"], st.angle_diff(r["phi"], r["phi_paper"]))
            worst["q"] = max(worst["q"], abs(r["q"] - r["q_paper"]))
            worst["leak"] = max(worst["leak"], r["leak"])
            worst["cov"] = max(worst["cov"], r["covariance"])
            worst["spread"] = max(worst["spread"], r["nontrivial_spread"])
            if p == 0:
                worst["unit"] = max(worst["unit"], r["unitarity"])
        lines.append(f"panel A/B exact channel, p={p}: max |diff| from Eqs. 14-17: prob {worst['prob']:.1e}, angle {worst['phi']:.1e}, dephasing {worst['q']:.1e}; "
                     f"leakage {worst['leak']:.1e}, non-Z-covariance {worst['cov']:.1e}, spread over the 7 nontrivial syndromes {worst['spread']:.1e}"
                     + (f", non-unitarity {worst['unit']:.1e}" if p == 0 else ""))
    for theta_pi in (0.0625, 0.125, 0.25, 0.5):
        i = int(np.argmin(np.abs(th - theta_pi)))
        pt = col(channel, "prob", p=0.0, cls="t")[i]
        ft = col(channel, "phi", p=0.0, cls="t")[i] / np.pi
        pn = col(channel, "prob", p=0.0, cls="n")[i]
        fn = f"{col(channel, 'phi', p=0.0, cls='n')[i] / np.pi:+.4f} pi" if pn > 1e-9 else "- (p_n = 0)"
        lines.append(f"  theta = {th[i]:.4f} pi, p = 0: p_t = {pt:.4f}, phi_t = {ft:+.4f} pi, phi_n = {fn}" + ("  (T^dagger up to Z_L: phi_n - phi_t = 1 pi)" if abs(th[i] - 0.25) < 1e-9 else ""))

    # (B) syndrome probabilities
    ax = axes[1]
    style(ax)
    for p, c in zip(meta["p_exact"], (INK, PURPLE, RED, ORANGE)):
        if p == 0.01:
            continue
        ax.plot(th, col(channel, "prob", p=p, cls="t"), color=c, lw=1.8 if p == 0 else 1.2, label=f"p_t exact, p = {p}")
        ax.plot(th, col(channel, "prob_paper", p=p, cls="t"), color=c, lw=1.0, ls=(0, (5, 3)), label="Eq. 16: 1/8 + 7λ⁴(3 + cos 4θ)/32" if p == 0 else None)
    for p, c, mk in zip(meta["p_ramsey"], (INK, RED), ("o", "s")):
        ax.errorbar(col(ramsey, "theta", p=p, cls="t") / np.pi, col(ramsey, "rate", p=p, cls="t"), yerr=2 * col(ramsey, "se_rate", p=p, cls="t"),
                    fmt=mk, ms=4, mfc=SURF, mec=c, ecolor=c, elinewidth=1, capsize=2, label=f"p_t from {meta['ramsey_shots']} shots, p = {p} (±2σ)")
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0.4, 1.02)
    ax.set_xlabel("physical rotation θ / π", color=INK2, fontsize=9.5)
    ax.set_ylabel("probability of the trivial syndrome", color=INK2, fontsize=9.5)
    ax.set_title("Probability of the trivial syndrome (the 7 others share 1 − p_t)", color=INK, fontsize=10.5, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.8, labelcolor=INK2, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
    rates = [r for r in ramsey if r["cls"] == "t" and r["se_rate"] > 0]
    worst = max(abs(r["rate"] - r["exact_rate"]) / r["se_rate"] for r in rates)
    lines.append(f"panel B syndrome rates from shots: worst deviation from the exact probability {worst:.1f} sigma over {len(rates)} trivial-syndrome rates "
                 f"(the nontrivial rate of a run is the complement)")

    # (C) Ramsey fringes
    ax = axes[2]
    style(ax)
    grid = np.linspace(0, np.pi / 2, 200)
    for p, ls, mk in zip(meta["p_ramsey"], ("-", (0, (4, 2))), ("o", "s")):
        for name, (c, lab) in CLS.items():
            s = 0 if name == "t" else 1
            chs = [st.exact_channel(g, p, s) for g in grid]
            exact = [(1 + (1 - 2 * ch.dephasing) * np.cos(ch.angle)) / 2 if ch.probability > 1e-9 else np.nan for ch in chs]
            ax.plot(grid / np.pi, exact, color=c, lw=1.6 if p == 0 else 1.2, ls=ls, label=f"{lab}, exact, p = {p}")
            ax.errorbar(col(ramsey, "theta", p=p, cls=name) / np.pi, col(ramsey, "p_plus", p=p, cls=name), yerr=2 * col(ramsey, "se", p=p, cls=name),
                        fmt=mk, ms=4, mfc=SURF if p == 0 else c, mec=c, ecolor=c, elinewidth=1, capsize=2, label=f"{lab}, shots, p = {p} (±2σ)")
    ax.set_xlim(0, 0.5)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("physical rotation θ / π", color=INK2, fontsize=9.5)
    ax.set_ylabel(r"P($\bar{X}$ = +1 | syndrome)", color=INK2, fontsize=9.5)
    ax.set_title(f"Logical Ramsey fringes, {PLUS_L} in, {meta['ramsey_shots']} shots per θ", color=INK, fontsize=10.5, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.8, labelcolor=INK2, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
    resolved = [r for r in ramsey if r["se"] > 0]
    empty = [r for r in ramsey if r["n"] == 0]
    saturated = [r for r in ramsey if r["n"] > 0 and not r["se"] > 0]
    worst = max(abs(r["p_plus"] - r["exact"]) / r["se"] for r in resolved)
    lines.append(f"panel C Ramsey fringes: worst deviation from the exact channel {worst:.1f} sigma over {len(resolved)} points with a resolvable error bar; "
                 f"{len(saturated)} points have every shot at the same outcome (exact P(X=+1) " + ", ".join(f"{r['exact']:.5f}" for r in saturated) + f"); "
                 f"{len(empty)} class/angle combinations have no shots (a syndrome class of probability 0)")

    fig.suptitle("Continuous-angle logical rotations in the Steane code (arXiv:2608.20676): transversal Rz(θ), X-syndrome extraction, correction; "
                 "the per-syndrome logical channel from exact Q# states and from shots", color=INK, fontsize=11.5, x=0.01, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(RESULTS, "rotations.png")
    fig.savefig(out, facecolor=SURF)
    print(f"wrote {out}")


def figure_noise(channel, tomography, tworounds, meta, lines):
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2), dpi=200)
    fig.patch.set_facecolor(SURF)

    # (A) logical dephasing against theta
    ax = axes[0]
    style(ax)
    th = col(channel, "theta", p=0.0, cls="t") / np.pi
    noisy = [p for p in meta["p_exact"] if p > 0]
    for k, p in enumerate(noisy):
        shade = (k + 1) / len(noisy)
        for name, (c, lab) in CLS.items():
            q = col(channel, "q", p=p, cls=name)
            q[col(channel, "prob", p=p, cls=name) < 1e-9] = np.nan
            ax.plot(th, q, color=c, lw=1.2 + 1.2 * shade, alpha=0.35 + 0.65 * shade,
                    label=f"{lab}, exact, p = {', '.join(str(x) for x in noisy)} (bottom to top)" if p == noisy[-1] else None)
            ax.plot(th, col(channel, "q_paper", p=p, cls=name), color=INK, lw=0.9, ls=(0, (5, 3)), label="Eqs. 14–17" if (k == 0 and name == "t") else None)
    for name, (c, lab) in CLS.items():
        upper_bound(ax, col(tomography, "theta", cls=name) / np.pi, col(tomography, "q", cls=name), col(tomography, "se_q", cls=name), c, "D",
                    f"{lab}, tomography from shots, p = {meta['p_tomo']} (±2σ)")
    ax.set_yscale("log")
    ax.set_xlim(0, 0.5)
    ax.set_ylim(1e-8, 1)
    ax.set_xlabel("physical rotation θ / π", color=INK2, fontsize=9.5)
    ax.set_ylabel("logical dephasing q_s", color=INK2, fontsize=9.5)
    ax.set_title("Logical dephasing q_s per syndrome", color=INK, fontsize=10.5, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.8, labelcolor=INK2, loc="lower left")
    lines.append("panel A tomography from shots (p = %s), fitted angle / residual infidelity / dephasing vs exact:" % meta["p_tomo"])
    for r in tomography:
        lines.append(f"  theta = {r['theta'] / np.pi:.2f} pi, {r['cls']}: rate {r['rate']:.4f} (exact {r['exact_rate']:.4f}), "
                     f"phi {r['phi'] / np.pi:+.4f} +- {r['se_phi'] / np.pi:.4f} pi (exact {r['exact_phi'] / np.pi:+.4f}), "
                     f"infidelity {r['infidelity']:.4f} +- {r['se_infidelity']:.4f} (exact {r['exact_infidelity']:.4f} = 2q/3), q {r['q']:.4f} +- {r['se_q']:.4f} (exact {r['exact_q']:.4f})")
    worst = max(max(st.angle_diff(r["phi"], r["exact_phi"]) / r["se_phi"], abs(r["q"] - r["exact_q"]) / r["se_q"]) for r in tomography)
    lines.append(f"  worst deviation of the fitted angle or dephasing from the exact channel: {worst:.1f} sigma over {2 * len(tomography)} comparisons")
    q0 = {p: (st.exact_channel(0.0, p, 0).dephasing, st.exact_channel(0.0, p, 1).dephasing) for p in (0.001, 0.01, 0.1)}
    lines.append("  theta = 0: " + ", ".join(f"p={p}: q_t = {a:.2e} (7p^3 = {7 * p**3:.2e}), q_n = {b:.4f} (3p = {3 * p:.4f})" for p, (a, b) in q0.items()))

    # (B) two rounds: angle per syndrome pair
    ax = axes[1]
    style(ax)
    p = meta["p_two"][-1]
    for name, (c, lab) in PAIR.items():
        tt = col(tworounds, "theta", p=p, pair=name) / np.pi
        exact = np.unwrap(col(tworounds, "exact_phi", p=p, pair=name))
        ideal = np.unwrap(col(tworounds, "ideal_phi", p=p, pair=name))
        shots = rewrap(col(tworounds, "phi", p=p, pair=name), exact)
        ax.plot(tt, exact / np.pi, color=c, lw=1.8, label=f"{lab}: exact, p = {p}")
        ax.plot(tt, ideal / np.pi, color=c, lw=0.9, ls=(0, (5, 3)), label=r"ideal $\varphi_{s_1}(\theta) - \varphi_{s_2}(\theta)$" if name == "tn" else None)
        ax.errorbar(tt, shots / np.pi, yerr=2 * col(tworounds, "se_phi", p=p, pair=name) / np.pi,
                    fmt="o", ms=4, mfc=SURF, mec=c, ecolor=c, elinewidth=1, capsize=2, label=f"shots, {meta['two_shots']} per readout (±2σ)" if name == "tt" else None)
    ax.axhline(0, color=INK3, lw=0.8)
    ax.set_xlim(0, 0.25)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlabel("physical rotation θ / π  (round 1: +θ, round 2: −θ)", color=INK2, fontsize=9.5)
    ax.set_ylabel("total logical angle / π", color=INK2, fontsize=9.5)
    ax.set_title("Two rounds, +θ then −θ: total logical angle per syndrome pair", color=INK, fontsize=10.5, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.2, labelcolor=INK2, loc="upper left")

    # (C) two rounds: dephasing and pair probabilities
    ax = axes[2]
    style(ax)
    for name, (c, lab) in PAIR.items():
        tt = col(tworounds, "theta", p=p, pair=name) / np.pi
        ax.plot(tt, col(tworounds, "exact_q", p=p, pair=name), color=c, lw=1.8, label=f"{lab}: exact")
        upper_bound(ax, tt, col(tworounds, "q", p=p, pair=name), col(tworounds, "se_q", p=p, pair=name), c, "o", "shots (±2σ)" if name == "tt" else "_nolegend_")
    ax.set_yscale("log")
    ax.set_xlim(0, 0.25)
    ax.set_ylim(1e-4, 1)
    ax.set_xlabel("physical rotation θ / π", color=INK2, fontsize=9.5)
    ax.set_ylabel("logical dephasing of the two-round channel", color=INK2, fontsize=9.5)
    ax.set_title(f"Two rounds at p = {p}: dephasing per syndrome pair", color=INK, fontsize=10.5, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.0, labelcolor=INK2, loc="upper left", bbox_to_anchor=(0.0, 0.66))
    ins = ax.inset_axes([0.6, 0.1, 0.37, 0.3])
    style(ins)
    for name, (c, lab) in PAIR.items():
        tt = col(tworounds, "theta", p=p, pair=name) / np.pi
        ins.plot(tt, col(tworounds, "exact_rate", p=p, pair=name), color=c, lw=1.2)
        ins.errorbar(tt, col(tworounds, "rate", p=p, pair=name), yerr=2 * col(tworounds, "se_rate", p=p, pair=name), fmt="o", ms=2.5, mfc=SURF, mec=c, ecolor=c, elinewidth=0.8)
    ins.set_title("pair probabilities", fontsize=7, color=INK2, loc="left")
    ins.tick_params(labelsize=6)
    ins.set_ylim(0, 1)
    for pp in meta["p_two"]:
        lines.append(f"panel B/C two rounds, p = {pp}, per syndrome pair (shots vs exact; the ideal angle is printed on the branch of the exact one):")
        worst, count = 0.0, 0
        for r in [r for r in tworounds if abs(r["p"] - pp) < 1e-9]:
            if r["n"] < 2:
                lines.append(f"  theta = {r['theta'] / np.pi:.3f} pi {r['pair']}: n = {int(r['n'])}, no estimate (exact rate {r['exact_rate']:.4f})")
                continue
            sgs = [st.angle_diff(r["phi"], r["exact_phi"]) / r["se_phi"] if r["se_phi"] > 0 else np.nan, abs(r["q"] - r["exact_q"]) / r["se_q"] if r["se_q"] > 0 else np.nan,
                   abs(r["rate"] - r["exact_rate"]) / r["se_rate"] if r["se_rate"] > 0 else np.nan]
            count += int(np.isfinite(sgs).sum())
            worst = max(worst, float(np.nanmax(sgs)))
            ideal = float(rewrap(r["ideal_phi"], r["exact_phi"]))
            lines.append(f"  theta = {r['theta'] / np.pi:.3f} pi {r['pair']}: n = {int(r['n'])}, phi {r['phi'] / np.pi:+.4f} +- {r['se_phi'] / np.pi:.4f} pi (exact {r['exact_phi'] / np.pi:+.4f}, "
                         f"ideal {ideal / np.pi:+.4f}), q {r['q']:.4f} +- {r['se_q']:.4f} (exact {r['exact_q']:.4f}), rate {r['rate']:.4f} (exact {r['exact_rate']:.4f}); {float(np.nanmax(sgs)):.1f} sigma")
        lines.append(f"  worst deviation over angle, dephasing and rate: {worst:.1f} sigma over {count} comparisons with a finite error bar (the four pair rates at one theta sum to 1)")

    fig.suptitle("Dephasing, tomography and the two-round protocol (arXiv:2608.20676)", color=INK, fontsize=11.5, x=0.01, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(RESULTS, "noise.png")
    fig.savefig(out, facecolor=SURF)
    print(f"wrote {out}")


def figure_bloch(tomography, tworounds, meta, lines):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), dpi=200)
    fig.patch.set_facecolor(SURF)
    circle = np.linspace(0, 2 * np.pi, 300)
    for ax, title in zip(axes, (f"One round at p = {meta['p_tomo']}: {PLUS_L} after Rz(θ)$^{{\\otimes 7}}$ and correction", f"Two rounds at p = {meta['p_two'][-1]}: +θ then −θ")):
        style(ax)
        ax.set_aspect("equal")
        ax.plot(np.cos(circle), np.sin(circle), color=INK3, lw=0.8)
        ax.plot([1], [0], "*", color=INK, ms=10, label=f"input {PLUS_L}")
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-1.15, 1.15)
        ax.set_xlabel(r"$\langle\bar{X}\rangle$", color=INK2, fontsize=9.5)
        ax.set_ylabel(r"$\langle\bar{Y}\rangle$", color=INK2, fontsize=9.5)
        ax.set_title(title, color=INK, fontsize=10.5, loc="left", fontweight="bold")
    ax = axes[0]
    for name, (c, lab) in CLS.items():
        rows = [r for r in tomography if r["cls"] == name]
        ax.errorbar([r["x_plus"] for r in rows], [r["y_plus"] for r in rows], xerr=[2 * r["se_x_plus"] for r in rows], yerr=[2 * r["se_y_plus"] for r in rows],
                    fmt="s", ms=5, color=c, ecolor=c, elinewidth=1, capsize=2, label=f"{lab}: shots (±2σ)")
        ax.plot([r["exact_x_plus"] for r in rows], [r["exact_y_plus"] for r in rows], "o", ms=7, mfc="none", mec=c, mew=1.2, label=f"{lab}: exact")
        for r in rows:
            if name == "n" or r["theta"] / np.pi > 0.12:
                ax.annotate(f"{r['theta'] / np.pi:.2f}π", (r["x_plus"], r["y_plus"]), textcoords="offset points", xytext=(5, 4), fontsize=6.5, color=c)
    ax.legend(frameon=False, fontsize=7.2, labelcolor=INK2, loc="lower left")
    ax = axes[1]
    p = meta["p_two"][-1]
    for name, (c, lab) in PAIR.items():
        rows = [r for r in tworounds if r["pair"] == name and abs(r["p"] - p) < 1e-9 and r["theta"] > 1e-9]
        ax.errorbar([r["x"] for r in rows], [r["y"] for r in rows], xerr=[2 * r["se_x"] for r in rows], yerr=[2 * r["se_y"] for r in rows],
                    fmt="s", ms=4, color=c, ecolor=c, elinewidth=0.8, capsize=1.5, label=f"{lab}: shots (±2σ)")
        ax.plot([r["exact_x"] for r in rows], [r["exact_y"] for r in rows], "o", ms=6, mfc="none", mec=c, mew=1.1, label=f"{lab}: exact")
    ax.legend(frameon=False, fontsize=6.8, labelcolor=INK2, loc="lower left")
    fig.tight_layout()
    out = os.path.join(RESULTS, "bloch.png")
    fig.savefig(out, facecolor=SURF)
    print(f"wrote {out}")


def main() -> None:
    with open(os.path.join(RESULTS, "meta.json")) as fh:
        meta = json.load(fh)
    channel = load("channel.csv")
    ramsey = load("ramsey.csv")
    tomography = load("tomography.csv")
    tworounds = load("tworounds.csv")
    lines: list[str] = []
    figure_rotations(channel, ramsey, meta, lines)
    figure_noise(channel, tomography, tworounds, meta, lines)
    figure_bloch(tomography, tworounds, meta, lines)
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(RESULTS, "summary.txt"), "w") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
