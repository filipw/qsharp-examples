"""Figures and results/summary.txt."""

from __future__ import annotations

import csv
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import magic

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
INK, INK2, INK3, SURF = "#0b0b0b", "#52514e", "#8a8880", "#fcfcfb"
BLUE, ORANGE, AQUA, RED, PINK = "#2a78d6", "#eb6834", "#1baf7a", "#d63a3a", "#e58fb3"


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


def blues(k: int):
    return [plt.cm.Blues(0.4 + 0.55 * i / max(k - 1, 1)) for i in range(k)]


def figure_rebirth(curves, thresholds, witness, meta, lines):
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2), dpi=200)
    fig.patch.set_facecolor(SURF)

    # (A) n = 2, alpha = 0.4: the three thresholds, exact curves and the shot witness
    ax = axes[0]
    style(ax)
    n, alpha = 2, 0.4
    g = col(curves, "gamma", n=n, alpha=alpha)
    R = col(curves, "R", n=n, alpha=alpha)
    neg = col(curves, "neg", n=n, alpha=alpha)
    th = [r for r in thresholds if r["n"] == n and abs(r["alpha"] - alpha) < 1e-9][0]
    ax.axvspan(th["g_minus"], th["g_plus"], color="#e9e8e3", lw=0, label="stabilizer window (magic = 0)")
    ax.plot(g, R, color=BLUE, lw=2.2, label="magic: robustness ℛ − 1, exact Q# state")
    ax.plot(g, neg, color=RED, lw=2.2, label="entanglement: negativity (exact)")
    wg = col(witness, "gamma", n=n, alpha=alpha)
    wv = col(witness, "witness", n=n, alpha=alpha)
    ws = col(witness, "se", n=n, alpha=alpha)
    ax.errorbar(wg, wv, yerr=2 * ws, fmt="o", ms=4, mfc=SURF, mec=BLUE, ecolor=BLUE, elinewidth=1, capsize=2,
                label=f"witness Tr(Wρ) − 1 from 2 × {meta['witness_shots']} shots (±2σ)")
    ax.axhline(0, color=INK3, lw=0.8)
    for x, y, lab, c in ((th["g_minus"], 0.34, "γ₋: magic dies", BLUE), (th["g_e"], 0.30, "γₑ: entanglement dies", RED), (th["g_plus"], 0.26, "γ₊: magic reborn", BLUE)):
        ax.axvline(x, color=c, lw=0.9, ls=":")
        ax.text(x + 0.008, y, f"{lab}   γ = {x:.3f}", color=c, fontsize=7.5, va="center", ha="left")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.17, 0.5)
    ax.set_xlabel("damping strength γ", color=INK2, fontsize=9.5)
    ax.set_ylabel("resource", color=INK2, fontsize=9.5)
    ax.set_title(f"n = 2, α = 0.4: magic dies, entanglement dies, magic returns", color=INK, fontsize=11, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc="upper right")
    lines.append(f"panel A (n=2, alpha=0.4): gamma_- = {th['g_minus']:.4f}, gamma_e = {th['g_e']:.4f}, gamma_+ = {th['g_plus']:.4f}, "
                 f"gamma_e + gamma_+ = {th['g_sum']:.6f}")
    for wn, wa in meta["witness"]:
        cert = [(int(r["certified"]), r["exact"] > 1e-12) for r in witness if r["n"] == wn and abs(r["alpha"] - wa) < 1e-9]
        lines.append(f"  witness from 2 x {meta['witness_shots']} shots, n={wn} alpha={wa}: certified magic (> 2 sigma) at "
                     f"{sum(c for c, _ in cert)} of {sum(e for _, e in cert)} magic points, "
                     f"{sum(1 for c, e in cert if c and not e)} false certifications at {sum(1 for _, e in cert if not e)} stabilizer points")

    # (B) thresholds against n
    ax = axes[1]
    style(ax)
    alpha = 0.3
    rows = sorted([r for r in thresholds if abs(r["alpha"] - alpha) < 1e-9], key=lambda r: r["n"])
    ns = np.array([r["n"] for r in rows])
    ax.fill_between(ns, [r["g_minus"] for r in rows], [r["g_plus"] for r in rows], color="#e9e8e3", lw=0, label="stabilizer window")
    nn = np.linspace(2, ns.max(), 200)
    beta = np.sqrt(1 - alpha**2)
    rr = alpha / beta
    ax.plot(nn, 1 - rr ** (2 / nn), color=BLUE, lw=1.2, ls=(0, (5, 3)), label="γ₊ = 1 − r^{2/n}  (Eq. 4)")
    ax.plot(nn, rr ** (2 / nn), color=RED, lw=1.2, ls=(0, (5, 3)), label="γₑ = r^{2/n}  (Prop. 2)")
    ax.plot(nn, [magic.thresholds(alpha, float(k))[0] for k in nn], color=BLUE, lw=1.2, ls=":", label="γ₋: root of P₀ = c")
    ax.plot(ns, [r["env_e"] for r in rows], "x", color=ORANGE, ms=9, mew=1.5, label="environment becomes entangled (measured)")
    ax.plot(ns, [r["g_plus"] for r in rows], "^", color=BLUE, ms=6, mec=SURF, mew=0.6, label="γ₊ measured (bisection on Q# states)")
    ax.plot(ns, [r["g_minus"] for r in rows], "v", color=BLUE, ms=6, mec=SURF, mew=0.6, label="γ₋ measured")
    ax.plot(ns, [r["g_e"] for r in rows], "o", color=RED, ms=6, mec=SURF, mew=0.6, label="γₑ measured")
    ax.plot(ns, [r["g_sum"] for r in rows], "s", color=INK2, ms=5, mfc=SURF, label="γₑ + γ₊ measured")
    ax.axhline(1, color=INK3, lw=0.8)
    ax.set_xlabel("qubits n", color=INK2, fontsize=9.5)
    ax.set_ylabel("damping strength γ", color=INK2, fontsize=9.5)
    ax.set_ylim(0, 1.5)
    ax.set_title(f"α = {alpha}: γₑ + γ₊ = 1 for every n", color=INK, fontsize=11, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.8, labelcolor=INK2, loc="upper left", ncol=2)
    lines.append("panel B thresholds (measured vs formula):")
    for a_ in meta["alphas"]:
        for r in sorted([r for r in thresholds if abs(r["alpha"] - a_) < 1e-9], key=lambda r: r["n"]):
            lines.append(f"  alpha={a_} n={int(r['n'])} ({r['regime']:>3}): gamma_- {r['g_minus']:.6f} ({r['th_minus']:.6f})  "
                         f"gamma_+ {r['g_plus']:.6f} ({r['th_plus']:.6f})  gamma_e {r['g_e']:.6f} ({r['th_e']:.6f})  "
                         f"sum {r['g_sum']:.6f}  env entangled from {r['env_e']:.6f}, env reborn at {r['env_plus']:.6f}")
    worst = max(max(abs(r["g_minus"] - r["th_minus"]), abs(r["g_plus"] - r["th_plus"]), abs(r["g_e"] - r["th_e"]), abs(r["g_sum"] - 1),
                    abs(r["env_e"] - r["g_plus"]), abs(r["env_plus"] - r["g_e"])) for r in thresholds)
    lines.append(f"  largest deviation from the formulas, from gamma_e + gamma_+ = 1, and from the environment mirror: {worst:.1e}")

    # (C) the Stinespring mirror
    ax = axes[2]
    style(ax)
    n, alpha = 3, 0.4
    g = col(curves, "gamma", n=n, alpha=alpha)
    R, Re = col(curves, "R", n=n, alpha=alpha), col(curves, "R_env", n=n, alpha=alpha)
    ng, nge = col(curves, "neg", n=n, alpha=alpha), col(curves, "neg_env", n=n, alpha=alpha)
    ax.plot(g, R, color=BLUE, lw=2.2, label="system: magic ℛ − 1")
    ax.plot(g, ng, color=RED, lw=2.2, label="system: negativity, qubit 1 | rest")
    ax.plot(g, Re, color=ORANGE, lw=2.2, label="environment: magic ℛ − 1")
    ax.plot(g, nge, color=PINK, lw=2.2, label="environment: negativity, qubit 1 | rest")
    ax.plot(g[::16], R[::-1][::16], "o", ms=4, mfc=SURF, mec=ORANGE, mew=1, label="system curves reflected, γ → 1 − γ")
    ax.plot(g[::16], ng[::-1][::16], "o", ms=4, mfc=SURF, mec=PINK, mew=1)
    ax.axvline(0.5, color=INK3, lw=0.8, ls=":")
    ax.set_xlim(0, 1)
    ax.set_xlabel("damping strength γ", color=INK2, fontsize=9.5)
    ax.set_ylabel("resource", color=INK2, fontsize=9.5)
    ax.set_title(f"n = {n}, α = {alpha}: environment = system at 1 − γ", color=INK, fontsize=11, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc="upper center")
    mirror = max(float(np.abs(Re - R[::-1]).max()), float(np.abs(nge - ng[::-1]).max()))
    pur = float(np.abs(col(curves, "purity", n=n, alpha=alpha) - col(curves, "purity_env", n=n, alpha=alpha)).max())
    lines.append(f"panel C mirror (n={n}, alpha={alpha}): max |R_env(gamma) - R_sys(1-gamma)| and negativity likewise = {mirror:.1e}; "
                 f"max |purity_sys - purity_env| = {pur:.1e}")

    fig.suptitle("Sudden death of entanglement, rebirth of magic (arXiv:2605.22603): amplitude damping of α|0ⁿ⟩ + β|1ⁿ⟩ "
                 "as a Q# Stinespring dilation, resources evaluated on the exact state",
                 color=INK, fontsize=11.5, x=0.01, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(RESULTS, "rebirth.png")
    fig.savefig(out, facecolor=SURF)
    print(f"wrote {out}")


def figure_extraction(curves, thresholds, extraction, stab, channel, chth, meta, lines):
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2), dpi=200)
    fig.patch.set_facecolor(SURF)

    # (A) the decoded qubit in the Bloch XZ plane
    ax = axes[0]
    style(ax)
    ax.set_aspect("equal")
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(t), np.sin(t), color=INK3, lw=0.8)
    for thr, c, lab in ((np.sqrt(2) * (1 - 2 * 0.1411), "#fbe3ec", "|H⟩-type distillable, |x|+|z| > 1.015"), (3 / np.sqrt(7), "#f3c6d8", "|T⟩-type distillable, |x|+|z| > 3/√7")):
        xs = np.linspace(0, 1, 300)
        zs_lo = np.clip(thr - xs, 0, None)
        zs_hi = np.sqrt(np.clip(1 - xs**2, 0, None))
        for sx in (1, -1):
            for sz in (1, -1):
                ax.fill_between(sx * xs, sz * zs_lo, sz * zs_hi, where=zs_hi > zs_lo, color=c, lw=0, label=lab if (sx, sz) == (1, 1) else None)
    ax.plot([1, 0, -1, 0, 1], [0, 1, 0, -1, 0], color=INK2, lw=1.2, label="stabilizer octahedron |x| + |z| ≤ 1")
    alpha = 0.2
    ns = [2, 4, 6]
    for n, c in zip(ns, blues(len(ns))):
        x = col(curves, "x", n=n, alpha=alpha)
        z = col(curves, "z", n=n, alpha=alpha)
        ax.plot(x, z, color=c, lw=1.8, label=f"n = {n}, exact")
        th = [r for r in thresholds if r["n"] == n and abs(r["alpha"] - alpha) < 1e-9][0]
        g = col(curves, "gamma", n=n, alpha=alpha)
        for gg in (th["g_minus"], th["g_plus"]):
            i = int(np.argmin(np.abs(g - gg)))
            ax.plot(x[i], z[i], "o", ms=4, color=c, mec=SURF)
        ex = col(extraction, "x", n=n, alpha=alpha)
        ez = col(extraction, "z", n=n, alpha=alpha)
        ax.errorbar(ex, ez, xerr=2 * col(extraction, "se_x", n=n, alpha=alpha), yerr=2 * col(extraction, "se_z", n=n, alpha=alpha),
                    fmt="s", ms=3.5, mfc=SURF, mec=c, ecolor=c, elinewidth=0.9, capsize=1.5, label=f"n = {n}, post-selected shots (±2σ)")
    ax.annotate("γ = 0", (col(curves, "x", n=2, alpha=alpha)[0], col(curves, "z", n=2, alpha=alpha)[0]), textcoords="offset points", xytext=(4, -10), fontsize=7.5, color=INK2)
    ax.annotate("γ = 1", (0, 1), textcoords="offset points", xytext=(4, -10), fontsize=7.5, color=INK2)
    ax.set_xlim(-0.6, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlabel("⟨X⟩", color=INK2, fontsize=9.5)
    ax.set_ylabel("⟨Z⟩", color=INK2, fontsize=9.5)
    ax.set_title(f"α = {alpha}: the parity-syndrome-extracted qubit", color=INK, fontsize=11, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.8, labelcolor=INK2, loc="center left")
    lines.append("panel A extraction (alpha=0.2), post-selected shots vs exact:")
    for r in extraction:
        lines.append(f"  n={int(r['n'])} gamma={r['gamma']:.4f}: p_success {r['p_success']:.4f} +- {r['se_p']:.4f} (exact {r['exact_p']:.4f}), "
                     f"x {r['x']:.3f} +- {r['se_x']:.3f} ({r['exact_x']:.3f}), z {r['z']:.3f} +- {r['se_z']:.3f} ({r['exact_z']:.3f}), "
                     f"quality |x|+|z| {r['quality']:.3f} ({r['exact_quality']:.3f}), yield {r['yield_']:.4f} = joint R-1 {r['R_joint']:.4f}")
    for n in ns:
        gp = [r for r in thresholds if r["n"] == n and abs(r["alpha"] - alpha) < 1e-9][0]["g_plus"]
        rows = [r for r in extraction if r["n"] == n and r["gamma"] > gp]
        lines.append(f"  n={n}, reborn branch (gamma > {gp:.3f}): best decoded quality |x|+|z| = {max(r['exact_quality'] for r in rows):.3f} "
                     f"(H-type distillable above 1.015, T-type above 1.134); best yield at the sampled points {max(r['yield_'] for r in rows):.4f}, "
                     f"the paper's maximum alpha^2/2 = {alpha**2 / 2:.4f}")

    # (B) stabilizer inputs: generators and insulators
    ax = axes[1]
    style(ax)
    cols = {"Phi+": BLUE, "Psi+": AQUA, "GHZ3": RED, "Phi+ x 0": INK3, "00+": ORANGE}
    labels = {"Phi+": "|Φ⁺⟩ = (|00⟩+|11⟩)/√2", "Psi+": "|Ψ⁺⟩ = (|01⟩+|10⟩)/√2", "GHZ3": "GHZ₃", "Phi+ x 0": "|Φ⁺⟩|0⟩", "00+": "|00+⟩"}
    gg = np.linspace(0, 1, 200)
    closed = {"Phi+": gg * (1 - gg), "Psi+": 0 * gg, "GHZ3": (1 - gg) ** 1.5 - (1 - gg) ** 3, "Phi+ x 0": gg * (1 - gg), "00+": np.sqrt(1 - gg) - (1 - gg)}
    for name in ("Phi+", "Psi+", "GHZ3", "00+"):
        g = col(stab, "gamma", state=name)
        R = col(stab, "R", state=name)
        ax.plot(gg, closed[name], color=cols[name], lw=1.2, ls=(0, (5, 3)))
        ax.plot(g, R, "o", ms=3.5, color=cols[name], mec=SURF, mew=0.5, label=f"{labels[name]}: robustness ℛ − 1 by LP")
    for name, ls in (("Phi+", "-"), ("Psi+", "-")):
        ax.plot(col(stab, "gamma", state=name), col(stab, "concurrence", state=name), color=cols[name], lw=1, ls=ls, alpha=0.45,
                label=f"{labels[name]}: concurrence")
    ax.set_xlim(0, 1)
    ax.set_xlabel("damping strength γ", color=INK2, fontsize=9.5)
    ax.set_ylabel("resource", color=INK2, fontsize=9.5)
    ax.set_title("Stabilizer inputs: magic-generators vs a magic-insulator", color=INK, fontsize=11, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7, labelcolor=INK2, loc="upper right")
    lines.append("panel B stabilizer inputs, LP robustness vs closed profiles (Eqs. 30, 33):")
    for name in cols:
        rows = [r for r in stab if r["state"] == name]
        lines.append(f"  {name:9s}: max |LP - closed| = {max(abs(r['R'] - r['closed']) for r in rows):.1e}, max R-1 = {max(r['R'] for r in rows):.4f}")

    # (C) channel conditions
    ax = axes[2]
    style(ax)
    n, alpha = meta["channel"]
    cs = [BLUE, AQUA, ORANGE, RED]
    for a, c in zip(meta["exponents"], cs):
        kt = col(channel, "kt", a=float(a))
        R = col(channel, "R", a=float(a))
        lab = "T₂/T₁ = " + {0.5: "2", 0.75: "4/3", 1.0: "1", 1.25: "4/5"}.get(a, f"{1 / a:.2f}") + ("  (pure amplitude damping)" if a == 0.5 else "")
        ax.plot(kt, R, color=c, lw=2, label=lab)
        th = [r for r in chth if abs(r["a"] - a) < 1e-9][0]
        if not np.isnan(th["kt_plus"]):
            ax.plot(th["kt_plus"], 0, "^", color=c, ms=7, mec=SURF)
    ax.axhline(0, color=INK3, lw=0.8)
    ax.set_xlim(0, meta_kt_max(channel))
    ax.set_xlabel("time κt  (γ = 1 − e^{−κt})", color=INK2, fontsize=9.5)
    ax.set_ylabel("magic ℛ − 1", color=INK2, fontsize=9.5)
    ax.set_title(f"n = {n}, α = {alpha}, damping + dephasing: rebirth iff T₂ > T₁", color=INK, fontsize=11, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc="upper right")
    lines.append(f"panel C channel conditions (n={n}, alpha={alpha}):")
    for r in chth:
        if np.isnan(r["kt_plus"]):
            lines.append(f"  a={r['a']:.2f} (T2/T1 = {r['T2_over_T1']:.2f}): no rebirth; max (c - P_n) after death = {r['max_margin_after_death']:.2e}")
        else:
            lines.append(f"  a={r['a']:.2f} (T2/T1 = {r['T2_over_T1']:.2f}): reborn at kappa t = {r['kt_plus']:.5f}, theory {r['kt_plus_theory']:.5f} "
                         f"(q_+ = r^(1/((1-a)n)) = {r['q_plus_theory']:.5f})")

    fig.suptitle("Extracting the reborn magic, which stabilizer inputs generate magic, and which channels revive it (arXiv:2605.22603)",
                 color=INK, fontsize=11.5, x=0.01, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(RESULTS, "extraction.png")
    fig.savefig(out, facecolor=SURF)
    print(f"wrote {out}")


def meta_kt_max(channel):
    return max(r["kt"] for r in channel)


def main() -> None:
    with open(os.path.join(RESULTS, "meta.json")) as fh:
        meta = json.load(fh)
    curves = load("curves.csv")
    thresholds = load("thresholds.csv")
    witness = load("witness.csv")
    extraction = load("extraction.csv")
    stab = load("stabilizer_inputs.csv")
    channel = load("channel.csv")
    chth = load("channel_thresholds.csv")
    lines: list[str] = []
    figure_rebirth(curves, thresholds, witness, meta, lines)
    figure_extraction(curves, thresholds, extraction, stab, channel, chth, meta, lines)
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(RESULTS, "summary.txt"), "w") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
