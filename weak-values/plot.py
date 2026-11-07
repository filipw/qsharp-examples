"""Figures and results/summary.txt."""

from __future__ import annotations

import csv
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize, stats

import weak

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
INK, INK2, INK3, SURF = "#0b0b0b", "#52514e", "#8a8880", "#fcfcfb"
BLUE, ORANGE, AQUA, RED, PURPLE = "#2a78d6", "#eb6834", "#1baf7a", "#d63a3a", "#8e5bd6"
W_COLORS = [BLUE, AQUA, ORANGE, PURPLE, RED]
BAND = "#ecebe6"
SZ = r"$\langle p\rangle_f\,/\,2\delta p$"


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
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(fh)]


def col(rows, key, **where):
    return np.array([r[key] for r in rows if all(abs(r[k] - v) < 1e-12 for k, v in where.items())])


def eigen_band(ax, label=True):
    ax.axhspan(-0.5, 0.5, color=BAND, zorder=0, label=r"eigenvalues of $S_z$: $\pm\frac{1}{2}$" if label else None)
    ax.axhline(0, color=INK3, lw=0.6, zorder=0)


def labels(ax, title, x, y):
    ax.set_title(title, color=INK, fontsize=10.5, loc="left", fontweight="bold")
    ax.set_xlabel(x, color=INK2, fontsize=9.5)
    ax.set_ylabel(y, color=INK2, fontsize=9.5)


def save(fig, name):
    fig.tight_layout()
    out = os.path.join(RESULTS, name)
    fig.savefig(out, facecolor=SURF)
    plt.close(fig)
    print(f"wrote {out}")


# --------------------------------------------------------------- figure 1: exact

def figure_exact(exact, strength, dists, meta, lines):
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4), dpi=200)
    d7 = meta["delta_seven"]
    sigma = meta["sigma"]

    ax = axes[0]
    style(ax)
    eigen_band(ax)
    d = col(exact, "delta", w=meta["w_curves"][0])
    wv = np.where(d != 0, 1 / np.tan(np.where(d != 0, d, 1) / 2), np.nan) / 2
    ax.plot(d, np.where(d > 0, wv, np.nan), "--", color=INK, lw=1.1, label=r"weak value $\frac{1}{2}\cot(\delta/2)$", zorder=6)
    ax.plot(d, np.where(d < 0, wv, np.nan), "--", color=INK, lw=1.1, zorder=6)
    for w, c in zip(meta["w_curves"], W_COLORS):
        ax.plot(col(exact, "delta", w=w), col(exact, "a_f", w=w) / 2, color=c, lw=1.8, label=f"w = δp·σ = {w:g}", zorder=3 - 0.1 * meta["w_curves"].index(w) + 0.5)
    ax.plot([d7], [weak.amplification(d7, meta["w_curves"][0]) / 2], "*", color=INK, ms=12, zorder=5,
            label=f"δ = 2·atan(1/14) = {d7:.4f}: reads 7")
    ax.set_xlim(-np.pi / 2, np.pi / 2)
    ax.set_ylim(-12, 12)
    ax.legend(frameon=False, fontsize=7.2, labelcolor=INK2, loc="lower right")
    labels(ax, "Post-selected pointer reading against the angle", "post-selection angle δ, |f⟩ = Ry(δ)|−x⟩  [rad]", f"$S_z$ read by the pointer, {SZ}")

    ax = axes[1]
    style(ax)
    ax.axhspan(1e-3, 0.5, color=BAND, zorder=0, label=r"inside the eigenvalue range, $\leq\frac{1}{2}$")
    w = col(strength, "w")
    ax.axhline(meta["sz_target"], color=INK, ls="--", lw=1.3, label="weak value at δ = 2·atan(1/14): 7")
    ax.plot(w, col(strength, "a_seven") / 2, color=BLUE, lw=2.2, label="at δ = 2·atan(1/14), Q# state")
    ax.plot(w, col(strength, "a_seven_closed") / 2, ":", color=INK, lw=1.0, label=r"$\frac{1}{2}\sin\delta / (1-\kappa\cos\delta)$, $\kappa = e^{-2w^2}$")
    ax.plot(w, col(strength, "a_max") / 2, color=ORANGE, lw=1.8, label=r"largest reading over δ, $\frac{1}{2}(1-\kappa^2)^{-1/2}$")
    for x, lab in ((meta["seven_w"], "1M-shot run"), (meta["fig2_w"], "Fig. 2 sweep")):
        ax.axvline(x, color=INK3, lw=0.8, ls=":")
        ax.annotate(lab, (x, 0.12), rotation=90, fontsize=7, color=INK2, xytext=(-9, 0), textcoords="offset points", ha="center", va="center")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(0.05, 300)
    ax.legend(frameon=False, fontsize=7.2, labelcolor=INK2, loc="upper right")
    labels(ax, "Weak enough? The reading against the kick strength", "w = δp·σ (kick against the coherence momentum 1/σ)", f"$S_z$ read by the pointer, {SZ}")

    ax = axes[2]
    style(ax)
    wd = meta["w_dist"]
    unit = 1 / (2 * sigma)
    p = col(dists, "p", w=wd) / unit
    ax.plot(p, col(dists, "rest", w=wd), color=INK3, lw=1.5, label="pointer at rest")
    d_star, a_max = weak.peak(wd)
    ax.plot(p, col(dists, "seven", w=wd), color=BLUE, lw=1.8, marker="o", ms=3, label=f"post-selected, δ = {d7:.3f} (reads {weak.amplification(d7, wd) / 2:.2f})")
    ax.plot(p, col(dists, "peak", w=wd), color=ORANGE, lw=1.8, marker="o", ms=3, label=f"post-selected, δ* = {d_star:.3f} (reads {a_max / 2:.1f})")
    kick = wd / sigma / unit
    ax.axvspan(-kick, kick, color=RED, alpha=0.35, lw=0, label=f"±δp: the kicks the two eigenstates get (±{kick:.2f})")
    ax.set_xlim(-4, 4)
    ax.set_ylim(0, 0.255)
    ax.legend(fontsize=7.2, labelcolor=INK2, loc="upper left", facecolor=SURF, edgecolor="none", framealpha=0.95)
    labels(ax, f"Momentum distributions at w = {wd:g}", "pointer momentum p·2σ (in units of its spread 1/2σ at rest)", "probability per read-out bin")

    fig.suptitle("A weak Stern-Gerlach measurement on qubits (arXiv:2609.37822): |+x⟩ in, kick exp(iδp σ_z z), post-select Ry(δ)|−x⟩, read the momentum",
                 color=INK, fontsize=11.5, x=0.01, ha="left", y=0.995)
    save(fig, "amplification.png")

    diff = max(abs(r["a_f"] - r["a_closed"]) for r in exact)
    rel = max(abs(r["a_f"] - r["a_closed"]) / max(1.0, abs(r["a_closed"])) for r in exact)
    pdiff = max(abs(r["p_f"] - r["p_closed"]) for r in exact)
    ndiff = max(abs(r["a_n"] - r["a_n_closed"]) for r in exact)
    lines.append(f"panel A exact: {len(exact)} (w, delta) points from Q# dumps; max |A_f - closed form| = {diff:.1e} (relative {rel:.1e}), "
                 f"|A_n - closed form| = {ndiff:.1e}, |P(f) - closed form| = {pdiff:.1e}")
    for w in meta["w_curves"]:
        ds, am = weak.peak(w)
        lines.append(f"  w = {w:<6g} peak at delta* = {ds:.4f}: largest reading S_z = {am / 2:8.3f}, post-selection probability {(1 - weak.kappa(w) ** 2) / 2:.2e}; "
                     f"at delta = {d7:.4f}: S_z = {weak.amplification(d7, w) / 2:.3f}")
    w1, w10 = (optimize.brentq(lambda w: weak.amplification(d7, w) / 2 - frac * meta["sz_target"], 1e-4, 1.0) for frac in (0.99, 0.9))
    lines.append(f"panel B: at delta = 2 atan(1/14) the reading is within 1% of 7 for w <= {w1:.4f} "
                 f"and within 10% for w <= {w10:.4f} (closed form, solved for w); max |Q# - closed form| = "
                 f"{max(abs(r['a_seven'] - r['a_seven_closed']) for r in strength):.1e}, at the peak {max(abs(r['a_max'] - r['a_max_closed']) for r in strength):.1e}")
    shift = [r["a_max"] * (r["w"] / sigma) / unit for r in strength]
    lines.append(f"  largest shift in units of 1/(2 sigma) (the Gaussian pointer's momentum width): {shift[0]:.4f} at w = {strength[0]['w']:g}, {shift[-1]:.4f} at w = {strength[-1]['w']:g} "
                 f"(tends to 1 as w -> 0)")
    for name in ("seven", "peak"):
        dist = col(dists, name, w=wd)
        m = float((col(dists, "p", w=wd) * dist).sum())
        lines.append(f"panel C w = {wd}: {name} distribution mean {m / unit:.4f} pointer widths = {m / (wd / sigma):.3f} dp")


# --------------------------------------------------------------- figure 2: shots

def figure_shots(fig2, seven, counts, dists, meta, lines):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4), dpi=200)
    sigma = meta["sigma"]
    w2 = meta["fig2_w"]

    ax = axes[0]
    style(ax)
    eigen_band(ax)
    dd = np.linspace(-np.pi / 2, np.pi / 2, 2001)
    wv = 0.5 / np.tan(np.where(np.abs(dd) > 0.004, dd, 1.0) / 2)
    ax.plot(dd, np.where(dd > 0.004, wv, np.nan), "--", color=INK, lw=1.3, label=r"weak value $\frac{1}{2}\cot(\delta/2)$")
    ax.plot(dd, np.where(dd < -0.004, wv, np.nan), "--", color=INK, lw=1.3)
    ax.plot(dd, weak.amplification(dd, w2) / 2, color=AQUA, lw=1.8, label=f"exact at w = {w2:g}")
    ax.plot(dd, weak.amplification_orthogonal(dd, w2) / 2, color=INK3, lw=1.0, label="exact, discarded |n⟩")
    d = col(fig2, "delta")
    ax.errorbar(d, col(fig2, "a_f") / 2, yerr=2 * col(fig2, "se_f") / 2, fmt="o", ms=5, color=BLUE, ecolor=BLUE, elinewidth=1.1, capsize=2.5,
                label="shots, post-selected |f⟩ (±2σ)")
    ax.errorbar(d, col(fig2, "a_n") / 2, yerr=2 * col(fig2, "se_n") / 2, fmt="s", ms=3.5, mfc=SURF, color=INK3, ecolor=INK3, elinewidth=0.8, capsize=1.5,
                label="shots, discarded |n⟩")
    ax.set_xlim(-np.pi / 2 - 0.05, np.pi / 2 + 0.05)
    ax.set_ylim(-7.5, 7.5)
    ax.legend(fontsize=6.8, labelcolor=INK2, loc="lower right", facecolor=SURF, edgecolor="none", framealpha=0.95)
    labels(ax, f"Fig. 2 from shots: w = {w2:g}, {int(meta['fig2_shots']) // 1000}k shots per angle", "post-selection angle δ  [rad]", f"$S_z$ read by the pointer, {SZ}")

    ax = axes[1]
    style(ax)
    r = seven[0]
    unit = 1 / (2 * sigma)
    dp = r["w"] / sigma
    p = col(counts, "p")
    nf = col(counts, "f")
    step = p[1] - p[0]
    assert abs(meta["w_dist"] - r["w"]) < 1e-12 and np.allclose(col(dists, "p"), p)
    exact_f = col(dists, "seven")
    rest = col(dists, "rest")
    ax.bar(p / unit, nf / nf.sum(), width=step / unit * 0.9, color=BLUE, alpha=0.55, label=f"shots: {int(r['n_f'])} post-selected of {int(r['shots']):,}")
    ax.plot(p / unit, exact_f, "o-", color=BLUE, ms=3, lw=1.2, label="exact post-selected distribution (Q# state)")
    ax.plot(p / unit, rest, color=INK3, lw=1.2, label="pointer at rest")
    mean = r["a_f"] * dp / unit
    ax.axvline(mean, color=BLUE, lw=1.4, label=f"shot mean: S_z = {r['sz']:.2f} ± {r['se_sz']:.2f} (±1σ; exact {r['exact_a_f'] / 2:.2f}, weak value 7)")
    ax.axvspan(-dp / unit, dp / unit, color=RED, alpha=0.4, lw=0, label=f"±δp, where S_z = ±½ would put it (±{dp / unit:.2f})")
    ax.set_xlim(-3.5, 3.5)
    ax.set_ylim(0, 0.235)
    ax.legend(fontsize=7.0, labelcolor=INK2, loc="upper left", ncol=2, facecolor=SURF, edgecolor="none", framealpha=0.95)
    labels(ax, f"Getting 7: δ = 2·atan(1/14), w = {r['w']:g}", "pointer momentum p·2σ (in units of its spread 1/2σ at rest)", "fraction of post-selected shots")

    fig.suptitle("Shots (arXiv:2609.37822): every shot records the spin outcome and the pointer momentum; the post-selected ones are averaged",
                 color=INK, fontsize=11.5, x=0.01, ha="left", y=0.995)
    save(fig, "shots.png")

    z_f = [(x["a_f"] - x["exact_a_f"]) / x["se_f"] for x in fig2]
    z_n = [(x["a_n"] - x["exact_a_n"]) / x["se_n"] for x in fig2]
    z_r = [(x["rate"] - x["exact_rate"]) / x["se_rate"] for x in fig2]
    lines.append(f"Fig. 2 sweep, w = {w2}, {int(meta['fig2_shots'])} shots per angle: worst |z| of the post-selected mean {max(map(abs, z_f)):.2f}, "
                 f"discarded mean {max(map(abs, z_n)):.2f}, post-selection rate {max(map(abs, z_r)):.2f} over {len(fig2)} angles")
    for x in fig2:
        lines.append(f"  delta = {x['delta']:+.3f}: {int(x['n_f']):6d} post-selected (rate {x['rate']:.4f}, exact {x['exact_rate']:.4f}), "
                     f"S_z = {x['a_f'] / 2:+7.3f} +- {x['se_f'] / 2:.3f} (exact {x['exact_a_f'] / 2:+7.3f}, weak value "
                     f"{'   -  ' if np.isnan(x['weak_value']) else format(x['weak_value'] / 2, '+7.3f')}); discarded S_z = {x['a_n'] / 2:+.3f} +- {x['se_n'] / 2:.3f}")
    lines.append(f"the 7: delta = {r['delta']:.5f}, w = {r['w']}, {int(r['shots'])} shots, {int(r['n_f'])} post-selected (rate {r['rate']:.5f} +- {r['se_rate']:.5f}, "
                 f"exact {r['exact_rate']:.5f}); S_z = {r['sz']:.3f} +- {r['se_sz']:.3f}, exact {r['exact_a_f'] / 2:.3f}, weak value {r['weak_value'] / 2:.3f}; "
                 f"{(r['a_f'] - r['exact_a_f']) / r['se_f']:+.2f} sigma from exact")
    lines.append(f"  discarded outcome: S_z = {r['a_n'] / 2:+.4f} +- {r['se_n'] / 2:.4f} (exact {r['exact_a_n'] / 2:+.4f})")
    expect = exact_f * nf.sum()
    keep = expect >= 5
    obs = np.append(nf[keep], nf[~keep].sum())
    exp = np.append(expect[keep], expect[~keep].sum())
    chi = stats.chisquare(obs, exp * obs.sum() / exp.sum())
    lines.append(f"  post-selected histogram against the exact distribution: chi^2 = {chi.statistic:.1f} over {len(obs) - 1} degrees of freedom, p = {chi.pvalue:.3f}")


def main() -> None:
    with open(os.path.join(RESULTS, "meta.json")) as fh:
        meta = json.load(fh)
    lines: list[str] = []
    figure_exact(load("exact.csv"), load("strength.csv"), load("distributions.csv"), meta, lines)
    figure_shots(load("fig2.csv"), load("seven.csv"), load("seven_counts.csv"), load("distributions.csv"), meta, lines)
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(RESULTS, "summary.txt"), "w") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
