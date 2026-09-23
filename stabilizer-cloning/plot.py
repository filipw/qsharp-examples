"""Figures: learning curve collapse, copies vs n, and the cloning lower bound's mechanism."""

from __future__ import annotations

import csv
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import learner

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
INK, INK2, INK3, SURF = "#0b0b0b", "#52514e", "#8a8880", "#fcfcfb"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def style(ax):
    ax.set_facecolor(SURF)
    ax.grid(True, color="#e6e5e0", lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#d8d7d1")
    ax.tick_params(colors=INK2, labelsize=8.5)


def load():
    with open(os.path.join(RESULTS, "meta.json")) as fh:
        meta = json.load(fh)
    with open(os.path.join(RESULTS, "learning.csv")) as fh:
        rows = [{k: int(v) for k, v in r.items()} for r in csv.DictReader(fh)]
    ranks = {n: np.load(os.path.join(RESULTS, f"ranks_n{n}.npy")) for n in meta["n_values"]}
    return meta, rows, ranks


def blues(k: int):
    """Sequential ramp (one hue, light -> dark) for the ordinal variable n."""
    return [plt.cm.Blues(0.35 + 0.6 * i / max(k - 1, 1)) for i in range(k)]


def main() -> None:
    meta, rows, ranks = load()
    ns = meta["n_values"]
    extra = meta["extra"]
    lines = []

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2), dpi=200)
    fig.patch.set_facecolor(SURF)
    cols = blues(len(ns))

    # (A) learning curve: P[group learned] vs samples beyond n
    ax = axes[0]
    style(ax)
    offs = np.arange(-4, extra + 1)
    for n, c in zip(ns, cols):
        R = ranks[n]
        p = [float(np.mean(R[:, n + o - 1] == n)) if 1 <= n + o <= R.shape[1] else np.nan for o in offs]
        ax.plot(offs, p, color=c, lw=1.6, marker="o", ms=3.5, mec=SURF, mew=0.6, label=f"n = {n}")
    theory = [learner.p_learned(128 + o, 128) if 128 + o >= 1 else 0.0 for o in offs]
    ax.plot(offs, theory, color=INK3, lw=1.6, ls=(0, (5, 3)), label="Lemma 5, n → ∞")
    ax.axhline(learner.p_span_limit(), color=INK3, lw=0.9, ls=":")
    ax.text(offs[-1], learner.p_span_limit() + 0.02, "0.289 at k = n", color=INK2, fontsize=8, ha="right")
    ax.set_xlabel("Bell difference samples beyond n  (k − n)", color=INK2, fontsize=9.5)
    ax.set_ylabel("P[stabilizer group fully learned]", color=INK2, fontsize=9.5)
    ax.set_title("Learning needs n + O(1) samples, for every n", color=INK, fontsize=11, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK2, ncol=2, loc="lower right")

    # (B) copies consumed vs n, linear fit, and the paper's lower bound
    ax = axes[1]
    style(ax)
    means, medians, q10, q90 = [], [], [], []
    for n in ns:
        c = np.array([r["copies"] for r in rows if r["n"] == n and r["learned"]])
        means.append(c.mean())
        medians.append(np.median(c))
        q10.append(np.quantile(c, 0.1))
        q90.append(np.quantile(c, 0.9))
    slope, intercept = np.polyfit(ns, means, 1)
    ax.fill_between(ns, q10, q90, color=BLUE, alpha=0.12, lw=0, label="10–90% of states")
    ax.plot(ns, means, color=BLUE, lw=2, marker="o", ms=4.5, mec=SURF, mew=0.8, label="measured mean copies to learn")
    xs = np.array([0, max(ns)])
    ax.plot(xs, slope * xs + intercept, color=BLUE, lw=1, ls=(0, (5, 3)), label=f"fit: {slope:.2f} n + {intercept:.1f}")
    ax.plot(xs, 4 * (xs + 1.6067) + 1, color=INK3, lw=1.2, ls=":", label="4·E[k] + 1 = 4n + 7.4 (Lemma 5)")
    ax.plot(xs, xs / 4, color=ORANGE, lw=2, label="paper: no cloner below ⌊n/4⌋ copies (Thm 11)")
    ax.set_xlabel("qubits n", color=INK2, fontsize=9.5)
    ax.set_ylabel("copies of the state", color=INK2, fontsize=9.5)
    ax.set_title("Both thresholds are linear in n", color=INK, fontsize=11, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc="upper left")
    lines.append(f"copies-to-learn fit: {slope:.3f} n + {intercept:.2f}   (Lemma 5 predicts 4 n + {4 * 1.6067 + 1:.2f})")
    for n, m, md in zip(ns, means, medians):
        lines.append(f"  n={n:4d}: mean copies {m:7.1f}  median {md:6.0f}  predicted {4 * learner.expected_samples(n) + 1:7.1f}")

    # (C) the mechanism: distinguisher advantage vs t
    ax = axes[2]
    style(ax)
    for n, c in zip(ns, cols):
        R = ranks[n]
        adv = []
        for o in offs:
            t = n + o
            if t < 1 or t > R.shape[1]:
                adv.append(np.nan)
                continue
            r_t = R[:, t - 1]
            r_tm1 = R[:, t - 2] if t >= 2 else np.zeros(R.shape[0], dtype=int)
            adv.append(learner.advantage_from_ranks(r_t, r_tm1, n))
        ax.plot(offs, adv, color=c, lw=1.6, marker="o", ms=3.5, mec=SURF, mew=0.6)
    th = [learner.amplification_advantage(128 + o, 128) if 128 + o >= 1 else 0.0 for o in offs]
    ax.plot(offs, th, color=INK3, lw=1.6, ls=(0, (5, 3)), label="Lemma 5, n → ∞")
    ax.axhline(learner.p_span_limit() / 2, color=ORANGE, lw=1.6, label="p/2 = 0.144: the paper's cloning error bound")
    ax.set_xlabel("Bell difference samples t, relative to n  (t − n)", color=INK2, fontsize=9.5)
    ax.set_ylabel("distinguisher advantage against the best Bell-sampling cloner", color=INK2, fontsize=9)
    ax.set_title("Cloning one extra sample at t = n fails by a constant", color=INK, fontsize=11, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc="upper right")
    for n in ns:
        R = ranks[n]
        a_n = learner.advantage_from_ranks(R[:, n - 1], R[:, n - 2], n)
        lines.append(f"  n={n:4d}: advantage at t=n measured {a_n:.3f}, Lemma 5 {learner.amplification_advantage(n, n):.3f}")

    fig.suptitle("Cloning is as hard as learning for stabilizer states (arXiv:2604.15269): "
                 f"{meta['states']} random states per n, Q# Bell difference sampling on the QDK stabilizer simulator",
                 color=INK, fontsize=11.5, x=0.01, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(RESULTS, "cloning.png")
    fig.savefig(out, facecolor=SURF)
    print(f"wrote {out}")
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(RESULTS, "summary.txt"), "w") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
