"""Figures: measured cost frontiers against the paper's objective, per alpha."""

from __future__ import annotations

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import predict as P

RESULTS = P.RESULTS
# Validated categorical slots (all-pairs safe); gray for the theory overlay.
C_COMP, C_QD, C_TR, C_TH = "#2a78d6", "#eb6834", "#1baf7a", "#6b6a66"
INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
EPS_GRID = np.logspace(-1, -7, 25)


def style(ax):
    ax.set_facecolor(SURF)
    ax.grid(True, which="major", color="#e6e5e0", lw=0.8)
    ax.grid(True, which="minor", color="#f0efea", lw=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#d8d7d1")
    ax.tick_params(colors=INK2, labelsize=8.5)


def analyse(alpha: float, t: float) -> dict:
    rows = P.load_sweep(alpha)
    a = P.load_coeffs(alpha)
    L = len(a)
    fam = lambda f: [r for r in rows if r["family"] == f]
    tro, qd = P.pareto(fam("trotter")), P.pareto(fam("qdrift"))
    comp_rows = fam("composite")
    env = P.pareto(rows)                              # min over everything, endpoints included
    out = dict(alpha=alpha, L=L, a=a, rows=rows, trotter=tro, qdrift=qd, composite=comp_rows, envelope=env)

    # Per-eps: measured cost of each family, best split, predicted K*. Only where both
    # the Trotter and the composite grids reach: past the end of either, "best split"
    # would report which grid happened to be longer, not which method is cheaper.
    reach = max(min(r["trace_distance"] for r in tro), min(r["trace_distance"] for r in comp_rows))
    eps, g_env, g_tr, g_qd, k_meas, k_pred, g_pred, k_fin, g_fin = [], [], [], [], [], [], [], [], []
    for e in EPS_GRID:
        ge = P.cost_at(env, e)
        if ge is None or e < reach:
            continue
        eps.append(e); g_env.append(ge)
        g_tr.append(P.cost_at(tro, e)); g_qd.append(P.cost_at(qd, e))
        k_meas.append(P.best_split_at(rows, e)[0]); k_pred.append(P.k_star(a, t, e)); g_pred.append(P.predicted_cost(a, t, e))
        k_fin.append(P.k_star_finite(a, t, e, 2)); g_fin.append(P.predicted_cost_finite(a, t, e, 2))
    out.update(eps=np.array(eps), g_env=np.array(g_env, float), g_tr=np.array([x or np.nan for x in g_tr], float),
               g_qd=np.array([x or np.nan for x in g_qd], float), k_meas=np.array(k_meas, float),
               k_pred=np.array(k_pred, float), g_pred=np.array(g_pred, float),
               k_fin=np.array(k_fin, float), g_fin=np.array(g_fin, float))

    # Exponent fit inside the window where the split is genuinely interior (2 <= K* <= L/4).
    win = (out["k_pred"] >= 2) & (out["k_pred"] <= L / 4)
    out["window"] = win
    out["slope"] = P.fit_exponent(out["eps"][win], out["g_env"][win]) if win.sum() >= 3 else np.nan
    out["slope_pred_curve"] = P.fit_exponent(out["eps"][win], out["g_pred"][win]) if win.sum() >= 3 else np.nan
    out["slope_fin_curve"] = P.fit_exponent(out["eps"][win], out["g_fin"][win]) if win.sum() >= 3 else np.nan
    # How fast the cut moves: K ~ eps^(-q). Measured on interior best-K points only.
    kw = win & (out["k_meas"] >= 1) & (out["k_meas"] < L)
    out["kslope_meas"] = P.fit_exponent(out["eps"][kw], out["k_meas"][kw]) if kw.sum() >= 3 else np.nan
    out["kslope_pred"] = P.fit_exponent(out["eps"][win], out["k_pred"][win]) if win.sum() >= 3 else np.nan
    out["kslope_fin"] = P.fit_exponent(out["eps"][win], np.maximum(out["k_fin"][win], 1)) if win.sum() >= 3 else np.nan
    return out


def fig_frontiers(res: list[dict], t: float, meta: dict) -> None:
    fig, axes = plt.subplots(1, len(res), figsize=(5.4 * len(res), 5.6), dpi=200, sharey=True)
    fig.patch.set_facecolor(SURF)
    for ax, R in zip(np.atleast_1d(axes), res):
        style(ax)
        a = R["alpha"]
        ax.scatter([r["twoq"] for r in R["composite"]], [r["trace_distance"] for r in R["composite"]],
                   s=7, color=C_COMP, alpha=0.18, linewidths=0, zorder=2)
        for front, col, lab in ((R["trotter"], C_TR, "Trotter, all L terms (K = L)"),
                                (R["qdrift"], C_QD, "qDRIFT (K = 0)"),
                                (R["envelope"], C_COMP, "best split K (envelope)")):
            ax.plot([p["twoq"] for p in front], [p["trace_distance"] for p in front],
                    color=col, lw=2, marker="o", ms=4, mec=SURF, mew=0.8, label=lab, zorder=4)
        # Paper's objective, scaled by one constant chosen at the middle of the window.
        w = R["window"]
        if w.sum() >= 2:
            mid = np.nonzero(w)[0][w.sum() // 2]
            scale = R["g_env"][mid] / R["g_pred"][mid]
            ax.plot(R["g_pred"] * scale, R["eps"], color=C_TH, lw=1.6, ls=(0, (5, 3)), zorder=3,
                    label="paper's objective, 1 fitted constant")
            scale_f = R["g_env"][mid] / R["g_fin"][mid]
            ax.plot(R["g_fin"] * scale_f, R["eps"], color=C_TH, lw=1.4, ls=(0, (1.5, 2.5)), zorder=3,
                    label="Fact 3 with order 2 explicit, 1 fitted constant")
            ax.axhspan(R["eps"][w].min(), R["eps"][w].max(), color=C_COMP, alpha=0.05, lw=0)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(f"α = {a}", color=INK, fontsize=12, loc="left", fontweight="bold")
        txt = (f"slope in the window,  G ∝ ε$^{{-p}}$\n"
               f"  measured envelope   p = {R['slope']:.2f}\n"
               f"  paper's objective   p = {R['slope_pred_curve']:.2f}   (p→∞ limit {P.ideal_exponent(a):.2f})\n"
               f"  order-2 explicit    p = {R['slope_fin_curve']:.2f}   (limit {P.finite_order_exponent(a, 2):.2f})")
        ax.text(0.03, 0.04, txt, transform=ax.transAxes, fontsize=8.5, color=INK2, va="bottom", linespacing=1.6)
        ax.set_xlabel("two-qubit gates (single-qubit gates free)", color=INK2, fontsize=9.5)
    np.atleast_1d(axes)[0].set_ylabel("trace distance to exact evolution", color=INK2, fontsize=9.5)
    handles, labels = np.atleast_1d(axes)[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5 if len(res) >= 3 else 2, frameon=False,
               fontsize=8.5, labelcolor=INK2, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(f"Cost vs accuracy for $a_j \\propto j^{{-\\alpha}}$ on {meta['n_qubits']} qubits, "
                 f"L = {meta['n_terms']}, t = {t}  —  measured frontiers (exact trace distance) vs the paper's objective",
                 color=INK, fontsize=12, x=0.01, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0.05 if len(res) >= 3 else 0.09, 1, 0.96))
    fig.savefig(os.path.join(RESULTS, "frontiers.png"), facecolor=SURF)
    print("wrote results/frontiers.png")


def fig_split(res: list[dict], t: float) -> None:
    fig, axes = plt.subplots(1, len(res), figsize=(5.0 * len(res), 4.4), dpi=200, sharey=True)
    fig.patch.set_facecolor(SURF)
    for ax, R in zip(np.atleast_1d(axes), res):
        style(ax)
        ax.plot(R["eps"], R["k_pred"], color=C_TH, lw=1.8, ls=(0, (5, 3)), label="paper's K*(ε)")
        ax.plot(R["eps"], np.maximum(R["k_fin"], 0.5), color=C_TH, lw=1.4, ls=(0, (1.5, 2.5)), label="K*(ε) with order 2 explicit (Fact 3)")
        ax.text(0.03, 0.96, f"K ∝ ε$^{{-q}}$ in the window:  measured q = {R['kslope_meas']:.2f}\n"
                f"paper q = {R['kslope_pred']:.2f}  (limit {P.ideal_exponent(R['alpha']):.2f}),  "
                f"order-2 explicit q = {R['kslope_fin']:.2f}",
                transform=ax.transAxes, fontsize=8.3, color=INK2, va="top", linespacing=1.5)
        ax.step(R["eps"], np.maximum(R["k_meas"], 0.5), where="mid", color=C_COMP, lw=2,
                label="cheapest measured split (0.5 = pure qDRIFT)")
        ax.axhline(R["L"], color=C_TR, lw=1.2, ls=":", label="K = L (pure Trotter)")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_xaxis()
        ax.set_title(f"α = {R['alpha']}", color=INK, fontsize=12, loc="left", fontweight="bold")
        ax.set_xlabel("target trace distance ε", color=INK2, fontsize=9.5)
    np.atleast_1d(axes)[0].set_ylabel("split K (terms Trotterized)", color=INK2, fontsize=9.5)
    handles, labels = np.atleast_1d(axes)[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4 if len(res) >= 3 else 2, frameon=False,
               fontsize=8.5, labelcolor=INK2, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Optimal split K: measured vs the paper's K*(ε) = argmin$_K$[Kt + t²λ$_K$²/ε]",
                 color=INK, fontsize=12, x=0.01, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0.07 if len(res) >= 3 else 0.14, 1, 0.95))
    fig.savefig(os.path.join(RESULTS, "split.png"), facecolor=SURF)
    print("wrote results/split.png")


def report(res: list[dict]) -> None:
    lines = []
    for R in res:
        a = R["alpha"]
        lines.append(f"\nalpha = {a}   (L = {R['L']}, window with 2 <= K* <= L/4: {R['window'].sum()} eps points)")
        lines.append(f"  slope of G(eps) in the window: measured envelope {R['slope']:.2f} | paper's objective "
                     f"{R['slope_pred_curve']:.2f} (p->inf limit {P.ideal_exponent(a):.2f}) | order-2-explicit objective "
                     f"{R['slope_fin_curve']:.2f} (limit {P.finite_order_exponent(a, 2):.2f})")
        lines.append(f"  cut exponent K ~ eps^-q in the window: measured {R['kslope_meas']:.2f} | paper {R['kslope_pred']:.2f} "
                     f"(limit {P.ideal_exponent(a):.2f}) | order-2 explicit {R['kslope_fin']:.2f}")
        lines.append(f"  {'eps':>9} {'best K':>7} {'K* paper':>9} {'K* ord2':>8} {'envelope':>10} {'Trotter':>10} {'qDRIFT':>11} {'vs Trotter':>11} {'vs qDRIFT':>10}")
        for e, km, kp, kf, ge, gt, gq in zip(R["eps"], R["k_meas"], R["k_pred"], R["k_fin"], R["g_env"], R["g_tr"], R["g_qd"]):
            f = lambda v: f"{v:10.0f}" if np.isfinite(v) else f"{'--':>10}"
            r = lambda v: f"{v / ge:10.1f}x" if np.isfinite(v) else f"{'--':>11}"
            lines.append(f"  {e:9.1e} {km:7.0f} {kp:9.0f} {kf:8.0f} {ge:10.0f} {f(gt)} {f(gq):>11} {r(gt)} {r(gq)}")
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(RESULTS, "summary.txt"), "w") as fh:
        fh.write(text + "\n")


def main() -> None:
    with open(os.path.join(RESULTS, "meta.json")) as fh:
        meta = json.load(fh)
    alphas = [float(a) for a in sys.argv[1:]] or meta["alphas"]
    res = [analyse(a, meta["time"]) for a in alphas]
    fig_frontiers(res, meta["time"], meta)
    fig_split(res, meta["time"])
    report(res)


if __name__ == "__main__":
    main()
