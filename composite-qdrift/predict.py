"""What arXiv:2607.19852 predicts for the sampled coefficient vectors.

The paper's cost objective, for sum_j a_j = 1 and trace-distance error eps, is

    cost(K, eps) = K t + t^2 lambda_K^2 / eps,     lambda_K = sum_{j>K} a_j,

up to constants, minimised over the split K. Everything here is derived from that
line: the predicted optimal split K*(eps), the predicted gate count up to one
overall constant, and for a_j ~ j^(-alpha) the predicted exponent of 1/eps.
"""

from __future__ import annotations

import csv
import os

import numpy as np

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")

# D = sqrt(1 - |<psi|phi>|^2) for a pure output saturates double precision around
# 1e-12 in the infidelity, i.e. ~1e-6 in D; below that the number measures the
# simulator's arithmetic, not the algorithm.
FLOOR = 3e-7


def load_coeffs(alpha: float, results: str | None = None) -> np.ndarray:
    a = np.load(os.path.join(results or RESULTS, f"coeffs_alpha{alpha}.npy"))
    return np.sort(a)[::-1]


def tails(a: np.ndarray) -> np.ndarray:
    """lambda_K for K = 0..L (lambda_0 = 1, lambda_L = 0)."""
    return np.concatenate([np.cumsum(a[::-1])[::-1], [0.0]])


def objective(a: np.ndarray, t: float, eps: float) -> np.ndarray:
    """cost(K, eps) for every K = 0..L."""
    K = np.arange(len(a) + 1)
    return K * t + (t * tails(a)) ** 2 / eps


def k_star(a: np.ndarray, t: float, eps: float) -> int:
    return int(np.argmin(objective(a, t, eps)))


def predicted_cost(a: np.ndarray, t: float, eps: float) -> float:
    return float(objective(a, t, eps).min())


def stages(order: int) -> int:
    return 1 if order == 1 else 2 * 5 ** (order // 2 - 1)


def objective_finite(a: np.ndarray, t: float, eps: float, order: int) -> np.ndarray:
    """The paper's Fact 3 (Hagan-Wiebe Thm 2.1) with the Trotter order kept explicit.

    Per outer step the nested channel costs Upsilon (Upsilon K + N_B) and needs
    r ~ (t/eps)^(1/order) steps, while qDRIFT needs 4 t^2 lambda_K^2 / eps samples in
    total (HW Thm 4). Dropping the o(1) idealisation therefore turns K t into
    Upsilon^2 K t (t/eps)^(1/order); the tail term keeps its constant 4.
    """
    K = np.arange(len(a) + 1)
    ups = stages(order)
    return ups * ups * K * t * (t / eps) ** (1.0 / order) + 4.0 * (t * tails(a)) ** 2 / eps


def k_star_finite(a: np.ndarray, t: float, eps: float, order: int) -> int:
    return int(np.argmin(objective_finite(a, t, eps, order)))


def predicted_cost_finite(a: np.ndarray, t: float, eps: float, order: int) -> float:
    return float(objective_finite(a, t, eps, order).min())


def ideal_exponent(alpha: float) -> float:
    """G ~ eps^(-p): the paper's p = 1/(2 alpha - 1) for lambda_K ~ K^(1-alpha)."""
    return 1.0 / (2.0 * alpha - 1.0)


def finite_order_exponent(alpha: float, order: int) -> float:
    """Same minimisation with the Trotter term carrying its real (t/eps)^(1/order)
    instead of the paper's (t/eps)^o(1): p = (1 - 1/order)/(2 alpha - 1) + 1/order."""
    q = 1.0 / order
    return (1.0 - q) / (2.0 * alpha - 1.0) + q


def load_sweep(alpha: float, results: str | None = None) -> list[dict]:
    with open(os.path.join(results or RESULTS, f"sweep_alpha{alpha}.csv")) as fh:
        rows = []
        for r in csv.DictReader(fh):
            rows.append({k: (float(v) if k in ("alpha", "twoq", "trace_distance", "infidelity", "seconds")
                             else int(v) if k not in ("family",) else v) for k, v in r.items()})
    kept = [r for r in rows if r["trace_distance"] > FLOOR]
    if len(kept) < len(rows):
        print(f"alpha={alpha}: dropped {len(rows) - len(kept)} configuration(s) at the precision floor (D <= {FLOOR:g})")
    return kept


def pareto(points: list[dict], cost="twoq", err="trace_distance") -> list[dict]:
    """Non-dominated configurations: nothing cheaper is also more accurate."""
    best, out = np.inf, []
    for p in sorted(points, key=lambda r: (r[cost], r[err])):
        if p[err] < best:
            out.append(p)
            best = p[err]
    return out


def cost_at(front: list[dict], eps: float, cost="twoq", err="trace_distance") -> float | None:
    """Cost the frontier needs for exactly this accuracy, interpolated in log-log."""
    pts = sorted(front, key=lambda p: p[err])
    y = np.log10([p[err] for p in pts])
    x = np.log10([p[cost] for p in pts])
    if len(pts) < 2 or not (y[0] <= np.log10(eps) <= y[-1]):
        return None
    return float(10 ** np.interp(np.log10(eps), y, x))


def best_split_at(rows: list[dict], eps: float) -> tuple[int | None, float | None]:
    """Cheapest measured configuration family (by K) at accuracy eps."""
    best = (None, None)
    for K in sorted({r["K"] for r in rows}):
        c = cost_at(pareto([r for r in rows if r["K"] == K]), eps)
        if c is not None and (best[1] is None or c < best[1]):
            best = (K, c)
    return best


def fit_exponent(eps_grid, costs) -> float:
    """Slope of log cost vs log(1/eps)."""
    e = np.array(eps_grid)
    g = np.array(costs)
    m = np.isfinite(g)
    slope, _ = np.polyfit(np.log(1 / e[m]), np.log(g[m]), 1)
    return float(slope)
