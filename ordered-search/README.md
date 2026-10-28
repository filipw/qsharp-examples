# Optimal quantum ordered search — in Q#

A Q# demo of *Optimal Quantum Algorithms for Ordered Search* (Joseph Carolan, Andrew Childs —
[arXiv:2609.34293](https://arxiv.org/abs/2609.34293)). It runs interactively, with the
figures, on Copenhagen's Papers page: https://copenhagen.strathweb.com/papers/ordered-search

## The claim

Finding the first 1 in a sorted bit string of length n takes ⌈log₂ n⌉ queries classically.
Høyer, Neerbek and Shi proved that any exact quantum algorithm needs at least
(1/π)(ln n − 1) ≈ 0.221 log₂ n, and algorithms had closed in from above for twenty years. The
paper reaches (1/π) ln n + o(log n) queries, optimal to leading order and asymptotically 4.5×
fewer than binary search. It gives two algorithms. The first, a zero-error algorithm from a
continuum relaxation, was discovered by Claude Fable 5. The second, an exact algorithm from
the Farhi–Goldstone–Gutmann–Sipser polynomial program, was discovered by GPT-5.6-Sol. This
demo is the first.

One query is U_s = F† O₀ F O_s: the comparison oracle O_s (a sign flip for x ≥ s), then the
sign function in momentum. Its eigenphase on the Mellin modes is θ(ω) = 2·atan(e^(−πω)), slope
−π at ω = 0, so each query moves a packet that is Gaussian in log-position about π closer to
the step in log-distance. The algorithm starts from such a packet at ln n + ln^(3/4) n with
width ln^(1/2) n, makes Q = ⌈(ln n + 2 ln^(3/4) n)/π⌉ queries, measures the position and
checks the answer with two classical comparisons, starting again if the check fails.

## The Q#

`src/OrderedSearch.qs`, with no entry point. The paper works on L²(ℝ) and gives no
discretisation; the one here is a ring of 2^m sites, two per list position, with the padded
list on one half and its complement on the other (the string x x̄ of the FGGS
translation-invariant algorithms). The header comment of the file describes it in full.

| | |
|---|---|
| `FlipAtOrAbove` | a phase of −1 on every site j ≥ c, with one multi-controlled Z per zero bit of c and one for j = c |
| `ComparisonOracle` | O_s on the ring: one query to the list |
| `SignInFourierBasis` | F† O₀ F: a single Z on the top frequency bit, between `ApplyQFT` and its adjoint |
| `SearchStep` | one query, U_s |
| `PrepareInitialState` | the paper's log-Gaussian, from the amplitudes of one half of the ring (`PreparePureStateD`), mirrored onto the other |

## Using it

From Python, with the `qdk` package (`pip install qdk numpy`), in this directory. This
searches a list of 4 on 9 qubits for the target 2:

```python
import numpy as np
from qdk import code, qsharp

qsharp.init(project_root=".")
qsharp.eval("""
operation Search(qubits : Int, halfAmplitudes : Double[], step : Int, queries : Int) : Int {
    use register = Qubit[qubits];
    OrderedSearch.PrepareInitialState(halfAmplitudes, register);
    for _ in 1..queries {
        OrderedSearch.SearchStep(step, register);
    }
    Std.Measurement.MeasureInteger(Std.Arrays.Reversed(register))
}
""")

n, m, target = 4, 9, 2
M, K = 2 ** m, 2                  # sites, sites per list position
N = M // (2 * K)                  # list positions in half the ring
L = np.log(n)
w, r = L ** 0.5, L ** 0.75
Q = int(np.ceil((L + 2 * r) / np.pi))
x = (np.arange(M) + 0.5) / K - N  # site positions, from the middle of the list
amps = np.abs(x) ** -0.5 * np.exp(-((np.log(np.abs(x)) - (L + r)) ** 2) / (4 * w * w))
amps /= np.linalg.norm(amps)
step = (N + target - n // 2) * K  # the first site of the target's cell

def candidate(j):                 # the list index a measured site stands for
    xj = (j + 0.5) / K - N
    xj = xj - N if xj >= N / 2 else xj + N if xj < -N / 2 else xj
    i = round(xj) + n // 2
    return i if 0 <= i < n else None

sites = qsharp.run(code.Search, 200, m, amps[M // 2:].tolist(), step, Q)
print(f"Q = {Q}, landed on {target}: {sum(candidate(j) == target for j in sites) / len(sites):.2f}")
```

After Q = 2 queries about two shots in three land on the target. A full zero-error search
checks each candidate with two classical comparisons and repeats on a failure. Lists of 8 and
16 need 11 and 12 qubits.

## Differences from the paper

- The discretisation: the paper's continuum algorithm, run on a finite ring. At n = 16 the
  ring holds 89% of the initial state, and the success after Q = 3 queries is 0.69 on the
  ring against 0.74 for the paper's operator on the line.
- Small n only. At these sizes a packet moves about 2.3 rather than π per query, the success
  probability is far from its asymptotic 1 − O(ln^(−1/2) n), and with checks and restarts the
  search needs fewer queries on average than binary search only from about n = 2,000.
- Neither the proofs nor the exact algorithm of Section 4.
