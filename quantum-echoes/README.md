# Quantum echoes (OTOC and OTOC(2)) — in Q#

A Q# demo of *Constructive interference at the edge of quantum ergodic dynamics* (Google
Quantum AI and collaborators: Abanin, Acharya, Aghababaie-Beni et al. —
[arXiv:2506.10191](https://arxiv.org/abs/2506.10191); in Nature, October 2025, as
"Observation of constructive interference at the edge of quantum ergodicity"). It runs
interactively, with the figures, on Copenhagen's Papers page:
https://copenhagen.strathweb.com/papers/quantum-echoes

## The claim

Run a random circuit U forward, flip one qubit (the *butterfly* q_b, B = X), run U backward,
and measure M = Z on another qubit q_m. In the Heisenberg picture that is
B(t) = U(t)† B U(t), and the out-of-time-order correlators are

  C^(2k) = ⟨ψ| (B(t) M)^(2k) |ψ⟩

C^(2) is the OTOC and C^(4) the second-order OTOC(2). Written as a sum over Pauli strings,
C^(2) is dominated by diagonal terms, while C^(4) has off-diagonal terms that pair strings
along large loops and add up constructively. Classical approximations that leave out those
loops match the processor on C^(2) but not on C^(4). On the 65-qubit, 23-cycle circuits of
the 103-qubit Willow processor the paper estimates exact simulation at about 13,000 times
the processor's data-taking time.

## The Q#

`src/Echoes.qs` is the protocol, with no entry point:

| | |
|---|---|
| `XYRotation` | the random single-qubit gate: a rotation by θ about an axis at angle φ in the XY plane (θ/π ∈ {¼, ½, ¾}, φ/π uniform in [−1, 1)) |
| `ISwapLike` | the paper's two-qubit gate, an iSWAP then a conditional phase of about 0.35 rad, merged into a SWAP and three diagonal rotations |
| `CouplerPatterns`, `ApplyCycle`, `ApplyCycles` | U(t) on a rows × cols grid: rotations on every qubit, then one coupler pattern per cycle |
| `ApplyCyclesWithPaulis` | U with a Pauli string inserted before one cycle (the paper's test of the interference) |
| `ApplyEchoedButterfly` | B(t), as `within { U } apply { X(q_b) }`: Q# derives U† from the gates |
| `ApplyNestedEcho` | B(t) [M B(t)]^(k−1). From \|0…0⟩, ⟨Z⟩ of q_m afterwards is C^(2k), since M\|0⟩ = \|0⟩ |
| `MeasureOtoc` | one shot of C^(2k): the nested echo, then q_m measured. C^(2k) = P(Zero) − P(One) |

Gate parameters are passed in as flat arrays (`cycles · rows · cols` values each, cycle by
cycle), so a seed on the caller's side names one circuit.

## Using it

From Python, with the `qdk` package (`pip install qdk numpy`), in this directory:

```python
import numpy as np
from qdk import code, qsharp

qsharp.init(project_root=".")
rows, cols, t = 3, 4, 6          # 12 qubits, 6 cycles
n = rows * cols
rng = np.random.default_rng(7)
thetas = rng.choice([0.25, 0.5, 0.75], size=t * n).tolist()
phis = rng.uniform(-1, 1, size=t * n).tolist()

for k in (1, 2):                 # C^(2), then C^(4)
    results = qsharp.run(code.Echoes.MeasureOtoc, 400, rows, cols, t, thetas, phis, 7, 5, k)
    c = sum(1 if r == qsharp.Result.Zero else -1 for r in results) / len(results)
    print(f"C^({2 * k}) ~ {c:+.3f}")
```

With q_m = 5 and the butterfly two sites along the same row (q_b = 7), this gives a C^(2)
of about 0.8 and a C^(4) of about 0.15 for this seed. The 800 shots take under a minute. Before the
butterfly's light cone reaches q_m, both correlators are exactly 1.

## Differences from the paper

- 10 and 12 qubits instead of 40 to 103, with no noise and no error mitigation.
- The couplers switch on in a simple fixed order (horizontal pairs from even columns, then
  odd columns, then vertical pairs from even rows, then odd rows), not the paper's
  slow-scrambling pattern, so a 12-qubit patch scrambles within about ten cycles.
- Every qubit starts in |0⟩; the paper starts from a product state with q_m in |0⟩.
- At these sizes Pauli insertion moves C^(4) more than C^(2), the direction the paper
  reports, but C^(4) is still mostly diagonal. The paper's result that large-loop
  interference dominates C^(4) is about 40-qubit circuits.
