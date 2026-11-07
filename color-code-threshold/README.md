# Colour codes below threshold — in Q#

A Q# demo of *Quantum error correction well below the color code threshold* (Quantinuum:
Lavasani, Reed, Dreiling et al. — [arXiv:2610.24077](https://arxiv.org/abs/2610.24077)). It
runs interactively, with the figures, on Copenhagen's Papers page:
https://copenhagen.strathweb.com/papers/color-code-threshold

## The claim

Below a threshold error rate, growing a code from distance d to d + 2 divides the logical
error by a constant factor Λ, so logical errors fall exponentially with distance:

  p_L(d) ≈ A · (p / p_th)^((d+1)/2),   Λ = p_L(d) / p_L(d+2)

On the trapped-ion processor Helios, triangular colour-code memories (20 rounds of syndrome
extraction, all six logical Pauli eigenstates) reach a logical error per round of 5.5×10⁻⁴
at d = 3 (13 qubits), 3.8×10⁻⁵ at d = 5 (37 qubits) and 5×10⁻⁶ at d = 7 (73 qubits): a fitted
Λ = 12 (+3/−2), decoded with Tesseract.

## The Q#

`src/ColorCode.qs` is the code as a quantum memory at d = 3 (Steane's [[7,1,3]]) and d = 5
([[19,1,5]]), with no entry point. Checks are measured with one ancilla, reset after each
measurement, so d = 5 needs 20 qubits.

| | |
|---|---|
| `Faces`, `DataQubitCount`, `LogicalSupport` | the code: faces from a triangular grid, (3d² + 1)/4 data qubits, one side of the triangle as the logical operator |
| `ApplyDepolarizing` | the noise: X, Y or Z, each with probability p/3 |
| `MeasureZCheck`, `MeasureXCheck`, `MeasureChecks` | the parity of a face in the Z or X basis, and every check of one type as a bitmask |
| `PrepareLogical` | logical \|0⟩ or \|+⟩, as the paper prepares it: one round of the other type of check projects the product state onto the code |
| `ReadOut` | transversal read-out as a bitmask over the data qubits |
| `MemoryShot` | one memory shot: prepare, `rounds` rounds of noise and both types of checks, read out. Returns one syndrome bitmask per round, then the read-out |
| `ZSectorShot` | the Z sector of a \|0⟩ memory alone. Under Pauli noise and perfect checks it gives the same decoded statistics as `MemoryShot` with `PauliZ`, at a fraction of the cost |

The noise is code-capacity depolarizing noise drawn in Q# itself: before each round every
data qubit suffers X, Y or Z with probability p/3 each, and the checks are perfect.

## Using it

From Python, with the `qdk` package (`pip install qdk`), in this directory. The decoder is a
lookup table: for each syndrome, the lowest-weight error that gives it. Each round the change
in syndrome is looked up and kept as a Pauli frame, and the read-out gives a last, perfect
syndrome.

```python
from itertools import combinations
from qdk import code, qsharp

qsharp.init(project_root=".")
d, rounds, p = 3, 10, 0.05
faces = qsharp.eval(f"ColorCode.Faces({d})")
logical = qsharp.eval(f"ColorCode.LogicalSupport({d})")
n = (3 * d * d + 1) // 4

def syndrome(flips):
    return sum(1 << i for i, f in enumerate(faces) if sum(flips >> q & 1 for q in f) % 2)

table = {}
for w in range(n + 1):
    for qs in combinations(range(n), w):
        table.setdefault(syndrome(sum(1 << q for q in qs)), sum(1 << q for q in qs))
    if len(table) == 1 << len(faces):
        break

def failed(record):
    frame, last = 0, 0
    for s in record[:-1]:
        frame ^= table[s ^ last]
        last = s
    readout = record[-1]
    frame ^= table[syndrome(readout) ^ last]
    return sum((readout ^ frame) >> q & 1 for q in logical) % 2

shots = qsharp.run(code.ColorCode.ZSectorShot, 1000, d, rounds, p)
print("logical failures:", sum(failed(r) for r in shots), "of", len(shots))
```

At p = 5% this gives about 170 failures in 1,000 ten-round shots, close to the exact rate of
this decoder, 2.0% per round, against 3.3% for a bare qubit. Use `d = 5` for the larger
code, and `MemoryShot(d, PauliZ, …)` or `MemoryShot(d, PauliX, …)` for the full memory.

## Differences from the paper

- Code-capacity noise with perfect checks, not Helios's circuit-level noise (two-qubit
  gates, transport, crosstalk, leakage). The threshold and Λ are this model's: the d = 3 and
  d = 5 curves cross at p ≈ 16%, and Λ reaches the paper's 12 only near p ≈ 1%.
- One bare ancilla per check, instead of the paper's SWAP-superdense extraction circuits.
  Under circuit-level noise a single ancilla fault can spread to two data qubits, and this
  d = 3 memory then fails linearly in p.
- A lookup-table decoder instead of Tesseract, and no d = 7.
