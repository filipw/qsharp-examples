# Unclonable Cryptography in Linear Quantum Memory

A runnable Q# demonstration of the core ideas in:

> **"Unclonable Cryptography in Linear Quantum Memory"**  
> Omri Shmueli and Mark Zhandry, arXiv:2511.04633

## What the demo shows

The paper constructs **one-shot signatures (OSS)**: a quantum signing key that
can be used to sign *exactly one* message (no-cloning prevents copying), while
achieving **O(λ) quantum memory** - linear in the security parameter - down
from the Ω(λ²)–Ω(λ³) of prior work.

The demo is structured as five Steps, one per core idea:

| Step | File | Paper idea |
|-----|------|------------|
| 1 | `CosetKey.qs` | Signing key = **coset state** `|sk⟩ = Σ_{u∈S} |u⟩` |
| 2 | `CosetKey.qs` | Measuring collapses to **one** element → one-shot |
| 3 | `Signing.qs` | **Measure-and-correct** signs a *chosen* message in O(log λ) rounds |
| 4 | `ClawFree.qs` | **Claw-free permutations** make forging a second signature hard |
| 5 | `Folding.qs` | **Folding** compresses λ² bits → O(λ) bits for the quantum key |

### Concrete instantiation

The demo uses a tiny 5-qubit coset state over `Z₂⁵`:

```
S = span{10100, 01010, 00001} = { (a, b, a, b, c) }
```

- Bits 0–1: message `(a, b)` - λ = 2 in the demo
- Bits 2–3: locked copies of the message (a forger must break this coupling)
- Bit 4: internal entropy (`c` free)

`|S|` = 8; all amplitudes are `1/√8 ≈ 0.3536`.

## Running the demo

```bash
pip install qdk
python run_demo.py
```

The narrated Python driver prints a guided walkthrough with live Q# output,
amplitude dumps, circuit diagrams, and statistics for all five acts.

---

## Sample output (abridged)

```
STEP 1 - coset key prepared; 8 basis states each with amplitude 0.3536
STEP 2 - 500 measurements: uniform over all 4 messages; all valid in S
STEP 3 - chosen-message signing: 500/500 valid, avg ~2.6 rounds per message
STEP 4 - claw-free example: claw at y=4 is (b=0,x=21) and (b=1,x=3)
STEP 5 - folding: λ=16 compresses 256-bit input → 31-bit key (8× shrink)
          roundtrip 50/50, trapdoor-free simulation 50/50
```


## Key Q# operations

| Operation | Signature | Description |
|-----------|-----------|-------------|
| `PrepareCosetKey` | `(qs : Qubit[]) : Unit is Adj+Ctl` | H on message + entropy qubits, CNOT to copy bits |
| `InSubspace` | `(bits : Bool[]) : Bool` | Membership test: bits[2]==bits[0] ∧ bits[3]==bits[1] |
| `SignChosenMessage` | `(m0 : Bool, m1 : Bool) : (Int, Bool[])` | Repeat-until: measure, re-coherise mismatches |
| `FindClaw` | `(lambda : Int) : (Int, Int, Int, Bool)` | Sample claw `(y, x0, x1, ok)` with `Π₀(x0)=Π₁(x1)=y`, using the trapdoor |
| `GenFoldingKey` | `(lam : Int) : FoldingKey` | Sample λ random 2-to-1 maps |
| `Fold` / `Unfold` | `(key, bs, xs)` | Round-trip through the folding transform |
| `UnfoldWithoutTrapdoor` | `(key, bs, ys, wsum, istar)` | Recover `xs[istar]` without knowing the istar-th trapdoor |
