# Unconditional unclonable encryption — in Q#

A Q# demo of *Unconditional Unclonable Encryption* (Prabhanjan Ananth, Seyoon Ragavan, Amit
Sahai — [arXiv:2607.21551](https://arxiv.org/abs/2607.21551), v3). Version 3 merges Ananth and
Sahai's first version with Ragavan's concurrent *Efficient unclonable encryption from Pauli
eigenstates* ([arXiv:2607.21811](https://arxiv.org/abs/2607.21811)), the original form of the
same scheme. It runs interactively, with the figures, on Copenhagen's Papers page:
https://copenhagen.strathweb.com/papers/unclonable-encryption

## The claim

Alice encrypts a random bit m. A cloner who doesn't know the key splits the ciphertext between
Bob and Charlie, who can't talk to each other. Then the key is revealed to both, and they win
if both output m. A shared coin wins half the time. The scheme: the key is a uniformly random
non-identity Pauli string P on n qubits (one of 4ⁿ − 1), and m is encrypted as a random
(−1)ᵐ eigenstate of P, the state (I + (−1)ᵐ P)/2ⁿ. Encryption needs only single-qubit
Cliffords, and decryption measures P (Lemma 3.1: it is perfectly correct). Theorem 4.1 proves
that every splitting attack has both receivers right with probability at most
½ + ½·√(2ⁿ/(4ⁿ − 1)), about ½ + ½·2^(−n/2), with no assumptions. The authors state (§1.2)
that the construction and proof were found by GPT-5.6 Sol Ultra and verified by them, and
their Lean 4 code (https://github.com/sragavan99/lean-unclonable-enc-paulis) proves the
correctness and security claims of Theorem 4.1.

Restricting the keys to X and Z (BB84 states) breaks it: an attack of Coladangelo, Liu and Xie
wins with cos²(π/8) ≈ 0.854 at every n. Using Y and I as well is what makes the bound work.

## The Q#

Two files, with no entry point.

| | |
|---|---|
| `Scheme.KeyGen` | Gen (§3): each position drawn from I, X, Y, Z, the all-I string redrawn |
| `Scheme.Encrypt` | Enc (§3.1, Proposition 3.2): bits r with parity m on the key's non-identity positions, each \|r_j⟩ turned into an eigenstate of its Pauli |
| `Scheme.FromZ` | the U_j of §3.1: nothing for I and Z, H for X, H then S for Y (adjointable) |
| `Scheme.Decrypt` | Dec (§3): one `Measure(key, qubits)`, a joint Pauli measurement; `One` is the bit 1 |
| `Attacks.MeasureDiagonal` | measures a qubit along (x + y + z)/√3, inside `within … apply` |
| `Attacks.UniversalClone` | Bužek and Hillery's symmetric universal cloner, as the network of Bužek, Braunstein, Hillery and Bruß: two copies, each with fidelity 5/6 |

**Both attacks are our own constructions, not the paper's.** Their success probabilities are
our derivations, each checked by an exact computation over every key at n = 1 to 4:

- *Measure and broadcast*: the cloner measures every qubit along the diagonal and sends the
  same signs to Bob and Charlie, who take their parity on the key's non-identity positions.
  P_A(n) = ½ + ((1 + √3)ⁿ − 1)/(2(4ⁿ − 1)): 0.789, 0.715, 0.654, 0.607 at n = 1 to 4. It equals
  the "identity attack" of Botteron, Designolle and Fawzi
  ([arXiv:2607.28602](https://arxiv.org/abs/2607.28602), Table 1). At n = 1 it is optimal.
- *Universal copier*: the cloner copies every qubit, and Bob and Charlie each decrypt their own
  copy. P_C(n) = [4ⁿ − 1 + 2(3ⁿ − 1) + (2ⁿ − 1)] / (4(4ⁿ − 1)): 0.667, 0.567, 0.484, 0.422 at
  n = 1 to 4, worse than a coin from n = 3.

## Using it

From Python, with the `qdk` package (`pip install qdk`), in this directory. This runs an
honest round and both attacks, 1,000 rounds each, and compares them with the bound:

```python
from math import sqrt

from qdk import code, qsharp

qsharp.init(project_root=".")
qsharp.eval("""
import Std.Random.*;
import Scheme.*;
import Attacks.*;

/// An honest round: encrypt a random bit under a fresh key, then decrypt it.
operation Honest(n : Int) : Bool {
    let key = KeyGen(n);
    let m = DrawRandomInt(0, 1) == 1;
    use qubits = Qubit[n];
    Encrypt(key, m, qubits);
    let decrypted = Decrypt(key, qubits);
    ResetAll(qubits);
    decrypted == m
}

/// Measure and broadcast: the cloner measures every qubit along the diagonal and
/// sends the signs to both receivers, who take their parity on the key's support.
operation Broadcast(n : Int) : Bool {
    let key = KeyGen(n);
    let m = DrawRandomInt(0, 1) == 1;
    use qubits = Qubit[n];
    Encrypt(key, m, qubits);
    mutable guess = false;
    for j in 0..n - 1 {
        let sign = MeasureDiagonal(qubits[j]);
        if key[j] != PauliI {
            guess = guess != sign;
        }
    }
    ResetAll(qubits);
    guess == m
}

/// The universal copier: every qubit is cloned, and Bob and Charlie each decrypt
/// their own copy, one position at a time on two reused qubits.
operation Copier(n : Int) : Bool {
    let key = KeyGen(n);
    let m = DrawRandomInt(0, 1) == 1;
    use qubits = Qubit[n];
    use (blank, machine) = (Qubit(), Qubit());
    Encrypt(key, m, qubits);
    mutable bob = false;
    mutable charlie = false;
    for j in 0..n - 1 {
        UniversalClone(qubits[j], blank, machine);
        if key[j] != PauliI {
            bob = bob != (Measure([key[j]], [qubits[j]]) == One);
            charlie = charlie != (Measure([key[j]], [blank]) == One);
        }
        Reset(blank);
        Reset(machine);
    }
    ResetAll(qubits);
    bob == m and charlie == m
}
""")

def bound(n):      # Theorem 4.1, Eq. 17
    return 0.5 + 0.5 * sqrt(2**n / (4**n - 1))

def broadcast(n):  # measure and broadcast, exact
    return 0.5 + ((1 + sqrt(3)) ** n - 1) / (2 * (4**n - 1))

def copier(n):     # universal copier, exact
    return (4**n - 1 + 2 * (3**n - 1) + 2**n - 1) / (4 * (4**n - 1))

shots = 1000
print(f"honest, n = 4: {sum(qsharp.run(code.Honest, shots, 4)) / shots:.3f} decrypted")
for n in [1, 2, 3, 4]:
    b = sum(qsharp.run(code.Broadcast, shots, n)) / shots
    c = sum(qsharp.run(code.Copier, shots, n)) / shots
    print(f"n = {n}: bound {bound(n):.3f}  broadcast {b:.3f} (exact {broadcast(n):.3f})"
          f"  copier {c:.3f} (exact {copier(n):.3f})")
```

Every honest round decrypts. The attacks land within shot noise (about ±0.015 at 1,000
rounds) of their exact values, which sit under the bound. From n = 3 the gap to the bound is
smaller than that noise, so a single estimate can even land above the bound: the theorem
limits the probability, not one estimate of it. Separating them takes thousands of rounds
per n.

## Differences from the paper

- The paper has no code and no attacks; both attacks here, and their formulas, are ours.
- The copier clones one position at a time onto two reused qubits (n + 2 qubits rather than
  3n). Nothing acts on two positions together, so this is the same as cloning everything first.
- Neither the proof of Theorem 4.1 nor the optimal attack: for n ≥ 2 that is open.
