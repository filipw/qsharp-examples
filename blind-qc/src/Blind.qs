/// Blind quantum computation by measurement, after Song, Knörzer, Dalton, Wallraff and Besse
/// (arXiv:2605.14656), whose protocol is Morimae and Fujii's (arXiv:1201.3966). A server makes a
/// cluster state and hands it over, one node at a time; a client measures it at angles it never
/// announces. Classical bits stay with the client: every `if` on an outcome and every controlled
/// correction below acts on a client qubit, and no outcome ever reaches a server qubit.
///
/// The basis B(theta) is { |0> + e^(i theta)|1>, |0> - e^(i theta)|1> } / sqrt 2, outcome 0 for
/// the + state. Measuring it is Rz(-theta), H, then a Z readout. In every cell below a client pulse
/// Rz(x), H and a readout s apply X^s H Rz(x) to the server's logical qubit. (In maths s is always
/// an outcome; in code `s` is the server's qubit.)
import Std.Diagnostics.DumpMachine;
import Std.Math.*;
import Std.Random.DrawRandomDouble;

/// Error rates, one number per kind of operation, as the paper's Appendix G quotes them.
struct Rates {
    OneQubit : Double,
    TwoQubit : Double,
    Readout : Double,
    Idle : Double,
}

function NoNoise() : Rates {
    new Rates { OneQubit = 0.0, TwoQubit = 0.0, Readout = 0.0, Idle = 0.0 }
}

/// X, Y or Z, each with probability rate / 2: average infidelity `rate`, and a bit flip with
/// probability exactly `rate` on a Z readout.
operation Noise1(rate : Double, q : Qubit) : Unit {
    if rate > 0.0 {
        let r = DrawRandomDouble(0.0, 1.0);
        if r < 0.5 * rate {
            X(q);
        } elif r < rate {
            Y(q);
        } elif r < 1.5 * rate {
            Z(q);
        }
    }
}

/// Each of the 15 non-identity Pauli pairs with probability rate / 12: average infidelity `rate`.
operation Noise2(rate : Double, a : Qubit, b : Qubit) : Unit {
    if rate > 0.0 {
        let k = Floor(DrawRandomDouble(0.0, 1.0) * 12.0 / rate);
        if k < 15 {
            let paulis = [PauliI, PauliX, PauliY, PauliZ];
            ApplyP(paulis[(k + 1) / 4], a);
            ApplyP(paulis[(k + 1) % 4], b);
        }
    }
}

/// H and its noise. Rz has none: it is a virtual gate on the hardware.
operation NoisyH(q : Qubit, noise : Rates) : Unit {
    H(q);
    Noise1(noise.OneQubit, q);
}

/// CZ and its noise.
operation NoisyCZ(a : Qubit, b : Qubit, noise : Rates) : Unit {
    CZ(a, b);
    Noise2(noise.TwoQubit, a, b);
}

/// CNOT as the hardware does it: a CZ between two H on the target.
operation NoisyCNOT(control : Qubit, target : Qubit, noise : Rates) : Unit {
    NoisyH(target, noise);
    NoisyCZ(control, target, noise);
    NoisyH(target, noise);
}

/// A qubit waiting.
operation Idle(q : Qubit, noise : Rates) : Unit {
    Noise1(noise.Idle, q);
}

/// A Z readout with its error, leaving the qubit in |0>.
operation Read(q : Qubit, noise : Rates) : Result {
    Noise1(noise.Readout, q);
    MResetZ(q)
}

/// One measurement cell: the server's node goes to the client qubit `c`, which is read after a
/// pulse Rz(x), and the server's qubit becomes the next node. Measuring B(theta) means x = -theta.
operation Cell(s : Qubit, c : Qubit, x : Double, noise : Rates) : Result {
    NoisyCNOT(s, c, noise);
    Rz(x, c);
    NoisyH(c, noise);
    let r = Read(c, noise);
    NoisyH(s, noise);
    r
}

/// The same cell for a node measured in Z: no pulse, and no H on the client.
operation ZCell(s : Qubit, c : Qubit, noise : Rates) : Result {
    NoisyCNOT(s, c, noise);
    let r = Read(c, noise);
    NoisyH(s, noise);
    r
}

/// Moves the server's qubit to the client's, which is in |0>: two CNOTs do what a SWAP would.
/// `c` holds nothing of the client's here: it was reset, so this gate carries no angle back.
operation HandOver(s : Qubit, c : Qubit, noise : Rates) : Unit {
    NoisyCNOT(s, c, noise);
    NoisyCNOT(c, s, noise);
}

/// Three cells and the hand-over, measured for real: H Rz(a3) Rx(a2) Rz(a1) on the server's qubit,
/// which ends on the client's. The outcomes only ever choose the client's next pulse and its
/// corrections. Returns (s1, s2, s3).
operation ThreeCells(s : Qubit, c : Qubit, a1 : Double, a2 : Double, a3 : Double, noise : Rates) : (Result, Result, Result) {
    let s1 = Cell(s, c, a1, noise);
    Idle(s, noise);
    let x2 = Sign(s1) * a2;
    let s2 = Cell(s, c, x2, noise);
    Idle(s, noise);
    let x3 = Sign(s2) * a3 + (s1 == One ? PI() | 0.0);
    let s3 = Cell(s, c, x3, noise);
    Idle(s, noise);
    HandOver(s, c, noise);
    if s3 == One {
        X(c);
        Noise1(noise.OneQubit, c);
    }
    if s2 == One {
        Z(c);
    }
    (s1, s2, s3)
}

/// (-1)^s.
function Sign(s : Result) : Double {
    s == One ? -1.0 | 1.0
}

/// The same three cells with the measurements deferred, so a dump after each stage is the state
/// averaged over the outcomes, which is what a server that never hears them holds. A client
/// qubit's Z value stands for its outcome, and each adaptive pulse becomes a controlled one.
/// Dumps five stages over (s, c[0], c[1], c[2], out).
operation DeferredCells(s : Qubit, c : Qubit[], out : Qubit, a1 : Double, a2 : Double, a3 : Double) : Unit {
    H(s);
    DumpMachine();

    CNOT(s, c[0]);
    Rz(a1, c[0]);
    H(c[0]);
    H(s);
    DumpMachine();

    CNOT(s, c[1]);
    Rz(a2, c[1]);
    Controlled Rz([c[0]], (-2.0 * a2, c[1]));
    H(c[1]);
    H(s);
    DumpMachine();

    CNOT(s, c[2]);
    Rz(a3, c[2]);
    Controlled Rz([c[1]], (-2.0 * a3, c[2]));
    Controlled Z([c[0]], c[2]);
    H(c[2]);
    H(s);
    DumpMachine();

    SWAP(s, out);
    Controlled X([c[2]], out);
    Controlled Z([c[1]], out);
    DumpMachine();
}

/// A node's measurement in B(theta), before its readout.
operation BasisPulse(theta : Double, q : Qubit) : Unit {
    Rz(-theta, q);
    H(q);
}

/// The 7-node Deutsch-Jozsa cluster of the paper's Fig. 4, every node at once. Rows D1, A and D2
/// run left to right; the program is only the three labels on the A row. Dumps twice over (a1, d1,
/// a2, d2, o1, a3, o2), with all measurements left unread: after the client's pulses on columns 1
/// and 2, and after column 3.
operation ClusterPattern(balanced : Bool) : Unit {
    use nodes = Qubit[7];
    let (a1, d1, a2, d2, o1, a3, o2) = (nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6]);
    ApplyToEach(H, nodes);
    CZ(a1, a2);
    CZ(a2, a3);
    CZ(d1, a2);
    CZ(a2, d2);
    CZ(d1, o1);
    CZ(d2, o2);

    if balanced {
        BasisPulse(PI(), a1);
        BasisPulse(0.0, a2);
    }
    BasisPulse(0.0, d1);
    BasisPulse(0.0, d2);
    DumpMachine();

    BasisPulse(balanced ? 0.0 | PI(), a3);
    DumpMachine();

    ResetAll(nodes);
}

/// The same pattern on the hardware's six qubits, rows s[0..2] and c[0..2] for D1, A and D2, each
/// server qubit reused as its row's next node. Returns the decoded string as 4 D1 + 2 A + D2.
/// Idling: the paper counts 12 for this size and doesn't say where, so 2 on each row's data after
/// cycle 2 and 2 after cycle 3 is this code's choice.
operation DjReplay(balanced : Bool, noise : Rates) : Int {
    use s = Qubit[3];
    use c = Qubit[3];

    // Cycle 1: the A row's first node.
    NoisyH(s[1], noise);
    let p = balanced ? Cell(s[1], c[1], -PI(), noise) | ZCell(s[1], c[1], noise);

    // Cycle 2: the data nodes join, and all three rows are measured.
    NoisyH(s[0], noise);
    NoisyH(s[2], noise);
    NoisyCZ(s[0], s[1], noise);
    NoisyCZ(s[1], s[2], noise);
    let u1 = Cell(s[0], c[0], 0.0, noise);
    let u2 = Cell(s[2], c[2], 0.0, noise);
    let q = balanced ? Cell(s[1], c[1], 0.0, noise) | ZCell(s[1], c[1], noise);
    for row in 0..2 {
        Idle(s[row], noise);
        Idle(s[row], noise);
    }

    // Cycle 3: the last nodes go to the client, which reads them.
    for row in 0..2 {
        HandOver(s[row], c[row], noise);
        Idle(c[row], noise);
        Idle(c[row], noise);
    }
    let m1 = Read(c[0], noise);
    let m2 = Read(c[2], noise);
    Rz(balanced ? 0.0 | -PI(), c[1]);
    NoisyH(c[1], noise);
    let mA = Read(c[1], noise);
    ResetAll(s);

    // The client's decoding: XORs, no angle depends on an outcome here.
    let b = Bit(balanced ? p | q);
    let d1 = Bit(m1) ^^^ Bit(u1) ^^^ b;
    let a = Bit(mA) ^^^ b;
    let d2 = Bit(m2) ^^^ Bit(u2) ^^^ b;
    4 * d1 + 2 * a + d2
}

function Bit(r : Result) : Int {
    r == One ? 1 | 0
}
