/// Bell difference sampling and generator measurement for stabilizer states.
///
/// The unknown state is prepared by a Clifford circuit passed in as gate codes
/// (the learner never sees it; it only ever gets copies). Every operation here is
/// Base-profile QIR, so it runs on the QDK's stabilizer simulator at any n.
import Std.Measurement.MeasureEachZ;

/// Gate codes: 0 = H(a), 1 = S(a), 2 = CNOT(a, b).
operation ApplyCircuit(gates : Int[], a : Int[], b : Int[], qs : Qubit[]) : Unit {
    for i in 0..Length(gates) - 1 {
        if gates[i] == 0 {
            H(qs[a[i]]);
        } elif gates[i] == 1 {
            S(qs[a[i]]);
        } else {
            CNOT(qs[a[i]], qs[b[i]]);
        }
    }
}

/// Measures every pair (x[i], y[i]) in the Bell basis. The first n results are the
/// Z-part bits, the last n the X-part bits of the Bell label y = (a, b) with
/// |Phi_y> = (Z^a X^b (x) I)|Phi_0>.
operation BellMeasure(x : Qubit[], y : Qubit[]) : Result[] {
    for i in 0..Length(x) - 1 {
        CNOT(x[i], y[i]);
        H(x[i]);
    }
    MeasureEachZ(x) + MeasureEachZ(y)
}

/// One Bell difference sample: four copies, two Bell measurements. Returns the
/// two Bell labels (2n bits each); the learner XORs them.
operation BellDifferenceSample(n : Int, gates : Int[], a : Int[], b : Int[]) : Result[] {
    use qs = Qubit[4 * n];
    for c in 0..3 {
        ApplyCircuit(gates, a, b, qs[c * n..c * n + n - 1]);
    }
    let first = BellMeasure(qs[0..n - 1], qs[n..2 * n - 1]);
    let second = BellMeasure(qs[2 * n..3 * n - 1], qs[3 * n..4 * n - 1]);
    first + second
}

/// One copy, and a joint measurement of every Pauli string in `paulis`. For the
/// generators of the state's own stabilizer group these commute and every outcome
/// is deterministic, which is how the learner reads off the signs.
operation MeasureGenerators(n : Int, gates : Int[], a : Int[], b : Int[], paulis : Pauli[][]) : Result[] {
    use qs = Qubit[n];
    ApplyCircuit(gates, a, b, qs);
    mutable out = [];
    for p in paulis {
        out += [Measure(p, qs)];
    }
    out
}
