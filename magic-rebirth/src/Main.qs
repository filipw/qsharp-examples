/// Entry points driven from Python.
///
///   Sim*      -- dump the joint system + environment state vector; Python traces out
///                whichever half it wants and evaluates magic and entanglement exactly.
///   Measure*  -- shot-based readout: the dual witness from two Pauli settings, and the
///                parity-syndrome extraction of a single magic qubit.
import Std.Diagnostics.DumpMachine;
import Std.Measurement.MResetEachZ;
import Std.Measurement.MResetX;
import Std.Measurement.MResetY;
import Std.Measurement.MResetZ;
import Damping.*;

/// The amplitude-damped cat, 2n qubits: system first, environment second.
operation SimDampedCat(n : Int, alpha : Double, gamma : Double) : Unit {
    use sys = Qubit[n];
    use env = Qubit[n];
    PrepareCat(alpha, sys);
    DampAll(gamma, sys, env);
    DumpMachine();
    ResetAll(sys + env);
}

/// Amplitude damping followed by a phase flip of strength p on every qubit, 3n qubits:
/// system, damping environment, dephasing environment.
operation SimDampedDephasedCat(n : Int, alpha : Double, gamma : Double, p : Double) : Unit {
    use sys = Qubit[n];
    use env = Qubit[n];
    use deph = Qubit[n];
    PrepareCat(alpha, sys);
    DampAll(gamma, sys, env);
    for i in 0..n - 1 {
        PhaseFlipDilate(p, sys[i], deph[i]);
    }
    DumpMachine();
    ResetAll(sys + env + deph);
}

/// An arbitrary input prepared from gate codes, then damped. Used for the stabilizer
/// inputs (Bell states, GHZ_3, |00+>) of the generator/insulator classification.
operation SimDampedCircuit(n : Int, gates : Int[], a : Int[], b : Int[], gamma : Double) : Unit {
    use sys = Qubit[n];
    use env = Qubit[n];
    ApplyCircuit(gates, a, b, sys);
    DampAll(gamma, sys, env);
    DumpMachine();
    ResetAll(sys + env);
}

/// Dumps the input state alone, so Python can check its mirror of PrepareCat.
operation SimCat(n : Int, alpha : Double) : Unit {
    use sys = Qubit[n];
    PrepareCat(alpha, sys);
    DumpMachine();
    ResetAll(sys);
}

/// Witness readout of the damped cat: all system qubits in the Z basis (basis = 0), which
/// gives the endpoint populations P_0 and P_n, or in the X basis (basis = 1), whose parity
/// gives <X^n> = 2c. Tr(W rho) = 1 + 2c - 2 min(P_0, P_n) > 1 certifies magic (Eq. 60-64).
operation MeasureCat(n : Int, alpha : Double, gamma : Double, basis : Int) : Result[] {
    use sys = Qubit[n];
    use env = Qubit[n];
    PrepareCat(alpha, sys);
    DampAll(gamma, sys, env);
    if basis == 1 {
        for q in sys {
            H(q);
        }
    }
    let r = MResetEachZ(sys);
    ResetAll(env);
    r
}

/// Parity-syndrome extraction (Sec. V): measure the n - 1 stabilizers Z_i Z_{i+1}, decode
/// with a CNOT cascade from qubit 0, and read qubit 0 in the Z (0), X (1) or Y (2) basis.
/// Returns the n - 1 syndrome bits followed by the readout; Python keeps the shots whose
/// syndrome is trivial (all Zero), which happens with probability P_0 + P_n.
operation ExtractMagic(n : Int, alpha : Double, gamma : Double, basis : Int) : Result[] {
    use sys = Qubit[n];
    use env = Qubit[n];
    PrepareCat(alpha, sys);
    DampAll(gamma, sys, env);
    mutable syndrome = [];
    for i in 0..n - 2 {
        syndrome += [Measure([PauliZ, PauliZ], [sys[i], sys[i + 1]])];
    }
    for i in 1..n - 1 {
        CNOT(sys[0], sys[i]);
    }
    mutable r = Zero;
    if basis == 0 {
        r = MResetZ(sys[0]);
    } elif basis == 1 {
        r = MResetX(sys[0]);
    } else {
        r = MResetY(sys[0]);
    }
    ResetAll(sys + env);
    syndrome + [r]
}
