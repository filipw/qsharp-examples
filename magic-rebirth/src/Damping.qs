/// State preparation and the noise channels, written as unitaries on the system plus
/// an explicit environment register (Stinespring dilations). Tracing the environment
/// out in Python gives the channel's output exactly; keeping it gives the paper's
/// complementary channel for free, which is where the gamma <-> 1 - gamma mirror lives.
import Std.Math.ArcCos;
import Std.Math.ArcSin;
import Std.Math.Sqrt;

/// alpha |0^n> + beta |1^n> with beta = sqrt(1 - alpha^2), alpha in (0, 1].
operation PrepareCat(alpha : Double, qs : Qubit[]) : Unit {
    Ry(2.0 * ArcCos(alpha), qs[0]);
    for i in 1..Length(qs) - 1 {
        CNOT(qs[0], qs[i]);
    }
}

/// Canonical dilation of amplitude damping of strength gamma (Eq. 83 of the paper):
///   |0>|0>_E -> |0>|0>_E
///   |1>|0>_E -> sqrt(1 - gamma) |1>|0>_E + sqrt(gamma) |0>|1>_E
/// so the environment qubit records whether the excitation leaked. Its own reduced
/// state is amplitude damping of the input at strength 1 - gamma (Appendix C).
operation AmplitudeDamp(gamma : Double, q : Qubit, e : Qubit) : Unit {
    Controlled Ry([q], (2.0 * ArcSin(Sqrt(gamma)), e));
    CNOT(e, q);
}

/// Dilation of the phase flip D_p(rho) = (1 - p) rho + p Z rho Z:
///   |psi>|0>_E -> sqrt(1 - p) |psi>|0>_E + sqrt(p) Z|psi> |1>_E.
/// Multiplies every coherence of q by 1 - 2p and leaves its populations alone.
operation PhaseFlipDilate(p : Double, q : Qubit, e : Qubit) : Unit {
    Ry(2.0 * ArcSin(Sqrt(p)), e);
    CZ(e, q);
}

/// Local amplitude damping of every qubit in `qs`, one fresh environment qubit each.
operation DampAll(gamma : Double, qs : Qubit[], env : Qubit[]) : Unit {
    for i in 0..Length(qs) - 1 {
        AmplitudeDamp(gamma, qs[i], env[i]);
    }
}

/// Gate codes for arbitrary (stabilizer) inputs: 0 = H(a), 1 = S(a), 2 = CNOT(a, b), 3 = X(a).
operation ApplyCircuit(gates : Int[], a : Int[], b : Int[], qs : Qubit[]) : Unit {
    for i in 0..Length(gates) - 1 {
        if gates[i] == 0 {
            H(qs[a[i]]);
        } elif gates[i] == 1 {
            S(qs[a[i]]);
        } elif gates[i] == 2 {
            CNOT(qs[a[i]], qs[b[i]]);
        } else {
            X(qs[a[i]]);
        }
    }
}
