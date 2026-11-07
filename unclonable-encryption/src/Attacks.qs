/// Two splitting attacks on the scheme in Scheme.qs. Neither circuit is from the
/// paper: they are Copenhagen's, to show the bound of Theorem 4.1 from below.
/// The cloner never sees the key; Bob and Charlie get it only after the split.
import Std.Math.*;

/// The diagonal n = (x + y + z) / sqrt(3): measure n . sigma, as a sign. Inside
/// `within`, Rz(-pi/4) then Ry(-theta) turn n onto the Z axis, with
/// cos(theta) = 1 / sqrt(3); `apply` measures Z there, and the rotations are undone.
/// On any eigenstate of X, Y or Z with sign s, the sign comes out s with
/// probability (1 + 1/sqrt(3)) / 2.
operation MeasureDiagonal(q : Qubit) : Bool {
    let theta = ArcCos(1.0 / Sqrt(3.0));
    within {
        Rz(-PI() / 4.0, q);
        Ry(-theta, q);
    } apply {
        M(q) == One
    }
}

/// Bužek and Hillery's symmetric universal cloner, as the network of Bužek,
/// Braunstein, Hillery and Bruß (quant-ph/9703046). The preparation puts the blank
/// and the machine in (2|00> + |01> + |11>) / sqrt(6), with cos(2 theta_1) =
/// 1/sqrt(5), cos(2 theta_2) = sqrt(5)/3 and cos(2 theta_3) = 2/sqrt(5) (Q#'s
/// Ry(phi) turns by phi, so phi = 2 theta). Four CNOTs then copy: the two copies
/// come out on `input` and `blank`, each with fidelity 5/6 for every input state,
/// and the machine is thrown away.
operation UniversalClone(input : Qubit, blank : Qubit, machine : Qubit) : Unit {
    // Preparation
    Ry(ArcCos(1.0 / Sqrt(5.0)), blank);
    CNOT(blank, machine);
    Ry(ArcCos(Sqrt(5.0) / 3.0), machine);
    CNOT(machine, blank);
    Ry(ArcCos(2.0 / Sqrt(5.0)), blank);
    // Copying
    CNOT(input, blank);
    CNOT(input, machine);
    CNOT(blank, input);
    CNOT(machine, input);
}
