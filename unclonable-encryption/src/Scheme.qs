/// Unclonable encryption of one bit: Ananth, Ragavan and Sahai, arXiv:2607.21551 v3,
/// §3 and its explicit instantiation in §3.1. The key is a uniformly random
/// non-identity Pauli string on n qubits, and the bit m is encrypted as a random
/// (-1)^m eigenstate of it. Written for Copenhagen; the paper has no code.
import Std.Arrays.*;
import Std.Logical.*;
import Std.Random.*;

/// Gen (§3): a Pauli string drawn uniformly from the 4^n - 1 that are not all I.
/// Position j is P_{x_j, z_j}: I for (0, 0), Z for (0, 1), X for (1, 0) and Y for
/// (1, 1). Drawing each position from the four and redrawing the all-I string is
/// uniform on the rest.
operation KeyGen(n : Int) : Pauli[] {
    mutable key = [];
    repeat {
        key = DrawMany(() => [PauliI, PauliX, PauliY, PauliZ][DrawRandomInt(0, 3)], n, ());
    } until Any(p -> p != PauliI, key);
    key
}

/// Enc (§3.1, Proposition 3.2): prepare U_k|r>, where r is uniform subject to the
/// parity of r over the key's non-identity positions being m. Qubit j starts in
/// |r_j>, and U_j (FromZ) makes it the (-1)^(r_j) eigenstate of P_j. Positions
/// where the key is I keep a random |r_j>: averaged over r, they are maximally
/// mixed.
operation Encrypt(key : Pauli[], m : Bool, qubits : Qubit[]) : Unit {
    let n = Length(key);
    mutable r = DrawMany(() => DrawRandomInt(0, 1) == 1, n, ());
    // j*, the first non-identity position, takes whatever bit makes the parity m.
    let first = IndexOf(p -> p != PauliI, key);
    mutable parity = false;
    for j in 0..n - 1 {
        if key[j] != PauliI {
            parity = Xor(parity, r[j]);
        }
    }
    r w/= first <- Xor(r[first], Xor(parity, m));
    for j in 0..n - 1 {
        if r[j] {
            X(qubits[j]);
        }
        FromZ(key[j], qubits[j]);
    }
}

/// U_j of §3.1: turns the Z axis onto P's, so |0> and |1> become P's +1 and -1
/// eigenstates. Nothing for I and Z, H for X, H then S for Y. Its adjoint turns
/// P back onto Z, so that measuring Z reads P.
operation FromZ(p : Pauli, q : Qubit) : Unit is Adj {
    if p == PauliX {
        H(q);
    } elif p == PauliY {
        H(q);
        S(q);
    }
}

/// Dec (§3): measure the key's Pauli string. Its eigenvalue +1 (Zero) decrypts to
/// 0 and -1 (One) to 1. One joint measurement; the paper's per-qubit form (§3.1)
/// measures each P_j and multiplies the signs, with the same outcome.
operation Decrypt(key : Pauli[], qubits : Qubit[]) : Bool {
    Measure(key, qubits) == One
}
