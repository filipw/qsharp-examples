/// The triangular colour code of arXiv:2610.24077 at distances 3 and 5, as a quantum memory.
///
/// Data qubits sit on the vertices of a honeycomb cut into a triangle; every face (plaquette)
/// carries one X check and one Z check on the same qubits, so the code is self-dual and the
/// X and Z decoding problems are the same. d = 3 is Steane's [[7,1,3]] code (three faces of
/// weight 4); d = 5 is [[19,1,5]] (six boundary faces of weight 4, three hexagons of weight 6).
///
/// Checks are measured with one ancilla, reused after a mid-circuit MResetZ. The noise is
/// part of the circuit: before each round of checks, every data qubit suffers X, Y or Z with
/// probability p/3 each (code-capacity depolarizing noise); the checks themselves are perfect.
import Std.Arrays.Subarray;
import Std.Random.DrawRandomDouble;

/// The faces of the distance-d code, as data-qubit indices. Laid out on a triangular grid
/// (r, c), 0 <= c <= r <= 3(d - 1)/2: a point with (r + c) mod 3 = 1 is a face, every other
/// point is a data qubit, and a face acts on its (up to six) grid neighbours.
function Faces(d : Int) : Int[][] {
    if d == 3 {
        [[0, 1, 2, 3], [1, 3, 5, 6], [2, 3, 4, 5]]
    } elif d == 5 {
        [
            [0, 1, 2, 3],
            [1, 3, 5, 6],
            [2, 3, 4, 5, 7, 8],
            [4, 7, 10, 11],
            [5, 6, 8, 9, 12, 13],
            [7, 8, 11, 12, 15, 16],
            [9, 13, 17, 18],
            [10, 11, 14, 15],
            [12, 13, 16, 17]
        ]
    } else {
        fail "Only d = 3 and d = 5 are built here."
    }
}

/// Data qubits of the distance-d triangular colour code: (3d^2 + 1)/4.
function DataQubitCount(d : Int) : Int {
    (3 * d * d + 1) / 4
}

/// One side of the triangle: d qubits that carry a logical Z (and, the code being
/// self-dual, a logical X). It meets every face in an even number of qubits.
function LogicalSupport(d : Int) : Int[] {
    if d == 3 { [0, 2, 4] } else { [0, 2, 4, 10, 14] }
}

/// The noise: X, Y or Z, each with probability p/3.
operation ApplyDepolarizing(p : Double, q : Qubit) : Unit {
    let r = DrawRandomDouble(0.0, 1.0);
    if r < p / 3.0 {
        X(q);
    } elif r < 2.0 * p / 3.0 {
        Y(q);
    } elif r < p {
        Z(q);
    }
}

/// The Z check of one face: the parity of its data qubits, copied onto the ancilla by CNOTs.
/// One means the face has an odd number of X (or Y) errors.
operation MeasureZCheck(face : Int[], data : Qubit[], ancilla : Qubit) : Result {
    for i in face {
        CNOT(data[i], ancilla);
    }
    MResetZ(ancilla)
}

/// The X check of one face: the ancilla in |+> controls an X on each of its data qubits,
/// which kicks the product of their X eigenvalues back onto it. One means an odd number of
/// Z (or Y) errors.
operation MeasureXCheck(face : Int[], data : Qubit[], ancilla : Qubit) : Result {
    H(ancilla);
    for i in face {
        CNOT(ancilla, data[i]);
    }
    H(ancilla);
    MResetZ(ancilla)
}

/// Every check of one type, as a bitmask: bit i is set when face i reads One.
/// `basis` is the check's Pauli: PauliZ checks catch X errors, PauliX checks catch Z errors.
operation MeasureChecks(basis : Pauli, faces : Int[][], data : Qubit[], ancilla : Qubit) : Int {
    mutable syndrome = 0;
    for i in 0..Length(faces) - 1 {
        let r = basis == PauliZ ? MeasureZCheck(faces[i], data, ancilla) | MeasureXCheck(faces[i], data, ancilla);
        if r == One {
            syndrome += 1 <<< i;
        }
    }
    syndrome
}

/// Logical |0> (basis PauliZ) or |+> (basis PauliX), as the paper prepares it. |0...0> is
/// already +1 for every Z check and for logical Z; measuring the X checks once projects it
/// onto the code, with random X-check signs that do not touch the logical state. |+> is the
/// same with the roles of X and Z swapped.
operation PrepareLogical(basis : Pauli, faces : Int[][], data : Qubit[], ancilla : Qubit) : Unit {
    if basis == PauliX {
        for q in data {
            H(q);
        }
    }
    let _ = MeasureChecks(basis == PauliZ ? PauliX | PauliZ, faces, data, ancilla);
}

/// Transversal read-out in `basis`, as a bitmask over the data qubits (and resets them).
operation ReadOut(basis : Pauli, data : Qubit[]) : Int {
    mutable bits = 0;
    for i in 0..Length(data) - 1 {
        let r = basis == PauliZ ? MResetZ(data[i]) | MResetX(data[i]);
        if r == One {
            bits += 1 <<< i;
        }
    }
    bits
}

/// One memory shot: prepare logical |0> or |+>, run `rounds` rounds of noise followed by
/// both types of checks, then read every data qubit out in the same basis. Returns the
/// syndromes of the checks that catch the errors this basis is sensitive to (Z checks for
/// |0>, X checks for |+>), one bitmask per round, then the read-out bitmask.
operation MemoryShot(d : Int, basis : Pauli, rounds : Int, p : Double) : Int[] {
    let faces = Faces(d);
    let other = basis == PauliZ ? PauliX | PauliZ;
    use data = Qubit[DataQubitCount(d)];
    use ancilla = Qubit();
    PrepareLogical(basis, faces, data, ancilla);
    mutable record = [];
    for _ in 1..rounds {
        for q in data {
            ApplyDepolarizing(p, q);
        }
        record += [MeasureChecks(basis, faces, data, ancilla)];
        let _ = MeasureChecks(other, faces, data, ancilla);
    }
    record + [ReadOut(basis, data)]
}

/// The Z sector of a |0> memory on its own: no preparation and no X checks, so the data stay
/// in a computational basis state. Under Pauli noise and perfect checks the X checks and the
/// code-state superposition change neither the Z checks nor the syndrome and logical parity
/// of the Z read-out, so this gives the same decoded statistics as MemoryShot(d, PauliZ, ...),
/// at a fraction of the simulator's cost.
operation ZSectorShot(d : Int, rounds : Int, p : Double) : Int[] {
    let faces = Faces(d);
    use data = Qubit[DataQubitCount(d)];
    use ancilla = Qubit();
    mutable record = [];
    for _ in 1..rounds {
        for q in data {
            ApplyDepolarizing(p, q);
        }
        record += [MeasureChecks(PauliZ, faces, data, ancilla)];
    }
    record + [ReadOut(PauliZ, data)]
}
