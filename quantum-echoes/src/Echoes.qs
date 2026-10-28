import Std.Math.PI;
import Std.Arrays.Filtered;
import Std.Canon.ApplyPauli;

// Quantum echoes (Google Quantum AI, arXiv:2506.10191): out-of-time-order
// correlators of a random circuit, measured by running the circuit forward,
// flipping one qubit, and running it backward.

/// The conditional phase of the paper's iSWAP-like gate, in radians.
function ConditionalPhase() : Double {
    0.35
}

/// The paper's random single-qubit gate: a rotation by θ about an axis in the
/// XY plane at angle φ from X, exp(−iθ/2 (cos φ X + sin φ Y)). θ and φ are in
/// units of π: θ ∈ {1/4, 1/2, 3/4}, φ ∈ [−1, 1).
operation XYRotation(theta : Double, phi : Double, q : Qubit) : Unit is Adj {
    within {
        Rz(-phi * PI(), q);
    } apply {
        Rx(theta * PI(), q);
    }
}

/// The paper's iSWAP-like gate: an iSWAP (|01⟩ → i|10⟩, |10⟩ → i|01⟩), then a
/// conditional phase e^(−iϕ) on |11⟩. Both commute with SWAP, so together they
/// are a SWAP after diag(1, i, i, e^(−iϕ)), which up to a global phase is
/// exp(iϕ/4 (Z ⊗ 1 + 1 ⊗ Z)) · exp(−i(π + ϕ)/4 Z ⊗ Z).
operation ISwapLike(a : Qubit, b : Qubit) : Unit is Adj {
    let phi = ConditionalPhase();
    Rz(-phi / 2.0, a);
    Rz(-phi / 2.0, b);
    Rzz((PI() + phi) / 2.0, a, b);
    SWAP(a, b);
}

/// The coupler patterns of a rows × cols grid, qubit r · cols + c: horizontal
/// pairs from even columns, from odd columns, then vertical pairs from even
/// rows, from odd rows. Empty patterns (the odd rows of a two-row ladder) are
/// left out. Cycle k switches on pattern k mod (number of patterns).
function CouplerPatterns(rows : Int, cols : Int) : (Int, Int)[][] {
    mutable patterns = [];
    for start in 0..1 {
        mutable pairs = [];
        for r in 0..rows - 1 {
            for c in start..2..cols - 2 {
                set pairs += [(r * cols + c, r * cols + c + 1)];
            }
        }
        set patterns += [pairs];
    }
    for start in 0..1 {
        mutable pairs = [];
        for r in start..2..rows - 2 {
            for c in 0..cols - 1 {
                set pairs += [(r * cols + c, (r + 1) * cols + c)];
            }
        }
        set patterns += [pairs];
    }
    Filtered(pairs -> Length(pairs) > 0, patterns)
}

/// One cycle: a random XY-plane rotation on every qubit, then iSWAP-like gates
/// on one coupler pattern. Gate parameters are flattened cycle by cycle.
operation ApplyCycle(qs : Qubit[], pairs : (Int, Int)[], cycle : Int, thetas : Double[], phis : Double[]) : Unit is Adj {
    let n = Length(qs);
    for i in 0..n - 1 {
        XYRotation(thetas[cycle * n + i], phis[cycle * n + i], qs[i]);
    }
    for (a, b) in pairs {
        ISwapLike(qs[a], qs[b]);
    }
}

/// The random circuit U: the given cycles, in order.
operation ApplyCycles(qs : Qubit[], rows : Int, cols : Int, cycles : Range, thetas : Double[], phis : Double[]) : Unit is Adj {
    let patterns = CouplerPatterns(rows, cols);
    for cycle in cycles {
        ApplyCycle(qs, patterns[cycle % Length(patterns)], cycle, thetas, phis);
    }
}

/// U with a Pauli string inserted before cycle `at`. In an echo the adjoint
/// applies the same string at the mirrored point of U†.
operation ApplyCyclesWithPaulis(qs : Qubit[], rows : Int, cols : Int, cycles : Int, thetas : Double[], phis : Double[], at : Int, paulis : Pauli[]) : Unit is Adj {
    ApplyCycles(qs, rows, cols, 0..at - 1, thetas, phis);
    ApplyPauli(paulis, qs);
    ApplyCycles(qs, rows, cols, at..cycles - 1, thetas, phis);
}

/// B(t) = U† B U, the butterfly B = X on qubit qb in the Heisenberg picture:
/// run U, flip qb, undo U.
operation ApplyEchoedButterfly(u : (Qubit[] => Unit is Adj), qs : Qubit[], qb : Int) : Unit is Adj {
    within {
        u(qs);
    } apply {
        X(qs[qb]);
    }
}

/// The nested echo U_k = B(t) [M B(t)]^(k−1), with M = Z on qubit qm. From
/// |0…0⟩, ⟨Z⟩ of qm afterwards is C^(2k): k = 1 is the OTOC, k = 2 OTOC(2).
operation ApplyNestedEcho(u : (Qubit[] => Unit is Adj), qs : Qubit[], qb : Int, qm : Int, k : Int) : Unit {
    ApplyEchoedButterfly(u, qs, qb);
    for _ in 2..k {
        Z(qs[qm]);
        ApplyEchoedButterfly(u, qs, qb);
    }
}

/// One shot of C^(2k) as the processor takes it: the nested echo of a t-cycle
/// circuit on |0…0⟩, then qm measured. C^(2k) = P(Zero) − P(One).
operation MeasureOtoc(rows : Int, cols : Int, cycles : Int, thetas : Double[], phis : Double[], qb : Int, qm : Int, k : Int) : Result {
    use qs = Qubit[rows * cols];
    ApplyNestedEcho(ApplyCycles(_, rows, cols, 0..cycles - 1, thetas, phis), qs, qb, qm, k);
    let result = M(qs[qm]);
    ResetAll(qs);
    result
}
