/// The zero-error ordered-search algorithm of arXiv:2609.34293 (Carolan, Childs, Section 3) on a
/// register of m qubits. The discretisation is Copenhagen's, not the paper's: the paper works on
/// L2(R) and gives none.
///
/// The register holds a ring of M = 2^m sites, K per cell of a ring of 2N cells (M = 2NK). Site j
/// sits at x_j = (j + 1/2)/K - N, in cells from the middle of the list, so the paper's x = 0 is
/// the boundary between sites M/2 - 1 and M/2 and no site is on it. The n-element list fills the
/// cells around x = 0; the rest of that half-ring is padding whose answers are known (0 to the
/// left of the list, 1 to the right), and the other half-ring is a copy of the list shifted by half
/// the ring with its bits complemented: the string x x-bar of the translation-invariant
/// algorithms of Farhi, Goldstone, Gutmann and Sipser. The register is big-endian
/// (register[0] is the most significant bit of j), the order PreparePureStateD and
/// DumpMachine use.
import Std.Canon.ApplyControlledOnBitString;
import Std.Canon.ApplyQFT;
import Std.Convert.IntAsBoolArray;
import Std.Arrays.Reversed;
import Std.StatePreparation.PreparePureStateD;

/// A phase of -1 on every site j >= c (0 <= c < 2^m). Reading from the top bit down, j >= c
/// exactly when j equals c, or j first differs from c by a 1 where c has a 0: one
/// multi-controlled Z for each 0 bit of c, and one for j = c.
operation FlipAtOrAbove(c : Int, register : Qubit[]) : Unit is Adj + Ctl {
    let m = Length(register);
    let bits = Reversed(IntAsBoolArray(c, m));
    for k in 0..m - 1 {
        if not bits[k] {
            if k == 0 {
                Z(register[0]);
            } else {
                ApplyControlledOnBitString(bits[...k - 1], Z, register[...k - 1], register[k]);
            }
        }
    }
    // j = c itself. With m = 1 there are no controls.
    let last = register[m - 1];
    within {
        if not bits[m - 1] {
            X(last);
        }
    } apply {
        if m == 1 {
            Z(last);
        } else {
            ApplyControlledOnBitString(bits[...m - 2], Z, register[...m - 2], last);
        }
    }
}

/// The comparison oracle O_s, (O_s psi)(x) = (-1)^[x >= s] psi(x), on the ring. `step` is the
/// first site of cell s, so on the list's half-ring a site is flipped when its cell is at or
/// past s: one query to the list's bit x_i = [i >= s] (the padding's bits are known). On the
/// copy's half-ring the bits are the complements, so the step repeats half a ring away, and
/// the sites between the two steps are flipped: up to a global sign, which nothing measures.
operation ComparisonOracle(step : Int, register : Qubit[]) : Unit is Adj + Ctl {
    let half = 1 <<< (Length(register) - 1);
    FlipAtOrAbove(step, register);
    FlipAtOrAbove((step + half) % (2 * half), register);
}

/// O_0 in the Fourier basis: the sign function. With the paper's F,
/// (F psi)(w) = (2 pi)^(-1/2) int e^(-i w x) psi(x) dx, Adjoint ApplyQFT on this big-endian
/// register leaves the frequency k little-endian (register[m - 1] its most significant bit),
/// k standing for w = 2 pi k / M per site and k >= M/2 for the negative frequencies. One Z on
/// that bit negates them, and ApplyQFT goes back to position.
operation SignInFourierBasis(register : Qubit[]) : Unit is Adj + Ctl {
    within {
        Adjoint ApplyQFT(register);
    } apply {
        Z(register[Length(register) - 1]);
    }
}

/// One query of the algorithm: U_s = F^dagger O_0 F O_s.
operation SearchStep(step : Int, register : Qubit[]) : Unit is Adj + Ctl {
    ComparisonOracle(step, register);
    SignInFourierBasis(register);
}

/// The paper's initial state psi_in(x) = |x|^(-1/2) h(ln|x|), a Gaussian h in log-position,
/// centred on the middle of the list, from amplitudes computed in TypeScript. It is even about
/// x = 0, the boundary between the halves of the register, so PreparePureStateD prepares only
/// the right half, on the low qubits (site M/2 + d for d = 0, 1, ...): half the amplitudes.
/// A Hadamard on the top qubit then makes two copies, and flipping the low qubits when the top
/// one is 0 sends site M/2 + d to its mirror image M/2 - 1 - d.
operation PrepareInitialState(halfAmplitudes : Double[], register : Qubit[]) : Unit {
    let top = register[0];
    PreparePureStateD(halfAmplitudes, register[1...]);
    H(top);
    within {
        X(top);
    } apply {
        for q in register[1...] {
            CNOT(top, q);
        }
    }
}
