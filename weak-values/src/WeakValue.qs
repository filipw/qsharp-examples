/// The weak Stern-Gerlach experiment of arXiv:2609.37822 on qubits.
///
/// The spin is one qubit, with |0> = |+z> and |1> = |-z>. The pointer, the atom's momentum along
/// the field gradient, is an n-qubit register holding a wave function on the position grid
///   z_j = j - (N - 1)/2,  j = 0..N - 1,  N = 2^n,
/// in units of the grid step, so momenta are in radians per grid step. The register is
/// big-endian (pointer[0] is the most significant bit of j), which is the order
/// PreparePureStateD uses.
import Std.Canon.ApplyQFT;
import Std.Convert.IntAsDouble;
import Std.Math.E;
import Std.Math.Sqrt;
import Std.StatePreparation.PreparePureStateD;

/// A Gaussian of position standard deviation sigma centred on the middle of the grid:
/// psi(z) proportional to exp(-z^2 / (4 sigma^2)), so |psi|^2 has standard deviation sigma and
/// the momentum distribution has standard deviation 1 / (2 sigma). Normalised.
function GaussianAmplitudes(n : Int, sigma : Double) : Double[] {
    let size = 1 <<< n;
    let centre = IntAsDouble(size - 1) / 2.0;
    mutable amps = [];
    mutable norm = 0.0;
    for j in 0..size - 1 {
        let z = IntAsDouble(j) - centre;
        let a = E()^(-z * z / (4.0 * sigma * sigma));
        amps += [a];
        norm += a * a;
    }
    mutable out = [];
    for a in amps {
        out += [a / Sqrt(norm)];
    }
    out
}

/// The pointer at rest: the Gaussian wave packet, zero mean momentum.
operation PreparePointer(sigma : Double, pointer : Qubit[]) : Unit {
    PreparePureStateD(GaussianAmplitudes(Length(pointer), sigma), pointer);
}

/// The pre-selected spin state |i> = |+x>.
operation PrepareSpin(spin : Qubit) : Unit {
    H(spin);
}

/// The weak Stern-Gerlach coupling exp(i dp sigma_z z): the |+z> component of the spin gets a
/// momentum kick +dp, the |-z> component -dp (Eqs. 3 and 8 of the paper).
///
/// With b_q the bits of j (big-endian), z = j - (N - 1)/2 = -sum_q 2^(n - 2 - q) Z_q, so the
/// coupling is a product of commuting ZZ rotations, one per pointer qubit:
///   exp(i dp Z_s z) = prod_q exp(-i dp 2^(n - 2 - q) Z_s Z_q) = prod_q Rzz(dp 2^(n - 1 - q)).
operation WeakKick(dp : Double, spin : Qubit, pointer : Qubit[]) : Unit {
    let n = Length(pointer);
    for q in 0..n - 1 {
        Rzz(dp * IntAsDouble(1 <<< (n - 1 - q)), spin, pointer[q]);
    }
}

/// Rotates the spin so that a Z measurement distinguishes the post-selected state
/// |f(delta)> = Ry(delta)|-x> (outcome One) from the orthogonal |n(delta)> = Ry(delta)|+x>
/// (outcome Zero): Ry(-delta) takes |f> to |-x>, and H takes |-x> to |1>.
operation RotateToPostSelection(delta : Double, spin : Qubit) : Unit {
    Ry(-delta, spin);
    H(spin);
}

/// Maps the pointer's momentum onto the computational basis, the job the quarter period in the
/// harmonic waveguide does in the experiment. ApplyQFT takes a little-endian input and leaves
/// a big-endian output, so on this big-endian register its adjoint leaves
///   sum_k (1/sqrt(N)) sum_j exp(-2 pi i j k / N) psi_j |k>
/// with k little-endian: a plane wave exp(i p z) with p = 2 pi k0 / N is read as k0, and
/// MeasureInteger returns k mod N. Python maps k >= N/2 to k - N.
operation MomentumToRegister(pointer : Qubit[]) : Unit {
    Adjoint ApplyQFT(pointer);
}

/// The whole experiment up to (not including) the two measurements: pre-selection, pointer,
/// kick, rotation to the post-selection basis, momentum read-out transform.
operation Experiment(sigma : Double, dp : Double, delta : Double, spin : Qubit, pointer : Qubit[]) : Unit {
    PrepareSpin(spin);
    PreparePointer(sigma, pointer);
    WeakKick(dp, spin, pointer);
    RotateToPostSelection(delta, spin);
    MomentumToRegister(pointer);
}
