/// Campbell's qDRIFT: a randomly compiled first-order channel whose cost depends
/// on the 1-norm lambda = sum_j |coeffs[j]| but *not* on the number of terms.
import Std.Convert.*;
import Std.Random.*;
import Hamiltonian.*;

/// Draws one term index from p_j = |coeffs[j]| / lambda.
operation SampleTermIndex(cdf : Double[]) : Int {
    SearchCdf(cdf, DrawRandomDouble(0.0, 1.0))
}

/// One qDRIFT channel realization for exp(-i H time) using `samples` rotations.
///
/// Each sampled term is applied for the *same* angle tau = lambda * time / samples,
/// carrying only the sign of its coefficient; that is what makes the gate count
/// independent of how many small terms the Hamiltonian has.
operation QDriftEvolve(
    terms : Pauli[][],
    coeffs : Double[],
    cdf : Double[],
    time : Double,
    samples : Int,
    qs : Qubit[]
) : Unit {
    if samples > 0 {
        let tau = OneNorm(coeffs) * time / IntAsDouble(samples);
        for _ in 1..samples {
            let j = SampleTermIndex(cdf);
            let sign = coeffs[j] >= 0.0 ? 1.0 | -1.0;
            ApplyTerm(terms[j], sign, tau, qs);
        }
    }
}
