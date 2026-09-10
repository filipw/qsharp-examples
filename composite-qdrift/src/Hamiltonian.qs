/// Shared plumbing for a Hamiltonian given as a weighted sum of Pauli strings,
///     H = sum_j coeffs[j] * terms[j]
/// where terms[j] is a full-width Pauli string (identities included).
import Std.Math.*;
import Std.Diagnostics.Fact;

/// Applies exp(-i * coeff * time * P) for the Pauli string P.
///
/// Q#'s `Exp` implements exp(+i * theta * P), so the angle is negated here to
/// match the physics convention exp(-iHt) used everywhere else in this demo.
operation ApplyTerm(term : Pauli[], coeff : Double, time : Double, qs : Qubit[]) : Unit is Adj + Ctl {
    Exp(term, -coeff * time, qs);
}

/// The 1-norm lambda = sum_j |coeffs[j]|, which is the quantity qDRIFT pays for.
function OneNorm(coeffs : Double[]) : Double {
    mutable total = 0.0;
    for c in coeffs {
        total += AbsD(c);
    }
    total
}

/// Cumulative distribution of the qDRIFT importance distribution p_j = |coeffs[j]| / lambda.
function ImportanceCdf(coeffs : Double[]) : Double[] {
    let n = Length(coeffs);
    Fact(n > 0, "ImportanceCdf needs at least one term");
    let lambda = OneNorm(coeffs);
    Fact(lambda > 0.0, "ImportanceCdf needs a non-zero 1-norm");
    mutable cdf = [0.0, size = n];
    mutable acc = 0.0;
    for j in 0..n - 1 {
        acc += AbsD(coeffs[j]) / lambda;
        cdf[j] = acc;
    }
    // Guard against the last bucket falling a few ulps short of 1.
    cdf[n - 1] = 1.0;
    cdf
}

/// Smallest index j with u <= cdf[j]. Binary search keeps qDRIFT sampling
/// O(log #terms) per gate, which matters when there are hundreds of small terms.
function SearchCdf(cdf : Double[], u : Double) : Int {
    mutable lo = 0;
    mutable hi = Length(cdf) - 1;
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if u <= cdf[mid] {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    lo
}

/// Reproducible, mildly entangled input state. `angles` has 2 * Length(qs) entries.
operation PrepareInput(angles : Double[], qs : Qubit[]) : Unit is Adj + Ctl {
    let n = Length(qs);
    Fact(Length(angles) >= 2 * n, "PrepareInput needs 2 * Length(qs) angles");
    for i in 0..n - 1 {
        Ry(angles[i], qs[i]);
    }
    for i in 0..n - 2 {
        CNOT(qs[i], qs[i + 1]);
    }
    for i in 0..n - 1 {
        Rz(angles[n + i], qs[i]);
    }
}
