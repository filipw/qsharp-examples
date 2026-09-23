/// Composite simulation in the Hagan-Wiebe form used by the paper's upper bound
/// (arXiv:2206.06409 Thm 2.1; arXiv:2607.19852 Fact 3).
///
/// Split H = A + B. Treat A and B as the two "terms" of an outer Suzuki formula of
/// order 2k. Every A-factor of duration s is implemented by an inner order-2k
/// product formula on A's terms; every B-factor is implemented by a qDRIFT
/// segment with a fixed number of samples. The base case is the symmetric
/// second-order split A(s/2) B(s) A(s/2); higher orders come from Suzuki's
/// five-fold recursion applied to that base.
import Std.Convert.*;
import Std.Diagnostics.Fact;
import ProductFormula.*;
import QDrift.QDriftEvolve;

/// One outer step of duration `time`.
operation CompositeStep(
    outerOrder : Int,
    innerOrder : Int,
    bigTerms : Pauli[][],
    bigCoeffs : Double[],
    smallTerms : Pauli[][],
    smallCoeffs : Double[],
    smallCdf : Double[],
    time : Double,
    samplesPerSegment : Int,
    qs : Qubit[]
) : Unit {
    if outerOrder == 2 {
        SuzukiStep(innerOrder, bigTerms, bigCoeffs, time / 2.0, qs);
        QDriftEvolve(smallTerms, smallCoeffs, smallCdf, time, samplesPerSegment, qs);
        SuzukiStep(innerOrder, bigTerms, bigCoeffs, time / 2.0, qs);
    } else {
        let u = SuzukiU(outerOrder);
        for slice in [u, u, 1.0 - 4.0 * u, u, u] {
            CompositeStep(
                outerOrder - 2, innerOrder,
                bigTerms, bigCoeffs, smallTerms, smallCoeffs, smallCdf,
                slice * time, samplesPerSegment, qs
            );
        }
    }
}

/// exp(-i (A + B) time) via `steps` outer steps of order `order`, with inner
/// product formulas of the same order on A and `samplesPerSegment` qDRIFT
/// samples in every B-segment (the paper's N_B).
operation CompositeEvolve(
    order : Int,
    bigTerms : Pauli[][],
    bigCoeffs : Double[],
    smallTerms : Pauli[][],
    smallCoeffs : Double[],
    smallCdf : Double[],
    time : Double,
    steps : Int,
    samplesPerSegment : Int,
    qs : Qubit[]
) : Unit {
    Fact(order > 0 and order % 2 == 0, "composite order must be a positive even integer");
    Fact(steps > 0, "steps must be positive");
    Fact(Length(bigTerms) > 0 and Length(smallTerms) > 0,
        "composite needs a non-empty A and B; use TrotterEvolve or QDriftEvolve for the endpoints");
    let dt = time / IntAsDouble(steps);
    for _ in 1..steps {
        CompositeStep(order, order, bigTerms, bigCoeffs, smallTerms, smallCoeffs, smallCdf,
                      dt, samplesPerSegment, qs);
    }
}
