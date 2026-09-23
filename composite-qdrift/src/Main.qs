/// Entry points driven from `experiment.py`.
///
/// Each simulator comes in two flavours:
///   Sim*  -- runs on the sparse simulator and dumps the final state vector, so
///            Python can score it against an exact classical reference.
///   Est*  -- identical circuit with no diagnostics, for `qsharp.logical_counts`
///            and the resource estimator.
import Std.Diagnostics.DumpMachine;
import Hamiltonian.*;
import ProductFormula.TrotterEvolve;
import QDrift.QDriftEvolve;
import Composite.CompositeEvolve;

operation SimTrotter(
    nQubits : Int,
    angles : Double[],
    terms : Pauli[][],
    coeffs : Double[],
    time : Double,
    order : Int,
    steps : Int
) : Unit {
    use qs = Qubit[nQubits];
    PrepareInput(angles, qs);
    TrotterEvolve(order, terms, coeffs, time, steps, qs);
    DumpMachine();
    ResetAll(qs);
}

operation SimQDrift(
    nQubits : Int,
    angles : Double[],
    terms : Pauli[][],
    coeffs : Double[],
    time : Double,
    samples : Int
) : Unit {
    use qs = Qubit[nQubits];
    PrepareInput(angles, qs);
    QDriftEvolve(terms, coeffs, ImportanceCdf(coeffs), time, samples, qs);
    DumpMachine();
    ResetAll(qs);
}

operation SimComposite(
    nQubits : Int,
    angles : Double[],
    bigTerms : Pauli[][],
    bigCoeffs : Double[],
    smallTerms : Pauli[][],
    smallCoeffs : Double[],
    time : Double,
    order : Int,
    steps : Int,
    samplesPerSegment : Int
) : Unit {
    use qs = Qubit[nQubits];
    PrepareInput(angles, qs);
    CompositeEvolve(
        order,
        bigTerms,
        bigCoeffs,
        smallTerms,
        smallCoeffs,
        ImportanceCdf(smallCoeffs),
        time,
        steps,
        samplesPerSegment,
        qs
    );
    DumpMachine();
    ResetAll(qs);
}

operation EstTrotter(
    nQubits : Int,
    angles : Double[],
    terms : Pauli[][],
    coeffs : Double[],
    time : Double,
    order : Int,
    steps : Int
) : Unit {
    use qs = Qubit[nQubits];
    PrepareInput(angles, qs);
    TrotterEvolve(order, terms, coeffs, time, steps, qs);
    ResetAll(qs);
}

operation EstQDrift(
    nQubits : Int,
    angles : Double[],
    terms : Pauli[][],
    coeffs : Double[],
    time : Double,
    samples : Int
) : Unit {
    use qs = Qubit[nQubits];
    PrepareInput(angles, qs);
    QDriftEvolve(terms, coeffs, ImportanceCdf(coeffs), time, samples, qs);
    ResetAll(qs);
}

operation EstComposite(
    nQubits : Int,
    angles : Double[],
    bigTerms : Pauli[][],
    bigCoeffs : Double[],
    smallTerms : Pauli[][],
    smallCoeffs : Double[],
    time : Double,
    order : Int,
    steps : Int,
    samplesPerSegment : Int
) : Unit {
    use qs = Qubit[nQubits];
    PrepareInput(angles, qs);
    CompositeEvolve(
        order,
        bigTerms,
        bigCoeffs,
        smallTerms,
        smallCoeffs,
        ImportanceCdf(smallCoeffs),
        time,
        steps,
        samplesPerSegment,
        qs
    );
    ResetAll(qs);
}

/// Dumps just the prepared input state, so Python can check that its classical
/// mirror of `PrepareInput` agrees with Q# before any of the numbers are trusted.
operation SimInputState(nQubits : Int, angles : Double[]) : Unit {
    use qs = Qubit[nQubits];
    PrepareInput(angles, qs);
    DumpMachine();
    ResetAll(qs);
}

/// State preparation only, so the rotation count it costs can be subtracted from
/// every simulator's total and the x-axis measures simulation work alone.
operation EstInputState(nQubits : Int, angles : Double[]) : Unit {
    use qs = Qubit[nQubits];
    PrepareInput(angles, qs);
    ResetAll(qs);
}
