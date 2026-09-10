/// Entry points driven from Python.
///
///   Sim*  -- dump the state vector; Python traces out the dephasing environment, projects
///            onto syndromes and extracts the exact logical channel.
///   Run*  -- shots: syndrome bits followed by the 7-bit transversal readout.
import Std.Diagnostics.DumpMachine;
import Std.Measurement.MResetEachZ;
import Std.Random.DrawRandomDouble;
import Steane.*;

/// The encoded input alone, 7 qubits.
operation SimEncoded(input : Int) : Unit {
    use data = Qubit[7];
    PrepareLogical(input, data);
    DumpMachine();
    ResetAll(data);
}

/// One round up to (not including) the syndrome measurement: transversal Rz(theta) and
/// dephasing of strength p, 14 qubits (data, then environment). Python projects the data onto
/// each of the 8 X-syndromes and applies the correction, which gives every branch exactly.
operation SimRotated(input : Int, theta : Double, p : Double) : Unit {
    use data = Qubit[7];
    use env = Qubit[7];
    PrepareLogical(input, data);
    Rotate(theta, data);
    DephaseDilated(p, data, env);
    DumpMachine();
    ResetAll(data + env);
}

/// One full round: rotation, dephasing, Steane syndrome extraction, correction; dumped after
/// the correction (14 qubits) and the measured syndrome returned.
operation SimOneRound(input : Int, theta : Double, p : Double) : Result[] {
    use data = Qubit[7];
    use env = Qubit[7];
    PrepareLogical(input, data);
    Rotate(theta, data);
    DephaseDilated(p, data, env);
    let s = MeasureXSyndrome(data);
    Correct(s, data);
    DumpMachine();
    ResetAll(data + env);
    s
}

/// Two rounds on |+_L>, +theta then -theta, each with its own environment; dumped before the
/// second syndrome measurement (21 qubits: data, env1, env2) with the first syndrome returned.
operation SimTwoRounds(theta : Double, p : Double) : Result[] {
    use data = Qubit[7];
    use env1 = Qubit[7];
    use env2 = Qubit[7];
    PrepareLogical(2, data);
    Rotate(theta, data);
    DephaseDilated(p, data, env1);
    let s1 = MeasureXSyndrome(data);
    Correct(s1, data);
    Rotate(-theta, data);
    DephaseDilated(p, data, env2);
    DumpMachine();
    ResetAll(data + env1 + env2);
    s1
}

/// One round with sampled dephasing and transversal readout in basis 0 = Z, 1 = X, 2 = Y.
/// Returns the 3 syndrome bits followed by the 7 readout bits.
operation RunOneRound(input : Int, theta : Double, p : Double, basis : Int) : Result[] {
    use data = Qubit[7];
    PrepareLogical(input, data);
    Rotate(theta, data);
    DephaseSampled(p, data);
    let s = MeasureXSyndrome(data);
    Correct(s, data);
    s + MeasureLogical(basis, data)
}

/// As RunOneRound with the flips drawn by DrawRandomBool, for the seeding check.
operation RunOneRoundDrawn(input : Int, theta : Double, p : Double, basis : Int) : Result[] {
    use data = Qubit[7];
    PrepareLogical(input, data);
    Rotate(theta, data);
    DephaseDrawn(p, data);
    let s = MeasureXSyndrome(data);
    Correct(s, data);
    s + MeasureLogical(basis, data)
}

/// Two classical coins and two quantum coins in one shot, for the random-stream check:
/// [u1 < 1/2, u2 < 1/2, M(|+>), M(|+>)] with u1, u2 drawn by DrawRandomDouble.
operation DrawAndMeasure() : Result[] {
    let u1 = DrawRandomDouble(0.0, 1.0);
    let u2 = DrawRandomDouble(0.0, 1.0);
    use qs = Qubit[2];
    for q in qs { H(q); }
    [u1 < 0.5 ? One | Zero, u2 < 0.5 ? One | Zero] + MResetEachZ(qs)
}

/// As RunOneRound, with the syndrome from joint measurements instead of the Steane block.
operation RunOneRoundDirect(input : Int, theta : Double, p : Double, basis : Int) : Result[] {
    use data = Qubit[7];
    PrepareLogical(input, data);
    Rotate(theta, data);
    DephaseSampled(p, data);
    let s = MeasureXSyndromeDirect(data);
    Correct(s, data);
    s + MeasureLogical(basis, data)
}

/// Two rounds (+theta, -theta) on |+_L> with sampled dephasing: s1 (3 bits), s2 (3 bits),
/// then the 7 readout bits.
operation RunTwoRounds(theta : Double, p : Double, basis : Int) : Result[] {
    use data = Qubit[7];
    PrepareLogical(2, data);
    Rotate(theta, data);
    DephaseSampled(p, data);
    let s1 = MeasureXSyndrome(data);
    Correct(s1, data);
    Rotate(-theta, data);
    DephaseSampled(p, data);
    let s2 = MeasureXSyndrome(data);
    Correct(s2, data);
    s1 + s2 + MeasureLogical(basis, data)
}
