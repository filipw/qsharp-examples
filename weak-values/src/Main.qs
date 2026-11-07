/// Entry points driven from Python.
///
///   Sim*  -- dump the state vector; Python reads the exact probabilities off it.
///   Run*  -- shots: the spin outcome (One = post-selected on |f>) and the pointer momentum.
///
/// The spin is allocated first, so it is the most significant qubit of every dump.
import Std.Convert.IntAsDouble;
import Std.Diagnostics.DumpMachine;
import Std.Math.PI;
import Std.Measurement.MeasureInteger;
import WeakValue.*;

/// The pointer alone, before anything happens to it.
operation SimPointer(n : Int, sigma : Double) : Unit {
    use pointer = Qubit[n];
    PreparePointer(sigma, pointer);
    DumpMachine();
    ResetAll(pointer);
}

/// Spin and pointer right after the weak kick, in position representation.
operation SimKicked(n : Int, sigma : Double, dp : Double) : Unit {
    use spin = Qubit();
    use pointer = Qubit[n];
    PrepareSpin(spin);
    PreparePointer(sigma, pointer);
    WeakKick(dp, spin, pointer);
    DumpMachine();
    ResetAll([spin] + pointer);
}

/// The experiment just before the two measurements. The dump holds the joint amplitude of every
/// (spin outcome, momentum) pair.
operation SimExperiment(n : Int, sigma : Double, dp : Double, delta : Double) : Unit {
    use spin = Qubit();
    use pointer = Qubit[n];
    Experiment(sigma, dp, delta, spin, pointer);
    DumpMachine();
    ResetAll([spin] + pointer);
}

/// One shot of the experiment: (spin outcome, momentum k mod N).
operation RunExperiment(n : Int, sigma : Double, dp : Double, delta : Double) : (Result, Int) {
    use spin = Qubit();
    use pointer = Qubit[n];
    Experiment(sigma, dp, delta, spin, pointer);
    let s = MResetZ(spin);
    let k = MeasureInteger(pointer);
    (s, k)
}

/// A plane wave exp(2 pi i k0 j / N) on the big-endian pointer, then the read-out transform:
/// checks that a momentum of 2 pi k0 / N is read as k0.
operation RunPlaneWave(n : Int, k0 : Int) : Int {
    use pointer = Qubit[n];
    for q in 0..n - 1 {
        H(pointer[q]);
        R1(2.0 * PI() * IntAsDouble(k0 * (1 <<< (n - 1 - q))) / IntAsDouble(1 <<< n), pointer[q]);
    }
    MomentumToRegister(pointer);
    MeasureInteger(pointer)
}
