import Std.Diagnostics.*;

// to see the effects, run with a larger sample size e.g. 1000 shots
operation Main() : Result[] {
    // create a GHZ state (cat state) with 8 qubits
    use qubits = Qubit[8];

    // apply Hadamard gate to the first qubit
    H(qubits[0]);

    // apply CNOT gates to create entanglement
    for qubit in qubits[1..Length(qubits) - 1] {
        CNOT(qubits[0], qubit);
    }

    // configure noise model
    // parameters represent probabilities of applying X, Y, and Z gates
    ConfigurePauliNoise(0.01, 0.01, 0.01);
    // same: DepolarizingNoise(0.03)

    // configure loss
    // parameter represents probability of losing a qubit
    // the returned Result is then "Loss"
    ConfigureQubitLoss(0.01);

    // dump the current state of the machine
    DumpMachine();

    // return measurement results
    MResetEachZ(qubits)
}