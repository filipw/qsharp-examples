import Std.Convert.*;
import Std.Math.*;

// A toy oracle for Simon's problem.
//
// This oracle acts on two registers:
//   - An "input" register of n qubits
//   - An "output" register of n qubits (initially |0...0⟩)
//
// It implements a function f such that for a given secret bitstring b:
//   f(x) = f(x ⊕ b)
//
// One simple (toy) approach is to first copy x into the output register,
// then conditionally (if x[0] is 0) flip the output bits for which b is true.
// This yields:
//   - If x[0] = 1, then f(x) = x.
//   - If x[0] = 0, then f(x) = x ⊕ b.
//
// With the assumption that the secret b is non-zero and (for simplicity)
// that its first entry is true, one can check that:
//   f(x) = f(x ⊕ b)
operation SimonOracle(secret : Bool[], qubits : (Qubit[], Qubit[])) : Unit {
    let (x_register, y_register) = qubits;
    let n = Length(x_register);

    // copy x into y via CNOTs.
    for i in 0..n - 1 {
        CNOT(x_register[i], y_register[i]);
    }

    if secret == [false, size = n] {
        fail ("This program assumes the secret is non-zero");
    }

    if (not secret[0]) {
        fail("This program assumes that secret[0] is true.");
    }

    // To condition on x[0] being |0⟩, we flip it (so |0⟩ becomes |1⟩),
    // then use it as a control.
    X(x_register[0]);
    for i in 0..n - 1 {
        if (secret[i]) {
            CNOT(x_register[0], y_register[i]);
        }
    }
    X(x_register[0]); // restore the state of x[0]
}

// Simon's algorithm
//
// This operation prepares 2n qubits (n for the input and n for the output),
// applies Hadamard to the input register, calls the oracle, applies Hadamard
// again to the input register, and finally measures the input register.
//
// In a full algorithm you would repeat this process (collecting a system of
// linear equations z·b = 0) until you have enough equations to solve for b.
// For demonstration, we simply run one iteration.
operation RunSimonsAlgorithm(secret : Bool[]) : Result[] {
    let n = Length(secret);
    use (input, output) = (Qubit[n], Qubit[n]);

    // apply Hadamard to all qubits in the input register.
    for i in 0..n - 1 {
        H(input[i]);
    }

    // call the oracle f: |x>|0> → |x>|f(x)>
    SimonOracle(secret, (input, output));

    // apply Hadamard to the input register again.
    for i in 0..n - 1 {
        H(input[i]);
    }

    // measure the input register.
    let results = MeasureEachZ(input);

    // reset all qubits.
    ResetAll(input + output);
    results
}

operation Main() : Bool {
    // secret is "110"
    let secret = [true, true, false];
    Message("Running Simon's Algorithm with secret b = 110");

    let measurement = RunSimonsAlgorithm(secret);
    Message($"Measured result from input register: {measurement}");

    let verification = DotProductMod2(secret, measurement);

    if verification {
        Message("The condition b ⋅ z = 0 (mod 2) is satisfied.");
    } else {
        Message("UH-OH! The condition b ⋅ z = 0 (mod 2) is NOT satisfied.");
    }

    verification
}

function DotProductMod2(secret : Bool[], z : Result[]) : Bool {
    let b = BoolArrayAsInt(secret);
    let z = ResultArrayAsInt(z);

    let and_result = b &&& z;

    // the dot product is 0 mod 2 if the number of 1 bits is even.
    HammingWeightI(and_result) % 2 == 0
}