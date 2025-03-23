/// # Hidden Shift Algorithm
///
/// The Hidden Shift Problem is one of the known problems whose quantum algorithm
/// solution shows exponential speedup over classical computing. Part of the
/// advantage lies on the ability to perform Fourier transforms efficiently. This
/// can be used to extract correlations between certain functions, as demonstrated here.
///
/// Let f and g be two functions {0,1}^N -> {0,1} which are the same
/// up to a hidden bit string s:
///
/// g(x) = f(x ⨁ s), for all x in {0,1}^N
///
/// The implementation considers the following (so-called "bent") functions:
///
/// f(x) = Σ_i x_(2i) x_(2i+1)
///
/// where x_i is the i-th bit of x and i runs from 0 to N/2 - 1.
///
/// Algorithm steps:
/// 1. Prepare the state |0⟩^N
/// 2. Apply Hadamard gates to create superposition of all inputs |x⟩
/// 3. Compute the shifted function g(x) = f(x ⨁ s) into the phase
/// 4. Apply Hadamard gates (Fourier transform)
/// 5. Query oracle f into the phase
/// 6. Apply Hadamard gates (inverse Fourier transform) to get state |s⟩
/// 7. Measure to get the hidden shift s
///
/// based on: 
///  - "Quantum Algorithms for some Hidden Shift Problems", Wim van Dam and Sean Hallgren and Lawrence Ip, 2002, https://arxiv.org/abs/quant-ph/0211140
///  - Cirq implementation: https://github.com/quantumlib/Cirq/blob/main/examples/hidden_shift_algorithm.py

import Std.Diagnostics.*;
import Std.Random.*;
import Std.Arrays.*;
import Std.Convert.*;

@EntryPoint()
operation Main() : Unit {
    let qubitCount = 6;
    let sampleCount = 100;
    
    // define random secret shift
    let shift = GenerateRandomBitString(qubitCount);
    Message($"Secret shift sequence: {shift}");
    
    // run the Hidden Shift algorithm
    mutable results = [];
    for _ in 1..sampleCount {
        set results += [HiddenShiftAlgorithm(qubitCount, shift)];
    }
    
    // count occurrences of each result
    let frequencies = CountOccurrences(results);
    Message($"Sampled results: {frequencies}");
    
    // check if we found the correct shift
    let (mostCommonResult, _) = FindMostCommon(frequencies);
    
    Message($"Most common bitstring: {BoolArrayAsBinaryString(mostCommonResult)}");
    Message($"Expected shift: {BoolArrayAsBinaryString(shift)}");
    Message($"Found a match: {mostCommonResult == shift}");
}

operation GenerateRandomBitString(length : Int) : Bool[] {
    mutable bits = [];
    for _ in 1..length {
        set bits += [DrawRandomBool(0.5)];
    }
    bits
}

function BoolArrayAsBinaryString(arr : Bool[]) : String {
    mutable output = "";
    for entry in arr {
        output += entry ? "1" | "0";
    }
    output
}

function CountOccurrences(results : Bool[][]) : (Bool[], Int)[] {
    mutable uniqueResults = [];
    mutable counts = [];
    
    for result in results {
        let idx = IndexOf(r -> r == result, uniqueResults);
        if idx == -1 {
            set uniqueResults += [result];
            set counts += [1];
        } else {
            set counts w/= idx <- counts[idx] + 1;
        }
    }
    
    Zipped(uniqueResults, counts)
}

function FindMostCommon(frequencies : (Bool[], Int)[]) : (Bool[], Int) {
    mutable maxCount = 0;
    mutable mostCommon = [];
    
    for (result, count) in frequencies {
        if count > maxCount {
            set maxCount = count;
            set mostCommon = result;
        }
    }
    
    (mostCommon, maxCount)
}

/// oracle f(x) = Σ_i x_(2i) x_(2i+1)
operation OracleF(qubits : Qubit[]) : Unit is Adj + Ctl {
    for i in 0..Length(qubits)/2 - 1 {
        CZ(qubits[2*i], qubits[2*i + 1]);
    }
}

/// Hidden Shift algorithm implementation
operation HiddenShiftAlgorithm(qubitCount : Int, shift : Bool[]) : Bool[] {
    use qubits = Qubit[qubitCount];
    
    // 1 & 2. initialize qubits in superposition
    ApplyToEach(H, qubits);
    
    // 3. query oracle g (which is f with a shift)
    // Apply shift
    for i in 0..qubitCount - 1 {
        if shift[i] {
            X(qubits[i]);
        }
    }
    
    // apply oracle f
    OracleF(qubits);
    
    // undo shift
    for i in 0..qubitCount - 1 {
        if shift[i] {
            X(qubits[i]);
        }
    }
    
    // 4. apply Hadamard gates
    ApplyToEach(H, qubits);
    
    // 5. query oracle f
    OracleF(qubits);
    
    // 6. apply Hadamard gates (inverse Fourier)
    ApplyToEach(H, qubits);
    
    // 7. measure the qubits to get the shift
    let result = MeasureEachZ(qubits);
    let boolResults = ResultArrayAsBoolArray(result);

    ResetAll(qubits);
    
    boolResults
}