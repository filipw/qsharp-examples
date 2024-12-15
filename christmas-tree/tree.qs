namespace QuantumDecoration {
    import Std.Diagnostics.DumpRegister;
    import Std.Diagnostics.DumpMachine;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;

    operation CreateQuantumDecoration(numPositions : Int) : Result[] {
        // Calculate total number of qubits needed
        let qubitsNeeded = numPositions * 3;
        // Allocate qubits
        use qubits = Qubit[qubitsNeeded] {
            // First Hadamard layer
            for i in 0..qubitsNeeded-1 {
                H(qubits[i]);
            }
            
            // Add phase gates for interference
            for i in 0..3..qubitsNeeded-1 {
                // Phase rotation on presence qubit
                S(qubits[i]);
                // Smaller phase rotations on type qubits
                T(qubits[i+1]);
                T(qubits[i+2]);
            }
            
            // Second Hadamard layer
            for i in 0..qubitsNeeded-1 {
                H(qubits[i]);
            }

            // Local entanglement
            for i in 0..3..qubitsNeeded-4 {
                CNOT(qubits[i], qubits[i+3]); // Connect presence qubits
            }

            // Connect presence to type within each position
            for i in 0..3..qubitsNeeded-1 {
                CNOT(qubits[i], qubits[i+1]);
                CNOT(qubits[i], qubits[i+2]);
            }

            // Final Hadamard layer
            for i in 0..qubitsNeeded-1 {
                H(qubits[i]);
            }
            
            // Measure all qubits
            let results = MeasureEachZ(qubits);
            ResetAll(qubits);
            
            return results;
        }
    }
}