namespace QuantumDecoration {
    operation CreateQuantumDecoration(numPositions : Int) : Result[] {
        let qubitsNeeded = numPositions * 3;
        use qubits = Qubit[qubitsNeeded];

        // first Hadamard layer
        for i in 0..qubitsNeeded-1 {
            H(qubits[i]);
        }
        
        // add phase gates for interference
        for i in 0..3..qubitsNeeded-1 {
            // phase rotation on presence qubit
            S(qubits[i]);
            // smaller phase rotations on type qubits
            T(qubits[i+1]);
            T(qubits[i+2]);
        }
        
        // second Hadamard layer
        for i in 0..qubitsNeeded-1 {
            H(qubits[i]);
        }

        // entanglement - connecting presence qubits
        for i in 0..3..qubitsNeeded-4 {
            CNOT(qubits[i], qubits[i+3]);
        }

        // connect presence to type within each position
        for i in 0..3..qubitsNeeded-1 {
            CNOT(qubits[i], qubits[i+1]);
            CNOT(qubits[i], qubits[i+2]);
        }

        // final Hadamard layer
        for i in 0..qubitsNeeded-1 {
            H(qubits[i]);
        }
        
        // measure all qubits
        let results = MeasureEachZ(qubits);
        ResetAll(qubits);
        
        return results;
    }
}