namespace QuantumDecoration {
    operation CreateQuantumDecoration(numPositions : Int) : Result[] {
        let qubitsNeeded = numPositions * 3;
        use qubits = Qubit[qubitsNeeded];

        for i in 0..qubitsNeeded-1 {
            H(qubits[i]);
        }

        for i in 0..3..qubitsNeeded-1 {
            S(qubits[i]);
            T(qubits[i+1]);
            T(qubits[i+2]);
        }
        
        for i in 0..qubitsNeeded-1 {
            H(qubits[i]);
        }

        for i in 0..3..qubitsNeeded-4 {
            CNOT(qubits[i], qubits[i+3]);
        }

        for i in 0..3..qubitsNeeded-1 {
            CNOT(qubits[i], qubits[i+1]);
            CNOT(qubits[i], qubits[i+2]);
        }

        for i in 0..qubitsNeeded-1 {
            H(qubits[i]);
        }
        
        let results = MeasureEachZ(qubits);
        ResetAll(qubits);
        
        return results;
    }
}