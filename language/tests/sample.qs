import Std.Diagnostics.Fact;
import Std.Math.*;
import Std.Convert.*;

/// sample Q# tests, showing how to do testing with Q#

// high-level integration test that combines multiple testing approaches
@Test()
operation BellStateEntanglementTests() : Unit {
    // test 1: perfect correlation in computational basis
    TestBellStateCorrelation(PrepareBellState00, "Bell |00⟩ + |11⟩");
    
    // test 2: perfect anti-correlation 
    TestBellStateAntiCorrelation(PrepareBellState01, "Bell |01⟩ + |10⟩");
    
    // test 3: superposition properties in X-basis
    TestBellStateSuperposition(PrepareBellState00, "Bell state X-basis measurement");
}

/// tests deterministic properties of Bell states using exact quantum simulation
@Test()
operation BellStateDeterministicTests() : Unit {
    // test perfect correlation: measuring both qubits in Z-basis should always give same result
    for _ in 1..10 {
        use (q1, q2) = (Qubit(), Qubit());
        PrepareBellState00(q1, q2);
        
        let result1 = M(q1);
        let result2 = M(q2);
        
        Fact(result1 == result2, "Bell state |00⟩ + |11⟩ should show perfect Z-correlation");
        
        Reset(q1);
        Reset(q2);
    }
}

/// tests statistical properties using multiple measurements
@Test()
operation BellStateProbabilisticTests() : Unit {
    let tolerance = 0.1; // 10% tolerance for statistical tests
    
    // test that we get roughly 50/50 distribution of |00⟩ and |11⟩ outcomes
    let (zeros, ones) = MeasureBellStateDistribution(PrepareBellState00, 1000);
    let zeroRatio = IntAsDouble(zeros) / 1000.0;
    let oneRatio = IntAsDouble(ones) / 1000.0;
    
    Fact(AbsD(zeroRatio - 0.5) < tolerance, 
         $"Expected ~50% |00⟩ outcomes, got {zeroRatio * 100.0}%");
    Fact(AbsD(oneRatio - 0.5) < tolerance, 
         $"Expected ~50% |11⟩ outcomes, got {oneRatio * 100.0}%");
}

/// tests violation of Bell's inequality (quantum advantage)
@Test()
operation BellInequalityViolationTest() : Unit {
    // measure correlations for different measurement settings based on your reference
    let pAB = RunBellTest(Uab);
    let pBC = RunBellTest(Ubc);
    let pAC = RunBellTest(Uac);
    
    // Bell's inequality: P(AB) + P(BC) >= P(AC)
    let bellInequalityHolds = pAB + pBC >= pAC;
    
    // for entangled states, this inequality should be violated
    Fact(not bellInequalityHolds, 
         $"Bell inequality should be violated: P(AB)={pAB} + P(BC)={pBC} = {pAB + pBC} should be < P(AC)={pAC}");
}

operation PrepareBellState00(q1 : Qubit, q2 : Qubit) : Unit {
    H(q1);
    CNOT(q1, q2);
}

operation PrepareBellState01(q1 : Qubit, q2 : Qubit) : Unit {
    H(q1);
    X(q2);
    CNOT(q1, q2);
}

/// tests perfect correlation in computational basis
operation TestBellStateCorrelation(prepareState : (Qubit, Qubit) => Unit, description : String) : Unit {
    mutable correlationCount = 0;
    let trials = 100;
    
    for _ in 1..trials {
        use (q1, q2) = (Qubit(), Qubit());
        prepareState(q1, q2);
        
        let result1 = M(q1);
        let result2 = M(q2);
        
        if result1 == result2 {
            set correlationCount += 1;
        }
        
        Reset(q1);
        Reset(q2);
    }
    
    Fact(correlationCount == trials, 
         $"{description}: Expected perfect correlation, got {correlationCount}/{trials}");
}

/// tests perfect anti-correlation in computational basis
operation TestBellStateAntiCorrelation(prepareState : (Qubit, Qubit) => Unit, description : String) : Unit {
    mutable antiCorrelationCount = 0;
    let trials = 100;
    
    for _ in 1..trials {
        use (q1, q2) = (Qubit(), Qubit());
        prepareState(q1, q2);
        
        let result1 = M(q1);
        let result2 = M(q2);
        
        if result1 != result2 {
            set antiCorrelationCount += 1;
        }
        
        Reset(q1);
        Reset(q2);
    }
    
    Fact(antiCorrelationCount == trials, 
         $"{description}: Expected perfect anti-correlation, got {antiCorrelationCount}/{trials}");
}

/// tests superposition properties by measuring in X-basis
operation TestBellStateSuperposition(prepareState : (Qubit, Qubit) => Unit, description : String) : Unit {
    mutable randomResults = 0;
    let trials = 100;
    let tolerance = 20; // allow some statistical variation
    
    for _ in 1..trials {
        use (q1, q2) = (Qubit(), Qubit());
        prepareState(q1, q2);
        
        // measure in X-basis
        let result1 = MResetX(q1);
        let result2 = MResetX(q2);
        
        // in X-basis, Bell state should give random-looking results
        if result1 == One or result2 == One {
            set randomResults += 1;
        }
    }
    
    // should see significant variation (not all zeros)
    Fact(randomResults > tolerance, 
         $"{description}: Expected randomness in X-basis, got {randomResults}/{trials} non-zero");
}

/// measures the distribution of Bell state outcomes
operation MeasureBellStateDistribution(prepareState : (Qubit, Qubit) => Unit, trials : Int) : (Int, Int) {
    mutable zeroZeroCount = 0;
    mutable oneOneCount = 0;
    
    for _ in 1..trials {
        use (q1, q2) = (Qubit(), Qubit());
        prepareState(q1, q2);
        
        let result1 = M(q1);
        let result2 = M(q2);
        
        if result1 == Zero and result2 == Zero {
            set zeroZeroCount += 1;
        }
        elif result1 == One and result2 == One {
            set oneOneCount += 1;
        }
        
        Reset(q1);
        Reset(q2);
    }
    
    return (zeroZeroCount, oneOneCount);
}

/// measures correlation between two qubits using specified measurement operations
operation MeasureCorrelation(
    prepareState : (Qubit, Qubit) => Unit,
    measureBoth : (Qubit, Qubit) => (Result, Result),
    trials : Int
) : Double {
    mutable correlation = 0;
    
    for _ in 1..trials {
        use (q1, q2) = (Qubit(), Qubit());
        prepareState(q1, q2);
        
        let (result1, result2) = measureBoth(q1, q2);
        
        // convert results to ±1 and calculate correlation
        let value1 = result1 == Zero ? 1 | -1;
        let value2 = result2 == Zero ? 1 | -1;
        set correlation += value1 * value2;
        
        Reset(q1);
        Reset(q2);
    }
    
    return IntAsDouble(correlation) / IntAsDouble(trials);
}

// === measurement Operations for CHSH Bell Test ===

/// runs Bell test with specific measurement operation
operation RunBellTest(op : ((Qubit, Qubit) => Unit)) : Double {
    mutable res = [0, 0, 0, 0];

    let runs = 1000;
    for i in 0..runs {
        use (q1, q2) = (Qubit(), Qubit());
        PrepareBellState(q1, q2);
        op(q1, q2);
        let (r1, r2) = (MResetX(q1), MResetX(q2));

        if (r1 == Zero and r2 == Zero) {
            set res[0] = res[0] + 1;
        }
        if (r1 == Zero and r2 == One) {
            set res[1] = res[1] + 1;
        }
        if (r1 == One and r2 == Zero) {
            set res[2] = res[2] + 1;
        }
        if (r1 == One and r2 == One) {
            set res[3] = res[3] + 1;
        }
    }
    let p = IntAsDouble(res[0]) / IntAsDouble(runs);
    return p;
}

/// prepares Bell state for the inequality test
operation PrepareBellState(q1 : Qubit, q2 : Qubit) : Unit {
    X(q1);
    X(q2);
    H(q1);
    CNOT(q1, q2);
}

// bell measurement operations
operation Uab(q1 : Qubit, q2 : Qubit) : Unit {
    R1(PI() / 3.0, q2);
}

operation Uac(q1 : Qubit, q2 : Qubit) : Unit {
    R1(2.0 * PI() / 3.0, q2);
}

operation Ubc(q1 : Qubit, q2 : Qubit) : Unit {
    R1(PI() / 3.0, q1);
    R1(2.0 * PI() / 3.0, q2);
}