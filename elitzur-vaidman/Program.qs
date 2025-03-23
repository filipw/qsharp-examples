import Std.Math.*;
import Std.Convert.*;
import Std.Arrays.*;

operation Main() : Unit {
    // use various values of n to see how the success rate changes
    let n = 16; // try also n = 16
    let iterations = 1000;

    Message("Simple case:");
    SimpleVariant(iterations);

    Message("");
    Message("******************************");
    Message("");
    Message("Advanced case:");
    AdvancedVariant(n, iterations);
}

operation SimpleVariant(iterations : Int) : Unit {
    for working in [false, true] {
        mutable results = [0, size = 4];

        // for not working alarms we should see inconclusive resullt in 100% of cases
        // for working alarms we should see alarm in 50% of cases
        // and in the remianing 50%
        //  - 25% we will see incoclusive (same result as for not working alarm)
        //  - 25% we will detect a working alarm (different result as for not working alarm)
        Message("");
        Message($"Alarm working? {working}");
        for i in 1..iterations {
            use (q_tester, q_alarm) = (Qubit(), Qubit());

            H(q_tester);

            // if alarm is working, there is entanglemet between testing photon and the alarm triggering
            if (working) {
                CNOT(q_tester, q_alarm);
            } else {
                I(q_tester);
            }

            H(q_tester);

            let detector = MResetZ(q_tester) == One;
            let alarmTriggered = MResetZ(q_alarm) == One;

            // |00⟩ - inconclusive, Santa did not catch us
            if (not detector and not alarmTriggered) {
                set results w/= 0 <- results[0] + 1;
            }

            // |01⟩ - inconclusive, Santa alerted
            if (not detector and alarmTriggered) {
                set results w/= 1 <- results[1] + 1;
            }

            // |10⟩ - working alarm detected, Santa did not catch us
            if (detector and not alarmTriggered) {
                set results w/= 2 <- results[2] + 1;
            }

            // |11⟩ - working alarm detected, Santa alerted
            if (detector and alarmTriggered) {
                set results w/= 3 <- results[3] + 1;
            }
        }
        Message($"Alarm triggered: {results[1] + results[3]}");
        Message($"Detector 1: {results[0]}");
        Message($"Detector 2: {results[2]} - working alarm detected safely");
        Message($"{results}");
    }
}

operation AdvancedVariant(n : Int, iterations : Int) : Unit {
    for working in [false, true] {
        mutable successes = 0;
        Message("");
        Message($"Alarm working? {working}");

        for i in 1..iterations {
            use (q_tester, q_alarm) = (Qubit(), Qubit());

            mutable alarmTriggered = [Zero, size = n];
            for j in 0..n-1 {

                // rotate by π/n
                Ry(PI() / IntAsDouble(n), q_tester);

                // if alarm is working, there is entanglemet between testing photon and the alarm triggering
                if (working) {
                    CNOT(q_tester, q_alarm);
                } else {
                    I(q_tester);
                }

                // we can now measure |00⟩ or |11⟩ only
                // the probability amplitude of triggering the alarm (second qubit measures to |1⟩ is now sin(π/n)
                set alarmTriggered w/= j <- MResetZ(q_alarm);
                Reset(q_alarm);
            }

            let workingAlarmIdentified = MResetZ(q_tester) == Zero;

            let success = All(r -> r == Zero, alarmTriggered) and workingAlarmIdentified;
            if (success) {
                successes += 1;
            }
        }

        Message($"Success rate: {IntAsDouble(successes) * 100.0 / IntAsDouble(iterations)}%");
    }
}