// =============================================================================
// Signing.qs
//
// STEP 3: signing a CHOSEN message with a single small coset state, by
// "measure and correct", in a constant number of parallel rounds.
//
// "Measure and correct" is the paper's own terminology (Sec. 3.2): their
// algorithm (Construction 27) "is a direct generalization of the 'measure and
// correct' signing algorithm of [Shm22]". The basic key only signs a *random*
// message (CosetKey.SignRandomMessage). To sign a message we *choose*, the
// signer measures ALL message coordinates at once, keeps the positions that
// already match the target, and "corrects" (re-entropizes) the positions that
// don't -- then re-measures. Two ideas make this efficient:
//
//   * Each coordinate independently collapses to the right value with prob
//     1/2, so each round halves the mismatches (1/2 -> 3/4 -> 7/8 of bits
//     match, Sec. 1.1) -- and corrections for different coordinates happen IN
//     THE SAME round. This is the paper's "parallel signing": signing bit j+1
//     does not wait for bit j.
//   * The correction returns mismatched coordinates to uniform superposition
//     "with a random choice of phase; the phase turns out not to matter, since
//     the next step is to just measure the bits anyway" (Sec. 1.1). The
//     within/apply block below realizes exactly that: it re-coheres the
//     coordinate up to an irrelevant sign.
//
// SIMPLIFICATIONS RELATIVE TO THE PAPER (flagged):
//   * Construction 27 runs a FIXED 3 iterations and signs a CODEWORD of an
//     error-correcting code (min distance > lambda/3); Ver then accepts any u
//     within Hamming distance lambda/6 of the codeword, so 3 rounds suffice
//     with overwhelming probability. The demo has no ECC: it repeats until the
//     message matches EXACTLY, giving O(log l) expected rounds (~2.7 observed
//     for l = 2). Unforgeability survives the ECC relaxation because at most
//     one codeword lies within distance lambda/6 of any u.
//   * In the paper the signer cannot see or rebuild the state; the correction
//     is one query to the dual oracle D, sandwiched between H^(x)k layers, and
//     D is augmented with l membership checks for subspaces of S_{y,0^l}-perp
//     (each one dimension smaller) so that ARBITRARY mismatched positions can
//     be fixed. In our tiny S = {(a,b,a,b,c)} the dual structure is explicit:
//     coordinate a lives on the locked pair (qubit 0, qubit 2) and coordinate
//     b on (qubit 1, qubit 3), so re-entropizing is the local circuit
//          within { CNOT(q0, q2); } apply { H(q0); }
//     which sends |a*>|a*> -> (|00> +/- |11>)/sqrt2. The +/- sign is the
//     "random choice of phase" the paper argues is harmless.
// =============================================================================

import CosetKey.PrepareCosetKey;

/// Re-entropizes one "locked coordinate pair" (seed qubit + its copy) back to
/// a fresh superposition, up to an irrelevant global sign. This is the
/// "correct" step of measure-and-correct, played by the explicit structure of
/// S where the paper's augmented dual oracle D would act.
operation RecohereCoordinate(seed : Qubit, copy : Qubit) : Unit {
    within {
        CNOT(seed, copy);   // unlink: copy <- a* XOR a* = 0, seed still holds a*
    } apply {
        H(seed);            // seed -> |+> (a*=0) or |-> (a*=1); sign irrelevant
    }                       // re-link on exit -> (|00> +/- |11>)/sqrt2
}

/// Signs a chosen 2-bit message (m0,m1) from ONE coset key by measure-and-
/// correct. Returns (number of rounds used, the signature bits). The signature
/// always verifies (CosetKey.VerifySignature); the harness gathers round
/// statistics over many shots.
operation SignChosenMessage(m0 : Bool, m1 : Bool) : (Int, Bool[]) {
    use qs = Qubit[5];
    PrepareCosetKey(qs);

    mutable rounds = 0;
    mutable a = false;
    mutable b = false;

    repeat {
        set rounds += 1;
        // Measure the two "message" coordinates of u in parallel.
        set a = M(qs[0]) == One;
        set b = M(qs[1]) == One;
        // Keep matches; correct mismatches -- both in the SAME round (parallel).
        if a != m0 { RecohereCoordinate(qs[0], qs[2]); }
        if b != m1 { RecohereCoordinate(qs[1], qs[3]); }
    } until (a == m0) and (b == m1);

    // Message coordinates now equal (m0,m1); read out the rest of u.
    let sig = [a, b, M(qs[2]) == One, M(qs[3]) == One, M(qs[4]) == One];
    ResetAll(qs);
    (rounds, sig)
}

/// Harness convenience: sign and immediately verify with the scheme's own
/// verifier. Returns (rounds, signature, verified) -- still pure returned
/// results, so the Python harness only aggregates.
operation SignAndVerify(m0 : Bool, m1 : Bool) : (Int, Bool[], Bool) {
    let (rounds, sig) = SignChosenMessage(m0, m1);
    (rounds, sig, CosetKey.VerifySignature(m0, m1, sig))
}
