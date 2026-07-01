// =============================================================================
// CosetKey.qs
//
// STEPS 1 & 2 of the demo: the unclonable signing key is a COSET STATE.
//
// Paper: "Unclonable Cryptography in Linear Quantum Memory" (Shmueli & Zhandry,
// arXiv:2511.04633). In the one-shot-signature (OSS) scheme (Construction 26)
// the quantum secret key is
//                  |sk> = (1/sqrt|S_y|) * SUM_{u in S_y} |u>
// where S_y = ColSpan(A_y) + b_y is a secret affine subspace tied to the
// classical public key y. A signature on a message is a vector u in S_y whose
// message coordinates open to that message. The no-cloning theorem makes the
// key unclonable; collision-resistance of the underlying hash makes it
// "one-shot" (see ClawFree.qs / Folding.qs).
//
// To make everything concrete and *runnable*, we instantiate one tiny coset.
// We work in Z_2^5 with the rank-3 subspace
//
//      S = span{ 10100, 01010, 00001 } = { (a, b, a, b, c) : a,b,c in {0,1} }.
//
// Reading a vector u = (u0,u1,u2,u3,u4) in S:
//   * the MESSAGE is the leading 2 bits (u0,u1) = (a,b);
//   * bits u2,u3 are "locked" copies of the message (u2=a, u3=b) -- this
//     coupling is what a forger would have to break;
//   * bit u4 = c is free "internal" entropy of the key.
//
// SIMPLIFICATIONS RELATIVE TO THE PAPER (flagged; see README for the full list):
//   * In Construction 26 the key pair is *derived*: SigGen evaluates the oracle
//     P on a uniform superposition, un-computes via P^{-1}, and MEASURES the
//     hash register to obtain pk := y; the state collapses to |sk> over S_y,
//     which stays hidden behind the oracles (P, P^{-1}, D). Here we fix one
//     small S and prepare |sk> directly from its known generators.
//   * Verification in the paper checks P^{-1}(y, u) != bottom (membership in
//     S_y via the inverse oracle) plus a Hamming-distance condition against an
//     error-correcting code (see Signing.qs). Here the verifier tests
//     membership in S explicitly (u2==u0 AND u3==u1) and exact message match.
//   * The paper takes the message coordinates to be the LAST lambda bits of a
//     k = 2*lambda-bit vector (Constructions 26/27); the demo puts the 2
//     message bits FIRST -- a pure relabeling, as the intro (Sec. 1.1) does too.
// =============================================================================

import Std.Diagnostics.*;

/// Prepares the coset signing key |sk> = SUM_{u in S} |u> on a 5-qubit register.
///
/// Construction: put the three "coordinate" qubits (0,1,4) into uniform
/// superposition, then realize the basis of S in place with two CNOTs. The
/// result is the equal superposition over all 8 elements of S -- exactly the
/// shape of the coset state that sits in long-term quantum memory in the paper.
operation PrepareCosetKey(qs : Qubit[]) : Unit is Adj + Ctl {
    Fact(Length(qs) == 5, "This demo coset lives in Z_2^5.");
    // Coordinates a, b, c become qubits 0, 1, 4.
    H(qs[0]);
    H(qs[1]);
    H(qs[4]);
    // Lock the redundant copies: u2 = a, u3 = b. (Column span of the basis.)
    CNOT(qs[0], qs[2]);
    CNOT(qs[1], qs[3]);
}

/// Classical membership test: is u in S?  (u2==u0 and u3==u1; u4 free.)
/// Stands in for the paper's oracle membership check P^{-1}(y, u) != bottom.
function InSubspace(u : Bool[]) : Bool {
    (u[2] == u[0]) and (u[3] == u[1])
}

/// Verifies a candidate signature `sig` for a claimed message (m0,m1): the
/// signature must be a genuine element of S AND open to that message exactly.
/// (Paper's Ver instead allows Hamming distance <= lambda/6 to the codeword of
/// the message -- the ECC relaxation the demo omits; see Signing.qs.)
function VerifySignature(m0 : Bool, m1 : Bool, sig : Bool[]) : Bool {
    InSubspace(sig) and (sig[0] == m0) and (sig[1] == m1)
}

/// STEP 2 -- "one key, one signature." Measuring the coset key collapses it to
/// a single random element of S, i.e. a valid signature on a *random* message
/// (its leading two bits). A single copy of the key yields exactly one such
/// classical string; the superposition is gone afterward (no-cloning), so you
/// cannot draw a second, independent signature from the same key.
operation SignRandomMessage() : Bool[] {
    use qs = Qubit[5];
    PrepareCosetKey(qs);
    let sig = [
        M(qs[0]) == One, M(qs[1]) == One, M(qs[2]) == One,
        M(qs[3]) == One, M(qs[4]) == One
    ];
    ResetAll(qs);
    sig
}
