// =============================================================================
// ClawFree.qs
//
// STEP 4: WHY you cannot sign two messages -- the security foundation.
//
// A signature is an element u in S that opens to a message. Two signatures on
// two DIFFERENT messages are two distinct elements of S that collide under the
// scheme's hash H -- a "collision". So one-shot unforgeability reduces to the
// collision-resistance of H. The paper's oracle-model construction builds that
// hardness from CLAW-FREE PERMUTATIONS:
//
//      H*(b, x) = Pi_b(x),   for two random permutations Pi_0, Pi_1.
//
// (Definition 33 states this as H*: {0,1}^{lambda+1} -> {0,1}^lambda; the
// folding section equivalently writes L: {0,1}^lambda -> {0,1}^{lambda-1}.
// This file follows the latter convention, since Folding.qs builds on it.)
//
// A collision -- a "claw" -- is a pair (0,x0),(1,x1) with Pi_0(x0)=Pi_1(x1).
// Two structural facts the paper uses heavily:
//   * collisions ALWAYS differ in the first bit b (each output has exactly one
//     preimage under each Pi_b), so H* is a 1-dimensional coset partition
//     function (Definition 32) whose "coordinate" is exactly that first bit;
//   * finding a claw is hard -- Lemma 35: any (unbounded) quantum algorithm
//     making q queries finds a collision with probability only O(q^3/2^lambda).
//     Hence producing two signatures is infeasible.
//
// This file is fully classical Q# (no qubits): we *build* the primitive and
// exhibit a claw so the structure is concrete. NOTE (simplification, flagged):
// FindClaw uses the inverse permutations -- the trapdoor -- to exhibit a claw
// instantly. The whole point of Lemma 35 is that WITHOUT the inverses, a
// q-query attacker needs ~2^{lambda/3} queries. The demo plays the role of the
// trapdoor holder, not the attacker.
// =============================================================================

import Std.Random.*;

/// A uniformly random permutation table of [0, 2^n) via Fisher-Yates.
operation RandomPermutation(n : Int) : Int[] {
    let size = 1 <<< n;
    mutable p = [0, size = size];
    for i in 0..size - 1 {
        set p[i] = i;
    }
    for i in size - 1..-1..1 {
        let j = DrawRandomInt(0, i);
        let tmp = p[i];
        set p[i] = p[j];
        set p[j] = tmp;
    }
    p
}

/// Inverse of a permutation table.
function InvertPermutation(p : Int[]) : Int[] {
    mutable inv = [0, size = Length(p)];
    for i in 0..Length(p) - 1 {
        set inv[p[i]] = i;
    }
    inv
}

/// Evaluate the claw-free permutation H*(b, x) = Pi_b(x).
function EvalClawFree(pi0 : Int[], pi1 : Int[], b : Bool, x : Int) : Int {
    (b ? pi1 | pi0)[x]
}

/// Builds a random claw-free permutation pair on {0,1}^{lambda-1} (so that
/// H*: {0,1}^lambda -> {0,1}^{lambda-1}) and exhibits its unique claw at a
/// random output y, using the trapdoor (the inverse tables):
///   (y, x0, x1, ok) with H*(0,x0) = H*(1,x1) = y; `ok` re-checks the claw in
///   the FORWARD direction so the harness can verify without the tables.
/// By construction the two colliding inputs differ in their first bit (b=0 vs
/// b=1), illustrating the 1-dimensional coset/CPF structure folding builds on.
operation FindClaw(lambda : Int) : (Int, Int, Int, Bool) {
    let pi0 = RandomPermutation(lambda - 1);
    let pi1 = RandomPermutation(lambda - 1);
    let inv0 = InvertPermutation(pi0);
    let inv1 = InvertPermutation(pi1);
    let range = 1 <<< (lambda - 1);
    let y = DrawRandomInt(0, range - 1);
    let x0 = inv0[y];
    let x1 = inv1[y];
    let ok = EvalClawFree(pi0, pi1, false, x0) == y
         and EvalClawFree(pi0, pi1, true, x1) == y;
    (y, x0, x1, ok)
}

/// log2 of the (inverse) success bound of Lemma 35: a q-query quantum attacker
/// finds a claw with probability only O(q^3 / 2^lambda). We return
/// lambda - 3*log2(q), the exponent of 1/Pr -- bigger means harder.
function ClawSecurityBitsExponent(lambda : Int, log2q : Int) : Int {
    lambda - 3 * log2q
}
