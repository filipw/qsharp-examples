// =============================================================================
// Folding.qs
//
// STEP 5: the headline result's engine -- FOLDING COSET PARTITION FUNCTIONS,
// which shrinks the quantum signing key from quadratic to LINEAR size.
//
// Background (Sec. 1.1, "Folding Coset Partition Functions"; formalized via
// Theorem 34). Prior work needs a coset partition function (CPF, Definition
// 32) with 2^lambda-size preimage cosets that is collision-resistant; the only
// known way is stacking lambda independent 2-to-1 maps in parallel, giving a
// lambda^2-bit input. In the security reduction that input size lower-bounds
// the signing-key size. The fix: keep the long lambda^2-bit input where it is
// needed, but expose a short FOLDED representation of every input, and use the
// folded part as the vectorial part of the oracles -- which is what actually
// sits in the quantum key.
//
// The folding map, exactly as in the paper (claw-free permutations as the
// per-block 2-to-1 maps):
//   Input  w in Z_2^{lambda^2}, read as lambda blocks (b_j, x_j),
//          b_j in {0,1} (the block's public coordinate), x_j in Z_2^{lambda-1}.
//   Q(w) = ( (y_1,...,y_lambda),              <- lambda*(lambda-1) bits
//            (b_1,...,b_lambda, wbar) )       <- 2*lambda - 1 bits = O(lambda)
//   where  y_j = Pi_{j,b_j}(x_j)   and   wbar = XOR_j x_j.
//
// The folded part (b_1..b_lambda, wbar) is only 2*lambda-1 bits -- linear --
// yet the map (input <-> (long part, folded part)) stays invertible, and the
// paper further shows the folded preimage sets are themselves affine cosets,
// which is why the folded part can serve as the signing key's vector space.
//
// Q^{-1} normally uses every block's trapdoor (Pi^{-1}) and does NOT need
// wbar. The reduction's key trick: Q and Q^{-1} can be perfectly simulated
// even while MISSING one block's trapdoor i* (the block being attacked) --
// invert every other block, then recover the missing x_{i*} from the stored
// sum alone:  x_{i*} = wbar XOR (XOR_{j != i*} x_j).
// =============================================================================

import Std.Random.*;
import ClawFree.RandomPermutation;
import ClawFree.InvertPermutation;

/// A folding CPF instance: for each of the lambda blocks, the two permutations
/// (b=0 / b=1) and their inverses (the per-block "trapdoor").
struct FoldingKey {
    Lambda : Int,
    Pi0 : Int[][],
    Pi1 : Int[][],
    Inv0 : Int[][],
    Inv1 : Int[][]
}

/// Sample a fresh folding CPF for parameter lambda.
operation GenFoldingKey(lambda : Int) : FoldingKey {
    mutable pi0 = [[0], size = lambda];
    mutable pi1 = [[0], size = lambda];
    mutable iv0 = [[0], size = lambda];
    mutable iv1 = [[0], size = lambda];
    for j in 0..lambda - 1 {
        let p0 = RandomPermutation(lambda - 1);
        let p1 = RandomPermutation(lambda - 1);
        set pi0[j] = p0;
        set pi1[j] = p1;
        set iv0[j] = InvertPermutation(p0);
        set iv1[j] = InvertPermutation(p1);
    }
    new FoldingKey { Lambda = lambda, Pi0 = pi0, Pi1 = pi1, Inv0 = iv0, Inv1 = iv1 }
}

/// Forward folding map Q. Input is given block-wise as bs[j] (coordinate bit)
/// and xs[j] (value in Z_2^{lambda-1}, packed as an Int). Returns the long part
/// `ys` and the sum `wbar` completing the short folded part (bs, wbar).
function Fold(key : FoldingKey, bs : Bool[], xs : Int[]) : (Int[], Int) {
    mutable ys = [0, size = key.Lambda];
    mutable wbar = 0;
    for j in 0..key.Lambda - 1 {
        let pj = bs[j] ? key.Pi1[j] | key.Pi0[j];
        set ys[j] = pj[xs[j]];
        set wbar = wbar ^^^ xs[j];
    }
    (ys, wbar)
}

/// Standard inverse Q^{-1} using all trapdoors. Returns the recovered xs.
/// (As in the paper: the stored sum wbar is NOT needed when all trapdoors are
/// known -- each block inverts from (b_j, y_j) alone.)
function Unfold(key : FoldingKey, bs : Bool[], ys : Int[]) : Int[] {
    mutable xs = [0, size = key.Lambda];
    for j in 0..key.Lambda - 1 {
        let ip = bs[j] ? key.Inv1[j] | key.Inv0[j];
        set xs[j] = ip[ys[j]];
    }
    xs
}

/// The reduction's trapdoor-free simulation of Q^{-1} for the targeted block
/// `istar`: invert every OTHER block, then recover x_{istar} from the folded
/// sum alone. Returns the recovered value -- no Pi_{istar}^{-1} is consulted.
function UnfoldWithoutTrapdoor(key : FoldingKey, bs : Bool[], ys : Int[], wbar : Int, istar : Int) : Int {
    mutable acc = 0;
    for j in 0..key.Lambda - 1 {
        if j != istar {
            let ip = bs[j] ? key.Inv1[j] | key.Inv0[j];
            set acc = acc ^^^ ip[ys[j]];
        }
    }
    wbar ^^^ acc
}

/// Bits in the original input w (= lambda^2).
function InputSizeBits(lambda : Int) : Int { lambda * lambda }

/// Bits in the short folded coordinate part (= 2*lambda - 1 = O(lambda)).
/// This is the part that ends up as the linear-size quantum signing key.
function FoldedSizeBits(lambda : Int) : Int { 2 * lambda - 1 }

/// One self-checking trial for the harness: sample a fresh folding key and a
/// random input, then return
///   (roundtrip ok, trapdoor-free simulation ok)
/// where the first checks Unfold(Fold(w)) == w over all blocks, and the second
/// checks that a uniformly chosen block i* is recovered from the folded sum
/// without its trapdoor.
operation FoldTrial(lambda : Int) : (Bool, Bool) {
    let key = GenFoldingKey(lambda);
    let range = 1 <<< (lambda - 1);
    mutable bs = [false, size = lambda];
    mutable xs = [0, size = lambda];
    for j in 0..lambda - 1 {
        set bs[j] = DrawRandomInt(0, 1) == 1;
        set xs[j] = DrawRandomInt(0, range - 1);
    }
    let (ys, wbar) = Fold(key, bs, xs);
    let back = Unfold(key, bs, ys);
    mutable roundtrip = true;
    for j in 0..lambda - 1 {
        if back[j] != xs[j] {
            set roundtrip = false;
        }
    }
    let istar = DrawRandomInt(0, lambda - 1);
    let recovered = UnfoldWithoutTrapdoor(key, bs, ys, wbar, istar);
    (roundtrip, recovered == xs[istar])
}
