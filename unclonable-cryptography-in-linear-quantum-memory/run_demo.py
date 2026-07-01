#!/usr/bin/env python3
"""
run_demo.py  --  Narrated walkthrough of the Q# demo for

    "Unclonable Cryptography in Linear Quantum Memory"
    Omri Shmueli and Mark Zhandry, arXiv:2511.04633.

The Q# project (./src) is pure quantum/classical code that only RETURNS
results; this harness orchestrates the shots, gathers the returned values,
and processes/presents them in five steps. Run:

    pip install qdk
    python run_demo.py
"""

from collections import Counter
from pathlib import Path
from qdk import qsharp

PROJECT = str(Path(__file__).resolve().parent)


def rule(title):
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74)


def bits(bools):
    return "".join("1" if b else "0" for b in bools)


def main():
    qsharp.init(project_root=PROJECT)

    # ----------------------------------------------------------------- Step 1
    rule("STEP 1  -  The signing key is a COSET STATE  |sk> = SUM_{u in S} |u>")
    print(
        "S = span{10100, 01010, 00001} = {(a,b,a,b,c)} <= Z_2^5.\n"
        "The leading 2 bits (a,b) are the MESSAGE; bits 2,3 are locked copies\n"
        "(a forger must break that coupling); bit 4 is internal key entropy.\n"
        "(Paper: |sk> = SUM over S_y = ColSpan(A_y)+b_y, derived by hashing a\n"
        "uniform superposition through the oracle P and measuring pk := y; the\n"
        "demo fixes one tiny known S instead -- see README.)\n"
        "Preparing the key and inspecting its amplitudes (all equal => uniform\n"
        "superposition over the 8 elements of S):\n"
    )
    qsharp.eval("use qs = Qubit[5]; CosetKey.PrepareCosetKey(qs);")
    print(qsharp.dump_machine())
    qsharp.eval("ResetAll(qs);")  # tidy the persisted simulator state

    print("\nCircuit that prepares the coset key:")
    print(qsharp.circuit("{ use qs = Qubit[5]; CosetKey.PrepareCosetKey(qs); }"))

    # ----------------------------------------------------------------- Step 2
    rule("STEP 2  -  One key  ->  one signature on one (random) message")
    shots = 500
    sigs = qsharp.run("CosetKey.SignRandomMessage()", shots)
    dist = Counter(bits(s[:2]) for s in sigs)
    # Membership test mirrors CosetKey.InSubspace: u2 == u0 and u3 == u1.
    allvalid = all((s[2] == s[0] and s[3] == s[1]) for s in sigs)
    print(
        f"Measuring a fresh key collapses it to ONE element of S, i.e. a signature\n"
        f"on a random message. Over {shots} fresh keys the message is ~uniform:\n"
    )
    for m in ("00", "10", "01", "11"):
        print(f"    message {m}: {dist[m]:5d}")
    print(f"\nEvery measured string is a valid element of S: {allvalid}")
    print(
        "A single copy yields exactly ONE classical signature; no-cloning forbids\n"
        "copying |sk> to draw a second independent one."
    )

    # ----------------------------------------------------------------- Step 3
    rule("STEP 3  -  MEASURE-AND-CORRECT: signing a CHOSEN message in parallel")
    print(
        "Measure ALL message coordinates at once, KEEP the ones already matching\n"
        "the target, and re-entropize the mismatched ones in the SAME round\n"
        "(the correction's random phase is irrelevant -- the next step measures).\n"
        "Each round halves the mismatches. Statistics over 500 signings each;\n"
        "signatures are checked by the project's own CosetKey.VerifySignature:\n"
    )
    print(f"    {'message':>8} | {'valid':>7} | {'avg rounds':>10} | {'max rounds':>10}")
    print(f"    {'-'*8}-+-{'-'*7}-+-{'-'*10}-+-{'-'*10}")
    for m0, m1 in [(False, False), (True, False), (False, True), (True, True)]:
        n = 500
        lit0, lit1 = str(m0).lower(), str(m1).lower()
        res = qsharp.run(f"Signing.SignAndVerify({lit0}, {lit1})", n)
        rounds = [r for (r, _, _) in res]
        valid = sum(1 for (_, _, ok) in res if ok)
        label = ("1" if m0 else "0") + ("1" if m1 else "0")
        print(f"    {label:>8} | {valid:>4}/{n} | {sum(rounds)/n:>10.3f} | {max(rounds):>10}")
    print(
        "\nAlways valid, in a small expected number of rounds. Flagged difference:\n"
        "the paper's Construction 27 runs a FIXED 3 rounds by signing a codeword\n"
        "of an error-correcting code (verification accepts Hamming distance up to\n"
        "lambda/6); the demo omits the ECC and repeats until exact match, which\n"
        "costs O(log l) expected rounds instead."
    )

    # ----------------------------------------------------------------- Step 4
    rule("STEP 4  -  CLAW-FREE PERMUTATIONS: why you cannot sign two messages")
    print(
        "Two signatures on two messages = two preimages of the public key = a\n"
        "COLLISION in the scheme's hash. The oracle construction makes that hard\n"
        "via claw-free permutations H*(b,x)=Pi_b(x). Sampling claws (each from a\n"
        "fresh permutation pair, exhibited USING the trapdoor -- the attacker has\n"
        "no trapdoor, which is the whole point):\n"
    )
    samples = 50
    claws = qsharp.run("ClawFree.FindClaw(6)", samples)
    all_genuine = all(ok for (_, _, _, ok) in claws)
    y, x0, x1, _ = claws[0]
    print(f"    example (lambda=6): claw at y={y} is (b=0, x={x0}) and (b=1, x={x1})")
    print(f"    {samples}/{samples} sampled claws verified forward (Pi_0(x0)=Pi_1(x1)=y): {all_genuine}")
    print(
        "    Colliding inputs ALWAYS differ in the first bit -> H* is a 1-dimensional\n"
        "    coset partition function (its coordinate is exactly that first bit).\n"
    )
    print("Security scaling (Lemma 35: a q-query attacker wins w.p. O(q^3 / 2^lambda)):")
    print(f"    {'lambda':>7} | {'log2 q':>7} | {'1/Pr >= 2^?':>12}")
    print(f"    {'-'*7}-+-{'-'*7}-+-{'-'*12}")
    for lam in (64, 128, 256):
        exp = qsharp.run(f"ClawFree.ClawSecurityBitsExponent({lam}, 20)", 1)[0]
        print(f"    {lam:>7} | {20:>7} | {('2^'+str(exp)):>12}")

    # ----------------------------------------------------------------- Step 5
    rule("STEP 5  -  FOLDING: quadratic quantum memory -> LINEAR")
    print(
        "Prior keys are lambda^2 bits (lambda stacked 2-to-1 maps). Folding exposes a\n"
        "short O(lambda) part per input: each block's coordinate bit b_j plus the\n"
        "XOR-sum wbar of the rest. Invertible, and - the reduction's key trick -\n"
        "Q^{-1} can be simulated even while MISSING one block's trapdoor, recovering\n"
        "that block from the stored sum alone. Each trial below is self-checking\n"
        "Q# (Folding.FoldTrial) returning (roundtrip ok, trapdoor-free sim ok):\n"
    )
    print(f"    {'lambda':>7} | {'input |w|':>10} | {'folded part':>11} | {'shrink':>7} | {'roundtrip':>9} | {'sim-no-td':>9}")
    print(f"    {'-'*7}-+-{'-'*10}-+-{'-'*11}-+-{'-'*7}-+-{'-'*9}-+-{'-'*9}")
    for lam in (4, 8, 16):
        trials = 50
        results = qsharp.run(f"Folding.FoldTrial({lam})", trials)
        ok_rt = sum(1 for (rt, _) in results if rt)
        ok_sim = sum(1 for (_, sim) in results if sim)
        inp = qsharp.run(f"Folding.InputSizeBits({lam})", 1)[0]
        fold = qsharp.run(f"Folding.FoldedSizeBits({lam})", 1)[0]
        print(
            f"    {lam:>7} | {inp:>10} | {fold:>11} | {inp//fold:>6}x | "
            f"{ok_rt:>4}/{trials} | {ok_sim:>4}/{trials}"
        )
    print(
        "\nThe folded part (2*lambda-1 bits) is what becomes the vectorial part of\n"
        "the oracles -- and hence the signing key: linear in lambda. That is the\n"
        "paper's main result: unclonable cryptography in LINEAR quantum memory."
    )


if __name__ == "__main__":
    main()