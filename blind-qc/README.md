# Blind quantum computation on a modular processor — in Q#

A Q# demo of *Blind Quantum Computation on a Modular Superconducting Processor* (Yongxin Song,
Johannes Knörzer, Kieran Dalton, Andreas Wallraff, Jean-Claude Besse —
[arXiv:2605.14656](https://arxiv.org/abs/2605.14656)). The protocol is Morimae and Fujii's
([arXiv:1201.3966](https://arxiv.org/abs/1201.3966)). It runs interactively, with the figures, on
Copenhagen's Papers page: https://copenhagen.strathweb.com/papers/blind-qc

## The claim

A client computes by measuring qubits a server sends it, at angles it never announces. All
communication runs from server to client, so the server learns nothing, under an explicit
assumption of no side channel. The paper runs this on two flip-chip-bonded superconducting modules,
including a three-qubit Deutsch–Jozsa that lands on the ideal string in 0.797 and 0.798 of its
shots, and measures the server's qubit at purity 0.501 to 0.503 mid-sequence.

## The Q#

One file, with no entry point.

| | |
|---|---|
| `Blind.Rates`, `Noise1`, `Noise2`, `NoisyH`, `NoisyCZ`, `NoisyCNOT`, `Idle`, `Read` | the noise: X, Y or Z with probability ε/2 each, the 15 Pauli pairs with ε/12 each, in the circuit |
| `Blind.Cell`, `Blind.ZCell` | one measurement cell: a client pulse Rz(x), H and a readout s apply X^s H Rz(x) to the server's logical qubit |
| `Blind.HandOver` | moves the server's qubit to the client's, which is in \|0⟩, by two CNOTs |
| `Blind.ThreeCells` | three cells with adaptive angles and the corrections, measured for real |
| `Blind.DeferredCells` | the same with measurement deferred, dumping after each stage: the server's reduced state is then the average over the outcomes |
| `Blind.ClusterPattern` | the 7-node Deutsch–Jozsa cluster of the paper's Fig. 4, with both programs |
| `Blind.DjReplay` | the same pattern on six reused qubits, as the hardware runs it, with decoding |

## Differences from the paper

- The paper has no code. These kernels are written from its figures and text.
- A statevector has no crosstalk, so a flat server state here is a theorem being displayed. It tests
  nothing about hardware.
- The noise is a model at the paper's rates, not a reproduction of the device: no ZZ crosstalk, no
  spread between qubits, no thermal relaxation, no reset error, and none of the paper's CPMG pulses
  (reset here is a perfect `MResetZ`). Where the idling goes is this code's choice, because the
  paper doesn't say.
- The sign of a cell is the one in the paper's main text and Fig. 1c. Its Fig. 6a and Methods B write
  the opposite one.
