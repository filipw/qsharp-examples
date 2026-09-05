# Q# Examples

A collection of various interesting Q# examples.

## Algorithms

| Sample Name | Description | Blog Post | Type |
|-------------|-------------|-----------|------|
| 💻 [Christmas Tree](./christmas-tree) | 📝 Decorating a Christmas tree with a quantum computer. Uses both Q# and Qiskit (Aer) simulators. | [Link](https://www.strathweb.com/2024/12/decorating-a-quantum-christmas-tree-with-qsharp-and-qiskit/) | Executable |
| 💻 [Quokka](./quokka) | 📝 A Jupyter notebook showing how to run Q# code on a [Quokka](https://www.quokkacomputing.com/) quantum simulator. | [Link](https://strathweb.com/2025/03/running-qsharp-code-on-quokka/) | Jupyter Notebook |
| 💻 [Certified deletion](./certified-deletion) | 📝 Q# sample illustrating [Quantum encryption with certified deletion](https://link.springer.com/chapter/10.1007/978-3-030-64381-2_4) from Anne Broadbent and Rabib Islam. | [Link](https://www.strathweb.com/2023/12/exploring-quantum-encryption-and-certified-deletion-with-qsharp/) | Executable |
| 💻 [Elitzur-Vaidman tester](./elitzur-vaidman) | 📝 [Elitzur-Vaidman](https://arxiv.org/abs/hep-th/9305002) bomb tester thought experiment, as a playful Q# Santa-themed sample | [Link](https://www.strathweb.com/2022/12/q-holiday-calendar-2022-peeking-into-santas-gifts-with-q/) | Executable
| 💻 [Hidden shift](./hidden-shift) | 📝 [Hidden shift](https://arxiv.org/abs/quant-ph/0211140) problem | N/A | Executable
| 💻 [Simon's algorithm](./simons-algorithm/) | 📝 A toy demonstration of [Simon's algorithm](https://epubs.siam.org/doi/10.1137/S0097539796298637) | N/A | Executable
| 💻 [Certified randomness amplification](./certified-randomness) | 📝 Q# sample illustrating [Certified randomness amplification](https://arxiv.org/abs/2511.03686), Liu et al, arXiv:2511.03686, 2025. | [Link](https://strathweb.com/2025/12/certified-randomness-amplification-with-qsharp/) | Jupyter Notebook
| 💻 [Composite qDRIFT](./composite-qdrift) | 📝 Q# implementation of the Hagan-Wiebe composite Trotter/qDRIFT channel, used to measure the cost objective of [Optimal Lower Bounds for Hamiltonian Simulation](https://arxiv.org/abs/2607.19852), Zlokapa, Allen and Harrow, arXiv:2607.19852, 2026, on power-law Hamiltonians with an exact channel evaluation. | [Link](https://www.strathweb.com/2026/09/where-to-cut-a-hamiltonian-composite-qdrift-with-qsharp/) | Executable

## Errors

| Sample Name | Description | Blog Post | Type |
|-------------|-------------|-----------|------|
| 💻 [Cat state with noise](./cat-state-with-noise/) | 📝 A simple 8-qubit GHZ state, with a basic noise model and, optionally, qubit loss, illustrating their impact on measurement. Run with 1000 shots. | N/A | Executable
| 💻 [Bit flip error correction](./error-correction/bitflip) | 📝 Bit flip error correction samples (auxiliary qubit-based syndrome extraction, automatic correction with auxiliary qubits and direct parity measurement) | N/A | Executable
| 💻 [Phase flip error correction](./error-correction/phaseflip) | 📝 Phase flip error correction samples (auxiliary qubit-based syndrome extraction, automatic correction with auxiliary qubits and direct parity measurement) | N/A | Executable

## Language features

| Sample Name | Description | Blog Post | Type |
|-------------|-------------|-----------|------|
| 💻 [Array shuffle](./language/shuffle) | 📝 Shuffling an array using Q# (classical logic). | [Link](https://www.strathweb.com/2023/12/shuffling-an-array-in-qsharp/) | Executable |
| 💻 [New array update syntax](./language/arrays/) | 📝 The improved array update syntax introduced in Q# 1.17 | [Link](https://www.strathweb.com/2025/06/a-cat-jumped-on-a-keyboard-and-fixed-qsharp-array-syntax/) | Executable |
| 💻 [Q# Testing (Bell state)](./language/tests/sample.qs) | 📝 Bell state tests demonstrating different Q# testing approaches | N/A | Tests |

## Further samples

You can find more Q# samples in my [Introduction to Quantum Computing with Q# and QDK](https://github.com/filipw/intro-to-qc-with-qsharp-book?tab=readme-ov-file#list-of-examples) book repository.