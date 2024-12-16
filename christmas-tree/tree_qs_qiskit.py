from abc import ABC, abstractmethod
from typing import Tuple, List
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator
import qsharp

class Colors:
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    RESET = '\033[0m'
    BOLD = '\033[1m'

class QuantumSimulator(ABC):
    @abstractmethod
    def run_simulation(self, num_positions: int) -> str:
        """Run quantum simulation and return binary result string"""
        pass

class QiskitSimulator(QuantumSimulator):
    def run_simulation(self, num_positions: int) -> str:
        qubits_needed = num_positions * 3
        qr = QuantumRegister(qubits_needed, 'q')
        cr = ClassicalRegister(qubits_needed, 'c')
        circuit = QuantumCircuit(qr, cr)
        
        for i in range(qubits_needed):
            circuit.h(qr[i])
        
        for i in range(0, qubits_needed, 3):
            circuit.s(qr[i])
            circuit.t(qr[i+1])
            circuit.t(qr[i+2])
        
        for i in range(qubits_needed):
            circuit.h(qr[i])
        
        for i in range(0, qubits_needed-3, 3):
            circuit.cx(qr[i], qr[i+3])
            
        for i in range(0, qubits_needed, 3):
            circuit.cx(qr[i], qr[i+1])
            circuit.cx(qr[i], qr[i+2])
        
        for i in range(qubits_needed):
            circuit.h(qr[i])
        
        circuit.measure(qr, cr)
        
        simulator = AerSimulator()
        result = simulator.run(circuit).result()
        counts = result.get_counts(circuit)
        return list(counts.keys())[0]

class QSharpSimulator(QuantumSimulator):
    def __init__(self, source_code: str):
        self.source_code = source_code
        qsharp.eval(self.source_code)

    def run_simulation(self, num_positions: int) -> str:
        result = qsharp.run(f"QuantumDecoration.CreateQuantumDecoration({num_positions})", shots=1)
        return ''.join(['1' if x == qsharp.Result.One else '0' for x in result[0]])

def parse_quantum_results(binary_result: str, num_positions: int) -> Tuple[List[int], List[int]]:
    """Parse binary quantum results into decorations and types."""
    decorations = []
    types = []
    
    for i in range(0, num_positions * 3, 3):
        presence = int(binary_result[i])
        if presence:
            type_bits = binary_result[i+1:i+3]
            decoration_type = int(type_bits, 2)
        else:
            decoration_type = 0
        decorations.append(presence)
        types.append(decoration_type)
    
    return decorations, types

def create_quantum_decorations(width: int, simulator: QuantumSimulator) -> Tuple[List[int], List[int]]:
    """Create quantum decorations using the provided simulator."""
    max_positions = 5
    results_decorations = []
    results_types = []
    
    for start_pos in range(0, width, max_positions):
        positions = min(max_positions, width - start_pos)
        binary = simulator.run_simulation(positions)
        decorations, types = parse_quantum_results(binary, positions)
        results_decorations.extend(decorations)
        results_types.extend(types)
    
    return results_decorations[:width], results_types[:width]

def get_decoration_char(type_num: int) -> str:
    """Return a festive decoration character based on the type"""
    chars = ['●', '★', '♦', '✶']
    return chars[type_num]

def get_color_code(type_num: int) -> str:
    """Return a color code based on the type"""
    colors = [Colors.RED, Colors.YELLOW, Colors.BLUE, Colors.MAGENTA]
    return colors[type_num]

def draw_christmas_tree(height: int, simulator: QuantumSimulator) -> None:
    """Draw a colorful tree with quantum decorations"""
    print(f"\n{Colors.BOLD}🎄 Quantum Christmas Tree! 🎄\n{Colors.RESET}")
    
    # draw star on top
    padding = height - 1
    print(" " * padding + f"{Colors.YELLOW}★{Colors.RESET}")
    
    # draw the tree from top to bottom
    for i in range(height):
        width = 2 * i + 1
        padding = height - i - 1
        
        # get quantum decorations and their types
        decorations, types = create_quantum_decorations(width, simulator)
        
        # create the row string
        row = " " * padding
        for j in range(width):
            if decorations[j]:
                color = get_color_code(types[j])
                decoration = get_decoration_char(types[j])
                row += f"{color}{decoration}{Colors.RESET}"
            else:
                row += f"{Colors.GREEN}*{Colors.RESET}"
        
        print(row)
    
    # Draw the trunk
    trunk_height = height // 3
    trunk_width = height // 2
    for i in range(trunk_height):
        trunk_padding = height - trunk_width//2 - 1
        print(" " * trunk_padding + f"{Colors.MAGENTA}#{Colors.RESET}" * trunk_width)

    # draw base decorations
    base_width = height * 2 - 1
    print(" " * 0 + f"{Colors.GREEN}~{Colors.RESET}" * base_width)
    print(f"\n{Colors.BOLD}🎁 Happy Quantum Holidays! 🎁{Colors.RESET}\n")

if __name__ == "__main__":
    tree_height = 12
    
    # try Q# first, fall back to Qiskit if the file is not there
    try:
        with open('tree.qs', 'r') as file:
            qs_content = file.read()
        simulator = QSharpSimulator(qs_content)
    except FileNotFoundError:
        print("Q# source not found, using Qiskit simulator")
        simulator = QiskitSimulator()
    
    draw_christmas_tree(tree_height, simulator)