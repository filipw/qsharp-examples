/// Deterministic Suzuki-Trotter product formulas of arbitrary even order.
import Std.Convert.*;
import Std.Math.*;
import Std.Diagnostics.Fact;
import Hamiltonian.ApplyTerm;

/// Suzuki's recursion parameter u_k for order = 2k:  u = 1 / (4 - 4^(1/(2k-1))).
function SuzukiU(order : Int) : Double {
    1.0 / (4.0 - (4.0 ^ (1.0 / IntAsDouble(order - 1))))
}

/// A single product-formula step approximating exp(-i H time).
///
/// order = 1 is the plain forward sweep; order = 2 is the symmetric (Strang)
/// sweep; higher even orders come from Suzuki's five-fold recursion.
operation SuzukiStep(
    order : Int,
    terms : Pauli[][],
    coeffs : Double[],
    time : Double,
    qs : Qubit[]
) : Unit is Adj + Ctl {
    Fact(order == 1 or (order > 0 and order % 2 == 0), "product formula order must be 1 or a positive even integer");
    Fact(Length(terms) == Length(coeffs), "terms and coeffs must have the same length");
    if order == 1 {
        for j in 0..Length(terms) - 1 {
            ApplyTerm(terms[j], coeffs[j], time, qs);
        }
    } elif order == 2 {
        let half = time / 2.0;
        for j in 0..Length(terms) - 1 {
            ApplyTerm(terms[j], coeffs[j], half, qs);
        }
        for j in Length(terms) - 1..-1..0 {
            ApplyTerm(terms[j], coeffs[j], half, qs);
        }
    } else {
        let u = SuzukiU(order);
        SuzukiStep(order - 2, terms, coeffs, u * time, qs);
        SuzukiStep(order - 2, terms, coeffs, u * time, qs);
        SuzukiStep(order - 2, terms, coeffs, (1.0 - 4.0 * u) * time, qs);
        SuzukiStep(order - 2, terms, coeffs, u * time, qs);
        SuzukiStep(order - 2, terms, coeffs, u * time, qs);
    }
}

/// exp(-i H time) approximated by `steps` repetitions of an order-`order` step.
operation TrotterEvolve(
    order : Int,
    terms : Pauli[][],
    coeffs : Double[],
    time : Double,
    steps : Int,
    qs : Qubit[]
) : Unit is Adj + Ctl {
    let dt = time / IntAsDouble(steps);
    for _ in 1..steps {
        SuzukiStep(order, terms, coeffs, dt, qs);
    }
}
