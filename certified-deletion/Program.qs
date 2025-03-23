import Std.Arrays.*;
import Std.Diagnostics.*;
import Std.Random.*;
import Std.Convert.*;
import Std.Logical.*;
import Std.Math.*;

operation Main() : Unit {
    // in QDK 1.x the entrypoint cannot have parameters
    // we have to fake it by having hardcoded data that changes from run to run
    let message = 4;
    let decryptFlow = false; // switch to true if necessary

    // we can reasonably only simulate this in one shot for messages up to 8 bits
    let message_bitsize = BitSizeI(message);
    Fact(message_bitsize < 8, "We can't simulate for messages of more than 8 bits");
    let algorithm_bitsize = 16;

    Message($"Running the algorithm for message: {message}");
    Message("");

    let theta = CreateRandomBoolArrayWithEqualDistribution(algorithm_bitsize);
    let r = DrawMany(() => DrawRandomInt(0, 1) == 1, algorithm_bitsize, ());

    use qubits = Qubit[algorithm_bitsize];

    mutable r_z = [];
    mutable r_x = [];

    for i in 0..algorithm_bitsize-1 {
        // X basis
        if theta[i] {
            if r[i] {
                // 1 is |-⟩
                X(qubits[i]);
                H(qubits[i]);
            } else {
                //0 is |+⟩
                H(qubits[i]);
            }

            // save the r value to r_x
            r_x += [r[i]];
        } else {
            // Z basis
            if r[i] {
                // 1 is |1⟩
                X(qubits[i]);
            } else {
                // 0 is |0⟩
                I(qubits[i]);
            }

            // save the r value to r_z
            r_z += [r[i]];
        }
    }

    Message($"θ: {BoolArrayAsBinaryString(theta)}");
    Message($"R: {BoolArrayAsBinaryString(r)}");
    Message($"R_z: {BoolArrayAsBinaryString(r_z)}");
    Message($"R_x: {BoolArrayAsBinaryString(r_x)}");
    Message("");

    let binaryMessage = IntAsBoolArray(message, algorithm_bitsize / 2);
    Message($"Binary message: {BoolArrayAsBinaryString(binaryMessage)}");

    // encrypt by doing XOR between the message and r_z
    let encrypted = MappedByIndex((i, x) -> Xor(binaryMessage[i], x), r_z);
    Message($"Encrypted message: {BoolArrayAsBinaryString(encrypted)}");
    Message("");

    if decryptFlow {
        // decrypt flow
        let decrypted = Decrypt(qubits, theta, encrypted);
        Fact(decrypted == binaryMessage, "The decrypted message is different than the original one!");
    } else {
        // delete and verify flow
        let deletion_proof = Delete(qubits);
        VerifyDeletion(theta, r_x, deletion_proof);
    }

    ResetAll(qubits);
}

operation Decrypt(qubits : Qubit[], theta : Bool[], encrypted : Bool[]) : Bool[] {
    // decrypt using theta as key
    // first obtain r_z by measuring only the qubits that were encoded in the Z basis
    mutable r_z_from_measurement = [];
    for i in 0..Length(theta)-1 {
        if not theta[i] {
            r_z_from_measurement += [M(qubits[i]) == One];
        }
    }
    Message($"R_z from qubits: {BoolArrayAsBinaryString(r_z_from_measurement)}");

    // now perform XOR between the encrypted data and the r_z
    let decrypted = MappedByIndex((i, x) -> Xor(encrypted[i], x), r_z_from_measurement);

    // the decrypted data should be identical to raw message
    Message($"Decrypted binary message: {BoolArrayAsBinaryString(decrypted)}");
    Message($"Decrypted message: {BoolArrayAsInt(decrypted)}");
    decrypted
}

operation Delete(qubits : Qubit[]) : Bool[] {
    mutable deletion_proof = [];
    for i in 0..Length(qubits)-1 {
        deletion_proof += [MResetX(qubits[i]) == One];
    }

    deletion_proof
}

function VerifyDeletion(theta : Bool[], r_x : Bool[], d : Bool[]) : Unit {
    mutable d_x = [];
    for i in 0..Length(theta)-1 {
        if theta[i] {
            d_x += [d[i]];
        }
    }

    // now verify the deletion by comparing d_x to r_x - they must be identical
    Message($"R_x from qubits: {BoolArrayAsBinaryString(d_x)}");
    Fact(d_x == r_x, "R_x obtained from measuring the qubits is different than the original one!");
}

function BoolArrayAsBinaryString(arr : Bool[]) : String {
    mutable output = "";
    for entry in arr {
        output += entry ? "1" | "0";
    }
    output
}

operation Shuffled<'T>(array : 'T[]) : 'T[] {
    let n = Length(array);
    mutable shuffled = array;

    for i in 0..n - 2 {
        let j = DrawRandomInt(i, n - 1);
        shuffled = Swapped(i, j, shuffled);
    }

    shuffled
}

operation CreateRandomBoolArrayWithEqualDistribution(size : Int) : Bool[] {
    Fact(size % 2 == 0, "Size must be divisble by 2");

    let array = [true, size = size / 2];
    return Shuffled(Padded(-size, false, array));
}