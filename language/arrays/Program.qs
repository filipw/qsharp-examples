operation Main() : Unit {
    // This sample shows the old array update syntax, compared to the new one, introduced in Q# 1.17
    // For more details, see blog post at https://www.strathweb.com/2025/06/a-cat-jumped-on-a-keyboard-and-fixed-qsharp-array-syntax/
    
    // [0, 0, 0]
    mutable array1_old = [0, size = 3];
    // [0, 5, 0]
    set array1_old w/= 1 <- 5;
    Message($"Old: {array1_old}");

    mutable array1_new = [0, size = 3];
    // [0, 5, 0]
    set array1_new[1] = 5;
    Message($"New: {array1_new}");

    mutable array2_old = [[["a", "b", "c"], ["d", "e", "f"]], [["g", "h", "i"], ["j", "k", "l"]]];
    // [[["a", "b", "c"], ["d", "e", "f"]], [["g", "h", "i"], ["d", "e", "x"]]
    set array2_old w/= 1 <- (array2_old[1] w/ 1 <- (array2_old[0][1] w/ 2 <- "x"));
    Message($"Old: {array2_old}");

    mutable array2_new = [[["a", "b", "c"], ["d", "e", "f"]], [["g", "h", "i"], ["j", "k", "l"]]];
    // [[["a", "b", "c"], ["d", "e", "f"]], [["g", "h", "i"], ["d", "e", "x"]]
    set array2_new[1][1][2] = "x";
    Message($"New: {array2_new}");
}