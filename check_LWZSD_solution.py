import os
import sys

def parse_input_file(file_name):
    """Parses the input file and extracts the necessary information, 
    then reconstructs H^T = [I | H'^T]."""

    with open(file_name, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    # The indices of the sections are identified.
    idx_n = lines.index("# n")
    idx_k = lines.index("# k")
    idx_w = lines.index("# w")
    idx_Ht = lines.index("# H^transpose (each line corresponds to a column of H, the identity part is omitted)")
    idx_s = lines.index("# s^transpose")

    # Extraction
    n = int(lines[idx_n + 1])
    k = int(lines[idx_k + 1])
    w = int(lines[idx_w + 1])

    # Part H'^T (each row is a column of H')
    Ht_partial = lines[idx_Ht + 1 : idx_s]
    r = n - k  # Number of lines of the syndrome (and size of identity)

    H_transpose = []
    for i in range(n):
        if i < r:
            # Column i of the identity
            identity_part = ['0'] * r
            identity_part[i] = '1'
            column = ''.join(identity_part)
        else:
            # Column of P at index i - r
            column = Ht_partial[i - r]
        H_transpose.append(column)
    
    # In the display, one row corresponds to one column, which is more practical for calculations.
    H_transpose = [''.join(row) for row in zip(*H_transpose)]

    # print("\nH^transpose ")
    # for col in H_transpose:
    #     print(col)

    # Syndrome
    s_transpose = lines[idx_s + 1]

    return n, k, w, H_transpose, s_transpose





def verify_solution(candidate, H_transpose, s, w, n, z=3):
    """
    Checks if a candidate is a correct solution to the LWSDP.

    Args:
        candidate (list[int] or str): Candidate vector of length n (each value in [0, q-1]).
        H_transpose (list[str]):Transposed H matrix (each element is a row of H^T).
        s (str): Expected syndrome (string of numbers).
        w (int): Maximum permitted weight.
        n (int): Vector length.
        q (int): Finite body used (default 3).

    Returns:
        (bool, str): True if valid + calculated syndrome, False otherwise + calculated syndrome.
    """
    
    # Ensure that the candidate is a list of integers
    if isinstance(candidate, str):
        candidate = [int(x) for x in candidate]

    if len(candidate) != n:
        raise ValueError(f"The candidate vector must be of length {n}, but it is of length {len(candidate)}.")

    # Calculating Hamming's weight
    weight = sum(1 for x in candidate if x != 0)
    print(f"--- Solution weight : {weight} ---")
    if weight < w:
        print(f"Failure: the candidate solution contains {weight} non-zero coefficients (min allowed: {w}).")
        return False, None

    # Convert H_transpose to an int and check the length
    H_int = [[int(x) for x in col] for col in H_transpose]

    # Syndrome calculation
    syndrome = []
    for col in H_int:
        if len(col) != n:
            raise ValueError(f"Column of H_transpose of length {len(col)} != n={n}")
        s_val = sum(c*v for c,v in zip(col, candidate)) % z
        syndrome.append(str(s_val))

    syndrome_str = "".join(syndrome)  # convert the list to a string

    if syndrome_str == s:
        return True, syndrome_str
    else:
        return False, syndrome_str



def main():
    if len(sys.argv) < 3:
        print("Usage : python script.py <fichier_entree> <solution_bianire>")
        sys.exit(1)

    input_file = sys.argv[1]  # file 
    candidate_arg = sys.argv[2]  # binary string
    
    # Parse the input file
    n, k, w, H, s = parse_input_file(input_file)

    candidate = candidate_arg
    print(f"Candidate solution: {candidate}")

    # Verification
    valid, computed_syndrome = verify_solution(candidate, H, s, w, n, z=3)
    if valid:
        print(f"The candidate solution is correct. {computed_syndrome}")
    else:
        print(f"Failure: The candidate solution is incorrect. Calculated syndrome: {computed_syndrome} | Expected syndrome: {s}")

if __name__ == "__main__":
    main()