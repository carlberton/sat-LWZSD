import os
import re
import subprocess

def parse_input_file(file_name):
    """
    Parse the input file and extract the values n, seed, k, w, H^T, and s^T.
    """

    with open(file_name, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    # On identifie les indices des sections
    idx_n = lines.index("# n")
    idx_seed = lines.index("# seed")
    idx_k = lines.index("# k")
    idx_w = lines.index("# w")
    idx_Ht = lines.index("# H^transpose (each line corresponds to a column of H, the identity part is omitted)")
    idx_s = lines.index("# s^transpose")

    # Extraction
    n = int(lines[idx_n + 1])
    print(f"n = {n}")
    seed = int(lines[idx_seed + 1])
    print(f"seed = {seed}")
    k = int(lines[idx_k + 1])
    print(f"k = {k}")
    w = int(lines[idx_w + 1])
    print(f"w = {w}")

    # Toutes les lignes entre H^T et s^T
    H_transpose = lines[idx_Ht + 1 : idx_s]

    # Transposition 
    H_transpose = [''.join(row) for row in zip(*H_transpose)]
    print("H^T (transposé) :")
    for row in H_transpose:
        print(row)

    # Syndrome
    s_transpose = lines[idx_s + 1]
    print("s^T :")
    print(s_transpose)

    return n, seed, w, k, H_transpose, s_transpose


def build_var_sets(H_transpose, s_transpose, n, k, t, z=3):
    """
    Construct the sets V and K.
    """

    m = n - k   # nombre d'équations
    
    # Build V_{E_i} sets for each equation
    V = []
    for i, row in enumerate(H_transpose):
        V_i = [i+1]  # Identity variable
        for j, val in enumerate(row):
            if val != '0':
                V_i.append(j + m + 1)
        V.append(V_i)

    # Build K_{E_i} sets (valid cardinalities mod z)
    K = []
    for i in range(m):
        s_i = int(s_transpose[i])
        # Récupérer les coefficients non nuls de l'équation E_i
        coeffs = []
        for j in V[i]:
            if j <= m:
                # Partie identité
                H_ij = 1
            else:
                # Partie H_transpose
                col = j - m - 1
                H_ij = int(H_transpose[i][col])

            assert H_ij != 0
            coeffs.append(H_ij)

        # Nombre minimal de variables non nulles dans E_i
        min_nonzero = max(0, len(V[i]) - (n - t))  # borne locale >= 0

        # Trier les coefficients pour construire la borne inférieure
        coeffs_sorted = sorted(coeffs)
        v_min = sum(coeffs_sorted[:min_nonzero])  # somme des plus petits coefficients
        while v_min % z != s_i:
            v_min += 1  # ajuster pour respecter la congruence mod z
        v_max = (z - 1) * sum(coeffs)             # borne sup. inchangée
        while v_max % z != s_i:
            v_max -= 1 # ajuster pour respecter la congruence mod z

        # Construire l’ensemble K_{E_i} en tenant compte de la borne inférieure
        K_i = K_i = list(range(v_min, v_max + 1, z))
        K.append(K_i)

    return V, K


def write_cnf_to_file(input_file, cc_encoding, pb_encoding, cnf, seed, variant):
    """
    Save CNF clauses to a file in a structured output directory.
    """

    # Generate the output file name based on input and encoding types
    file_name, _ = os.path.splitext(input_file)
    output_file = f"./Challenges/seed_{seed}/{variant}/PB_{pb_encoding}/encoding_{cc_encoding}/{file_name.split('/')[-1]}.cnf"
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)  # Create directories if they don't exist
    
    # Save CNF clauses to the file
    cnf.to_file(output_file)
    
    print(f"CNF clauses have been saved to {output_file}")


def extract_LW_n(filename):
    """
    Extrait le nombre 'n' des fichiers "LargeWeight_n_seed"
    Exemple : "LargeWeight_130_5" -> 130
    """
    # On cherche "LargeWeight_" suivi par un ou plusieurs chiffres
    match = re.search(r"LargeWeight_(\d+)", filename)
    if match:
        return int(match.group(1))
    return float('inf')


def verify_sol(input_file, solution) :
    """
    Verifies the correctness of a candidate solution for a Syndrome Decoding Problem (SDP)
    by calling an external checker script.

    Args:
        input_file (str): Path to the input file describing the SDP instance.
        solution (str): A binary string representing the candidate solution.

    Returns:
        bool: True if the solution is correct according to the checker, False otherwise.
    """

    if solution:
        try:
            # Call the external Python script to check the solution
            result = subprocess.run(
                ['python3', 'check_LWZSD_solution.py', input_file, solution],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True
            )
            output = result.stdout.strip()

            # Check if the checker script confirms correctness
            if "The candidate solution is correct." in output:
                return True
            else:
                return False
        except subprocess.CalledProcessError as e:
            # If the subprocess fails (e.g., script crashes), print the error
            print("Error during external verification:", e.stderr)
    else:
        print("No solution found.")