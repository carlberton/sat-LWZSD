import os
import sys

def parse_input_file(file_name):
    """Parse le fichier d'entrée et extrait les informations nécessaires, 
       puis reconstruit H^T = [I | H'^T]."""
    with open(file_name, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    # On identifie les indices des sections
    idx_n = lines.index("# n")
    idx_k = lines.index("# k")
    idx_w = lines.index("# w")
    idx_Ht = lines.index("# H^transpose (each line corresponds to a column of H, the identity part is omitted)")
    idx_s = lines.index("# s^transpose")

    # Extraction
    n = int(lines[idx_n + 1])
    k = int(lines[idx_k + 1])
    w = int(lines[idx_w + 1])

    # Partie H'^T (chaque ligne est une colonne de H')
    Ht_partial = lines[idx_Ht + 1 : idx_s]
    r = n - k  # Nombre de lignes du syndrome (et taille identité)

    H_transpose = []
    for i in range(n):
        if i < r:
            # Colonne i de l'identité
            identity_part = ['0'] * r
            identity_part[i] = '1'
            column = ''.join(identity_part)
        else:
            # Colonne de P à l’index i - r
            column = Ht_partial[i - r]
        H_transpose.append(column)
    
    # Dans l'affichage, une ligne correspond à une colonne, plus pratique pour les calculs
    H_transpose = [''.join(row) for row in zip(*H_transpose)]

    # print("\nH^transpose ")
    # for col in H_transpose:
    #     print(col)

    # Syndrome
    s_transpose = lines[idx_s + 1]

    return n, k, w, H_transpose, s_transpose





def verify_solution(candidate, H_transpose, s, w, n, z=3):
    """
    Vérifie si un candidat est une solution correcte au LWSDP.

    Args:
        candidate (list[int] or str): Vecteur candidat de longueur n (chaque valeur dans [0, q-1]).
        H_transpose (list[str]): Matrice H transposée (chaque élément est une ligne de H^T).
        s (str): Syndrome attendu (chaîne de chiffres).
        w (int): Poids maximal autorisé.
        n (int): Longueur du vecteur.
        q (int): Corps fini utilisé (par défaut 3).

    Returns:
        (bool, str): True si valide + syndrome calculé, False sinon + syndrome calculé.
    """
    
    # Assurer que le candidat est une liste d'entiers
    if isinstance(candidate, str):
        candidate = [int(x) for x in candidate]

    if len(candidate) != n:
        raise ValueError(f"Le vecteur candidat doit être de longueur {n}, or il est de longueur {len(candidate)}.")

    # Calcul du poids de Hamming
    weight = sum(1 for x in candidate if x != 0)
    print(f"--- Poids de la solution : {weight} ---")
    if weight < w:
        print(f"Échec : la solution candidate contient {weight} coefficients non nulles (min autorisé : {w}).")
        return False, None

    # Convertir H_transpose en int et vérifier la longueur
    H_int = [[int(x) for x in col] for col in H_transpose]

    # Calcul du syndrome
    syndrome = []
    for col in H_int:
        if len(col) != n:
            raise ValueError(f"Colonne de H_transpose de longueur {len(col)} != n={n}")
        s_val = sum(c*v for c,v in zip(col, candidate)) % z
        syndrome.append(str(s_val))

    syndrome_str = "".join(syndrome)  # convertir la liste en chaîne

    if syndrome_str == s:
        return True, syndrome_str
    else:
        return False, syndrome_str



def main():
    if len(sys.argv) < 3:
        print("Usage : python script.py <fichier_entree> <solution_bianire>")
        sys.exit(1)

    input_file = sys.argv[1]  # fichier 
    candidate_arg = sys.argv[2]  # une chaîne binaire
    
    # Parse le fichier d'entrée
    n, k, w, H, s = parse_input_file(input_file)

    candidate = candidate_arg
    print(f"Solution candidate : {candidate}")

    # Vérification
    valid, computed_syndrome = verify_solution(candidate, H, s, w, n, z=3)
    if valid:
        print(f"The candidate solution is correct. {computed_syndrome}")
    else:
        print(f"Échec : La solution candidate n'est pas correcte. Syndrome calculé : {computed_syndrome} | Syndrome attendu : {s}")

if __name__ == "__main__":
    main()