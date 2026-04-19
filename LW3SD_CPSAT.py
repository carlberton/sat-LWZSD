import time
import csv
import argparse
from ortools.sat.python import cp_model
from utils import *

def build_and_solve_CP1(n, w, k, H_transpose, s_transpose, Z=3, forward=False, equiv=False, timeout=3600):
    model = cp_model.CpModel()
    m = n - k

    # Variables d'erreur e_{j,c} 
    # e_vars[(j, c)] correspond à e_{j,c} dans le modèle
    # j de 1 à n, c de 1 à Z-1
    e_vars = {}
    for j in range(1, n + 1):
        for c in range(1, Z):
            e_vars[(j, c)] = model.NewBoolVar(f"e_{j}_{c}")

    # Chaîne d'implication (Unary encoding) 
    # -e_{j,v} V e_{j,v-1}  <=>  e_{j,v} => e_{j,v-1}
    for j in range(1, n + 1):
        for v in range(2, Z):
            model.AddImplication(e_vars[(j, v)], e_vars[(j, v-1)])

    # Build sets V and K
    V, K = build_var_sets(H_transpose, s_transpose, n, k, w, Z)

    # Pour chaque équation E_i
    for i in range(m):
        # On ne crée plus que les variables x_{i,v}
        x_vars_dict = {}
        for v in K[i]:
            x_vars_dict[v] = model.NewBoolVar(f"x_{i}_{v}")

        # Exactement une variable x_{i,v} est à 1
        model.AddExactlyOne(x_vars_dict.values())

        # Équation Pseudo-Booléenne
        lhs_terms = []
        
        # Partie Erreurs : H_ij * e_{j,c}
        for j_idx in V[i]:
            if j_idx <= m:
                H_ij = 1
            else:
                col = j_idx - m - 1
                H_ij = int(H_transpose[i][col])
            
            if H_ij != 0:
                for c in range(1, Z):
                    lhs_terms.append(e_vars[(j_idx, c)] * H_ij)

        # Partie Auxiliaire : v * Not(x_{i,v})
        for v, xv in x_vars_dict.items():
            if v != 0:
                lhs_terms.append(xv.Not() * int(v))

        rhs = sum(int(v) for v in K[i])
        model.Add(sum(lhs_terms) == rhs)

    # Poids de Hamming 
    t = w 
    model.Add(sum(e_vars[(j, 1)] for j in range(1, n + 1)) >= t)

    # Résolution
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = timeout
    
    start_time = time.time()
    status = solver.Solve(model)
    res_time = time.time() - start_time

    solution = None
    if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
        status_str = 'sat'
        # On reconstruit la valeur de l'erreur pour chaque j
        # La valeur de e_j est la somme des e_{j,c}
        sol_list = []
        for j in range(1, n + 1):
            val_j = sum(solver.Value(e_vars[(j, c)]) for c in range(1, Z))
            sol_list.append(val_j)
        solution = ''.join(map(str, sol_list))
    elif status == cp_model.INFEASIBLE:
        status_str = 'unsat'
    else:
        status_str = 'timeout'

    return status_str, f"{res_time:.5f}", solution


def build_and_solve_CP2(n, w, k, H_transpose, s_transpose, Z=3,forward=False, equiv=False, timeout=3600):
    model = cp_model.CpModel()
    m = n - k

    # Variables d'erreur e_{j,c}
    e_vars = {}
    for j in range(1, n + 1):
        for c in range(1, Z):
            e_vars[(j, c)] = model.NewBoolVar(f"e_{j}_{c}")

    # Chaîne d'implication pour e
    for j in range(1, n + 1):
        for v in range(2, Z):
            model.AddImplication(e_vars[(j, v)], e_vars[(j, v-1)])

    # Build sets V and K
    V, K = build_var_sets(H_transpose, s_transpose, n, k, w, Z)

    # Variables auxiliaires x_{i,v} (Encodage Unaire) 
    x_vars_dict = {}
    for i in range(m):
        for v in K[i]:
            x_vars_dict[(i, v)] = model.NewBoolVar(f"x_{i}_{v}")

        # Chaîne d'implication pour x (¬x_{i,v} ∨ x_{i,v-Z})
        # Équivalent à : x_{i,v} => x_{i,v-Z}
        for v in K[i]:
            prev_v = v - Z
            if prev_v in K[i]:
                model.AddImplication(x_vars_dict[(i, v)], x_vars_dict[(i, prev_v)])

    # Pour chaque équation E_i
    for i in range(m):
        v_min = min(K[i])
        v_max = max(K[i])
        
        lhs_terms = []

        # Contribution des erreurs (Double somme j et c)
        for j_idx in V[i]:
            if j_idx <= m:
                H_ij = 1
            else:
                col = j_idx - m - 1
                H_ij = int(H_transpose[i][col])
            
            if H_ij != 0:
                for c in range(1, Z):
                    lhs_terms.append(e_vars[(j_idx, c)] * H_ij)

        # Contribution des variables auxiliaires : Z * bar_x_{i,v}
        for v in K[i]:
            if v == v_min:
                continue
            
            xv = x_vars_dict[(i, v)]
            lhs_terms.append(xv.Not() * Z)

        model.Add(sum(lhs_terms) == int(v_max))

    # Poids de Hamming 
    t = w
    model.Add(sum(e_vars[(j, 1)] for j in range(1, n + 1)) >= t)

    # Résolution
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = timeout
    
    start_time = time.time()
    status = solver.Solve(model)
    res_time = time.time() - start_time

    solution = None
    if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
        status_str = 'sat'
        sol_list = []
        for j in range(1, n + 1):
            val_j = sum(solver.Value(e_vars[(j, c)]) for c in range(1, Z))
            sol_list.append(val_j)
        solution = ''.join(map(str, sol_list))
    elif status == cp_model.INFEASIBLE:
        status_str = 'unsat'
    else:
        status_str = 'timeout'

    return status_str, f"{res_time:.5f}", solution

def build_and_solve_CP3(n, w, k, H_transpose, s_transpose, Z=3, forward=False, equiv=False, timeout=3600):
    model = cp_model.CpModel()
    m = n - k

    # Variables d'erreur e_{j,c} 
    e_vars = {}
    for j in range(1, n + 1):
        for c in range(1, Z):
            e_vars[(j, c)] = model.NewBoolVar(f"e_{j}_{c}")

    for j in range(1, n + 1):
        for v in range(2, Z):
            model.AddImplication(e_vars[(j, v)], e_vars[(j, v-1)])

    # Build sets V and K
    V, K = build_var_sets(H_transpose, s_transpose, n, k, w, Z)

    # Variables auxiliaires x_{i,l} (Représentation binaire du quotient q) 
    x_vars_dict = {} 
    for i in range(m):
        v_max = max(K[i])
        q_max = (v_max - int(s_transpose[i])) // Z
        L_i = max(1, q_max.bit_length())
        for l in range(L_i):
            x_vars_dict[(i, l)] = model.NewBoolVar(f"x_{i}_{l}")

    # Pour chaque équation E_i
    for i in range(m):
        v_max = max(K[i])
        q_max = (v_max - int(s_transpose[i])) // Z
        L_i = max(1, q_max.bit_length())
        
        lhs_terms = []

        # Contribution des erreurs
        for j_idx in V[i]:
            H_ij = 1 if j_idx <= m else int(H_transpose[i][j_idx - m - 1])
            if H_ij != 0:
                for c in range(1, Z):
                    lhs_terms.append(e_vars[(j_idx, c)] * H_ij)

        # Contribution des bits du quotient 
        for l in range(L_i):
            xv = x_vars_dict[(i, l)]
            lhs_terms.append(xv.Not() * (Z * (2**l)))

        rhs = int(s_transpose[i]) + Z * (2**L_i - 1)
        model.Add(sum(lhs_terms) == rhs)

        # Filtrage Exhaustif (Exclusion des valeurs interdites) 
        # On identifie les valeurs de v (quotient) qui ne sont pas dans K[i]
        # s_i + Z*v doit être dans K[i]
        for v in range(2**L_i):
            if (int(s_transpose[i]) + Z * v) not in K[i]:
                # On ajoute une clause pour interdire la combinaison de bits correspondant à v
                # La clause doit être vraie si AU MOINS UN bit est différent de la config de v
                forbidden_clause = []
                for l in range(L_i):
                    bit = (v >> l) & 1 # Récupère le l-ième bit de v
                    if bit == 0:
                        # Si le bit de v est 0, on veut que l'un des bits soit 1 (x)
                        forbidden_clause.append(x_vars_dict[(i, l)])
                    else:
                        # Si le bit de v est 1, on veut que l'un des bits soit 0 (Not(x))
                        forbidden_clause.append(x_vars_dict[(i, l)].Not())
                
                # model.AddBoolOr([Lits]) force au moins un littéral à être vrai
                model.AddBoolOr(forbidden_clause)

    # Contrainte de poids de Hamming 
    model.Add(sum(e_vars[(j, 1)] for j in range(1, n + 1)) >= w)

    # Résolution
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = timeout
    
    start_time = time.time()
    status = solver.Solve(model)
    res_time = time.time() - start_time

    solution = None
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        status_str = 'sat'
        sol_list = [sum(solver.Value(e_vars[(j, c)]) for c in range(1, Z)) for j in range(1, n+1)]
        solution = ''.join(map(str, sol_list))
    else:
        status_str = 'unsat' if status == cp_model.INFEASIBLE else 'timeout'

    return status_str, f"{res_time:.5f}", solution


def build_and_solve_CP4(n, w, k, H_transpose, s_transpose, Z=3, forward=False, equiv=False, timeout=3600):
    model = cp_model.CpModel()
    m = n - k

    # Variables d'erreur e_{j,c}
    e_vars = {}
    for j in range(1, n + 1):
        for c in range(1, Z):
            e_vars[(j, c)] = model.NewBoolVar(f"e_{j}_{c}")

    # Encodage unaire 
    for j in range(1, n + 1):
        for v in range(2, Z):
            model.AddImplication(e_vars[(j, v)], e_vars[(j, v-1)])

    # Chargement des ensembles V et K via utils
    V, K = build_var_sets(H_transpose, s_transpose, n, k, w, Z)

    # Variables auxiliaires x_{i,l} (Bits du quotient)
    x_vars_dict = {}
    for i in range(m):
        v_max = max(K[i])
        q_max = (v_max - int(s_transpose[i])) // Z
        L_i = max(1, q_max.bit_length())
        for l in range(L_i):
            x_vars_dict[(i, l)] = model.NewBoolVar(f"x_{i}_{l}")

    # Pour chaque équation E_i
    for i in range(m):
        v_min, v_max = min(K[i]), max(K[i])
        s_i = int(s_transpose[i])
        q_min, q_max = (v_min - s_i) // Z, (v_max - s_i) // Z
        L_i = max(1, q_max.bit_length())

        # Équation Pseudo-Booléenne 
        lhs_terms = []
        for j_idx in V[i]:
            H_ij = 1 if j_idx <= m else int(H_transpose[i][j_idx - m - 1])
            if H_ij != 0:
                for c in range(1, Z):
                    lhs_terms.append(e_vars[(j_idx, c)] * H_ij)

        for l in range(L_i):
            lhs_terms.append(x_vars_dict[(i, l)].Not() * (Z * (2**l)))

        rhs = s_i + Z * (2**L_i - 1)
        model.Add(sum(lhs_terms) == rhs)

        # Filtrage Compact (Bornes q_max / q_min) 
        bits_max = [(q_max >> l) & 1 for l in range(L_i)]
        bits_min = [(q_min >> l) & 1 for l in range(L_i)]
        
        J_max = sorted([l for l in range(L_i) if bits_max[l] == 1], reverse=True)
        J_min = sorted([l for l in range(L_i) if bits_min[l] == 0], reverse=True)

        # Exclusion q > q_max 
        pmax = {}

        if J_max:
            for j in J_max:
                if j > 0 and bits_max[j - 1] == 0:
                    pmax[j] = model.NewBoolVar(f"pmax_{i}_{j}")
            
            block_pmax = sorted(pmax.keys(), reverse=True)

            # Implication de chaîne et segments
            if not forward or equiv:
                for idx, j in enumerate(block_pmax):
                    p_j = pmax[j]
                    if idx == 0:
                        for k in range(j, L_i):
                            model.AddImplication(p_j, x_vars_dict[(i, k)])
                    else:
                        j_prime = block_pmax[idx - 1]
                        model.AddImplication(p_j, pmax[j_prime])
                        for k in range(j, j_prime):
                            if bits_max[k] == 1:
                                model.AddImplication(p_j, x_vars_dict[(i, k)])
                            else:
                                model.AddImplication(p_j, x_vars_dict[(i, k)].Not())

            # Forward direction: (prefix match) => p_j
            if forward:
                for idx, j in enumerate(block_pmax):
                    p_j = pmax[j]
                    if idx == 0:
                        clause = [p_j]
                        for k in range(j, L_i):
                            xk = x_vars_dict[(i, k)]
                            if bits_max[k] == 1:
                                clause.append(xk.Not())
                            else:
                                clause.append(xk)
                        model.AddBoolOr(clause)
                    else:
                        j_prime = block_pmax[idx - 1]
                        clause = [p_j, pmax[j_prime].Not()]

                        for k in range(j, j_prime):
                            xk = x_vars_dict[(i, k)]
                            if bits_max[k] == 1:
                                clause.append(xk.Not())
                            else:
                                clause.append(xk)
                        model.AddBoolOr(clause)

            # Clauses d'élimination (7a et 7b du modèle mathématique)
            for idx, j in enumerate(block_pmax):
                p_j = pmax[j]
                if idx < len(block_pmax) - 1:
                    j2 = block_pmax[idx + 1]
                    for k in range(j2 + 1, j):
                        model.AddImplication(p_j, x_vars_dict[(i, k)].Not())
                else:
                    last = 0
                    if J_max and min(J_max) == 0:
                        for jj in sorted(J_max):
                            if bits_max[jj] == 0: break
                            last = jj
                    for k in range(last, j):
                        model.AddImplication(p_j, x_vars_dict[(i, k)].Not())

        # Exclusion q < q_min 
        pmin = {}
        if J_min:
            # (7c) Bits au dessus du MSB à 0
            j_msb0 = J_min[0]
            if j_msb0 != L_i - 1:
                for k in range(j_msb0 + 1, L_i):
                    model.Add(x_vars_dict[(i, k)] == 1)

            for j in J_min:
                if j > 0 and bits_min[j - 1] == 1:
                    pmin[j] = model.NewBoolVar(f"pmin_{i}_{j}")
            
            block_pmin = sorted(pmin.keys(), reverse=True)

            if not forward or equiv:
                for idx, j in enumerate(block_pmin):
                    p_j = pmin[j]
                    if idx == 0:
                        for k in range(j, L_i):
                            if bits_min[k] == 0:
                                model.AddImplication(p_j, x_vars_dict[(i, k)].Not())
                            else:
                                model.AddImplication(p_j, x_vars_dict[(i, k)])
                    else:
                        j_prime = block_pmin[idx - 1]
                        model.AddImplication(p_j, pmin[j_prime])
                        for k in range(j, j_prime):
                            if bits_min[k] == 0:
                                model.AddImplication(p_j, x_vars_dict[(i, k)].Not())
                            else:
                                model.AddImplication(p_j, x_vars_dict[(i, k)])

            if forward:
                for idx, j in enumerate(block_pmin):
                    p_j = pmin[j]
                    if idx == 0:
                        clause = [p_j]
                        for k in range(j, L_i):
                            xk = x_vars_dict[(i, k)]
                            if bits_min[k] == 0:
                                clause.append(xk)
                            else:
                                clause.append(xk.Not())
                        model.AddBoolOr(clause)
                    else:
                        j_prime = block_pmin[idx - 1]
                        clause = [p_j, pmin[j_prime].Not()]
                        for k in range(j, j_prime):
                            xk = x_vars_dict[(i, k)]
                            if bits_min[k] == 0:
                                clause.append(xk)
                            else:
                                clause.append(xk.Not())
                        model.AddBoolOr(clause)

            for idx, j in enumerate(block_pmin):
                p_j = pmin[j]
                if idx < len(block_pmin) - 1:
                    j2 = block_pmin[idx + 1]
                    for k in range(j2 + 1, j):
                        model.AddImplication(p_j, x_vars_dict[(i, k)])
                else:
                    last = 0
                    if J_min and min(J_min) == 0:
                        for jj in sorted(J_min):
                            if bits_min[jj] == 1: break
                            last = jj
                    for k in range(last, j):
                        model.AddImplication(p_j, x_vars_dict[(i, k)])
        else:
            # Si q_min n'a aucun bit à 0 (cas q_min = 2^L - 1)
            for k in range(L_i):
                model.Add(x_vars_dict[(i, k)] == 1)

    # Contrainte de Poids de Hamming
    model.Add(sum(e_vars[(j, 1)] for j in range(1, n + 1)) >= w)

    # Résolution 
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = timeout
    
    start_time = time.time()
    status = solver.Solve(model)
    res_time = time.time() - start_time

    solution = None
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        status_str = 'sat'
        sol_list = [sum(solver.Value(e_vars[(j, c)]) for c in range(1, Z)) for j in range(1, n+1)]
        solution = ''.join(map(str, sol_list))
    elif status == cp_model.INFEASIBLE:
        status_str = 'unsat'
    else:
        status_str = 'timeout'

    return status_str, f"{res_time:.5f}", solution


def process_file(file_path, solve_function, forward=False, equiv=False):
    """
    Process a single input file and solve the syndrome decoding problem.

    Args:
        file_path (str): Path to the input file containing problem parameters.
        solve_function (function): Function to solve the problem.
        forward (bool): Use forward filtering.
        equiv (bool): Use equivalence filtering.
    
    Returns:
        tuple: (file, status, res_time, sol)
    """

    # Parse the input file
    n, seed, w, k, H_transpose, s_transpose = parse_input_file(file_path)
    
    status, res_time, sol = solve_function(
        n, w, k,
        H_transpose,
        s_transpose,
        Z=3,
        forward=forward,
        equiv=equiv
    )
    
    file = os.path.basename(file_path)

    # Verification
    if status == 'sat' and sol is not None:
        is_valid = verify_sol(file_path, sol)
        if is_valid: 
            return file, status, res_time, sol
        else:
            return file, status, res_time, "Invalid solution"
    else:
        return file, status, res_time, "No solution"



def main():
    # Configuration du dictionnaire de méthodes
    methods = {
        'CNF1': build_and_solve_CP1,
        'CNF2': build_and_solve_CP2,
        'CNF3': build_and_solve_CP3,
        'CNF4': build_and_solve_CP4
    }

    # Configuration de l'argument parser
    parser = argparse.ArgumentParser(description="CPSAT solver for syndrome decoding.")
    parser.add_argument('-m', '--method', choices=methods.keys(), required=True, 
                        help="Resolution method to use (CNF1, CNF2, CNF3 or CNF4)")
    
    parser.add_argument(
        "--forward",
        action="store_true",
        help="Use forward implication only"
    )

    parser.add_argument(
        "--equiv",
        action="store_true",
        help="Use equivalence (forward + backward)"
    )
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', '--file', help='Path to an input file')
    group.add_argument('-d', '--dir', help='Path to a directory to process')
    args = parser.parse_args()

    forward = args.forward
    equiv = args.equiv

    # Déterminer suffixe pour CNF4
    suffix = ""
    if args.method == "CNF4":
        if args.equiv:
            suffix = "_equiv"
        elif args.forward:
            suffix = "_forward"

    # Logique des modes
    if equiv:
        forward = True

    # Sélection de la fonction de résolution
    solve_function = methods[args.method]

    if args.file:
        # Traitement d'un fichier unique
        file, status, res_time, sol = process_file(args.file, solve_function, forward=forward, equiv=equiv)
        directory = os.path.dirname(args.file)
        csv_filepath = os.path.join(directory, f"CPSAT_{args.method}{suffix}_{os.path.splitext(file)[0]}.csv")
        
        with open(csv_filepath, mode='w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(["File", "Result", "Time (s)", "Solution"])
            csv_writer.writerow([file, status, res_time, sol])
            csvfile.flush()
    else:
        # Traitement d'un répertoire complet
        csv_filepath = os.path.join(args.dir, f"CPSAT_{args.method}{suffix}.csv")
        
        # On ouvre le fichier CSV en mode append (au cas où tu relances) ou write
        with open(csv_filepath, mode='w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(["File", "Result", "Time (s)", "Solution"])
            csvfile.flush()

            # Lister et filtrer les fichiers d'instances
            all_files = os.listdir(args.dir)
            entries = [
                f for f in all_files 
                if os.path.isfile(os.path.join(args.dir, f)) 
                and f.startswith("LargeWeight_") 
                and not f.endswith(".csv")
            ]
            
            # Tri
            entries.sort(key=extract_LW_n)

            if not entries:
                print(f"Aucun fichier d'instance 'LargeWeight_' trouvé dans : {args.dir}")
                return

            print(f"Début du traitement de {len(entries)} fichiers dans l'ordre croissant de n.")

            # Boucle de traitement
            for entry in entries:
                path = os.path.join(args.dir, entry)
                n_val = extract_LW_n(entry)
                
                print(f"\n[n={n_val}] Traitement de : {entry} ...")
                
                # Exécution du solver CP-SAT sélectionné
                file_name, status, res_time, sol = process_file(path, solve_function, forward=forward, equiv=equiv)
                
                # Écriture immédiate dans le CSV (sécurité si crash)
                csv_writer.writerow([file_name, status, res_time, sol])
                csvfile.flush()
                
                print(f"Terminé : {status} en {res_time}s")

    print(f"\nTous les fichiers ont été traités. Résultats : {csv_filepath}")


if __name__ == "__main__":
    main()