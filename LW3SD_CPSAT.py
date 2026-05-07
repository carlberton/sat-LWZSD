import time
import csv
import argparse
from ortools.sat.python import cp_model
from utils import *

def build_and_solve_CP1(n, w, k, H_transpose, s_transpose, Z=3, forward=False, backward=False, timeout=3600):
    model = cp_model.CpModel()
    m = n - k

    # Error variables e_{j,c}
    # e_vars[(j, c)] corresponds to e_{j,c} in the model
    # j from 1 to n, c from 1 to Z-1
    e_vars = {}
    for j in range(1, n + 1):
        for c in range(1, Z):
            e_vars[(j, c)] = model.NewBoolVar(f"e_{j}_{c}")

    # Chain of implication (Unary encoding)
    for j in range(1, n + 1):
        for v in range(2, Z):
            model.AddImplication(e_vars[(j, v)], e_vars[(j, v-1)])

    # Build sets V and K
    V, K = build_var_sets(H_transpose, s_transpose, n, k, w, Z)

    # For each equation E_i
    for i in range(m):
        # We only create the variables x_{i,v}
        x_vars_dict = {}
        for v in K[i]:
            x_vars_dict[v] = model.NewBoolVar(f"x_{i}_{v}")

        # Exactly one variable x_{i,v} is at 1
        model.AddExactlyOne(x_vars_dict.values())

        # Pseudo-Boolean Equation
        lhs_terms = []
        
        # Error section: H_ij * e_{j,c}
        for j_idx in V[i]:
            if j_idx <= m:
                H_ij = 1
            else:
                col = j_idx - m - 1
                H_ij = int(H_transpose[i][col])
            
            if H_ij != 0:
                for c in range(1, Z):
                    lhs_terms.append(e_vars[(j_idx, c)] * H_ij)

        # Auxiliary Part: v * Not(x_{i,v})
        for v, xv in x_vars_dict.items():
            if v != 0:
                lhs_terms.append(xv.Not() * int(v))

        rhs = sum(int(v) for v in K[i])
        model.Add(sum(lhs_terms) == rhs)

    # Hamming weight
    t = w 
    model.Add(sum(e_vars[(j, 1)] for j in range(1, n + 1)) >= t)

    # Resolution
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = timeout
    
    start_time = time.time()
    status = solver.Solve(model)
    res_time = time.time() - start_time

    solution = None
    if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
        status_str = 'sat'
        # We reconstruct the error value for each j
        # The value of e_j is the sum of e_{j,c}
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


def build_and_solve_CP2(n, w, k, H_transpose, s_transpose, Z=3,forward=False, backward=False, timeout=3600):
    model = cp_model.CpModel()
    m = n - k

    # Error variables e_{j,c}
    e_vars = {}
    for j in range(1, n + 1):
        for c in range(1, Z):
            e_vars[(j, c)] = model.NewBoolVar(f"e_{j}_{c}")

    # Chain of involvement for e
    for j in range(1, n + 1):
        for v in range(2, Z):
            model.AddImplication(e_vars[(j, v)], e_vars[(j, v-1)])

    # Build sets V and K
    V, K = build_var_sets(H_transpose, s_transpose, n, k, w, Z)

    # Auxiliary variables x_{i,v} (Unary Encoding)
    x_vars_dict = {}
    for i in range(m):
        for v in K[i]:
            x_vars_dict[(i, v)] = model.NewBoolVar(f"x_{i}_{v}")

        # Chain of implication for x
        for v in K[i]:
            prev_v = v - Z
            if prev_v in K[i]:
                model.AddImplication(x_vars_dict[(i, v)], x_vars_dict[(i, prev_v)])

    # For each equation E_i
    for i in range(m):
        v_min = min(K[i])
        v_max = max(K[i])
        
        lhs_terms = []

        # Contribution of errors
        for j_idx in V[i]:
            if j_idx <= m:
                H_ij = 1
            else:
                col = j_idx - m - 1
                H_ij = int(H_transpose[i][col])
            
            if H_ij != 0:
                for c in range(1, Z):
                    lhs_terms.append(e_vars[(j_idx, c)] * H_ij)

        # Contribution of auxiliary variables
        for v in K[i]:
            if v == v_min:
                continue
            
            xv = x_vars_dict[(i, v)]
            lhs_terms.append(xv.Not() * Z)

        model.Add(sum(lhs_terms) == int(v_max))

    # Hamming weight
    t = w
    model.Add(sum(e_vars[(j, 1)] for j in range(1, n + 1)) >= t)

    # Resolution
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


def build_and_solve_CP3(n, w, k, H_transpose, s_transpose, Z=3, forward=False, backward=False, timeout=3600):
    model = cp_model.CpModel()
    m = n - k

    # Error variables e_{j,c}
    e_vars = {}
    for j in range(1, n + 1):
        for c in range(1, Z):
            e_vars[(j, c)] = model.NewBoolVar(f"e_{j}_{c}")

    for j in range(1, n + 1):
        for v in range(2, Z):
            model.AddImplication(e_vars[(j, v)], e_vars[(j, v-1)])

    # Build sets V and K
    V, K = build_var_sets(H_transpose, s_transpose, n, k, w, Z)

    # Auxiliary variables x_{i,l} (Binary representation of the quotient q)
    x_vars_dict = {} 
    for i in range(m):
        v_max = max(K[i])
        q_max = (v_max - int(s_transpose[i])) // Z
        L_i = max(1, q_max.bit_length())
        for l in range(L_i):
            x_vars_dict[(i, l)] = model.NewBoolVar(f"x_{i}_{l}")

    # For each equation E_i
    for i in range(m):
        v_max = max(K[i])
        q_max = (v_max - int(s_transpose[i])) // Z
        L_i = max(1, q_max.bit_length())
        
        lhs_terms = []

        # Contribution of errors
        for j_idx in V[i]:
            H_ij = 1 if j_idx <= m else int(H_transpose[i][j_idx - m - 1])
            if H_ij != 0:
                for c in range(1, Z):
                    lhs_terms.append(e_vars[(j_idx, c)] * H_ij)

        # Contribution of the bits of the quotient
        for l in range(L_i):
            xv = x_vars_dict[(i, l)]
            lhs_terms.append(xv.Not() * (Z * (2**l)))

        rhs = int(s_transpose[i]) + Z * (2**L_i - 1)
        model.Add(sum(lhs_terms) == rhs)

        # Exhaustive Filtering (Exclusion of forbidden values)
        # We identify the values ​​of v (quotient) that are not in K[i]
        # s_i + Z*v must be in K[i]
        for v in range(2**L_i):
            if (int(s_transpose[i]) + Z * v) not in K[i]:
                # We add a clause to prohibit the combination of bits corresponding to v
                # The clause must be true if AT LEAST ONE bit is different from the v configuration
                forbidden_clause = []
                for l in range(L_i):
                    bit = (v >> l) & 1 # Retrieves the lth bit of v
                    if bit == 0:
                        # If the bit of v is 0, we want one of the bits to be 1 (x)
                        forbidden_clause.append(x_vars_dict[(i, l)])
                    else:
                        # If the bit of v is 1, we want one of the bits to be 0 (Not(x))
                        forbidden_clause.append(x_vars_dict[(i, l)].Not())
                
                # model.AddBoolOr([Lits]) forces at least one literal to be true
                model.AddBoolOr(forbidden_clause)

    # Hamming weight constraint
    model.Add(sum(e_vars[(j, 1)] for j in range(1, n + 1)) >= w)

    # Resolution
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


def build_and_solve_CP4(n, w, k, H_transpose, s_transpose, Z=3, forward=False, backward=False, timeout=3600):
    model = cp_model.CpModel()
    m = n - k

    # Error variables e_{j,c}
    e_vars = {}
    for j in range(1, n + 1):
        for c in range(1, Z):
            e_vars[(j, c)] = model.NewBoolVar(f"e_{j}_{c}")

    # Unary encoding 
    for j in range(1, n + 1):
        for v in range(2, Z):
            model.AddImplication(e_vars[(j, v)], e_vars[(j, v-1)])

    # Build sets
    V, K = build_var_sets(H_transpose, s_transpose, n, k, w, Z)

    # Auxiliary variables
    x_vars_dict = {}
    for i in range(m):
        v_max = max(K[i])
        q_max = (v_max - int(s_transpose[i])) // Z
        L_i = max(1, q_max.bit_length())
        for l in range(L_i):
            x_vars_dict[(i, l)] = model.NewBoolVar(f"x_{i}_{l}")

    # For each equation E_i
    for i in range(m):
        v_min, v_max = min(K[i]), max(K[i])
        s_i = int(s_transpose[i])
        q_min, q_max = (v_min - s_i) // Z, (v_max - s_i) // Z
        L_i = max(1, q_max.bit_length())

        # Pseudo-Boolean constraints
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

        # Filtering of forbidden values (Exclusion of q > q_max and q < q_min)
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

            # Implication 
            if backward:
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

            # Forward clause
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

            # Elimination of forbidden values between blocks
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
            # Bits above MSB set to 0
            j_msb0 = J_min[0]
            if j_msb0 != L_i - 1:
                for k in range(j_msb0 + 1, L_i):
                    model.Add(x_vars_dict[(i, k)] == 1)

            for j in J_min:
                if j > 0 and bits_min[j - 1] == 1:
                    pmin[j] = model.NewBoolVar(f"pmin_{i}_{j}")
            
            block_pmin = sorted(pmin.keys(), reverse=True)

            if backward:
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
            # If q_min has no bits set to 0 (case q_min = 2^L - 1)
            for k in range(L_i):
                model.Add(x_vars_dict[(i, k)] == 1)

    # Hamming Weight Constraint
    model.Add(sum(e_vars[(j, 1)] for j in range(1, n + 1)) >= w)

    # Resolution 
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


def process_file(file_path, solve_function, forward=False, backward=True):
    """
    Process a single input file and solve the syndrome decoding problem.

    Args:
        file_path (str): Path to the input file containing problem parameters.
        solve_function (function): Function to solve the problem.
        forward (bool): Use forward filtering.
        backward (bool): Use backward filtering.
    
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
        backward=backward
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
    # Method dictionary configuration
    methods = {
        'CNF1': build_and_solve_CP1,
        'CNF2': build_and_solve_CP2,
        'CNF3': build_and_solve_CP3,
        'CNF4': build_and_solve_CP4
    }

    # Configuring the parser argument
    parser = argparse.ArgumentParser(description="CPSAT solver for syndrome decoding.")
    parser.add_argument('-m', '--method', choices=methods.keys(), required=True, 
                        help="Resolution method to use (CNF1, CNF2, CNF3 or CNF4)")
    
    parser.add_argument(
        "--forward",
        action="store_true",
        help="Use forward implication only"
    )

    parser.add_argument(
        "--backward",
        action="store_true",
        help="Use backward implication"
    )
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', '--file', help='Path to an input file')
    group.add_argument('-d', '--dir', help='Path to a directory to process')
    args = parser.parse_args()

    forward = args.forward
    backward = args.backward

    # Determine suffix for CNF4
    suffix = ""
    if args.method == "CNF4":
        if forward and backward:
            suffix = "_equiv"
        elif forward:
            suffix = "_forward"
        elif backward:
            suffix = "_backward"
        else:
            suffix = "_none"

    # Selecting the resolution function
    solve_function = methods[args.method]

    if args.file:
        # Processing a single file
        file, status, res_time, sol = process_file(args.file, solve_function, forward=forward, backward=backward)
        directory = os.path.dirname(args.file)
        csv_filepath = os.path.join(directory, f"CPSAT_{args.method}{suffix}_{os.path.splitext(file)[0]}.csv")
        
        with open(csv_filepath, mode='w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(["File", "Result", "Time (s)", "Solution"])
            csv_writer.writerow([file, status, res_time, sol])
            csvfile.flush()
    else:
        # Processing a complete directory
        csv_filepath = os.path.join(args.dir, f"CPSAT_{args.method}{suffix}.csv")
        
        # Open the CSV file in append mode (in case you need to restart) or write mode.
        with open(csv_filepath, mode='w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(["File", "Result", "Time (s)", "Solution"])
            csvfile.flush()

            # List and filter instance files
            all_files = os.listdir(args.dir)
            entries = [
                f for f in all_files 
                if os.path.isfile(os.path.join(args.dir, f)) 
                and f.startswith("LargeWeight_") 
                and not f.endswith(".csv")
            ]
            
            # Sort
            entries.sort(key=extract_LW_n)

            if not entries:
                print(f"No 'LargeWeight_' instance file found in: {args.dir}")
                return

            print(f"Starting to process {len(entries)} files in ascending order of n.")

            # Processing loop
            for entry in entries:
                path = os.path.join(args.dir, entry)
                n_val = extract_LW_n(entry)
                
                print(f"\n[n={n_val}] Treatment of: {entry} ...")
                
                # Running the selected CP-SAT solver
                file_name, status, res_time, sol = process_file(path, solve_function, forward=forward, backward=backward)
                
                # Immediate writing to the CSV (safety in case of crash)
                csv_writer.writerow([file_name, status, res_time, sol])
                csvfile.flush()
                
                print(f"Completed: {status} in {res_time}s")

    print(f"\nAll files have been processed. Results: {csv_filepath}")


if __name__ == "__main__":
    main()