from pysat.card import CardEnc
from pysat.pb import PBEnc
from pysat.formula import CNF
from utils import *

# one-hot encoding
def build_CNF1(n, w, k, H_transpose, s_transpose, cc_encoding, pb_encoding, Z=3):

    m = n - k  # Number of equations
    cnf = CNF()

    top_id = 0  # Highest variable index so far
    # Assign unique variable ids to error variables e_{j,c}
    e_var = {}
    for j in range(1, n+1):
        for c in range(1, Z):
            top_id += 1
            e_var[(j, c)] = top_id

    # -e_{j,v} V e_{j,v-1}
    for j in range(1, n+1):
        for c in range(2, Z):
            cnf.append([
                -e_var[(j, c)],
                e_var[(j, c-1)]
            ])

    # Build sets V and K
    V , K = build_var_sets(H_transpose, s_transpose, n, k, w, Z)
    
    # Introduce auxiliary variables x_{i,v} 
    x_vars_dict = {}  # Mapping (i, v) to the variable x_{i,v}
    for i in range(m):
        x_vars = []
        for v in K[i]:
            top_id += 1  # Allocate a new variable
            x_vars_dict[(i, v)] = top_id
            x_vars.append(top_id)

        # For each equation, enforce that exactly one auxiliary variable x_{i,v} is set to true 
        if x_vars:
            res = CardEnc.equals(lits=x_vars, bound=1, top_id=top_id, encoding=cc_encoding)
            cnf.extend(res.clauses)
            top_id = res.nv # Update top_id
        
    # Encode the pseudo-Boolean equality constraints for each equation E_i 
    for i in range(m):
        lits = []
        weights = []

        # contribution from e_{j,c}
        for j in V[i]:
            if j <= m:
                H_ij = 1
            else:
                col = j - m - 1
                H_ij = int(H_transpose[i][col])

            for c in range(1, Z):
                coeff = H_ij        
                lits.append(e_var[(j, c)])
                weights.append(coeff)

        # contribution from x_{i,v}
        rhs = sum(int(v) for v in K[i])
        for v in K[i]:
            vv = int(v)
            xv = x_vars_dict[(i, v)]
            lits.append(-xv)  
            weights.append(vv)

        res = PBEnc.equals(lits=lits, weights=weights, bound=rhs, top_id=top_id, encoding=pb_encoding)
        cnf.extend(res.clauses)
        top_id = res.nv

    
    # Encode the constraint on the total Hamming weight of e 
    nonzero_lits = [e_var[(j,1)] for j in range(1, n+1)]
    if nonzero_lits:
        res_w = CardEnc.atleast(lits=nonzero_lits, bound=w, top_id=top_id, encoding=cc_encoding)
        cnf.extend(res_w.clauses)

    return cnf


# unary encoding
def build_CNF2(n, w, k, H_transpose, s_transpose, cc_encoding, pb_encoding, Z=3):
    m = n - k  # Number of equations
    cnf = CNF()

    top_id = 0  # Highest variable index so far
    # Assign unique variable ids to error variables e_{j,c}
    e_var = {}
    for j in range(1, n+1):
        for c in range(1, Z):
            top_id += 1
            e_var[(j, c)] = top_id

    # -e_{j,v} V e_{j,v-1}
    for j in range(1, n+1):
        for c in range(2, Z):
            cnf.append([
                -e_var[(j, c)],
                e_var[(j, c-1)]
            ])

    # Build sets V and K
    V , K = build_var_sets(H_transpose, s_transpose, n, k, w, Z)

    # Introduce auxiliary variables x_{i,v} with chaining clauses
    x_vars_dict = {}
    for i in range(m):
        for v in K[i]:
            top_id += 1
            x_vars_dict[(i, v)] = top_id

            # Clause: ¬x_{i,v} ∨ x_{i,v-Z}
            prev_var = x_vars_dict.get((i, v - Z))
            if prev_var is not None:
                cnf.append([-top_id, prev_var])

    # Encode the pseudo-Boolean equality constraints for each equation E_i (3a dans le rapport)
    for i in range(m):
        lits = []
        weights = []

        # contribution from e_{j,c}
        for j in V[i]:
            if j <= m:  
                # identity side, H_ij = 1 if i == j, else 0
                H_ij = 1 
            else:
                # parity-check side, H_ij from H^T
                col = j - m - 1
                H_ij = int(H_transpose[i][col])

            for c in range(1, Z):
                coeff = H_ij        
                lits.append(e_var[(j, c)])
                weights.append(coeff)

        # contribution from x_{i,v} with v >= Z
        v_min = min(K[i]) 
        rhs = max(K[i])
        for v in K[i]:
            vv = int(v)
            if vv == v_min:
                continue

            xv = x_vars_dict[(i, v)]
            lits.append(-xv)
            weights.append(Z)

        # Encode the PB equality
        res = PBEnc.equals(lits=lits, weights=weights, bound=rhs, top_id=top_id, encoding=pb_encoding)
        cnf.extend(res.clauses)
        top_id = res.nv

    
    # Encode the constraint on the total Hamming weight of e 
    nonzero_lits = [e_var[(j,1)] for j in range(1, n+1)]
    if nonzero_lits:
        res_w = CardEnc.atleast(lits=nonzero_lits, bound=w, top_id=top_id, encoding=cc_encoding)
        cnf.extend(res_w.clauses)
    
    return cnf


# Filtrage exhaustif 
def build_CNF3(n, w, k, H_transpose, s_transpose, cc_encoding, pb_encoding, Z=3):
    m = n - k  # Number of equations
    cnf = CNF()

    top_id = 0  # Highest variable index so far
    # Assign unique variable ids to error variables e_{j,c}
    e_var = {}
    for j in range(1, n+1):
        for c in range(1, Z):
            top_id += 1
            e_var[(j, c)] = top_id

    # -e_{j,v} V e_{j,v-1}
    for j in range(1, n+1):
        for c in range(2, Z):
            cnf.append([
                -e_var[(j, c)],
                e_var[(j, c-1)]
            ])

    # Build sets V and K
    V , K = build_var_sets(H_transpose, s_transpose, n, k, w, Z)

    # Introduce auxiliary variables x_{i,l} for binary representation
    x_vars_dict = {}  
    for i in range(m):
        v_max = max(K[i]) # maximum value in K_i
        q_max = (v_max - int(s_transpose[i])) // Z # max quotient
        L_i = max(1, q_max.bit_length())  # Number of bits needed
        for l in range(L_i):
            top_id += 1
            x_vars_dict[(i, l)] = top_id


    # Encode the pseudo-Boolean equality constraints for each equation E_i
    for i in range(m):
        lits = []
        weights = []

        # contribution from e_{j,c}
        for j in V[i]:
            if j <= m:
                H_ij = 1
            else:
                col = j - m - 1
                H_ij = int(H_transpose[i][col])
            for c in range(1, Z):
                lits.append(e_var[(j, c)])
                weights.append(H_ij)

        # contribution from x_{i,l}
        v_max = max(K[i])
        q_max = (v_max - int(s_transpose[i])) // Z
        L_i = max(1, q_max.bit_length())
        for l in range(L_i):
            xv = x_vars_dict[(i,l)]
            lits.append(-xv) 
            weights.append(Z * 2**l)

        rhs = int(s_transpose[i]) + Z * (2**L_i - 1) 

        res = PBEnc.equals(lits=lits, weights=weights, bound=rhs, top_id=top_id, encoding=pb_encoding)
        cnf.extend(res.clauses)
        top_id = res.nv

        # Add clauses to exclude invalid v values not in K[i]
        F_i = [v for v in range(0, 2**L_i) if (int(s_transpose[i]) + Z*v) not in K[i]]
        for v in F_i:
            clause = []
            # Représentation binaire LSB vers MSB
            bin_v = list(reversed([int(b) for b in format(v, f'0{L_i}b')]))

            for l, bit in enumerate(bin_v):
                xv = x_vars_dict[(i, l)]
                if bit == 0:
                    clause.append(xv)
                else:
                    clause.append(-xv)

            cnf.append(clause)

    # Encode the constraint on the total Hamming weight of e 
    nonzero_lits = [e_var[(j,1)] for j in range(1, n+1)]
    if nonzero_lits:
        res_w = CardEnc.atleast(lits=nonzero_lits, bound=w, top_id=top_id, encoding=cc_encoding)
        cnf.extend(res_w.clauses)

    return cnf


# Filtrage compact
def build_CNF4(n, w, k, H_transpose, s_transpose, cc_encoding, pb_encoding, forward, equiv, Z=3):
    m = n - k  # Number of equations
    cnf = CNF()

    top_id = 0  # Highest variable index so far
    # Assign unique variable ids to error variables e_{j,c}
    e_var = {}
    for j in range(1, n+1):
        for c in range(1, Z):
            top_id += 1
            e_var[(j, c)] = top_id

    # -e_{j,v} V e_{j,v-1}
    for j in range(1, n+1):
        for c in range(2, Z):
            cnf.append([
                -e_var[(j, c)],
                e_var[(j, c-1)]
            ])

    # Build sets V and K
    V , K = build_var_sets(H_transpose, s_transpose, n, k, w, Z)

    # Introduce auxiliary variables x_{i,l} for binary representation
    x_vars_dict = {}  
    for i in range(m):
        v_max = max(K[i]) # maximum value in K_i
        q_max = (v_max - int(s_transpose[i])) // Z # max quotient
        L_i = max(1, q_max.bit_length())  # Number of bits needed
        for l in range(L_i):
            top_id += 1
            x_vars_dict[(i, l)] = top_id


    # Encode the pseudo-Boolean equality constraints for each equation E_i
    for i in range(m):
        lits = []
        weights = []

        # contribution from e_{j,c}
        for j in V[i]:
            if j <= m:
                H_ij = 1
            else:
                col = j - m - 1
                H_ij = int(H_transpose[i][col])
            for c in range(1, Z):
                lits.append(e_var[(j, c)])
                weights.append(H_ij)

        v_min = min(K[i])
        v_max = max(K[i])
        s_i = int(s_transpose[i])

        q_min = (v_min - s_i) // Z
        q_max = (v_max - s_i) // Z

        L_i = max(1, q_max.bit_length())

        # contribution from x_{i,l}
        for l in range(L_i):
            xv = x_vars_dict[(i,l)]
            lits.append(-xv) 
            weights.append(Z * 2**l)

        rhs = int(s_transpose[i]) + Z * (2**L_i - 1) 

        res = PBEnc.equals(lits=lits, weights=weights, bound=rhs, top_id=top_id, encoding=pb_encoding)
        cnf.extend(res.clauses)
        top_id = res.nv

        bits_max = [(q_max >> l) & 1 for l in range(L_i)] # bits of q_max
        bits_min = [(q_min >> l) & 1 for l in range(L_i)] # bits of q_min

        J_max = [l for l in range(L_i) if bits_max[l] == 1] # positions of bits 1 in q_max
        J_min = [l for l in range(L_i) if bits_min[l] == 0] # positions of bits 0 in q_min
        J_max.sort(reverse=True)
        J_min.sort(reverse=True)

        # Exclude invalid values above q_max
        pmax = {}

        if J_max:
            # Create one pmax variable per block-start position in q_max
            for j in J_max:
                if j > 0 and bits_max[j - 1] == 0:
                    top_id += 1
                    pmax[j] = top_id

            block_pmax = sorted(pmax.keys(), reverse=True)

            # backward
            if not forward or equiv:
                for idx, j in enumerate(block_pmax):
                    p_j = pmax[j]
                    if idx == 0:
                        # p_j => all bits from j to L_i-1 match q_max
                        for kk in range(j, L_i):
                            cnf.append([-p_j, x_vars_dict[(i, kk)]])
                    else:
                        j_prime = block_pmax[idx - 1]
                        p_jp = pmax[j_prime]
                        cnf.append([-p_j, p_jp])
                        for kk in range(j, j_prime):
                            xk = x_vars_dict[(i, kk)]
                            if bits_max[kk] == 1:
                                cnf.append([-p_j, xk])
                            else:
                                cnf.append([-p_j, -xk])

            # forward
            if forward:
                for idx, j in enumerate(block_pmax):
                    p_j = pmax[j]
                    if idx == 0:
                        # (all bits from j to L_i-1 match q_max) => p_j
                        clause = [p_j]
                        for kk in range(j, L_i):
                            xk = x_vars_dict[(i, kk)]
                            if bits_max[kk] == 1:
                                clause.append(-xk)   # differs if x_k = 0
                            else:
                                clause.append(xk)    # differs if x_k = 1
                        cnf.append(clause)
                    else:
                        j_prime = block_pmax[idx - 1]
                        p_jp = pmax[j_prime]
                        # (p_{j'} AND bits in [j, j'-1] match q_max) => p_j
                        clause = [p_j, -p_jp]
                        for kk in range(j, j_prime):
                            xk = x_vars_dict[(i, kk)]
                            if bits_max[kk] == 1:
                                clause.append(-xk)
                            else:
                                clause.append(xk)
                        cnf.append(clause)

            # Forbidding clauses (same structure, now triggered by p_j=True)
            for idx, j in enumerate(block_pmax):
                p_j = pmax[j]
                if idx < len(block_pmax) - 1:
                    j2 = block_pmax[idx + 1]
                    for kk in range(j2 + 1, j):
                        cnf.append([-p_j, -x_vars_dict[(i, kk)]])
                else:
                    last = 0
                    if min(J_max) == 0:
                        for jj in sorted(J_max):
                            if bits_max[jj] == 0:
                                break
                            last = jj
                    for kk in range(last, j):
                        cnf.append([-p_j, -x_vars_dict[(i, kk)]])

        # Exclude invalid values below q_min 
        pmin = {}
        if J_min:
            j = J_min[0]  # Most significant 0-bit of q_min
            if j != L_i - 1:
                for kk in range(j+1, L_i):
                    cnf.append([x_vars_dict[(i, kk)]])

            # Create pmin variables
            for j in J_min:
                if j > 0 and bits_min[j - 1] == 1:
                    top_id += 1
                    pmin[j] = top_id

            block_pmin = sorted(pmin.keys(), reverse=True)

            # backward
            if not forward or equiv:
                for idx, j in enumerate(block_pmin):
                    p_j = pmin[j]
                    if idx == 0:
                        for kk in range(j, L_i):
                            xk = x_vars_dict[(i, kk)]
                            if bits_min[kk] == 0:
                                cnf.append([-p_j, -xk])
                            else:
                                cnf.append([-p_j, xk])
                    else:
                        j_prime = block_pmin[idx - 1]
                        p_jp = pmin[j_prime]
                        cnf.append([-p_j, p_jp])
                        for kk in range(j, j_prime):
                            xk = x_vars_dict[(i, kk)]
                            if bits_min[kk] == 0:
                                cnf.append([-p_j, -xk])
                            else:
                                cnf.append([-p_j, xk])

            # forward
            if forward:
                for idx, j in enumerate(block_pmin):
                    p_j = pmin[j]
                    if idx == 0:
                        # (bits from j to L_i-1 match q_min) => p_j
                        clause = [p_j]
                        for kk in range(j, L_i):
                            xk = x_vars_dict[(i, kk)]
                            if bits_min[kk] == 0:
                                clause.append(xk)    # differs if x_k = 1
                            else:
                                clause.append(-xk)   # differs if x_k = 0
                        cnf.append(clause)
                    else:
                        j_prime = block_pmin[idx - 1]
                        p_jp = pmin[j_prime]
                        # (p_{j'} AND bits in [j, j'-1] match q_min) => p_j
                        clause = [p_j, -p_jp]
                        for kk in range(j, j_prime):
                            xk = x_vars_dict[(i, kk)]
                            if bits_min[kk] == 0:
                                clause.append(xk)
                            else:
                                clause.append(-xk)
                        cnf.append(clause)

            for idx, j in enumerate(block_pmin):
                p_j = pmin[j]
                if idx < len(block_pmin) - 1:
                    j2 = block_pmin[idx + 1]
                    for kk in range(j2 + 1, j):
                        cnf.append([-p_j, x_vars_dict[(i, kk)]])
                else:
                    last = 0
                    if min(J_min) == 0:
                        for jj in sorted(J_min):
                            if bits_min[jj] == 1:
                                break
                            last = jj
                    for kk in range(last, j):
                        cnf.append([-p_j, x_vars_dict[(i, kk)]])

        else:
            # All bits of q_min are 1: force all x bits to 1
            for kk, b in enumerate(bits_min):
                cnf.append([x_vars_dict[(i, kk)]])

    # Encode the constraint on the total Hamming weight of e 
    nonzero_lits = [e_var[(j,1)] for j in range(1, n+1)]
    if nonzero_lits:
        res_w = CardEnc.atleast(lits=nonzero_lits, bound=w, top_id=top_id, encoding=cc_encoding)
        cnf.extend(res_w.clauses)

    return cnf

# 0 Filtrage 
def build_CNF5(n, w, k, H_transpose, s_transpose, cc_encoding, pb_encoding, Z=3):
    m = n - k  # Number of equations
    cnf = CNF()

    top_id = 0  # Highest variable index so far
    # Assign unique variable ids to error variables e_{j,c}
    e_var = {}
    for j in range(1, n+1):
        for c in range(1, Z):
            top_id += 1
            e_var[(j, c)] = top_id

    # -e_{j,v} V e_{j,v-1}
    for j in range(1, n+1):
        for c in range(2, Z):
            cnf.append([
                -e_var[(j, c)],
                e_var[(j, c-1)]
            ])

    # Build sets V and K
    V , K = build_var_sets(H_transpose, s_transpose, n, k, w, Z)

    # Introduce auxiliary variables x_{i,l} for binary representation
    x_vars_dict = {}  
    for i in range(m):
        v_max = max(K[i]) # maximum value in K_i
        q_max = (v_max - int(s_transpose[i])) // Z # max quotient
        L_i = max(1, q_max.bit_length())  # Number of bits needed
        for l in range(L_i):
            top_id += 1
            x_vars_dict[(i, l)] = top_id


    # Encode the pseudo-Boolean equality constraints for each equation E_i
    for i in range(m):
        lits = []
        weights = []

        # contribution from e_{j,c}
        for j in V[i]:
            if j <= m:
                H_ij = 1
            else:
                col = j - m - 1
                H_ij = int(H_transpose[i][col])
            for c in range(1, Z):
                lits.append(e_var[(j, c)])
                weights.append(H_ij)

        # contribution from x_{i,l}
        v_max = max(K[i])
        q_max = (v_max - int(s_transpose[i])) // Z
        L_i = max(1, q_max.bit_length())
        for l in range(L_i):
            xv = x_vars_dict[(i,l)]
            lits.append(-xv) 
            weights.append(Z * 2**l)

        rhs = int(s_transpose[i]) + Z * (2**L_i - 1) 

        res = PBEnc.equals(lits=lits, weights=weights, bound=rhs, top_id=top_id, encoding=pb_encoding)
        cnf.extend(res.clauses)
        top_id = res.nv

    # Encode the constraint on the total Hamming weight of e 
    nonzero_lits = [e_var[(j,1)] for j in range(1, n+1)]
    if nonzero_lits:
        res_w = CardEnc.atleast(lits=nonzero_lits, bound=w, top_id=top_id, encoding=cc_encoding)
        cnf.extend(res_w.clauses)

    return cnf