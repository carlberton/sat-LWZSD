#!/usr/bin/env python3

import sys
import random
import os
import math

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def usage():
    eprint("ERROR the script expects 2 integer arguments: 'n' and 'seed'.")
    eprint("\b - 'n' is an integer corresponding to the size of the matrix: H will be of size (n-k) * n.")
    eprint("\b - 'seed' is an integer corresponding to the initial value of the random seed")
    eprint("This script generates an instance of the ternary syndrome decoding problem.")
    eprint("The length is n and the dimension is k = floor(R*n),  where R=1-log_3(2)=0.36907.")
    eprint("This matrix H is given in systematic form. The identity part is omitted.")
    eprint("The instance is stored in 'Challenges/seed_${seed}/LargeWeight/LargeWeight_n_seed'.")

def main(n, seed):
    R = 1 - math.log(2,3)
    w = math.floor(0.99 * n)
    k = math.floor(R*n)
    random.seed(seed)

    # Nouveau chemin avec seed
    prefix = os.path.join(os.getcwd(), f"Challenges/seed_{seed}/LargeWeight/")
    os.makedirs(prefix, exist_ok=True)  # créer les dossiers si nécessaire

    text = ""
    text += "# n\n"
    text += str(n) + "\n"
    text += "# seed\n"
    text += str(seed) + "\n"
    text += "# k\n"
    text += str(k) + "\n"
    text += "# w\n"
    text += str(w) + "\n"
    text += "# H^transpose (each line corresponds to a column of H, the identity part is omitted)\n"
    for i in range(k):
        line = ""
        for j in range(n-k):
            line += str(random.randint(0,2))
        line += "\n"
        text += line
    text += "# s^transpose\n"
    line = ""
    for j in range(n-k):
        line += str(random.randint(0,2))
    line += "\n"
    text += line

    filename = os.path.join(prefix, f"LargeWeight_{n}_{seed}")
    with open(filename, "w") as file:
        file.write(text)

if __name__ == "__main__":
    if len(sys.argv)!=3:
        usage()
        exit(1)
    try:
        n = int(sys.argv[1])
        seed = int(sys.argv[2])
    except:
        usage()
        exit(1)
    main(n, seed)
