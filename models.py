import sys
from pysat.card import CardEnc
from pysat.pb import PBEnc
from pysat.formula import CNF
from utils import *

from LW3SD_CNF import build_CNF1, build_CNF2, build_CNF3, build_CNF4, build_CNF5


def main():
    if len(sys.argv) < 7:
        print("Usage: python3 models.py <variant: CNF1 | CNF2 | CNF3 | CNF4 | CNF5> <chemin_instance> <cc_encoding> <pb_encoding> <forward:0|1> <equiv:0|1>")
        print("Exemple: python3 models.py CNF4 Challenges/LargeWeight/LargeWeight_10_0 3 5 1 0")
        sys.exit(1)

    variant = sys.argv[1]
    input_file = sys.argv[2]
    cc_encoding = int(sys.argv[3])
    pb_encoding = int(sys.argv[4])
    forward = bool(int(sys.argv[5]))
    equiv = bool(int(sys.argv[6]))
    suffix = ""
    if variant == "CNF4":
        if equiv:
            suffix = "_equiv"
        elif forward:
            suffix = "_forward"
    Z = 3

    # Vérification logique
    if equiv and not forward:
        print("Erreur : equiv nécessite forward")
        sys.exit(1)

    # Parse instance
    n, seed, w, k, H_transpose, s_transpose = parse_input_file(input_file)

    if variant == "CNF1":
        cnf = build_CNF1(n, w, k, H_transpose, s_transpose, cc_encoding, pb_encoding, Z)

    elif variant == "CNF2":
        cnf = build_CNF2(n, w, k, H_transpose, s_transpose, cc_encoding, pb_encoding, Z)

    elif variant == "CNF3":
        cnf = build_CNF3(n, w, k, H_transpose, s_transpose, cc_encoding, pb_encoding, Z)

    elif variant == "CNF4":
        cnf = build_CNF4(
            n, w, k,
            H_transpose,
            s_transpose,
            cc_encoding,
            pb_encoding,
            forward,
            equiv,
            Z
        )

    elif variant == "CNF5":
        cnf = build_CNF5(n, w, k, H_transpose, s_transpose, cc_encoding, pb_encoding, Z)

    else:
        print("variant doit être CNF1, CNF2, CNF3, CNF4 ou CNF5")
        sys.exit(1)

    print(f"{variant} construite : {len(cnf.clauses)} clauses, {cnf.nv} variables")

    variant_name = variant + suffix
    write_cnf_to_file(
        input_file,
        cc_encoding,
        pb_encoding,
        cnf,
        seed,
        variant_name
    )


if __name__ == "__main__":
    main()