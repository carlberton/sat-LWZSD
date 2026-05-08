import sys
from pysat.formula import CNF
from utils import *

# Import the specific encoding builders defined in LW3SD_CNF.py
from LW3SD_CNF import (
    build_CNF1, 
    build_CNF2, 
    build_CNF3, 
    build_CNF3_Exhaustive, 
    build_CNF3_Compact
)

def main():
    # Map command-line arguments to the corresponding builder functions
    # CNF3_E and CNF3_C are used as convenient aliases for filtering variants
    variants = {
        "CNF1": build_CNF1,
        "CNF2": build_CNF2,
        "CNF3": build_CNF3,
        "CNF3_E": build_CNF3_Exhaustive,
        "CNF3_C": build_CNF3_Compact
    }

    if len(sys.argv) < 5:
        print("Usage: python3 models.py <variant> <instance_path> <cc_encoding> <pb_encoding>")
        print("Variants: CNF1, CNF2, CNF3, CNF3_E, CNF3_C")
        print("Example: python3 models.py CNF3_C Challenges/LargeWeight/LargeWeight_10_0 3 5")
        sys.exit(1)

    variant_key = sys.argv[1]
    input_file = sys.argv[2]
    cc_encoding = int(sys.argv[3])
    pb_encoding = int(sys.argv[4])
    Z = 3

    # Check if the requested variant is supported
    if variant_key not in variants:
        print(f"Error: Variant '{variant_key}' is unknown.")
        print(f"Supported variants: {', '.join(variants.keys())}")
        sys.exit(1)

    # Parse the input instance file
    try:
        n, seed, w, k, H_transpose, s_transpose = parse_input_file(input_file)
    except Exception as e:
        print(f"Error parsing file {input_file}: {e}")
        sys.exit(1)

    # Select the appropriate builder function
    build_func = variants[variant_key]
    
    # Generate the CNF formula
    cnf = build_func(
        n, w, k, 
        H_transpose, 
        s_transpose, 
        cc_encoding, 
        pb_encoding, 
        Z
    )

    print(f"[{variant_key}] build complete: {len(cnf.clauses)} clauses, {cnf.nv} variables")

    # Save the generated CNF to a file using the variant name
    write_cnf_to_file(
        input_file,
        cc_encoding,
        pb_encoding,
        cnf,
        seed,
        variant_key
    )

if __name__ == "__main__":
    main()