import os
import subprocess
import argparse
import csv
import re
import time

def extract_n(filename):
    match = re.search(r"LargeWeight_(\d+)_\d+\.cnf", filename)
    return int(match.group(1)) if match else float('inf')

def test_instance(cnf_filepath, cryptominisat_dir, timeout=None):
    """
    Tests a CNF instance with cryptominisat and returns the result, execution time, and solution.
    """
    try:
        # Path to the cryptominisat executable
        cryptominisat_path = os.path.join(cryptominisat_dir, "cryptominisat5")

        # Extract seed and variant
        parts = cnf_filepath.strip("/").split(os.sep)
        seed = next(p for p in parts if p.startswith("seed_"))
        variant = next(p for p in parts if p.startswith("CNF"))

        # Full log path: logs/seed/variant
        logs_dir = os.path.join("logs", seed, variant)
        os.makedirs(logs_dir, exist_ok=True)

        # Log file name
        filename = os.path.basename(cnf_filepath)
        log_filename = os.path.splitext(filename)[0] + "_cryptominisat.log"
        log_path = os.path.join(logs_dir, log_filename)

        with open(log_path, "w") as log_file:
            # Measuring the resolution time
            start_time = time.time()

            run_args = {
                "stdout": log_file,
                "stderr": subprocess.PIPE,
            }
            if timeout is not None:
                run_args["timeout"] = timeout

            # Run cryptominisat
            subprocess.run(
                [cryptominisat_path, cnf_filepath],
                **run_args
            )

            end_time = time.time()
            exec_time = end_time - start_time  # Elapsed time in seconds

        # Read the output from the log file
        with open(log_path, "r") as f:
            output = f.read()

        # Check if the instance is SAT or UNSAT
        if "s UNSATISFIABLE" in output:
            return "unsat", exec_time, None
        elif "s SATISFIABLE" in output:
            result_status = "sat"
        else:
            return "ERREUR", exec_time, None

        # Extract the solution (all lines starting with 'v')
        solution = " ".join(re.findall(r"v ([\d\s-]+)", output)).replace("\n", " ")

        return result_status, exec_time, solution.strip()

    except subprocess.TimeoutExpired:
        return "timeout", None, None

    except Exception as e:
        print(f"Error executing cryptominisat : {e}")
        return "error", None, None
    


def process_single_file(cnf_filepath, cryptominisat_dir, timeout):
    """
    Test a single CNF instance and write the results to a CSV file.
    """
    if not os.path.exists(cnf_filepath):
        print(f"Error: The file {cnf_filepath} does not exist.")
        return

    directory = os.path.dirname(cnf_filepath)
    filename = os.path.basename(cnf_filepath)
    csv_filepath = os.path.join(directory, f"{os.path.splitext(filename)[0]}.csv")

    print(f"Test of {filename}...")
    result, exec_time, solution = test_instance(cnf_filepath, cryptominisat_dir, timeout)

    # Writing to the CSV file
    with open(csv_filepath, mode='w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Fichier", "Résultat", "Temps (s)", "Solution"])
        csv_writer.writerow([filename, result, exec_time, solution])
        csvfile.flush()

    print(f"Results written in {csv_filepath}")



def process_directory(input_dir, cryptominisat_dir, timeout):
    """
    Browses the CNF instance folder, sorts the files by n, tests each instance and writes the results to a CSV.
    """
    if not os.path.exists(input_dir):
        print(f"Error: The folder {input_dir} does not exist.")
        return

    csv_filepath = os.path.join(input_dir, "LWZSD_cryptominisat.csv")

    # Extract the .cnf files and sort by the value of n in the name
    cnf_files = sorted(
        [f for f in os.listdir(input_dir) if f.endswith(".cnf")],
        key=extract_n
    )

    with open(csv_filepath, mode='w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Fichier", "Résultat", "Temps (s)", "Solution"])
        csvfile.flush()

        for filename in cnf_files:
            cnf_filepath = os.path.join(input_dir, filename)
            print(f"Test of {filename}...")
            result, exec_time, solution = test_instance(cnf_filepath, cryptominisat_dir, timeout)

            if solution is not None:
                n = int(re.match(r"LargeWeight_(\d+)_", filename).group(1))

            csv_writer.writerow([os.path.splitext(filename)[0], result, exec_time, solution])
            csvfile.flush()

            if result == "unsat":
                print(f"Instance {filename} est UNSAT.")
            elif result == "sat":
                print(f"Instance {filename} est SAT en {exec_time} secondes.")
            elif result == "timeout":
                print(f"Instance {filename} : TIMEOUT.")
            else:
                print(f"Error with instance {filename}.")

    print(f"Results written in {csv_filepath}")

def main():
    # Parsing command-line arguments
    parser = argparse.ArgumentParser(description="Testing CNF instances with cryptominisat.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-d', '--directory', type=str, help="Directory containing the CNF instances to be tested.")
    group.add_argument('-f', '--file', type=str, help="CNF file to be tested.")

    parser.add_argument('--cryptominisat', type=str, required=True, help="Cryptominisat directory.")
    parser.add_argument('-t', '--timeout', type=int, default=None,
                        help="Timeout in seconds. If not specified, the execution time is unlimited.")

    args = parser.parse_args()

    if args.directory:
        process_directory(args.directory, args.cryptominisat, args.timeout)
    elif args.file:
        process_single_file(args.file, args.cryptominisat, args.timeout)

if __name__ == "__main__":
    main()
