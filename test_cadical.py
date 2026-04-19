import os
import subprocess
import argparse
import csv
import re
import time

def extract_n(filename):
    match = re.search(r"LargeWeight_(\d+)_\d+\.cnf", filename)
    return int(match.group(1)) if match else float('inf')

def test_instance(cnf_filepath, cadical_dir, timeout=None):
    """
    Teste une instance CNF avec cadical et renvoie le résultat, le temps d'exécution et la solution.
    Si timeout est None, le script attend indéfiniment.
    """
    try:
        # Chemin vers l'exécutable cadical
        cadical_path = os.path.join(cadical_dir, "cadical")

        # Extraire seed et variante
        parts = cnf_filepath.strip("/").split(os.sep)
        seed = next(p for p in parts if p.startswith("seed_"))
        variant = next(p for p in parts if p.startswith("CNF"))

        # Chemin complet des logs : logs/seed/variant
        logs_dir = os.path.join("logs", seed, variant)
        os.makedirs(logs_dir, exist_ok=True)  # création récursive si nécessaire

        # Nom du fichier de log
        filename = os.path.basename(cnf_filepath)
        log_filename = os.path.splitext(filename)[0] + "_cadical.log"
        log_path = os.path.join(logs_dir, log_filename)

        with open(log_path, "w") as log_file:
            # Mesurer le temps de résolution
            start_time = time.time()

            # Préparer les arguments pour subprocess.run
            run_args = {
                "stdout": log_file,
                "stderr": subprocess.PIPE,
            }
            if timeout is not None:
                run_args["timeout"] = timeout

            # Lancer cadical
            subprocess.run(
                [cadical_path, cnf_filepath],
                **run_args
            )

            end_time = time.time()
            exec_time = end_time - start_time  # Temps écoulé en secondes

        # Lire la sortie depuis le fichier log
        with open(log_path, "r") as f:
            output = f.read()

        # Vérifier si l'instance est SAT ou UNSAT
        if "s UNSATISFIABLE" in output:
            return "unsat", exec_time, None
        elif "s SATISFIABLE" in output:
            result_status = "sat"
        else:
            return "ERREUR", exec_time, None

        # Extraire la solution (toutes les lignes commençant par 'v')
        solution = " ".join(re.findall(r"v ([\d\s-]+)", output)).replace("\n", " ")

        return result_status, exec_time, solution.strip()

    except subprocess.TimeoutExpired:
        return "timeout", None, None

    except Exception as e:
        print(f"Erreur lors de l'exécution de cadical : {e}")
        return "error", None, None
    


def process_single_file(cnf_filepath, cadical_dir, timeout):
    """
    Teste une seule instance CNF et écrit les résultats dans un fichier CSV.
    """
    if not os.path.exists(cnf_filepath):
        print(f"Erreur : Le fichier {cnf_filepath} n'existe pas.")
        return

    directory = os.path.dirname(cnf_filepath)
    filename = os.path.basename(cnf_filepath)
    csv_filepath = os.path.join(directory, f"{os.path.splitext(filename)[0]}.csv")

    print(f"Test de {filename}...")
    result, exec_time, solution = test_instance(cnf_filepath, cadical_dir, timeout)

    # Écriture dans le fichier CSV
    with open(csv_filepath, mode='w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Fichier", "Résultat", "Temps (s)", "Solution"])
        csv_writer.writerow([filename, result, exec_time, solution])
        csvfile.flush()

    print(f"Résultats écrits dans {csv_filepath}")



def process_directory(input_dir, cadical_dir, timeout):
    """
    Parcourt le dossier d'instances CNF, trie les fichiers par n, teste chaque instance et écrit les résultats dans un CSV.
    """
    if not os.path.exists(input_dir):
        print(f"Erreur : Le dossier {input_dir} n'existe pas.")
        return

    csv_filepath = os.path.join(input_dir, "LWZSD_cadical.csv")

    # Extraire les fichiers .cnf et trier par la valeur de n dans le nom
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
            print(f"Test de {filename}...")
            result, exec_time, solution = test_instance(cnf_filepath, cadical_dir, timeout)

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
                print(f"Erreur avec l'instance {filename}.")

    print(f"Résultats écrits dans {csv_filepath}")

def main():
    # Parser des arguments en ligne de commande
    parser = argparse.ArgumentParser(description="Test des instances CNF avec cadical.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-d', '--directory', type=str, help="Répertoire contenant les instances CNF à tester.")
    group.add_argument('-f', '--file', type=str, help="Fichier CNF à tester.")

    parser.add_argument('--cadical', type=str, required=True, help="Répertoire de cadical.")
    parser.add_argument('-t', '--timeout', type=int, default=None,
                        help="Timeout en secondes. Si non précisé, le temps d'exécution est illimité.")

    args = parser.parse_args()

    if args.directory:
        process_directory(args.directory, args.cadical, args.timeout)
    elif args.file:
        process_single_file(args.file, args.cadical, args.timeout)

if __name__ == "__main__":
    main()
