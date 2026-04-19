# SAT-Based Large Weight Syndrome Decoding (LWSD)

This repository contains source code and instance generators for modelling the Large Weight Syndrome Decoding problem for $Z \ge 2$ using Satisfiability (SAT) and Constraint Programming (CP). 
---

## Environment Setup

It is recommended to use a Python virtual environment (Python 3.8+) to manage dependencies and avoid conflicts.

```bash
# Create the virtual environment
python3 -m venv venv

# Activate the environment
source venv/bin/activate

# Install required dependencies
pip install -r requirements.txt
```
## Large Weight Syndrome Decoding (LW3SD)

1. Generate original instance

```bash
# Generate instances for a specific size 'n' and seed 's'
python3 largeweight_generate.py ${n} ${s}
```

2. Generate CNF models

Transform the generated instances into SAT models (CNF) using different encoding variants and filtering techniques.

```bash
python3 models.py <variant> <instance_file> <card_enc> <pb_enc> <forward:0|1> <equiv:0|1>
```

3. Build and solve with CP-SAT
```bash
python3 LW3SD_CPSAT.py -m <variant> -f <instance_file>
```
