# SK1

# Repository for Molecular Dynamics and DEER Data

# "Structural dynamics of sphingosine kinase 1 regulation and inhibition" by Baharak Abd Emami, Ahmed Shubbar, Hope Woods, Mahmoud Moradi, and Reza Dastvan 

This repository contains double electron–electron resonance (DEER) spectroscopy data and molecular dynamics (MD) simulations data including input files, simulation outputs, and analysis scripts.  
Each system is organized with a consistent sub-structure for reproducibility.

---

## Directory Layout

Top-level system directories:

- `apo/` – Apo system
- `apo_S225/` – Apo system with Ser225 phosphorylation
- `sqs_atp/` – ATP-bound system
- `sqs_atp_S225/` – ATP-bound system with Ser225 phosphorylation

Each system directory contains three subfolders:
system_name/
│
├── input/      # Input files (.inp) used to start MD runs
├── data/       # Raw outputs (.dcd, .xst, .log, .out, restart files, etc.)
└── analysis/   # Post-processing scripts and processed results


## Running Simulations

1. Navigate to the `input/` folder of the system you want to run:

   ```bash
   cd apo/input

2.	Use the appropriate MD engine to launch the run with the .inp file:
	•	With NAMD:
namd production.inp > production.log

3.	Output trajectories (final coordinate), restart files, and logs will be saved under the corresponding input/



## Running Analysis Scripts

All post-processing is done inside the analysis/ subfolder of each system.

Shell Scripts (.sh)
	•	Used to automate multiple analyses or batch runs.
	•	./{script}.sh

Python Scripts (.py)
	•	Run with Python 3: python3 {script}.py

Tcl Scripts (.tcl)
	•	Run inside VMD in text mode:
 vmd -dispdev text -e {script}.tcl

	•	Inputs: trajectory files (final coordinate) and topology files (.psf / .pdb) from input/
	•	Outputs: text files, CSVs, and plots written to data/





 
