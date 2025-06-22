#!/usr/bin/env python
"""
MCCE4 Tool: Setup Reaction field energy training data
This script prepares the training data for reaction field energy calculations.
The set up includes:
1. Amino acid level reaction field energy training data
2. Small, medium, and large protein level reaction field energy training data

# The required files are:
- amino_acid.pdb: microstate PDB file of amino acid
- small_protein.pdb: microstate PDB file of small protein
- medium_protein.pdb: microstate PDB file of medium protein
- large_protein.pdb: microstate PDB file of large protein

The output files are:
- amino_acid.csv: training data for amino acid level reaction field energy
- small_protein.csv: training data for small protein level reaction field energy
- medium_protein.csv: training data for medium protein level reaction field energy
- large_protein.csv: training data for large protein level reaction field energy
- mixed_protein.csv: training data for mixed amino acid and protein level reaction field energy

# The output files contain the following columns:
DensityAverage_Near, DensityAverage_Mid, DensityAverage_Far, PBRXN
"""

import logging
import argparse
import os
import pandas as pd
import subprocess

TEMPLATE_PDBS = ["amino_acid.pdb", "small_protein.pdb", "medium_protein.pdb", "large_protein.pdb"]

def parse_arguments():
    helpmsg = "Setup reaction field energy training data."
    parser = argparse.ArgumentParser(description=helpmsg, formatter_class=argparse.RawTextHelpFormatter)
    return parser.parse_args()

def setup_residue():
    """
    Set up residue level reaction field energy training data.
    This function reads the amino acid PDB file and prepares the training data.
    we will place a charge on a single atom and go through all side chain atoms and run PB solver on each of them.
    The output is saved to amino_acid.csv.
    """
    rxn_lines = []  # List of lines to write out

    # Calculate local density scores for the amino acid PDB file
    logging.info("Calculating local density scores for amino acid PDB file...")
    result = subprocess.run(["local_density.py", "amino_acid.pdb"], capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"local_density.py failed with exit code {result.returncode}")
        logging.error(f"stderr: {result.stderr.strip()}")
        return
    
    # Read the local density scores from the output file
    density_data = {}  # Dictionary to store local density scores, keyed by atom ID, value is a list of [Near, Mid, Far] scores
    density_file = "amino_acid.density"
    density_lines = open(density_file, 'r').readlines()  # Read the local density scores
    for line in density_lines:
        atomname = line[12:16]
        resname = line[17:20]
        chainid = line[21]
        resseq = line[22:26]
        atom_id = (atomname, resname, chainid, resseq)
        fields = line[54:].strip().split()
        if len(fields) == 3:  # Expecting three fields: near, mid, far
            near = float(fields[0])
            mid = float(fields[1])
            far = float(fields[2])
            density_data[atom_id] = [near, mid, far]

    # For amino acid level reaction field energy,
    atom_lines = [line.strip() for line in open("amino_acid.pdb").readlines() if line.startswith("ATOM  ") or line.startswith("HETATM")]
    # set all atoms to 0.0 charge
    for i, line in enumerate(atom_lines):
        atom_lines[i] = line[:62] + "%12.3f" % (0.0) + line[74:]  # Set charge to 0.0
    # divide atoms into backbone atoms and side chain atoms
    backbone_atoms = []
    sidechain_atoms = []
    for line in atom_lines:
        if line[80:82] == "BK":
            backbone_atoms.append(line)
        else:
            sidechain_atoms.append(line)

    # Infer conformer name from the side chain atom and compose the output file name
    sidechain_atom = sidechain_atoms[0]
    conf_name = sidechain_atom[17:20] + sidechain_atom[80:82] + sidechain_atom[21:30]  
    raw_file = f"energies/{conf_name}.raw"

    # Loop through each side chain atom, and set the charge to 1.0 in step2_out.pdb file
    for i_chargedatom in range(len(sidechain_atoms)):
        new_lines = backbone_atoms.copy() + sidechain_atoms.copy()  # Create a new list to modify and make sure backbone atoms are on top
        # Set the charge of the i-th side chain atom to 1.0
        new_lines[len(backbone_atoms) + i_chargedatom] = sidechain_atoms[i_chargedatom][:62] + "%12.3f" % (1.0) + sidechain_atoms[i_chargedatom][74:]

        # Write the new PDB file to step2_out.pdb
        with open("step2_out.pdb", "w") as f:
            f.write("\n".join(new_lines) + "\n")

        # Run MCCE step3 to calculate the reaction field energy
        logging.info(f"Running MCCE step3 with charge on atom \"{sidechain_atoms[i_chargedatom][12:16]}\"...")
        result = subprocess.run(["step3.py", "-s", "delphi"], capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f"step3.py failed with exit code {result.returncode}")
            logging.error(f"stderr: {result.stderr.strip()}")

        # Read the output reaction field energies from energies/*.raw
        atom_id = (sidechain_atoms[i_chargedatom][12:16], sidechain_atoms[i_chargedatom][17:20], sidechain_atoms[i_chargedatom][21], sidechain_atoms[i_chargedatom][22:26])
        density = density_data.get(atom_id, [0.0, 0.0, 0.0])  # Get the local density scores for this atom
        raw_lines = open(raw_file, 'r').readlines()
        pbrxn = None
        for line in raw_lines:
            if line.startswith("[RXN"):
                fields = line.strip().split()
                if len(fields) >= 2:
                    pbrxn = float(fields[-1])
                    break

        if pbrxn is not None:
            out_line = f"{density[0]:.3f}, {density[1]:.3f}, {density[2]:.3f}, {pbrxn:.3f}\n"
            rxn_lines.append(out_line)
    # Write the reaction field energy training data to amino_acid.csv
    output_file = "amino_acid_rxn.csv"
    with open(output_file, "w") as f:
        f.write("DensityAverage_Near,DensityAverage_Mid,DensityAverage_Far,PBRXN\n")
        f.writelines(rxn_lines)




if __name__ == "__main__":
    args = parse_arguments()
    logging.basicConfig(level=logging.INFO)
    logging.info("Setting up reaction field energy training data...")
    # Check if the required PDB files exist
    for pdb_file in TEMPLATE_PDBS:
        if not os.path.isfile(pdb_file):
            logging.error(f"Required PDB file {pdb_file} does not exist. Please provide the file.")
            exit(1)

    # Residue level reaction field energy
    logging.info("Processing residue level reaction field energy...")
    setup_residue()
    logging.info("Residue level reaction field energy training data setup complete.")


