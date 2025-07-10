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

N_THREADS = 3  # Number of threads to use for MCCE calculations, set to 1 for now
AMINO_ACID_PDBS = ["ala.pdb", "arg.pdb", "asp.pdb", "cys.pdb", "glu.pdb", "his.pdb", "ile.pdb", "leu.pdb", "lys.pdb", "met.pdb", "phe.pdb", "pro.pdb", "ser.pdb", "thr.pdb", "trp.pdb", "tyr.pdb", "val.pdb"]
PROTEIN_PDBS = ["small_protein.pdb", "medium_protein.pdb", "large_protein.pdb"]

# Make sure this is consistent with the one in embedding_score.py
ATOM_RADII = {
    ' H': 1.2,
    ' C': 1.7,
    ' N': 1.55,
    ' O': 1.52,
    ' S': 1.8,
    ' P': 1.8,
    # Add more atom types and their radii as needed
}
ATOM_RADIUS_UNKNOWN = 1.5  # Default radius for unknown atom types


def parse_arguments():
    helpmsg = "Setup reaction field energy training data."
    parser = argparse.ArgumentParser(description=helpmsg, formatter_class=argparse.RawTextHelpFormatter)
    return parser.parse_args()

def get_local_density_scores(pdb_file):
    """
    Get local density scores for the atoms in the PDB file.
    This function runs the local_density.py script to calculate the local density scores.
    The output is saved to a file with the same name as the PDB file but with .density extension..
    This fuction also returns a dictionary with atom IDs as keys and their local density scores as values.
    """
    # Calculate local density scores for the amino acid PDB file
    result = subprocess.run(["local_density.py", pdb_file], capture_output=True, text=True)
    # result = subprocess.run(["local_embedding.py", pdb_file], capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"local_density.py failed with exit code {result.returncode}")
        logging.error(f"stderr: {result.stderr.strip()}")
        return
    
    # Read the local density scores from the output file
    density_data = {}  # Dictionary to store local density scores, keyed by atom ID, value is a list of [Near, Mid, Far] scores
    density_file = pdb_file.replace(".pdb", ".density")
    # density_file = pdb_file.replace(".pdb", ".embedding")  # Use the same naming convention as local_embedding.py
    density_lines = open(density_file, 'r').readlines()  # Read the local density scores
    for line in density_lines:
        atomname = line[12:16]
        resname = line[17:20]
        chainid = line[21]
        resseq = line[22:26]
        atom_id = (atomname, resname, chainid, resseq)
        fields = line[54:].strip().split()
        if len(fields) == 4:  # Expecting four fields: near, mid, far, d2surface
            near = float(fields[0])
            mid = float(fields[1])
            far = float(fields[2])
            d2surface = float(fields[3])
            density_data[atom_id] = [near, mid, far, d2surface]

    return density_data


def setup_residue(pdb_file):
    """
    Set up residue level reaction field energy training data.
    This function reads the amino acid PDB file and prepares the training data.
    we will place a charge on a single atom and go through all side chain atoms and run PB solver on each of them.
    The output is saved to amino_acid.csv.
    """
    density_data = get_local_density_scores(pdb_file)

    rxn_lines = []  # List of lines to write out

    # For amino acid level reaction field energy,
    atom_lines = [line.strip() for line in open(pdb_file).readlines() if line.startswith("ATOM  ") or line.startswith("HETATM")]
    # Assign atom radii to remove 0 radius atoms
    for i, line in enumerate(atom_lines):
        atomname = line[12:16]
        if len(atomname.strip()) == 4 and atomname[0] == "H":
            element = " H"
        else:
            element = atomname[:2]
        radius = ATOM_RADII.get(element, ATOM_RADIUS_UNKNOWN)  # Get the radius from the dictionary or use default
        atom_lines[i] = line[:54] + "%8.3f" % radius + line[62:]  # Set the radius in the PDB line
    
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
        result = subprocess.run(["step3.py", "-s", "delphi", "-p", f"{N_THREADS}"], capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f"step3.py failed with exit code {result.returncode}")
            logging.error(f"stderr: {result.stderr.strip()}")

        # Read the output reaction field energies from energies/*.raw
        atom_id = (sidechain_atoms[i_chargedatom][12:16], sidechain_atoms[i_chargedatom][17:20], sidechain_atoms[i_chargedatom][21], sidechain_atoms[i_chargedatom][22:26])
        density = density_data.get(atom_id, [0.0, 0.0, 0.0, 0.0])  # Get the local density scores for this atom
        raw_lines = open(raw_file, 'r').readlines()
        pbrxn = None
        for line in raw_lines:
            if line.startswith("[RXN"):
                fields = line.strip().split()
                if len(fields) >= 2:
                    pbrxn = float(fields[-1])
                    break

        if pbrxn is not None:
            out_line = f"{density[0]:.3f}, {density[1]:.3f}, {density[2]:.3f}, {density[3]:.3f}, {pbrxn:.3f}\n"
            rxn_lines.append(out_line)
    # Write the reaction field energy training data to amino_acid.csv
    output_file = f"{pdb_file.replace('.pdb', '')}_rxn.csv"
    with open(output_file, "w") as f:
        f.write("Density_Near,Density_Mid,Density_Far,D2surface,PBRXN\n")
        f.writelines(rxn_lines)


def setup_protein(pdb_file):
    """
    Set up protein level reaction field energy training data.
    This function takes a microstate PDB file and place a charge on a single atom of every side chain atom.
    To sample all atoms, we will place the charge on each side chain atom, by looping from 0 to max_side_chain_atoms.
    The output is saved to a csv file.
    """
    density_data = get_local_density_scores(pdb_file)
    rxn_lines = []  # List of lines to write out
    # Read the PDB file and assign atom radii
    atom_lines = [line.strip() for line in open(pdb_file).readlines() if line.startswith("ATOM  ") or line.startswith("HETATM")]
    for i, line in enumerate(atom_lines):
        atomname = line[12:16]
        if len(atomname.strip()) == 4 and atomname[0] == "H":
            element = " H"
        else:
            element = atomname[:2]
        radius = ATOM_RADII.get(element, ATOM_RADIUS_UNKNOWN)
        atom_lines[i] = line[:54] + "%8.3f" % radius + line[62:]  # Set the radius in the PDB line
    # Set all atoms to 0.0 charge
    for i, line in enumerate(atom_lines):
        atom_lines[i] = line[:62] + "%12.3f" % (0.0) + line[74:]

    # We have many residues in the protein, so we will loop through each residue and set the charge on each side chain atom
    all_residues = []  # List to store all residues
    residue_lines = []  # List to store lines for each residue
    for line in atom_lines:
        resid = line[17:20] + line[21] + line[22:26]  # Residue ID
        # Check if the resid is the same as the last one, if so, append the line to the residue_lines
        last_resid = residue_lines[-1][17:20] + residue_lines[-1][21] + residue_lines[-1][22:26] if residue_lines else None
        if resid == last_resid or last_resid is None: # If it's the same residue, append the line
            residue_lines.append(line)
        else:  # If it's a new residue, process the last residue and start a new one
            all_residues.append(residue_lines)  # Add the last residue to the list
            residue_lines = [line]

    # Process the last residue if it exists
    if residue_lines:
        all_residues.append(residue_lines)

    # Now we have all_residues populated, we can process each residue
    i_charge = 0  # Index for the charged atom in the residue side chain
    assigned_charges = 1  # initially set to 1 to enter the loop
    while assigned_charges > 0:  # Continue until no charges are assigned
        assigned_charges = 0  # Count of assigned charges
        map_conformer_to_atom = {}  # Dictionary to map conformer name to atom ID, one conformer has one charged atom
        new_pdb_lines = []  # List to store new PDB lines
        for residue in all_residues:
            # group atoms into backbone and sidechain in the residue
            backbone_atoms = []
            sidechain_atoms = []
            for line in residue:
                if line[80:82] == "BK":
                    backbone_atoms.append(line)
                else:
                    sidechain_atoms.append(line)
            # If there are no side chain atoms, skip this residue
            if not sidechain_atoms:
                continue
            # Get the conformer name from the first atom of the side chain, this is used to identify the conformer in energies folder
            sidechain_atom = sidechain_atoms[0]
            conf_name = sidechain_atom[17:20] + sidechain_atom[80:82] + sidechain_atom[21:30]  # Conformer name
            # Assign the charge to the i-th side chain atom
            new_lines = backbone_atoms.copy() + sidechain_atoms.copy()  # Create a new list to modify and make sure backbone atoms are on top
            if i_charge < len(sidechain_atoms):
                new_lines[len(backbone_atoms) + i_charge] = sidechain_atoms[i_charge][:62] + "%12.3f" % (1.0) + sidechain_atoms[i_charge][74:]  # Set charge to 1.0 on the i-th side chain atom
                # We need know which atom is placed with charge on which conformer, so we will make a dictionary to map conformer name to atom ID
                charged_atom = new_lines[len(backbone_atoms) + i_charge]  # Get the charged atom line
                atom_id = (charged_atom[12:16], charged_atom[17:20], charged_atom[21], charged_atom[22:26])  # Atom ID
                map_conformer_to_atom[conf_name] = atom_id  # Map conformer name to atom ID

                assigned_charges += 1  # Increment the count of assigned charges
            new_pdb_lines.extend(new_lines)  # Add the new lines to the new PDB lines
        # Write the new PDB file to step2_out.pdb
        logging.info(f"Running MCCE step3 with charge on atom {i_charge + 1}, charged assigned in this cycle is {assigned_charges}...")
        if assigned_charges > 0:
            with open("step2_out.pdb", "w") as f:
                f.write("\n".join(new_pdb_lines) + "\n")
            # Run MCCE step3 to calculate the reaction field energy
            result = subprocess.run(["step3.py", "-s", "delphi", "-p", f"{N_THREADS}"], capture_output=True, text=True)
            if result.returncode != 0:
                logging.error(f"step3.py failed with exit code {result.returncode}")
                logging.error(f"stderr: {result.stderr.strip()}")
                return
        # Read the output reaction field energies from energies/*.raw. Valid conformers are in the map_conformer_to_atom dictionary
        for conformer, atom_id in map_conformer_to_atom.items():
            # Read the raw file for this conformer
            raw_file = f"energies/{conformer}.raw"
            if not os.path.isfile(raw_file):
                logging.error(f"Required raw file is missing: {raw_file}")
                continue
            raw_lines = open(raw_file, 'r').readlines()
            pbrxn = None
            for line in raw_lines:
                if line.startswith("[RXN"):
                    fields = line.strip().split()
                    if len(fields) >= 2:
                        pbrxn = float(fields[-1])
                        break
            if pbrxn is not None:
                # Get the local density scores for the charged atom
                density = density_data.get(atom_id, [0.0, 0.0, 0.0, 0.0])
                out_line = f"{density[0]:.3f}, {density[1]:.3f}, {density[2]:.3f}, {density[3]:.3f}, {pbrxn:.3f}\n"
                rxn_lines.append(out_line)  # Append the line to the reaction field energy training data
        i_charge += 1  # Increment the index for the next charged atom

    # Write the reaction field energy training data to a csv file
    output_file = f"{pdb_file.replace('.pdb', '')}_rxn.csv"
    with open(output_file, "w") as f:
        f.write("Density_Near,Density_Mid,Density_Far,D2surface,PBRXN\n")
        f.writelines(rxn_lines)


if __name__ == "__main__":
    args = parse_arguments()
    logging.basicConfig(level=logging.INFO)
    logging.info("Setting up reaction field energy training data...")
    # Check if the required PDB files exist
    missing_files = []
    for pdb_file in AMINO_ACID_PDBS+PROTEIN_PDBS:
        if not os.path.isfile(pdb_file):
            missing_files.append(pdb_file)
    if missing_files:
        logging.error(f"Required PDB files are missing:\n{', '.join(missing_files)}")
        exit(1)

    # Residue level reaction field energy
    logging.info("Processing residue level reaction field energy...")
    for pdb_file in AMINO_ACID_PDBS:
        logging.info(f"Setting up reaction field energy for {pdb_file}...")
        setup_residue(pdb_file)
    logging.info("Residue level reaction field energy training data setup complete.")

    # Combine amino residue csv files into a single file
    logging.info("Combining amino acid level reaction field energy training data...")
    amino_acid_files = [f"{pdb_file.replace('.pdb', '')}_rxn.csv" for pdb_file in AMINO_ACID_PDBS]
    combined_amino_acid_file = "amino_acid_rxn.csv"
    with open(combined_amino_acid_file, "w") as outfile:
        outfile.write("Density_Near,Density_Mid,Density_Far,D2surface,PBRXN\n")  # Write header
        for amino_file in amino_acid_files:
            with open(amino_file, "r") as infile:
                next(infile)  # Skip header
                outfile.writelines(infile.readlines())  # Write the rest of the lines
    logging.info(f"Combined amino acid level reaction field energy training data saved to {combined_amino_acid_file}.")

    # Protein level reaction field energy
    logging.info("Processing protein level reaction field energy...")
    for pdb_file in PROTEIN_PDBS:
        logging.info(f"Setting up reaction field energy for {pdb_file}...")
        # Call the setup_protein function with the PDB file and output file
        setup_protein(pdb_file)
    logging.info("Protein level reaction field energy training data setup complete.")
