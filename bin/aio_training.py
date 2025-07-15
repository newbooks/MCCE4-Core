#!/usr/bin/env python3
"""
All in one training script for fast force field by Machine Learning.
"""

import logging
import argparse
import os
import shutil
import subprocess
import tempfile
import random
import time
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.neural_network import MLPRegressor
import joblib
from matplotlib import pyplot as plt
import seaborn as sns
from mcce.geom import Vector


# Atom radii are for step2_out.pdb
ATOM_RADII = {
    " H": 1.20,  # Hydrogen
    " C": 1.70,  # Carbon
    " N": 1.55,  # Nitrogen
    " O": 1.52,  # Oxygen
    " S": 1.80,  # Sulfur
    " P": 1.80,  # Phosphorus
}
ATOM_RADIUS_UNKNOWN = 1.80  # Default radius for unknown elements

class Atom:
    def __init__(self):
        self.line = ""  # Original line from PDB file
        self.element = ""
        self.radius = 0.0
        self.charge = 0.0


def read_pdb(pdb_file):
    """
    Read a PDB file and extract atom information.
    """
    atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                atom = Atom()
                atom.line = line.strip()
                atomname = line[12:16]
                if len(atomname.strip()) == 4 and atomname[0] == "H":
                    atom.element = " H"
                else:
                    atom.element = atomname[:2]
                atom.radius = ATOM_RADII.get(atom.element, ATOM_RADIUS_UNKNOWN)

                atoms.append(atom)
    return atoms

def assign_charges(atoms):
    """
    Assign charges to atoms of the protein. +1 charge is assigned to a random atom of each side chain conformer as long as the atom.
    The reason a random atom is chosen is that more atoms with different embedding depth will be sampled.
    """
    sidechain_atoms = {}
    for atom in atoms:
        # collect all atoms of the same side chain
        sidechain_id = atom.line[17:30]
        if sidechain_id[-3:] != "000":
            if sidechain_id not in sidechain_atoms:
                sidechain_atoms[sidechain_id] = []
            sidechain_atoms[sidechain_id].append(atom)
    for sidechain_id, sidechain_atoms in sidechain_atoms.items():
        random_atom = random.choice(sidechain_atoms)
        random_atom.charge = 1.0

def write_pdb(atoms, output_file):
    """
    Write the modified atoms to a new PDB file.
    """
    with open(output_file, 'w') as f:
        for atom in atoms:
            line = atom.line[:54] + f"{atom.radius:8.3f}" + f"{atom.charge:12.3f}" + atom.line[74:]
            f.write(line + '\n')

class AtomProperties:
    __slots__ = ("confid", "xyz", "radius", "charge", "density_near", "density_mid", "density_far", "d2surface")

    def __init__(self):
        self.confid = ""
        self.xyz = Vector()
        self.radius = 0.0
        self.charge = 0.0
        self.density_near = 0
        self.density_mid = 0
        self.density_far = 0
        self.d2surface = 0.0 

    def __repr__(self):
        return (f"{self.confid} {self.xyz.x:8.3f} {self.xyz.y:8.3f} {self.xyz.z:8.3f} "
                f"{self.radius:8.3f} {self.charge:8.3f} {self.density_near:6d} {self.density_mid:6d} "
                f"{self.density_far:6d} {self.d2surface:6.3f}")

def load_atoms():
    """
    Load atom properties from step2_out.pdb.
    Only charged atoms are considered.
    """
    pdb_file = "step2_out.pdb"
    if not os.path.isfile(pdb_file):
        logging.error(f"PDB file {pdb_file} not found")
        exit(1)
    atoms = {}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                atom = AtomProperties()
                atomname = line[12:16]
                resname = line[17:20]
                chainid = line[21]
                resseq = line[22:26]
                atom_id = (atomname, resname, chainid, resseq)
                confid = resname + line[80:82] + line[21:30]
                atom.xyz = Vector((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                atom.radius = float(line[54:62])
                atom.charge = float(line[62:74])
                if abs(atom.charge) > 0.001:  # Only consider charged atoms
                    atom.confid = confid
                    atoms[atom_id] = atom
    return atoms

def update_density_score(atoms, fname):
    """
    Update local density scores for atoms.
    Local density is calculated by step2_out.pdb and stored in a file with the same name as the state name.
    The file contains lines like:
    ATOM  1  N   ALA A   1      11.104  12.345  13.456  1.00  0.00           N
    The last columns are the local density scores.
    """
    if not os.path.isfile(fname):
        logging.error(f"Local density file {fname} not found")
        exit(1)
    with open(fname) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                atom_id = (line[12:16], line[17:20], line[21], line[22:26])
                local_density = [x for x in line[54:].strip().split()]
                if atom_id in atoms:
                    atoms[atom_id].density_near = int(local_density[0])
                    atoms[atom_id].density_mid = int(local_density[1])
                    atoms[atom_id].density_far = int(local_density[2])
                    atoms[atom_id].d2surface = float(local_density[3])


def get_electrostatic_energy(atoms):
    """
    Get electrostatic energy from opp files.
    Returns a dictionary mapping (atom_id1, atom_id2) to electrostatic energy.
    """
    # Map confid to atom_id for quick lookup
    conformers = {}
    for atom_id, atom in atoms.items():
        if atom.confid not in conformers:
            conformers[atom.confid] = atom_id
        else:
            logging.warning(f"Confid {atom.confid} already exists, skipping {atom_id}")

    pairwise_ele = {}
    opp_dir = "energies"
    for conf1, atom_id1 in conformers.items():
        fname = os.path.join(opp_dir, f"{conf1}.opp")
        if not os.path.isfile(fname):
            logging.error(f"Opp file {fname} not found")
            exit(1)
        with open(fname) as opp_file:
            for line in opp_file:
                fields = line.strip().split()
                if len(fields) < 6:
                    continue
                conf2 = fields[1]
                atom_id2 = conformers.get(conf2)
                if atom_id2 is not None:
                    try:  # This version averages the energy if both directions are recorded
                        ele = float(fields[5])
                        reverse_key = (atom_id2, atom_id1)
                        if pairwise_ele.get(reverse_key) is None:  # the other direction is not recorded
                            pairwise_ele[(atom_id1, atom_id2)] = ele
                        else:
                            pairwise_ele[reverse_key] = (pairwise_ele[reverse_key] + ele) / 2  # Average the energy if both directions are recorded
                    except ValueError:
                        continue
                    # try: # This version only records the energy in one direction
                    #     ele = float(fields[5])
                    #     pairwise_ele[(atom_id1, atom_id2)] = ele
                    # except ValueError:
                    #     logging.error(f"Invalid energy value in {fname} for {conf1} and {conf2}: {fields[5]}")
                    #     continue
    return pairwise_ele


def fit_ann(X_scaled, y, title):
    """
    Fit an ANN model to the scaled features and target variable.
    This function is a placeholder for the actual ANN fitting logic.
    """
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=int(time.time()) % 1000)

    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    model = MLPRegressor(hidden_layer_sizes=(50, 20), alpha=0.01, learning_rate_init=0.001, learning_rate='adaptive', max_iter=500, random_state=int(time.time()) % 1000, early_stopping=True)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    r2 = r2_score(y_test, y_pred)

    logging.info(f"Results for {title}:")
    logging.info(f"Root Mean Squared Error: {rmse}")
    logging.info(f"R^2 Score: {r2}")

    # Plot the predictions vs actual values
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_test, y=y_pred)
    plt.xlabel("Actual Values")
    plt.ylabel("Predicted Values")
    plt.title(f"Actual vs Predicted - {title}")
    # Plot the diagonal line
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)
    plt.grid(True)
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    plt.text(0.05, 0.95, f"R^2: {r2:.3f}\nRMSE: {rmse:.3f}", transform=plt.gca().transAxes)
    plt.tight_layout()
    plt.savefig(f"{title.replace(' ', '_')}_predictions.png")
    joblib.dump(model, f"{title.replace(' ', '_')}.pkl")
    logging.info(f"Model saved as {title.replace(' ', '_')}.pkl")



if __name__ == "__main__":
    # Set up command line arguments
    helpmsg = "Train a fast force field model from a pdb folder. This script attempts to reproduce PB solver delphi potential using machine learning techniques."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument('--fit_pdbs', metavar="pdbfolder", default=None, help='Fit by pdb files. The folder should contain proteins at different size and shape.')
    parser.add_argument('--fit_csv', metavar="csvfile", default=None, help='Fit by csv file. The file should contain features and target values.')
    parser.add_argument('--debug', default=False, action='store_true', help='If set, create and preserve temporary directory in current directory.')
    args = parser.parse_args()

    # Set up logging
    logging_format = "%(asctime)s %(levelname)s: %(message)s"
    logging_datefmt='%Y-%m-%d %H:%M:%S'
    logging.basicConfig(format=logging_format, datefmt=logging_datefmt, level=logging.INFO)

    # check which fitting method is used
    if args.fit_pdbs is None and args.fit_csv is None:
        logging.error("Please provide either --fit_pdbs or --fit_csv argument.")
        exit(1)
    
    if args.fit_pdbs is not None and args.fit_csv is not None:
        logging.error("Please provide only one of --fit_pdbs or --fit_csv argument.")
        exit(1)
    
    if args.fit_pdbs is not None:
        pdbfolder = args.fit_pdbs
        if not os.path.isdir(pdbfolder):
            logging.error(f"The provided path '{pdbfolder}' is not a valid directory.")
            exit(1)
        if not any(file.endswith('.pdb') for file in os.listdir(pdbfolder)):
            logging.error(f"The directory '{pdbfolder}' does not contain any PDB files.")
            exit(1)
        
        # get the absolute path to pdb files in the folder
        args.pdb_folder = os.path.abspath(pdbfolder)
        pdbs = [os.path.join(args.pdb_folder, file) for file in os.listdir(args.pdb_folder) if file.endswith('.pdb')]

        # save current working directory
        current_dir = os.getcwd()

        # Work under a temporary directory
        logging.info(f"Using temporary directory for training: {args.pdb_folder}")
        # Use TemporaryDirectory only if not debugging, otherwise use a persistent directory
        if args.debug:
            temp_dir = "aiofff_training_debug"
            os.makedirs(temp_dir, exist_ok=True)
            logging.info(f"Using persistent debug directory at {temp_dir}")
            cleanup_temp_dir = False
        else:
            temp_dir_obj = tempfile.TemporaryDirectory(prefix="aiofff_training_")
            temp_dir = temp_dir_obj.name
            logging.info(f"Created temporary directory at {temp_dir}")
            cleanup_temp_dir = True

        # go to the temporary directory
        os.chdir(temp_dir)
        # print current working directory
        logging.info(f"Current working directory: {os.getcwd()}")

        feature_lines = ["Distance,AverageDensity_Near,AverageDensity_Mid,AverageDensity_Far,AverageD2surface,PBPotential"]
        for pdb in pdbs:
            base_name = os.path.basename(pdb)
            logging.info(f"Processing PDB file: {base_name}")
            
            logging.info(f"Running step1.py")
            shutil.copy(pdb, base_name)
            result = subprocess.run([f"step1.py", f"{base_name}"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if result.returncode != 0:
                logging.error(f"step1.py failed on {base_name} with exit code {result.returncode}")
                exit(result.returncode)

            logging.info(f"Running step2.py")
            result = subprocess.run([f"step2.py", "--writepdb"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if result.returncode != 0:
                logging.error(f"step2.py failed on {base_name} with exit code {result.returncode}")
                exit(result.returncode)
            shutil.copy("ga_output/state_0001.pdb", "microstate.pdb")
            
            logging.info(f"Running local_density.py")
            result = subprocess.run([f"local_density.py", "microstate.pdb"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if result.returncode != 0:
                logging.error(f"local_density.py failed on microstate.pdb with exit code {result.returncode}")
                exit(result.returncode)

            logging.info(f"Setting charge and radius and write step2_out.pdb")
            atoms = read_pdb("microstate.pdb")
            assign_charges(atoms)
            write_pdb(atoms, "step2_out.pdb")

            logging.info(f"Running step3.py")
            result = subprocess.run([f"step3.py", "-s", "delphi", "-p", "3"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if result.returncode != 0:
                logging.error(f"step3.py failed with exit code {result.returncode}")
                exit(result.returncode)
            else:
                logging.info(f"step3.py completed successfully.")

            logging.info("Compile density and electrostatics energy to a csv file.")
            atoms = load_atoms()
            update_density_score(atoms, "microstate.density")
            pairwise_ele = get_electrostatic_energy(atoms)

            # Append to the output lines
            for (atom_id1, atom_id2), ele in pairwise_ele.items():
                atom1, atom2 = atoms[atom_id1], atoms[atom_id2]
                distance = atom1.xyz.distance(atom2.xyz)
                density1_near = atom1.density_near
                density1_mid = atom1.density_mid
                density1_far = atom1.density_far
                d2surface1 = atom1.d2surface
                density2_near = atom2.density_near
                density2_mid = atom2.density_mid
                density2_far = atom2.density_far
                d2surface2 = atom2.d2surface

                average_density_near = (density1_near + density2_near) / 2
                average_density_mid = (density1_mid + density2_mid) / 2
                average_density_far = (density1_far + density2_far) / 2
                average_d2surface = (d2surface1 + d2surface2) / 2

                feature_lines.append(f"{distance:.3f},{average_density_near:.1f},{average_density_mid:.1f},{average_density_far:.1f},{average_d2surface:.3f},{ele:.3f}")

            logging.info(f"Processed pdb {base_name} with {len(pairwise_ele)} pairs of atoms.")
        # Only cleanup if not debugging
        if cleanup_temp_dir:
            temp_dir_obj.cleanup()

        # Change back to the original directory
        os.chdir(current_dir)
        # Write the feature lines to a CSV file
        output_file = f"features_ele.csv"
        logging.info(f"Writing features and electrostatic potential to {output_file}")
        # print current working directory
        with open(output_file, 'w') as f:
            f.write("\n".join(feature_lines))
        logging.info(f"Feature writing completed. The energy unit is kcal/mol, distance unit is Angstrom.")
        csv_file = output_file
    else:  # it has to be fit by CSV file
        csv_file = args.fit_csv


    # Move on to training the model
    logging.info(f"Training the model with Neural Network from feature file {csv_file}")
    df = pd.read_csv(csv_file)
    df["iDistance"] = 1 / df["Distance"] # Inverse distance
    features = ['iDistance', 'AverageDensity_Near', 'AverageDensity_Mid', 'AverageDensity_Far', 'AverageD2surface']
    target = 'PBPotential'
    X = df[features]
    y = df[target]
    fit_ann(X, y, "ANN Model")
