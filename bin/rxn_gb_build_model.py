#!/usr/bin/env python
"""
Use Generalized Born (GB) model to calculate reaction energies.

When a single atom with a unit charge is placed in a multi-atom system, the relationship between the Born radius and reaction field can be described by the inverted linear relationship.
RXN = -0.5*(1 - 1/eps) * q^2 / r_Born  (q = 1.0)
This is to convert the reaction field to Born radius, where RXN is the reaction energy, eps is the dielectric constant, and r_Born is the Born radius.

Use Density_Near, Density_Mid, Density_Far and D2Surface as features, the Born radius as target, train machine learning model to calculate Born radius from local density terms.
The Born radius is then calculated for each atom.

When we need reaction field energy, use Born model to calculate the reaction field energy analytically.

Since we need reference reaction field energy, we need to have a training set of amino acids and cofactors.

Step 1: Generalized Born effective distance f_ij
The effective distance between two atoms i and j is calculated as:
f_ij = sqrt( (r_ij^2 + r_Born_i * r_Born_j * exp(-r_ij^2/(4*r_Born_i*r_Born_j))) )
This is the effective distance formula for two atoms i and j, where r_ij is the distance between the two atoms, r_Born_i and r_Born_j are the Born radii of the two atoms.

Note when r_ij is much larger than r_Born_i and r_Born_j, the effective distance approaches f_ij.

Step 2. Calculate reaction energies using the effective distance
RXN = -0.5 * (1 - 1/eps) * Sum(qi * qj / f_ij)

This script is going to read training structure files, set up delphi calculations, and train a machine model to predict Born radii 
for atoms from local density terms.

A separate script will be used to calculate reaction energies using the trained model.
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
from scipy.spatial import cKDTree, ConvexHull
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.neural_network import MLPRegressor
import joblib
from matplotlib import pyplot as plt
import seaborn as sns
from mcce.geom import Vector
import atexit

DEFAULT_TRAINING = "trainingset"

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

# Conversion factors
# RXN = 166 * (-0.5 * (1/eps_in - 1/eps_out) * qi * qj / r_Born)
# where 166 is the conversion factor to calculate RXN in kcal/mol from unit charge and distance in Angstroms.
K_rxn = 166.0 * (-0.5 * (1/4.0 - 1/80.0))  # eps_in = 4, eps_out = 80 
# RXN = K_rxn * qi * qj / r_Born

# Constants
Far_Radius = 15.0   # Far limit to count far local density
Mid_Radius = 6.0    # Mid limit to count mid local density
Near_Radius = 3.0   # Near limit to count near local density

CSV_OUTPUT_FILE = "local_density.csv"

class Atom:
    def __init__(self, line):
        self.line = line
        atom_name = line[12:16]
        self.name = atom_name
        if len(atom_name.strip()) == 4 and atom_name[0] == 'H':
            self.element = " H"
        else:
            self.element = atom_name[:2]
        self.xyz = Vector((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        self.charge = 0.0 # set all charges to 0
        self.radius = ATOM_RADII.get(self.element, ATOM_RADIUS_UNKNOWN)  # Get radius from ATOM_RADII to avoid or use default
        self.conf_id = line[17:20] + line[80:82] + line[21:30]

        self.density_near = 0
        self.density_mid = 0
        self.density_far = 0
        self.d2surface = 0.0  # Distance to the nearest surface
        self.rxn = None  # Placeholder for reaction field energy
        self.radius_born = 0.0  # Placeholder for Born radius

    def csvline(self):
        """
        Return a CSV line for this atom.
        The fields are:
        - Density_Near
        - Density_Mid
        - Density_Far
        - D2Surface
        - RXN
        - Radius_Born
        """
        if self.rxn is None:
            logging.error(f"RXN of side chain atom {self.name} on conformer {self.conf_id} was not calculated. This is not expected.")
            exit(1)

        # Calculate the Born radius from the reaction field energy
        # The formula is derived from the RXN = -0.5 * (1/eps_in - 1/eps_out) * q^2 / r_Born
        # Rearranging gives us r_Born = -0.5 * (1/eps_in - 1/eps_out) * q^2 / RXN
        # Assuming eps_in = 4 and eps_out = 80, and q = 1.0

        self.radius_born = K_rxn / self.rxn

        return f"{self.density_near},{self.density_mid},{self.density_far},{self.d2surface:.3f},{self.rxn:.3f},{self.radius_born:.3f}\n"


class Protein:
    def __init__(self):
        self.pdb_file = ""
        self.atoms = []  # List of Atom objects

    def load_pdb(self, pdb_file):
        """
        Load PDB file and initialize atoms.
        """
        self.pdb_file = pdb_file
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom = Atom(line)
                    self.atoms.append(atom)

    def calculate_local_density(self):
        """
        Calculate local density (number of atoms within Local_Density_Radius) for each atom
        using cKDTree for efficient nearest neighbor search.
        """
        if not self.atoms:
            logging.warning("No atoms loaded, skipping local density calculation.")
            return
        # Create a list of coordinates for cKDTree
        coords = np.array([atom.xyz.to_np() for atom in self.atoms])
        # Create a cKDTree for fast nearest neighbor search
        tree = cKDTree(coords)
        # Query the tree for neighbors
        indices_far = tree.query_ball_point(coords, r=Far_Radius)  # Get all neighbors within Far_Radius
        indices_mid = tree.query_ball_point(coords, r=Mid_Radius)  # Get all neighbors within Mid_Radius
        indices_near = tree.query_ball_point(coords, r=Near_Radius)  # Get all neighbors within Near_Radius
        # Calculate local density for each atom
        for i, atom in enumerate(self.atoms):
            # Get the indices of neighbors for this atom and exclude itself
            neighbor_count_far = len(indices_far[i]) - 1
            neighbor_count_mid = len(indices_mid[i]) - 1
            neighbor_count_near = len(indices_near[i]) - 1
            atom.density_near = neighbor_count_near
            atom.density_mid = neighbor_count_mid - neighbor_count_near
            atom.density_far = neighbor_count_far - neighbor_count_mid


    def calculate_distance_to_surface(self):
        """
        Calculate distance to the nearest surface for each atom.
        Step 1: Create a set pf mesh points from the atoms by using ConvexHull.
        Step 2: For each atom, find the nearest point on the surface mesh using cKDTree.
        """
        # Step 1: Create a set of mesh points from the atoms
        atom_coords = np.array([atom.xyz.to_np() for atom in self.atoms])
        hull = ConvexHull(atom_coords)
        surface_points = atom_coords[hull.vertices]

        # Step 2: For each atom, find the nearest point on the surface mesh
        tree = cKDTree(surface_points)
        for atom in self.atoms:
            dist, _ = tree.query(atom.xyz.to_np())
            atom.d2surface = dist

    def write_pqr(self, filename):
        """
        Write the protein atoms to a PQR file.
        The PQR format is similar to PDB but includes charge and radius.
        """
        with open(filename, 'w') as f:
            for atom in self.atoms:
                line = f"{atom.line[:54]}{atom.radius:8.3f}{atom.charge:12.3f}{atom.line[74:].rstrip()}\n"
                f.write(line)


def parse_args():
    parser = argparse.ArgumentParser(description="Model atom Born radii from local density. It can train from pdb files, or a precompiled csv file.")
    parser.add_argument("--from_pdb", type=str, default=None, help="The path to the folder with training set as pdb files")
    parser.add_argument("--from_csv", type=str, default=None, help="The path to the precompiled csv file")
    parser.add_argument("--debug", action='store_true', help="Enable debug mode to reserve temporary folder")
    return parser.parse_args()


def pdb2csv(pdb_file):
    """
    Setup delphi calculation to get reaction field energy and local density terms.
    In the resulted csv file, we will have the following columns:
    - Density_Near
    - Density_Mid
    - Density_Far
    - D2Surface
    - RXN
    - BORN_RADIUS
    """
    # Set up temporary directory to run step1.py and step2.py
    temp_dir = tempfile.mkdtemp(prefix="rxn_training_")

    logging.info(f"Using temporary directory: {temp_dir}")
    cwd = os.getcwd()
    shutil.copy(pdb_file, os.path.join(temp_dir, "prot.pdb"))

    # Change to the temporary directory
    os.chdir(temp_dir)

    # Run step1.py to generate a mcce protein object
    logging.info(f"Running step1.py to generate mcce protein object for {pdb_file} ...")
    step1_cmd = ["step1.py", "prot.pdb", "--no_hoh", "--no_ter"]
    subprocess.run(step1_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # Run step2.py to get single conormer microstates
    logging.info(f"Running step2.py to get single conformer microstates for {pdb_file} ...")
    step2_cmd = ["step2.py", "--writepdb"]
    subprocess.run(step2_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # copy ga_output/state_0001.pdb to microstate.pdb
    shutil.copy("ga_output/state_0001.pdb", "microstate.pdb")

    # Get local density terms
    logging.info(f"Getting local density terms ...")
    protein = Protein()
    protein.load_pdb("microstate.pdb")
    protein.calculate_local_density()
    protein.calculate_distance_to_surface()

    # Set up delphi calculation and get reaction field energy
    logging.info(f"Setting up delphi calculation for {pdb_file} ...")
    # Make a dictionary to hold side chain atoms, atom is referenced and we will alter its charge
    sidechain_atoms = {}
    for atom in protein.atoms:
        # collect all atoms of the same side chain
        sidechain_id = atom.line[17:30]
        if sidechain_id[-3:] != "000":
            if sidechain_id not in sidechain_atoms:
                sidechain_atoms[sidechain_id] = []
            sidechain_atoms[sidechain_id].append(atom)

    i_charge = 0  # Loop through side conformer atoms and set their charge to +1
    assigned_charges = 1 # counter of assigned charges, initially set to 1 so it enters the loop
    while assigned_charges > 0:
        assigned_charges = 0
        atom_lookup_by_conf_id = {}
        for sidechain_id, atoms in sidechain_atoms.items():
            if i_charge < len(atoms):
                atoms[i_charge].charge = 1.0  # set the charge of the atom to +1
                assigned_charges += 1
                atom_lookup_by_conf_id[atoms[i_charge].conf_id] = atoms[i_charge]
        # Write out protein as a step2_out.pdb file
        protein.write_pqr("step2_out.pdb")

        # Restore the original charges in protein object after writing out the pdb file
        for sidechain_id, atoms in sidechain_atoms.items():
            if i_charge < len(atoms):
                atoms[i_charge].charge = 0.0  # set the charge of the atom back to 0


        # Now that we have step2_out.pdb, call step3.py
        logging.info(f"Running step3.py on the {i_charge}th atom of {assigned_charges} side chains ...")
        result = subprocess.run(
            ["step3.py", "-s", "delphi", "-p", "6"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        if result.returncode != 0:
            logging.error("step3.py failed to run successfully.")
            exit(result.returncode)
        
        # Now collect rxn from raw files and assign back to atoms
        logging.info("Collecting reaction field energies from raw files ...")
        for conf_id, atom in atom_lookup_by_conf_id.items():
            raw_file = f"energies/{conf_id}.raw"
            if os.path.isfile(raw_file):
                with open(raw_file, 'r') as f:
                    pbrxn = None
                    for line in f:
                        if line.startswith("[RXN"):
                            fields = line.strip().split()
                            if len(fields) >= 2:
                                pbrxn = float(fields[-1])
                                break
                    if pbrxn is not None:
                        atom.rxn = pbrxn
                        logging.debug(f"Found RXN {pbrxn:.3f} for charged atom {atom.name} in {conf_id}.")
                    else:
                        logging.error(f"RXN value not found in {raw_file} for charged atom {atom.name} in {conf_id}.")
                        exit(1)

            else:
                logging.warning(f"Raw file {raw_file} not found for charged atom {atom.name} in {conf_id}.")

        i_charge += 1  # Move to the next atom in the side chain


    csv_lines = []
    for atom in protein.atoms:
        if atom.rxn is not None:
            csv_line = atom.csvline()
            csv_lines.append(csv_line)

    os.chdir(cwd)  # Change back to the original directory
    if not args.debug:
        # Clean up the temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    logging.info(f"Finished processing {pdb_file} and cleaned the temporary directory.")
    return csv_lines

if __name__ == "__main__":
    args = parse_args()
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.INFO)

    # Check which input to use, either from PDB files or from a precompiled CSV file
    csv_header = "Density_Near,Density_Mid,Density_Far,D2Surface,RXN,BORN_RADIUS\n"
    csv_all_lines = [csv_header]
    training_from_pdb = False
    if args.from_pdb is None and args.from_csv is None:
        # Use default training set
        logging.info("No input specified, using default training set.")
        pdb_files = [os.path.join(DEFAULT_TRAINING, f) for f in os.listdir(DEFAULT_TRAINING) if f.endswith('.pdb')]
        training_from_pdb = True

    elif args.from_pdb and args.from_csv is None:
        logging.info("Training from PDB files.")
        # get list of pdb files
        pdb_files = [os.path.join(args.from_pdb, f) for f in os.listdir(args.from_pdb) if f.endswith('.pdb')]
        training_from_pdb = True

    if training_from_pdb:
        if not pdb_files:
            logging.error("No PDB files found in the specified directory.")
            exit(1)

        for pdb_file in pdb_files:
            logging.info(f"Processing {pdb_file} ...")
            csv_lines = pdb2csv(pdb_file)
            logging.info(f"Processed {pdb_file} to get local density and rxn.")
            csv_all_lines.extend(csv_lines)
        # Write the CSV file
        open(CSV_OUTPUT_FILE, 'w').writelines(csv_all_lines)


    elif args.from_csv and args.from_pdb is None:
        logging.info("Training from CSV file.")

    else:
        logging.error("Please specify either --from_pdb or --from_csv, not both.")
        exit(1)

