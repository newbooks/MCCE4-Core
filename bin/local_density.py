#!/usr/bin/env python

"""
MCCE4 Tool: Standalone Program to Calculate Atom Local Density Score
Local density score measures how many atoms are within near, mid and far distances from the atom to account for the local environment.
"""

import logging
import argparse
import time
import numpy as np
from scipy.spatial import cKDTree
from mcce.geom import *

# Constants
Far_Radius = 15.0   # Far limit to count far local density
Mid_Radius = 6.0    # Mid limit to count mid local density

logging_format = "%(asctime)s %(levelname)s: %(message)s"
logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
logging_datefmt='%Y-%m-%d %H:%M:%S'

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
    helpmsg = "Calculate atom local density scores for a given protein structure (microstate)"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-d", "--debug", default=False, action="store_true", help="Print debug information")
    parser.add_argument("pdb_file", nargs='+', default=[], help="PDB file(s) to be processed")
    
    return parser.parse_args()

class Atom:
    def __init__(self):
        self.line = ""  # Original line from PDB file
        self.element = ""
        self.xyz = Vector()
        self.radius = 0.0
        self.local_density_far = 0   # Local density score, integer, count of atoms within Far_Radius
        self.local_density_mid = 0   # Local density score, integer, count of atoms within Mid_Radius

    def __repr__(self):
        return f"{self.line[:30]}{self.xyz.x:8.3f}{self.xyz.y:8.3f}{self.xyz.z:8.3f}{self.radius:8.3f}{self.local_density_mid:8d}{self.local_density_far:8d}"

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
                    atom = Atom()
                    atom.line = line.strip()
                    atom.element = line[76:78].strip()  # Extract element symbol
                    atom.xyz = Vector((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                    atom.radius = ATOM_RADII.get(atom.element, ATOM_RADIUS_UNKNOWN)  # Set radius based on element
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
        indices_far = tree.query_ball_tree(tree, r=Far_Radius)  # Get all neighbors within Far_Radius
        indices_mid = tree.query_ball_tree(tree, r=Mid_Radius)  # Get all neighbors within Mid_Radius
        for i, atom in enumerate(self.atoms):
            # Get the indices of neighbors for this atom
            neighbor_indices_far = indices_far[i]
            neighbor_indices_mid = indices_mid[i]
            # Count neighbors excluding itself
            atom.local_density_far = len(neighbor_indices_far) - 1
            atom.local_density_mid = len(neighbor_indices_mid) - 1  

    def write_local_density(self):
        """
        Write local density scores to a file.
        """
        output_file = f"{self.pdb_file.rsplit('.', 1)[0]}.density"
        with open(output_file, 'w') as f:
            for atom in self.atoms:
                f.write(f"{atom}\n")



if __name__ == "__main__":
    args = parse_arguments()
    logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)
    logging.info("Starting local density calculation...")
    if args.debug:
        run_times = []
        for i in range(5):
            logging.info(f"Run {i+1} of 5")
            # Start timing
            time_start = time.time()
            for pdb_file in args.pdb_file:
                logging.info(f"Processing PDB file: {pdb_file}")
                protein = Protein()
                protein.load_pdb(pdb_file)
                logging.info(f"Loaded {len(protein.atoms)} atoms from {pdb_file}")
                protein.calculate_local_density()
                protein.write_local_density()
                logging.info(f"Local density scores written to {pdb_file.rsplit('.', 1)[0]}.density")
            time_elapsed = time.time() - time_start
            run_times.append(time_elapsed)
        
        # report average run time and standard deviation
        avg_time = sum(run_times) / len(run_times)
        std_dev_time = (sum((x - avg_time) ** 2 for x in run_times) / len(run_times)) ** 0.5
        print(f"Average run time: {avg_time:.2f} seconds, Standard Deviation: {std_dev_time:.2f} seconds")
    else:
        for pdb_file in args.pdb_file:
            logging.info(f"Processing PDB file: {pdb_file}")
            protein = Protein()
            protein.load_pdb(pdb_file)
            logging.info(f"Loaded {len(protein.atoms)} atoms from {pdb_file}")
            protein.calculate_local_density()
            protein.write_local_density()
            logging.info(f"Local density scores written to {pdb_file.rsplit('.', 1)[0]}.density")


