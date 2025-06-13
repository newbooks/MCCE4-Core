#!/usr/bin/env python

"""
MCCE4 Tool: Standalone Program to Calculate Atom Local Density Score
Local density score measures how many atoms are within a certain distance from the atom.
The cut off distance is different than embedding depth in three ways:
- it is not a grid based calculation, but a distance based atom counts.
- The search sphere is bigger than embedding depth in order to capture more local environment.
- The search sphere is a constant radius, not a depth that depends on atom radius + probe radius.
"""

import logging
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time
from mcce.geom import *

# Constants
Local_Density_Radius = 10.0  # Radius in Angstroms to count local density, embedding depth box size is about 7.0 Angstroms

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
    helpmsg = "Calculate atom local density score for a given protein structure (microstate)"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("pdb_file", nargs='+', default=[], help="PDB file(s) to be processed")
    
    return parser.parse_args()

class Atom:
    def __init__(self):
        self.line = ""  # Original line from PDB file
        self.element = ""
        self.xyz = Vector()
        self.radius = 0.0
        self.local_density = 0   # Local density score, integer, count of atoms within Local_Density_Radius
        self.inrange_atoms = []  # Atoms within the local density radius

    def __repr__(self):
        return f"{self.line[:30]}{self.xyz.x:8.3f}{self.xyz.y:8.3f}{self.xyz.z:8.3f}{self.radius:8.3f}{self.local_density:8d}"
    
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
        Calculate local density for each atom using a 3D grid to accelerate neighbor search.
        Local density is the number of atoms within Local_Density_Radius from the atom.
        """
        # Build 3D grid
        grid_size = Local_Density_Radius
        grid = {}
        def grid_index(x, y, z):
            return (int(x // grid_size), int(y // grid_size), int(z // grid_size))

        # Assign atoms to grid cells
        for atom in self.atoms:
            idx = grid_index(atom.xyz.x, atom.xyz.y, atom.xyz.z)
            grid.setdefault(idx, []).append(atom)

        range_squared = Local_Density_Radius * Local_Density_Radius

        # For each atom, only check atoms in the same and neighboring grid cells
        for atom in self.atoms:
            idx = grid_index(atom.xyz.x, atom.xyz.y, atom.xyz.z)
            neighbors = []
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        neighbor_idx = (idx[0]+dx, idx[1]+dy, idx[2]+dz)
                        for other in grid.get(neighbor_idx, []):
                            # Avoid double counting by only considering atoms with higher index
                            if id(other) < id(atom):
                                continue
                            if abs(atom.xyz.x - other.xyz.x) < Local_Density_Radius and \
                               abs(atom.xyz.y - other.xyz.y) < Local_Density_Radius and \
                               abs(atom.xyz.z - other.xyz.z) < Local_Density_Radius:
                                distance_squared = atom.xyz.distance_squared(other.xyz)
                                if distance_squared < range_squared:
                                    atom.inrange_atoms.append(other)
                                    other.inrange_atoms.append(atom)
            # No need to set local_density here; do it after all pairs are processed

        # Set local density score as the number of in-range atoms
        for atom in self.atoms:
            atom.local_density = len(atom.inrange_atoms)



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
    # Timing of the execution
    # Original: 12.59+-0.14 (Office PC), 13.52+-0.25 (Home PC)
    # Improved version:


        # Compare with embedding scores if available
        # embedding_file = f"{pdb_file.rsplit('.', 1)[0]}.embedding"
        # if os.path.isfile(embedding_file):
        #     logging.info(f"Comparing with embedding scores from {embedding_file}")
        #     embedding_scores = {}
        #     with open(embedding_file, 'r') as ef:
        #         for line in ef:
        #             parts = line.strip().split()
        #             if len(parts) >= 6:
        #                 atom_line = line[:54]
        #                 embedding_score = float(parts[-1])
        #                 embedding_scores[atom_line] = embedding_score
        #     with open(f"{pdb_file.rsplit('.', 1)[0]}_density_vs_embedding.csv", 'w') as df:
        #         df.write("Atom,LocalDensity,EmbeddingScore\n")
        #         for atom in protein.atoms:
        #             atom_line = atom.line[:54]
        #             embedding_score = embedding_scores.get(atom_line, 0.0)
        #             df.write(f"{atom_line},{atom.local_density},{embedding_score}\n")
        #     logging.info(f"Density vs embedding scores written to {pdb_file.rsplit('.', 1)[0]}_density_vs_embedding.csv")
        #     df = pd.read_csv(f"{pdb_file.rsplit('.', 1)[0]}_density_vs_embedding.csv")
        #     plt.figure(figsize=(10, 6))
        #     sns.scatterplot(data=df, x="LocalDensity", y="EmbeddingScore")
        #     plt.title(f"Local Density vs Embedding Score for {pdb_file}")
        #     plt.xlabel("Local Density")
        #     plt.ylabel("Embedding Score")
        #     plt.grid()
        #     plt.savefig(f"{pdb_file.rsplit('.', 1)[0]}_density_vs_embedding.png")
        #     plt.show()
        #     logging.info(f"Plot saved as {pdb_file.rsplit('.', 1)[0]}_density_vs_embedding.png")
        #     plt.close()
        # else:
        #     logging.info(f"No embedding scores found for {pdb_file}, skipping comparison.")
