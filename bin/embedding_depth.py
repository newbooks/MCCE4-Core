#!/usr/bin/env python

"""
MCCE4 Tool: Standalone Program to Calculate Atom Embedding Score
"""

import logging
import argparse
from scipy.spatial import cKDTree
import time
import numpy as np
from mcce.geom import *

logging_format = "%(asctime)s %(levelname)s: %(message)s"
logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
logging_datefmt='%Y-%m-%d %H:%M:%S'

# Constants
PROBE_RAD = 1.4    # Probe radius in Angstrom
GRID_SIZE = 1.0    # size of the grid in Angstroms
GRID_DEPTH = 3.0   # depth to look into, how far to count the embedding depth in addition to atom radius + probe radius

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

Max_Radius = max(list(ATOM_RADII.values()) + [ATOM_RADIUS_UNKNOWN])  # Maximum radius of atoms in the system
Grid_Expansion = Max_Radius + PROBE_RAD + GRID_DEPTH + GRID_SIZE # Maximum expansion of the grid box to account for atom radii and probe radius

def parse_arguments():
    helpmsg = "Calculate atom embedding score for a given protein structure (microstate)"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("pdb_file", nargs='+', default=[], help="PDB file(s) to be processed")
    
    return parser.parse_args()


class Atom:
    def __init__(self):
        self.line = ""  # Original line from PDB file
        self.element = ""
        self.xyz = Vector()
        self.radius = 0.0
        self.embedding = 0.0

    def __repr__(self):
        return f"{self.line[:30]}{self.xyz.x:8.3f}{self.xyz.y:8.3f}{self.xyz.z:8.3f}{self.radius:8.3f}{self.embedding:8.3f}"


class Protein:
    def __init__(self):
        self.pdb_file = ""
        self.atoms = []  # List of Atom objects
        self.box = None  # Placeholder for the box object if needed
        self.grid = None  # Placeholder for the grid object if needed
        self.neighbor_candidates = []  # List of points in the grid that are occupied by atoms and their probe radius

    def __repr__(self):
        return f"Protein: {self.pdb_file}, Atoms: {len(self.atoms)}"


    def load_atoms_from_pdb(self, pdb_file):
        self.pdb_file = pdb_file
        lines = open(pdb_file).readlines()
        
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom = Atom()
                atom.line = line.strip()
                atomname = line[12:16]
                if len(atomname.strip()) == 4 and atomname[0] == "H":
                    atom.element = " H"
                else:
                    atom.element = atomname[:2]
                atom.xyz = Vector((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                atom.radius = ATOM_RADII.get(atom.element, ATOM_RADIUS_UNKNOWN)
                self.atoms.append(atom)

    def calculate_box(self):
        """
        Calculate the grid box for the protein based on atom coordinates.
        This grid bos is a list of coordinates that occupied by atoms and their probe radius.
        """

        min_x = min(atom.xyz.x for atom in self.atoms) - Grid_Expansion
        max_x = max(atom.xyz.x for atom in self.atoms) + Grid_Expansion
        min_y = min(atom.xyz.y for atom in self.atoms) - Grid_Expansion
        max_y = max(atom.xyz.y for atom in self.atoms) + Grid_Expansion
        min_z = min(atom.xyz.z for atom in self.atoms) - Grid_Expansion
        max_z = max(atom.xyz.z for atom in self.atoms) + Grid_Expansion


        # Create a box object (if needed)
        self.box = ((min_x, min_y, min_z), (max_x, max_y, max_z))

        # Initialize the grid object (if needed)
        self.grid = np.zeros((
            int(np.ceil((max_x - min_x) / GRID_SIZE)),
            int(np.ceil((max_y - min_y) / GRID_SIZE)),
            int(np.ceil((max_z - min_z) / GRID_SIZE))
        ), dtype=bool)
    
        # Fill the grid with atom spheres
        for atom in self.atoms:
            x_idx = int((atom.xyz.x - min_x) / GRID_SIZE)
            y_idx = int((atom.xyz.y - min_y) / GRID_SIZE)
            z_idx = int((atom.xyz.z - min_z) / GRID_SIZE)
            radius = int((atom.radius + PROBE_RAD) / GRID_SIZE)
            radius_sq = radius * radius

            # Mark the grid cells within the atom's radius
            for dx in range(-radius, radius + 1):
                for dy in range(-radius, radius + 1):
                    for dz in range(-radius, radius + 1):
                        if dx*dx + dy*dy + dz*dz <= radius_sq: 
                            x_cell = x_idx + dx
                            y_cell = y_idx + dy
                            z_cell = z_idx + dz
                            self.grid[x_cell, y_cell, z_cell] = True
        
        # Convert the grid to points list. This list will be the neighbor candidates
        self.neighbor_candidates = []
        for x in range(self.grid.shape[0]):
            for y in range(self.grid.shape[1]):
                for z in range(self.grid.shape[2]):
                    if self.grid[x, y, z]:
                        # Convert grid indices back to coordinates
                        coord = [
                            min_x + x * GRID_SIZE,
                            min_y + y * GRID_SIZE,
                            min_z + z * GRID_SIZE
                        ]
                        self.neighbor_candidates.append(coord)  


    def atom_embedding_score(self):
        """Calculate the embedding depth for each atom."""
        # The query set is the atoms, and the reference set is the neighbor_candidates
        # Use KD Tree to find all neighbors within the embedding depth radius, then subtract the neighbors within the atom radius + probe radius
        # Create a KDTree from the neighbor candidates
        tree = cKDTree(self.neighbor_candidates)
        for atom in self.atoms:
            # Calculate the embedding depth radius
            embedding_depth_radius = atom.radius + PROBE_RAD + GRID_DEPTH
            
            # Query the KDTree for neighbors within the embedding depth radius
            indices = tree.query_ball_point(atom.xyz.to_np(), r=embedding_depth_radius)
            
            # Count the number of neighbors found
            total_neighbors = len(indices)
            
            # Now query for neighbors within the atom radius + probe radius
            inner_radius = atom.radius + PROBE_RAD
            inner_indices = tree.query_ball_point(atom.xyz.to_np(), r=inner_radius)
            
            # Count the number of neighbors within the inner radius
            inner_neighbors = len(inner_indices)
            
            # Estimate the total number of grids between embedding depth radius and inner radius by the volume of the spherical shell divided by the grid size
            points_in_shell = (4/3 * np.pi * (embedding_depth_radius**3 - inner_radius**3))/ (GRID_SIZE**3)

            # Calculate the embedding score as the difference between total neighbors and inner neighbors, normalized by the estimated points in the shell
            if points_in_shell > 0:
                atom.embedding = (total_neighbors - inner_neighbors) / points_in_shell
                # If the embedding score is over 1, set it to 1
                # if atom.embedding > 1:
                #     atom.embedding = 1
            else:
                # If no points in shell, set embedding to 0
                atom.embedding = 0

    def write_embedding_scores(self):
        """Write the embedding scores to a file."""
        output_file = f"{self.pdb_file.rsplit('.', 1)[0]}.depth"
        with open(output_file, "w") as f:
            for atom in self.atoms:
                f.write(f"{atom}\n")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)
    args = parse_arguments()
    
    run_times = []
    N = 5  # Number of runs for averaging
    for i in range(N):
        logging.info(f"Run {i+1} of {N} ")
        # Start timing
        time_start = time.time()
        
        for pdb in args.pdb_file:
            logging.info(f"Processing {pdb}...")            
            prot = Protein()
            prot.load_atoms_from_pdb(pdb)
            logging.info(f"Loaded {len(prot.atoms)} atoms from {pdb}")
            prot.calculate_box()
            logging.info(f"Grid Box with dimension {prot.grid.shape} calculated: Memory usage = {prot.grid.nbytes / 1024/1024:.2f} MB")
            prot.atom_embedding_score()
            logging.info(f"Atom embedding scores calculated for {pdb}")
            # 157, 149, 154 s before optimization (100 runs)
            # 50, 53, 53 s after optimization (100 runs)

            prot.write_embedding_scores()
            output_file = f"{pdb.rsplit('.', 1)[0]}.embedding"
            logging.info(f"Scores written to {output_file}")
        time_elapsed = time.time() - time_start
        run_times.append(time_elapsed)
    # report average run time and standard deviation
    avg_time = sum(run_times) / len(run_times)
    std_dev_time = (sum((x - avg_time) ** 2 for x in run_times) / len(run_times)) ** 0.5
    print(f"Average run time: {avg_time:.2f} seconds, Standard Deviation: {std_dev_time:.2f} seconds")

    # Timing of the execution on large/state_0001.pdb:
    # embedding_score.py: 3.47 +- 0.04 (Home PC)
    # embedding_depth.py: 1.85 +- 0.03 (Home PC)

    # plot the embedding scores vs embedding depth
