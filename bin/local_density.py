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
from scipy.spatial import ConvexHull
from mcce.geom import *

# Constants
Far_Radius = 15.0   # Far limit to count far local density
Mid_Radius = 6.0    # Mid limit to count mid local density
Near_Radius = 3.0   # Near limit to count near local density

logging_format = "%(asctime)s %(levelname)s: %(message)s"
logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
logging_datefmt='%Y-%m-%d %H:%M:%S'

class LocalDensity:
    """
    Class to hold local density scores for an atom.
    The scores are stored as a list of integers:
    [Near_Count, Mid_Count, Far_Count, Quadrant_Counts...]
    where Quadrant_Counts is a list of 8 integers representing counts in each octant.
    """
    def __init__(self):
        self.near_count = 0
        self.mid_count = 0
        self.far_count = 0
        self.variance = 0.0


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
        self.local_density = LocalDensity()  # Local density scores
        self.d2surface = 0.0  # Distance to the nearest surface

    def __repr__(self):
        return f"{self.line[:30]}{self.xyz.x:8.3f}{self.xyz.y:8.3f}{self.xyz.z:8.3f}{self.local_density.near_count:6d}{self.local_density.mid_count:6d}{self.local_density.far_count:6d}{self.d2surface:8.3f}"

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
            atom.local_density.near_count = neighbor_count_near
            atom.local_density.mid_count = neighbor_count_mid - neighbor_count_near
            atom.local_density.far_count = neighbor_count_far - neighbor_count_mid


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
            protein.calculate_distance_to_surface()
            protein.write_local_density()
            logging.info(f"Local density scores written to {pdb_file.rsplit('.', 1)[0]}.density")


