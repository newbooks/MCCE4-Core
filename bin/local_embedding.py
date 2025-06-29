#!/usr/bin/env python

"""
MCCE4 Tool: Standalone Program to Calculate Atom Local Embedding Score
Local embedding score measures how many embedded grids are within near, mid and far distances from the atom to account for the local environment.
"""

import logging
import argparse
import time
import numpy as np
from scipy.spatial import cKDTree
from mcce.geom import *

# Constants
Far_Radius = 15.0   # Far limit
Mid_Radius = 6.0    # Mid limit
Near_Radius = 3.0   # Near limit

logging_format = "%(asctime)s %(levelname)s: %(message)s"
logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
logging_datefmt='%Y-%m-%d %H:%M:%S'


class Box:
    def __init__(self):
        self.min_x = 0.0
        self.max_x = 0.0
        self.min_y = 0.0
        self.max_y = 0.0
        self.min_z = 0.0
        self.max_z = 0.0

    def from_atoms(self, atoms):
        if not atoms:
            logging.warning("No atoms provided to define the box.")
            return
        coords = np.array([[atom.xyz.x, atom.xyz.y, atom.xyz.z] for atom in atoms])
        self.min_x, self.min_y, self.min_z = np.min(coords, axis=0)
        self.max_x, self.max_y, self.max_z = np.max(coords, axis=0)

    def size(self):
        return (self.max_x - self.min_x, self.max_y - self.min_y, self.max_z - self.min_z)
    
    def dimensions(self):
        return (self.min_x, self.max_x, self.min_y, self.max_y, self.min_z, self.max_z)
    
    def volume(self):
        return self.size()[0] * self.size()[1] * self.size()[2]
    
    def expand(self, delta):
        """
        Expand the box by a given delta in all directions.
        """
        self.min_x -= delta
        self.max_x += delta
        self.min_y -= delta
        self.max_y += delta
        self.min_z -= delta
        self.max_z += delta

    def __repr__(self):
        return f"Box(min_x={self.min_x:.3f}, max_x={self.max_x:.3f}, min_y={self.min_y:.3f}, max_y={self.max_y:.3f}, min_z={self.min_z:.3f}, max_z={self.max_z:.3f})"


def parse_arguments():
    helpmsg = "Calculate atom local embedding scores for a given protein structure (microstate)"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("pdb_file", nargs='+', default=[], help="PDB file(s) to be processed")
    
    return parser.parse_args()

class Atom:
    def __init__(self):
        self.line = ""  # Original line from PDB file
        self.element = ""
        self.xyz = Vector()
        self.local_embedding = [0, 0, 0]  # Local embedding score, integer, count of atoms within Near, Mid and Far radius

    def __repr__(self):
        return f"{self.line[:30]}{self.xyz.x:8.3f}{self.xyz.y:8.3f}{self.xyz.z:8.3f}{self.local_embedding[0]:8d}{self.local_embedding[1]:8d}{self.local_embedding[2]:8d}"

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

    def calculate_embedding(self):
        """
        Define the grid for the protein structure at resolution of 1 A.
        """
        if not self.atoms:
            logging.warning("No atoms loaded, skipping grid definition.")
            return
        # Create a grid for the protein structure
        logging.info("   Using meshed grid for defining the protein structure...")
        grid_size = 1.0  # 1 A resolution

        box = Box()
        box.from_atoms(self.atoms)
        box.expand(Far_Radius)  # Expand the box by Far_Radius to ensure the embedding search is far enough

        # Use numpy meshgrid for faster grid point generation
        xs = np.arange(box.min_x, box.max_x + grid_size, grid_size)
        ys = np.arange(box.min_y, box.max_y + grid_size, grid_size)
        zs = np.arange(box.min_z, box.max_z + grid_size, grid_size)
        grid_coords = np.stack(np.meshgrid(xs, ys, zs, indexing='ij'), axis=-1).reshape(-1, 3)
        # print the RAM usage of grid_coords
        ram_usage = np.array(grid_coords).nbytes / (1024 ** 2)  # Convert bytes to MB
        logging.info(f"Total Grid points defined: {len(grid_coords)}, RAM usage: {ram_usage:.2f} MB")

        atom_coords = np.array([[atom.xyz.x, atom.xyz.y, atom.xyz.z] for atom in self.atoms])
        
        # get a new set of grid_coords to only include points in grid_coords within 3 A distance from any atom in atom_coords
        tree = cKDTree(atom_coords)
        near_mask = tree.query_ball_point(grid_coords, 3.0)
        mask = [len(indices) > 0 for indices in near_mask]
        buried_grid_coords = [xyz for xyz in grid_coords[mask]] # replace grid_coords with only those within 3 A of any atom
        # delete grid_coords to free up memory
        del grid_coords
        logging.info(f"Buried grid points within 3 A of any atom: {len(buried_grid_coords)}, RAM usage: {np.array(buried_grid_coords).nbytes / (1024 ** 2):.2f} MB")

        # Create a KDTree for fast nearest neighbor search
        tree = cKDTree(buried_grid_coords)
        # Calculate local embedding scores for each atom
        near_neighbors = tree.query_ball_point(atom_coords, Near_Radius)
        mid_neighbors = tree.query_ball_point(atom_coords, Mid_Radius)
        far_neighbors = tree.query_ball_point(atom_coords, Far_Radius)
        for i, atom in enumerate(self.atoms):
            # Count the number of neighbors in each radius
            near_count = len(near_neighbors[i])
            mid_count = len(mid_neighbors[i]) - near_count
            far_count = len(far_neighbors[i]) - mid_count - near_count
            atom.local_embedding = [near_count, mid_count, far_count]

    def write_local_embedding(self):
        """
        Write local embedding scores to a file.
        """
        output_file = f"{self.pdb_file.rsplit('.', 1)[0]}.embedding"
        with open(output_file, 'w') as f:
            for atom in self.atoms:
                f.write(f"{atom}\n")


if __name__ == "__main__":
    args = parse_arguments()
    logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)
    logging.info("Starting local embedding calculation...")
    for pdb_file in args.pdb_file:
        logging.info(f"Processing PDB file: {pdb_file}")
        protein = Protein()
        protein.load_pdb(pdb_file)
        logging.info(f"Loaded {len(protein.atoms)} atoms from {pdb_file}")
        protein.calculate_embedding()
        protein.write_local_embedding()
        logging.info(f"Local embedding scores written to {pdb_file.rsplit('.', 1)[0]}.embedding")