#!/usr/bin/env python
"""
MCCE4 Tool: Find the box size for protein structure
"""

import argparse
import logging

class Box:
    def __init__(self, min_x, max_x, min_y, max_y, min_z, max_z):
        self.min_x = min_x
        self.max_x = max_x
        self.min_y = min_y
        self.max_y = max_y
        self.min_z = min_z
        self.max_z = max_z

    def size(self):
        return (self.max_x - self.min_x, self.max_y - self.min_y, self.max_z - self.min_z)
    
    def dimensions(self):
        return (self.min_x, self.max_x, self.min_y, self.max_y, self.min_z, self.max_z)
    
    def volume(self):
        return self.size()[0] * self.size()[1] * self.size()[2]

    def __repr__(self):
        return f"Box(min_x={self.min_x:.3f}, max_x={self.max_x:.3f}, min_y={self.min_y:.3f}, max_y={self.max_y:.3f}, min_z={self.min_z:.3f}, max_z={self.max_z:.3f})"


def parse_arguments():
    helpmsg = "Find the box size for protein structure"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("pdb_file", nargs='+', default=[], help="PDB file(s) to be processed")
    
    return parser.parse_args()

def find_box_size(pdb_file):
    """
    Find the box size for the given PDB file.
    """

    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    # Extract coordinates from ATOM/HETATM lines
    coords = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append((x, y, z))
            except ValueError:
                logging.warning(f"Invalid coordinate in line: {line.strip()}")
    if not coords:
        logging.error("No valid ATOM/HETATM lines found in the PDB file.")
        return None
    # Calculate the bounding box
    min_x = min(coord[0] for coord in coords)
    max_x = max(coord[0] for coord in coords)
    min_y = min(coord[1] for coord in coords)
    max_y = max(coord[1] for coord in coords)
    min_z = min(coord[2] for coord in coords)
    max_z = max(coord[2] for coord in coords)

    box = Box(min_x, max_x, min_y, max_y, min_z, max_z)
    return box

if __name__ == "__main__":
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Parse command line arguments
    args = parse_arguments()

    if not args.pdb_file:
        logging.error("No PDB files provided.")
        exit(1)

    for pdb_file in args.pdb_file:
        box = find_box_size(pdb_file)
        if box:
            print(f"Box size for {pdb_file}: ({box.size()[0]:.3f}, {box.size()[1]:.3f}, {box.size()[2]:.3f})")
            print(f"Box dimensions: {box}")
            print(f"Box volume: {box.volume():.3f} cubic Angstroms")
        else:
            print(f"Failed to find box size for {pdb_file}.")