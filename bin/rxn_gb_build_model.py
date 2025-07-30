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
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.neural_network import MLPRegressor
import joblib
from matplotlib import pyplot as plt
import seaborn as sns
from mcce.geom import Vector

DEFAULT_TRAINING = "trainingset"

class Atom:
    def __init__(self, line):
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.charge = float(line[54:60]) if len(line) > 60 else 0.0
        self.radius = float(line[60:66]) if len(line) > 66 else 0.0
        self.radius_born = 0.0  # Placeholder for Born radius

def parse_args():
    parser = argparse.ArgumentParser(description="Model atom Born radii from local density. It can train from pdb files, or a precompiled csv file.")
    parser.add_argument("--from_pdb", type=str, default=None, help="The path to the folder with training set as pdb files")
    parser.add_argument("--from_csv", type=str, default=None, help="The path to the precompiled csv file")
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
    pass


if __name__ == "__main__":
    args = parse_args()
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.INFO)

    # Check which input to use, either from PDB files or from a precompiled CSV file
    if args.from_pdb is None and args.from_csv is None:
        # Use default training set
        logging.info("No input specified, using default training set.")
        pdb_files = [f for f in os.listdir(DEFAULT_TRAINING) if f.endswith('.pdb')]
        if not pdb_files:
            logging.error("No PDB files found in the specified directory.")
            exit(1)

        for pdb_file in pdb_files:
            csv_lines = pdb2csv(pdb_file)
            logging.info(f"Processed {pdb_file} to get local density and rxn.")

    elif args.from_pdb and args.from_csv is None:
        logging.info("Training from PDB files.")
        # get list of pdb files
        pdb_files = [f for f in os.listdir(args.from_pdb) if f.endswith('.pdb')]
        if not pdb_files:
            logging.error("No PDB files found in the specified directory.")
            exit(1)

        for pdb_file in pdb_files:
            csv_lines = pdb2csv(pdb_file)
            logging.info(f"Processed {pdb_file} to get local density and rxn.")


    elif args.from_csv and args.from_pdb is None:
        logging.info("Training from CSV file.")
    else:
        logging.error("Please specify either --from_pdb or --from_csv, not both.")
        exit(1)

    
