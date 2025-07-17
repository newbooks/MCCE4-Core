#! /usr/bin/env python3
"""
Step 3 Using Fast Force Field (Machine Learning based)
"""


import argparse
import logging
import subprocess
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial import ConvexHull
import pandas as pd
import joblib
from mcce.pdbio import *
from mcce.main import *
from mcce.constants import *


# Constants
Far_Radius = 15.0   # Far limit to count far local density
Mid_Radius = 6.0    # Mid limit to count mid local density
Near_Radius = 3.0   # Near limit to count near local density

def train_inplace(protein):
    """
    Train a Fast Force Field model using the pre-calculated opp files.
    For now, this function mocks training and returns a placeholder model file path.
    """
    # Given a step2_out.pdb, we will run PB solver for one picked atom on the first side chain only.
    # This gives us an idea of the protein shape and how it affects the electrostatic energy.
    # For now, we just use the pre-calculated opp files.
    # In 4lzt, we have 389 opp files, but we only need 140 PB runs at most for 140 residues.
    pass

def ann_predict(model, df):
    # Standardize the features
    X = df[model['features']]
    #rint(X.head())
    scaler = model['scaler']
    X_scaled = scaler.transform(X)
    ele_modifiers = model['model'].predict(X_scaled)
    return ele_modifiers


def compose_opp(protein, model_file):
    """
    Compose the output files for the given protein using the specified model file.
    Optimized: Only consider atom pairs within 20A, and batch ann_predict calls.
    """
    output_folder = STEP3_LOOKUP
    native_pdb = "native.pdb"
    DIST_CUTOFF = 20.0  # Only consider atom pairs within 20A

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    else:  # Clear the folder if it exists
        for f in os.listdir(output_folder):
            os.remove(os.path.join(output_folder, f))

    # native protein to get local density reference
    reference_points = []
    for res in protein.residues:
        for conf in res.conformers[:2]:  # Use only the first two conformers
            for atom in conf.atoms:
                reference_points.append(atom.xyz.to_np())

    # make a cKDTree using the reference points
    reference_points = np.array(reference_points)
    density_tree = cKDTree(reference_points)
    hull = ConvexHull(reference_points)
    surface_points = reference_points[hull.vertices]
    surface_tree = cKDTree(surface_points)
    logging.info("   Reference points organized into a spatial index and hull vertices.")

    # For each atom in the protein, calculate the local density and distance to surface
    for res in protein.residues:
        for conf in res.conformers:
            for atom in conf.atoms:
                xyz = atom.xyz.to_np()
                # Query the tree for neighbors
                indices_far = density_tree.query_ball_point(xyz, r=Far_Radius)
                indices_mid = density_tree.query_ball_point(xyz, r=Mid_Radius)
                indices_near = density_tree.query_ball_point(xyz, r=Near_Radius)
                # Calculate local density for this atom
                atom.density_near = len(indices_near) - 1  # Exclude itself
                atom.density_mid = len(indices_mid) - len(indices_near)
                atom.density_far = len(indices_far) - len(indices_mid)
                # Query the surface tree for distance to surface
                atom.d2surface = surface_tree.query(xyz)[0]  # Get the distance
    logging.info("   Local density and distance to surface calculated for each atom.")

    # Load the model file and apply it to the protein
    model = joblib.load(model_file)

    # Loop over source conformers and calculate energies
    for res in protein.residues:
        if len(res.conformers) > 1:
            for conf in res.conformers[1:]:
                logging.info(f"Calculating energies for {conf.confid} ...")
                raw_lines = []
                raw_file = f"{conf.confid}.raw"
                source_atoms = conf.atoms
                source_coords = np.array([atom.xyz.to_np() for atom in source_atoms])

                # Precompute source atom features for efficiency
                source_features = [
                    (atom.density_near, atom.density_mid, atom.density_far, atom.d2surface, atom.charge)
                    for atom in source_atoms
                ]

                for target_res in protein.residues:
                    if target_res == res or len(target_res.conformers) < 2:
                        continue
                    for target_conf in target_res.conformers[1:]:  # Skip bkb
                        target_atoms = target_conf.atoms
                        target_coords = np.array([atom.xyz.to_np() for atom in target_atoms])
                        target_features = [
                            (atom.density_near, atom.density_mid, atom.density_far, atom.d2surface, atom.charge)
                            for atom in target_atoms
                        ]

                        # Use broadcasting to compute all pairwise distances
                        dists = np.linalg.norm(
                            source_coords[:, None, :] - target_coords[None, :, :], axis=2
                        )
                        mask = dists <= DIST_CUTOFF

                        # Prepare feature vectors for all pairs within cutoff
                        pairs = np.argwhere(mask)
                        if pairs.shape[0] == 0:
                            continue

                        # Build DataFrame in batch
                        data = []
                        for i, j in pairs:
                            Distance = dists[i, j]
                            iDistance = 1 / Distance if Distance != 0 else 0
                            s_dn, s_dm, s_df, s_d2s, s_q = source_features[i]
                            t_dn, t_dm, t_df, t_d2s, t_q = target_features[j]
                            AverageDensity_Near = (s_dn + t_dn) / 2
                            AverageDensity_Mid = (s_dm + t_dm) / 2
                            AverageDensity_Far = (s_df + t_df) / 2
                            AverageD2surface = (s_d2s + t_d2s) / 2
                            data.append([i, j, iDistance, AverageDensity_Near, AverageDensity_Mid, AverageDensity_Far, AverageD2surface, s_q, t_q])

                        df = pd.DataFrame(data, columns=['i', 'j', 'iDistance', 'AverageDensity_Near', 'AverageDensity_Mid', 'AverageDensity_Far', 'AverageD2surface', 's_q', 't_q'])

                        # Only pass the feature columns to ann_predict
                        features_df = df[['iDistance', 'AverageDensity_Near', 'AverageDensity_Mid', 'AverageDensity_Far', 'AverageD2surface']]
                        ele_modifiers = ann_predict(model, features_df)

                        # Calculate the energy in batch
                        ele = np.sum(df['s_q'].values * df['t_q'].values * ele_modifiers)
                        if abs(ele) > 0.001:  # Only write if the energy is significant
                            raw_lines.append(f"{target_conf.confid} {ele:8.3f}\n")
                # Write the raw file
                with open(os.path.join(output_folder, raw_file), 'w') as f:
                    f.writelines(raw_lines)

if __name__ == "__main__":
    # set up multiline help message
    helpmsg = "MCCE4 Step 3: Calculate Energy Lookup Table using Fast Force Field."
    parser = argparse.ArgumentParser(description=helpmsg, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", metavar="ftpl_folder", default="", help="Load from this ftpl folder instead of the default one")
    parser.add_argument("-m", metavar="model_file", default="", help="Load a force field model, DEFAULT for premade model or a user trained pkl file.")
    parser.add_argument("--debug", default=False, action="store_true", help="Print debug information")
    args = parser.parse_args()

    logging_format = "%(asctime)s %(levelname)s: %(message)s"
    logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
    logging_datefmt='%Y-%m-%d %H:%M:%S'

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format=logging_format_debug, datefmt=logging_datefmt)
    else:
        logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)

    logging.info("Step 3: Calculate Energy Lookup Table using Fast Force Field.")

    # Get prm
    prm = Runprm()                  # Create a Runprm instance and load the default runprm file
    prm.update_by_cmd(args)         # Update parameters using command line arguments
    # prm.update_by_files(args.r)     # Update parameters using additional runprm files
    prm.dump(comment="Step 3 uses these runprm parameters") # Save the parameters to run.prm.record

    # Get tpl
    tpl = Tpl()                     # Create a Tpl instance
    tpl.load_ftpl_folder(prm._FTPL_FOLDER.value) # Load the ftpl folder specified in runprm
    if os.path.isdir(USER_PARAM):
        tpl.load_ftpl_folder(USER_PARAM)  # Load user defined ftpl files
    if os.path.isfile(NEW_FTPL):
        tpl.load_ftpl_file(NEW_FTPL)    # Load new.ftpl
    tpl.dump(comment="Step 3 uses these tpl parameters") # Save the parameters to tpl.dat.record


    # Read step2_out.pdb and construct mcce protein object
    # Get protein from mccepdb step1_out.pdb
    mcce = MCCE(prm=prm, tpl=tpl)
    mcce.load_mccepdb(STEP2_OUT)    # load mccepdb within mcce object as tpl is required to load the mccepdb
    logging.info(f"   Protein loaded from {STEP2_OUT}")

    # Decide if we have a model file, if not, construct one from precaculated data
    if args.m:
        if args.m.upper() == "DEFAULT":
            model_file = prm._FASTFF_MODEL.value
        else:
            model_file = args.m
    else: # Train a new model in place
        model_file = train_inplace(mcce.protein)

    logging.info(f"   Using Fast Force Field model: {model_file}")


    # Calculate energy lookup table using the Fast Force Field model
    compose_opp(mcce.protein, model_file)