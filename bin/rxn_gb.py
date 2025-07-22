#!/usr/bin/env python
"""
Use Generalized Born (GB) model to calculate reaction energies and use ML to fit to Delphi rxn reaction energies.
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


def parse_args():
    parser = argparse.ArgumentParser(description="Calculate reaction energies using GB model and fit ML model to Delphi energies.")
    parser.add_argument("--rxn_source", type=str, default="./", help="Source working directory with delphi calculated rxn, default is the current working directory")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.INFO)

    # Load reaction data
    rxn_source = os.path.join(args.rxn_source, "energies")
    if not os.path.exists(rxn_source):
        logging.error(f"Reaction source directory {rxn_source} does not exist.")
        exit(1)

    # Collect reaction energies from Delphi
    delphi_rxn = {}
    for filename in os.listdir(rxn_source):
        if filename.endswith(".opp"):
            conf_id = filename[:14]
            lines = open(os.path.join(rxn_source, filename), 'r').readlines()
            for line in lines:
                fields = line.split()
                conf2_id = fields[1][:14]
                ele_average = float(fields[2])
                ele_raw = float(fields[5])
                delphi_energies[(conf_id, conf2_id)] = (ele_average, ele_raw)

    logging.info(f"Loaded Delphi reaction energies: {len(delphi_energies)} data points")

    # Prepare data for ML model
    X, y = [], []
    for (conf_id, conf2_id), (ele_average, ele_raw) in delphi_energies.items():
        X.append([ele_average])  # Features can be extended later
        y.append(ele_raw)

    X = np.array(X)
    y = np.array(y)

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Scale the features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Train a simple MLP regressor
    mlp_model = MLPRegressor(hidden_layer_sizes=(100,), max_iter=500, random_state=42)
    mlp_model.fit(X_train_scaled, y_train)

    # Save the model and scaler
    joblib.dump(mlp_model, 'mlp_model.pkl')
    joblib.dump(scaler, 'scaler.pkl')

    # Evaluate the model
    y_pred = mlp_model.predict(X_test_scaled)
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)