#!/usr/bin/env python
"""
Use Random Forest to predict electrostatic energy based on embedding scores and distances.
The input is a CSV file with the following columns:
- Conf1: Conformer ID for atom 1
- Conf2: Conformer ID for atom 2
- Distance: Distance between two atoms in Angstroms
- Radius1: Radius for atom 1
- Radius2: Radius for atom 2
- Embedding1: Embedding score for atom 1
- Embedding2: Embedding score for atom 2
- Density1: Density score for atom 1
- Density2: Density score for atom 2
- CoulombPotential: Coulomb potential between two atoms
- AdjustedCoulombPotential: Adjusted Coulomb potential based on embedding scores
- PBPotential: Electrostatic energy from Poisson-Boltzmann calculation
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np
from sklearn.ensemble import RandomForestRegressor
import joblib
import argparse
import logging
import time

def predict_rf(features, data, model, title):
    X = data[features]
    y = data['PBPotential']
    rf = model['model']
    scaler = model['scaler']
    # Standardize the features
    X_scaled = scaler.transform(X)
    # Predict using the model
    logging.info(f"Predicting with {title}...")
    y_pred = rf.predict(X_scaled)

    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    y_range = np.ptp(y)  # Range of true values
    normalized_rmse = rmse / y_range if y_range != 0 else 0
    r2 = r2_score(y, y_pred)
    logging.info(f"{title} - R2: {r2:.3f}, RMSE: {normalized_rmse:.3f}")
    # get feature importances
    feature_importances = rf.feature_importances_
    feature_names = X.columns
    # Log the feature importances
    for name, importance in zip(feature_names, feature_importances):
        logging.info(f"Feature: {name}, Importance: {importance:.4f}")
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y, y=y_pred)
    plt.xlabel("True PBPotential")
    plt.ylabel("Predicted PBPotential")
    plt.title(f"Prediction Results: {title}")
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'r--', lw=2)
    plt.grid(True)
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    # print the RMSE and R2 score on the plot
    plt.text(0.05, 0.95, f"R2: {r2:.3f}\nRMSE: {normalized_rmse:.3f}", transform=plt.gca().transAxes, fontsize=12,
             verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5, edgecolor='black'))
    # print the feature importances on the plot
    feature_importance_text = "\n".join([f"{name}: {importance:.4f}" for name, importance in zip(feature_names, feature_importances)])
    plt.text(0.05, 0.85, "Feature Importances:\n" + feature_importance_text, transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5, edgecolor='black'))
    plt.tight_layout()
    # replace the space with underscore in title for saving the figure
    title = title.replace(" ", "_")
    plt.savefig(f"{title}_prediction_results.png")


if __name__ == "__main__":
    # Set up logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s", datefmt='%Y-%m-%d %H:%M:%S')

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Predict electrostatic energy using Random Forest.")
    parser.add_argument("input_csv", help="Input CSV file")
    parser.add_argument("--model", default="local_density_with_scaler.pkl", help="Path to the saved model and scaler in pkl format")
    args = parser.parse_args()


    # Load the data
    logging.info(f"Loading data from {args.input_csv} ...")
    data = pd.read_csv(args.input_csv)
    # Prepare features and target variable
    data['DensityNearAverage'] = (data['Density1_Near'] + data['Density2_Near']) / 2
    data['DensityMidAverage'] = (data['Density1_Mid'] + data['Density2_Mid']) / 2
    data['DensityFarAverage'] = (data['Density1_Far'] + data['Density2_Far']) / 2

    # predict using the saved model
    features = ['DensityNearAverage', 'DensityMidAverage', 'DensityFarAverage', 'CoulombPotential']
    fname = args.model     # Load model and scaler from args.model
    title = "Predicted Electrostatic Energy with pre-trained model"
    logging.info("Loading model %s" % fname)
    try:
        model = joblib.load(fname)
    except FileNotFoundError as e:
        logging.error(f"Model not found: {e}")
        exit(1)
    predict_rf(features, data, model, title)


    plt.show()  # Show the plot interactively