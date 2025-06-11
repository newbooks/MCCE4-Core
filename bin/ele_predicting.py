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

def predict_rf(X, y, model, title):
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
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'k--', lw=2)
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
    parser.add_argument("input_csv", help="Input CSV file with columns: Conf1, Conf2, Distance, Radius1, Radius2, Embedding1, Embedding2, Density1, Density2, CoulombPotential, AdjustedCoulombPotential, PBPotential")
    args = parser.parse_args()


    # Load the data
    logging.info(f"Loading data from {args.input_csv} ...")
    data = pd.read_csv(args.input_csv)


    # predict using the model without local density
    X = data[['Distance', 'AdjustedCoulombPotential']]
    y = data['PBPotential']
    # Load model and scaler from without_local_density
    logging.info("Loading model and scaler from without_local_density...")
    try:
        model = joblib.load("without_local_density_model_with_scaler.pkl")
    except FileNotFoundError as e:
        logging.error(f"Error loading model or scaler: {e}")
        exit(1)
    predict_rf(X, y, model, "Without Local Density")


    # predict using the model with local density
    X = data[['Distance', 'Density1', 'Density2', 'AdjustedCoulombPotential']]
    X = X.copy()
    X['DensityAverage'] = (X['Density1'] + X['Density2']) / 2
    X.drop(columns=['Density1', 'Density2'], inplace=True)  # Drop individual density columns
    y = data['PBPotential']
    # Load model and scaler from without_local_density
    logging.info("Loading model and scaler from with_local_density...")
    try:
        model = joblib.load("with_local_density_model_with_scaler.pkl")
    except FileNotFoundError as e:
        logging.error(f"Error loading model or scaler: {e}")
        exit(1)
    predict_rf(X, y, model, "With Local Density")


    plt.show()  # Show the plot interactively