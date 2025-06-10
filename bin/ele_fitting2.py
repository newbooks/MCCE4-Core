#!/usr/bin/env python
"""
Use compiled electrostatic energy to fit a model to predict electrostatic energy based on embedding scores and distances.

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


# import the modules
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
import joblib
import argparse
import logging
import time

D_in = 4.0    # inner dielectric constant (Coulomb potential)
D_out = 80.0  # outter dielectric constant

def fit_rf(X, y, title):
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train a Random Forest model
    logging.info(f"Training with {title}...")
    rf_adjusted = RandomForestRegressor(n_estimators=100, random_state=int(time.time()))
    rf_adjusted.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = rf_adjusted.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2_adjusted = r2_score(y_val, y_pred_adjusted)
    logging.info(f"Adjusted Coulomb Potential - R2: {r2_adjusted:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
    # get feature importances
    feature_importances_adjusted = rf_adjusted.feature_importances_
    feature_names_adjusted = X.columns
    # Log the feature importances
    for name, importance in zip(feature_names_adjusted, feature_importances_adjusted):
        logging.info(f"Feature: {name}, Importance: {importance:.4f}")
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val, y=y_pred_adjusted, alpha=0.5)
    plt.grid(True)  # add grid lines
    plt.plot([y_val.min(), y_val.max()], [y_val.min(), y_val.max()], 'r--', lw=2)  # Diagonal line
    plt.xlabel("True PB Potential")
    plt.ylabel("Predicted PB Potential")
    plt.title(f"{title}")
    # print R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2_adjusted:.3f} Good if > 0.9\nRMSE: {normalized_rmse_adjusted:.3f} Good if <0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    # Sort and add feature importances to the plot
    feature_importance_adjusted_df = pd.DataFrame({'Feature': feature_names_adjusted, 'Importance': feature_importances_adjusted})
    feature_importance_adjusted_df = feature_importance_adjusted_df.sort_values(by='Importance', ascending=False)
    plt.text(0.05, 0.85, "\n".join([f"{row['Feature']}: {row['Importance']:.4f}" for _, row in feature_importance_adjusted_df.iterrows()]), transform=plt.gca().transAxes, fontsize=10, verticalalignment='top')
    plt.xlim(y_val.min(), y_val.max())
    plt.ylim(y_val.min(), y_val.max())
    plt.savefig(f"{title}.png")


if __name__ == "__main__":
    # Set up logging    
    logging_format = "%(asctime)s %(levelname)s: %(message)s"
    logging_datefmt='%Y-%m-%d %H:%M:%S'
    logging.basicConfig(format=logging_format, datefmt=logging_datefmt, level=logging.INFO)    

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Fit a model to predict electrostatic energy based on embedding scores and distances.")
    parser.add_argument("input_csv", help="Input CSV file with columns: conf1, conf2, distance, embedding1, embedding2, CoulombPotential, PBPotential")
    parser.add_argument("--output_model", default="", help="Output model file name")
    args = parser.parse_args()

    # Load the data
    logging.info(f"Loading data from {args.input_csv} ...")
    data = pd.read_csv(args.input_csv)
    # Check if required columns are present
    required_columns = ['Conf1', 'Conf2', 'Distance', 'Radius1', 'Radius2', 'Embedding1', 'Embedding2', 'Density1', 'Density2', 'CoulombPotential', 'AdjustedCoulombPotential', 'PBPotential']
    for col in required_columns:
        if col not in data.columns:
            logging.error(f"Missing required column: {col}")
            exit(1)
    
    # Evaluate all features as is
    X = data[['Distance', 'Radius1', 'Radius2', 'Embedding1', 'Embedding2', 'Density1', 'Density2', 'CoulombPotential']]
    y = data['PBPotential']
    title = "All_Features"
    logging.info("Features and target variable prepared.")
    fit_rf(X, y, title)

    # Use Distance and Adjusted Coulomb Potential
    X = data[['Distance', 'AdjustedCoulombPotential']]
    y = data['PBPotential']
    title = "Distance+AdjustedCoulomb"
    logging.info("Features and target variable prepared.")
    fit_rf(X, y, title)

    # Add density scores
    X = data[['Distance', 'Density1', 'Density2', 'AdjustedCoulombPotential']]
    y = data['PBPotential']
    title = "Distance+Densities+Coulomb"
    logging.info("Features and target variable prepared.")
    fit_rf(X, y, title)

    # Use average density scores
    X = data[['Distance', 'Density1', 'Density2', 'AdjustedCoulombPotential']]
    X = X.copy()
    X['DensityAverage'] = (X['Density1'] + X['Density2']) / 2
    y = data['PBPotential']
    X.drop(columns=['Density1', 'Density2'], inplace=True)  # Remove individual densities
    title = "Distance+DensitiesAverage+Coulomb"
    fit_rf(X, y, title)
    logging.info("Features and target variable prepared.")

    plt.show()
    plt.close()  # Close the plot to free memory

