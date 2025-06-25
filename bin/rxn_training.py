#!/usr/bin/env python
"""
Train a model to predict electrostatic energy based on embedding scores and distances.

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
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np
from sklearn.ensemble import RandomForestRegressor
import joblib
import argparse
import logging
import time

def fit_rf(features, data, title):
    X = data[features]
    y = data['PBRXN']
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train a Random Forest model
    logging.info(f"Training with {title}...")
    rf = RandomForestRegressor(n_estimators=100, random_state=int(time.time()))
    rf.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = rf.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2 = r2_score(y_val, y_pred_adjusted)
    logging.info(f"R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
    # get feature importances
    feature_importances_adjusted = rf.feature_importances_
    feature_names_adjusted = X.columns
    # Log the feature importances
    for name, importance in zip(feature_names_adjusted, feature_importances_adjusted):
        logging.info(f"Feature: {name}, Importance: {importance:.4f}")
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val, y=y_pred_adjusted, alpha=0.5)
    plt.grid(True)  # add grid lines
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)  # Diagonal line
    plt.xlabel("True PB RXN")
    plt.ylabel("Predicted PB RXN")
    plt.title(f"{title}")
    # print R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2:.3f} Good if > 0.9\nRMSE: {normalized_rmse_adjusted:.3E} Good if <0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    # Sort and add feature importances to the plot
    feature_importance_adjusted_df = pd.DataFrame({'Feature': feature_names_adjusted, 'Importance': feature_importances_adjusted})
    feature_importance_adjusted_df = feature_importance_adjusted_df.sort_values(by='Importance', ascending=False)
    plt.text(0.05, 0.85, "\n".join([f"{row['Feature']}: {row['Importance']:.4f}" for _, row in feature_importance_adjusted_df.iterrows()]), transform=plt.gca().transAxes, fontsize=10, verticalalignment='top')
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    plt.savefig(f"{title}.png")
    # Save the trained model
    logging.info(f"Saving the trained model and scaler for {title} ...")
    model_filename = f"{title.replace(' ', '_').lower()}_with_scaler.pkl"
    feature_names = [f.replace(' ', '_') for f in features]  # Replace spaces with underscores in feature names
    joblib.dump({'model': rf, 'scaler': scaler, 'features': feature_names}, model_filename)
    logging.info(f"Saved the trained model and scaler to {model_filename}.")

if __name__ == "__main__":
    # Set up logging    
    logging_format = "%(asctime)s %(levelname)s: %(message)s"
    logging_datefmt='%Y-%m-%d %H:%M:%S'
    logging.basicConfig(format=logging_format, datefmt=logging_datefmt, level=logging.INFO)    

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Train a model to predict electrostatic energy based on embedding scores and distances.")
    parser.add_argument("input_csv", help="Input CSV file")
    args = parser.parse_args()

    # Load the data
    logging.info(f"Loading data from {args.input_csv} ...")
    data = pd.read_csv(args.input_csv)


    # Train the model
    features = ['Density_Near', 'Density_Mid', 'Density_Far']
    title = "rxn"
    fit_rf(features, data, title)
    logging.info(f"Model trained by {title}.")


    plt.show()
    plt.close('all')
