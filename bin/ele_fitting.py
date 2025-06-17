#!/usr/bin/env python
"""
Use compiled electrostatic energy to fit a model to predict electrostatic energy based on embedding scores and distances.

The input is a CSV file with the following columns:
Columns in CSV file:
- Conf1
- Conf2
- Distance
- Radius1
- Radius2
- Density1_Near
- Density2_Near
- Density1_Mid
- Density2_Mid
- Density1_Near
- Density2_Near
- CoulombPotential
- PBPotential
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
import argparse
import logging
import time


def fit_rf(features, data, title):
    X = data[features]
    y = data['PBPotential']
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
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'r--', lw=2)  # Diagonal line
    plt.xlabel("True PB Potential")
    plt.ylabel("Predicted PB Potential")
    plt.title(f"{title}")
    # print R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2_adjusted:.3f} Good if > 0.9\nRMSE: {normalized_rmse_adjusted:.3f} Good if <0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    # Sort and add feature importances to the plot
    feature_importance_adjusted_df = pd.DataFrame({'Feature': feature_names_adjusted, 'Importance': feature_importances_adjusted})
    feature_importance_adjusted_df = feature_importance_adjusted_df.sort_values(by='Importance', ascending=False)
    plt.text(0.05, 0.85, "\n".join([f"{row['Feature']}: {row['Importance']:.4f}" for _, row in feature_importance_adjusted_df.iterrows()]), transform=plt.gca().transAxes, fontsize=10, verticalalignment='top')
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
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
    data['DensityNearAverage'] = (data['Density1_Near'] + data['Density2_Near']) / 2
    data['DensityMidAverage'] = (data['Density1_Mid'] + data['Density2_Mid']) / 2
    data['DensityFarAverage'] = (data['Density1_Far'] + data['Density2_Far']) / 2

    # Evaluate all features as is
    features = ['Distance', 'Radius1', 'Radius2', 'Density1_Near', 'Density2_Near', 'Density1_Mid', 'Density2_Mid', 'Density1_Far', 'Density2_Far', 'CoulombPotential']
    title = "All Native Features"
    fit_rf(features, data, title)

    # Reduce some features by averaging the symmetric ones and removing the distance
    logging.info("Reducing features by averaging symmetric ones...")
    title = "Independent features"
    features = ['DensityNearAverage', 'DensityMidAverage', 'DensityFarAverage', 'CoulombPotential']
    fit_rf(features, data, title)

    # Use the density mid average
    logging.info("Using density mid average...")
    features = ['DensityMidAverage', 'CoulombPotential']
    title = "Mid Range Density"
    fit_rf(features, data, title)

    # Use the density far average
    logging.info("Using density far average...")
    features = ['DensityFarAverage', 'CoulombPotential']
    title = "Density Far Average"
    fit_rf(features, data, title)


    plt.show()
    plt.close()  # Close the plot to free memory

