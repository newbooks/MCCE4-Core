#!/usr/bin/env python
"""
Use Random Forest to fit a model, using 
1. the original Coulomb potential 
2. distance
3. radius 
4. embedding score
to predict electrostatic energy.
"""

import logging
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np
from sklearn.ensemble import RandomForestRegressor
import joblib
import time


def parse_arguments():
    helpmsg = "Use Random Forest to fit a model, using the original Coulomb potential, distance, radius, and embedding score to predict electrostatic energy.\n"
    helpmsg += "The input is a CSV file with the following columns:\n"
    helpmsg += "- Conf1: Conformer ID for atom 1\n"
    helpmsg += "- Conf2: Conformer ID for atom 2\n"
    helpmsg += "- Distance: Distance between two atoms in Angstroms\n"
    helpmsg += "- Radius1: Radius for atom 1\n"
    helpmsg += "- Radius2: Radius for atom 2\n"
    helpmsg += "- Embedding1: Embedding score for atom 1\n"
    helpmsg += "- Embedding2: Embedding score for atom 2\n"
    helpmsg += "- CoulombPotential: Coulomb potential between two atoms\n"
    helpmsg += "- PBPotential: Electrostatic energy from Poisson-Boltzmann calculation\n"
    parser = argparse.ArgumentParser(description=helpmsg, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("csv_file", help="One CSV file containing embedding scores and electrostatic energies")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s", datefmt='%Y-%m-%d %H:%M:%S')
    logging.info("Starting evaluation of ML models...")

    # Use the first CSV file as the main dataset for training and evaluation
    main_csv_file = args.csv_file
    logging.info(f"Using {main_csv_file} as the main dataset.")

    # Load the main dataset
    df = pd.read_csv(main_csv_file)
    logging.info(f"Loaded {len(df)} rows from {main_csv_file}.")
    # Check if the required columns are present
    required_columns = ['Conf1', 'Conf2', 'Distance', 'Radius1', 'Radius2', 'Embedding1', 'Embedding2', 'CoulombPotential', 'PBPotential']
    if not all(col in df.columns for col in required_columns):
        logging.error(f"CSV file {main_csv_file} is missing required columns: {required_columns}")
        exit(1)
    # Split the data into features and target variable
    X = df[['Distance', 'Radius1', 'Radius2', 'Embedding1', 'Embedding2', 'CoulombPotential']]
    y = df['PBPotential']

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    logging.info(f"Split data into {len(X_train)} training and {len(X_test)} testing samples.")
    # Standardize the features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    logging.info("Standardized the features.")

    # Train a Random Forest Regressor
    rf_model = RandomForestRegressor(n_estimators=100, random_state=int(time.time()))
    rf_model.fit(X_train_scaled, y_train)
    logging.info("Trained Random Forest Regressor.")
    # Save the trained model
    joblib.dump(rf_model, 'rf_model.pkl')
    logging.info("Saved the trained Random Forest model to rf_model.pkl.")
    # Make predictions on the test set
    y_pred_rf = rf_model.predict(X_test_scaled)
    # Calculate R^2 score and RMSE for Random Forest
    r2_rf = r2_score(y_test, y_pred_rf)
    rmse_rf = np.sqrt(mean_squared_error(y_test, y_pred_rf))
    y_range = np.ptp(y_test)  # Range of actual values
    normalized_rmse_rf = rmse_rf / y_range if y_range != 0 else rmse_rf  # Normalize RMSE by range of actual values
    logging.info(f"Random Forest Regressor R^2: {r2_rf:.3f}, RMSE: {rmse_rf:.3f}, Normalized RMSE: {normalized_rmse_rf:.3f}")
    # Get the feature importances
    feature_importances = rf_model.feature_importances_
    feature_names = X.columns
    # Create a DataFrame for feature importances
    feature_importance_df = pd.DataFrame({'Feature': feature_names, 'Importance': feature_importances})
    # Sort the DataFrame by importance
    feature_importance_df = feature_importance_df.sort_values(by='Importance', ascending=False)
    # Print the feature importances
    logging.info("Feature Importances from Random Forest Regressor:")
    logging.info(feature_importance_df)
    # Plot the predictions vs actual values
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_test, y=y_pred_rf, alpha=0.6)
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'r--', lw=2)
    plt.xlabel("Actual PB Potential")
    plt.ylabel("Predicted PB Potential")
    plt.title("Random Forest Regressor: Actual vs Predicted")
    # print the R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2_rf:.3f} (Good if > 0.8)\nRMSE: {normalized_rmse_rf:.3f} (Good if < 0.05)", transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    # print the feature importances on the plot
    plt.text(0.05, 0.85, "Feature Importances:\n" + "\n".join([f"{row['Feature']}: {row['Importance']:.3f}" for _, row in feature_importance_df.iterrows()]), transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    plt.grid(True)
#    plt.show()  # Uncomment this line to display the plot interactively
    # save the plot
    png_fname = main_csv_file.rsplit('.', 1)[0] + '_rf_predictions_vs_actual.png'  # save as the same name as the main CSV file
    plt.savefig(png_fname)
    logging.info(f"Saved the plot of Random Forest predictions vs actual values to {png_fname}.")

    plt.show()
