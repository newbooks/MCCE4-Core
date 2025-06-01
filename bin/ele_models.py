#!/usr/bin/env python
"""
Evaulate ML models
- Randowm Forest Regressor
- Neural Network Regressor
- Polynomial Regression
for fitting Coulomb's potential energy with embedding score to delphi electrostatic energy

This script requires precompiled embedding score in the form of a CSV file. The input is a CSV file with the following columns:
- conf1: Conformer ID for atom 1
- conf2: Conformer ID for atom 2
- distance: Distance between two atoms in Angstroms
- embedding1: Embedding score for atom 1
- embedding2: Embedding score for atom 2
- CoulombPotential: Coulomb potential between two atoms
- PBPotential: Electrostatic energy from Poisson-Boltzmann calculation

The evaluation will cover two aspects:
1. The performance of of each model measured by R^2 score and RMSE
2. The capability of each model trained by one microstate to predict the electrostatic energy of other microstates
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
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import joblib
import time

# Constants
D_in = 4.0    # inner dielectric constant (Coulomb potential)
D_out = 80.0  # outer dielectric constant

def parse_arguments():
    helpmsg = "Evaluate ML models for fitting Coulomb's potential energy with embedding score to delphi electrostatic energy"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("csv_file", nargs='+', help="One or more CSV files containing embedding scores and electrostatic energies")
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s", datefmt='%Y-%m-%d %H:%M:%S')
    logging.info("Starting evaluation of ML models...")

    # Use the first CSV file as the main dataset for training and evaluation
    main_csv_file = args.csv_file[0]
    logging.info(f"Using {main_csv_file} as the main dataset.")

    # Load the main dataset
    df = pd.read_csv(main_csv_file)
    logging.info(f"Loaded {len(df)} rows from {main_csv_file}.")
    # Check if the required columns are present
    required_columns = ['Conf1', 'Conf2', 'Distance', 'Embedding1', 'Embedding2', 'CoulombPotential', 'PBPotential']
    if not all(col in df.columns for col in required_columns):
        logging.error(f"CSV file {main_csv_file} is missing required columns: {required_columns}")
        exit(1)
    # Split the data into features and target variable
    X = df[['Distance', 'Embedding1', 'Embedding2', 'CoulombPotential']]
    y = df['PBPotential']

    # Calculate the average embedding score between the two atoms
    X = X.copy()  # Create a copy to avoid SettingWithCopyWarning
    X["EmbeddingAverage"] = (X['Embedding1'] + X['Embedding2']) / 2.0

    # Adjust the Coulomb potential based on the average embedding and dielectric constants
    # If EmbeddingAverage=0: the pair is totally exposed, use Coulomb potential scaled by D_in/D_out
    # If EmbeddingAverage=1: the pair is totally buried, use Coulomb potential
    # Linearly interpolate between these two extremes
    effective_dielectric = D_out - (D_out - D_in) * X['EmbeddingAverage']
    # Adjust the Coulomb potential
    X['CoulombPotentialAdjusted'] = X['CoulombPotential'] * (D_in / effective_dielectric)

    # Since CoulombPotentialAdjusted is derived from CoulombPotential, EmbeddingAverage, Embedding1, and Embedding2, we will drop other columns
    X_dropped = X.drop(columns=['Embedding1', 'Embedding2', 'EmbeddingAverage', 'CoulombPotential'])

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X_dropped, y, test_size=0.2, random_state=int(time.time()))
    logging.info(f"Split data into {len(X_train)} training and {len(X_test)} testing samples.")
    # Standardize the features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    logging.info("Standardized the features.")

    # plot y_test vs x_train_scaled
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=X_test_scaled[:, 1], y=y_test, alpha=0.6)
    plt.xlabel("Coulomb Potential Adjusted")
    plt.ylabel("PB Potential")
    plt.title("PB Potential vs Coulomb Potential Adjusted")
    plt.grid(True)
    # plt.show()  # Uncomment this line to display the plot interactively
    # save the plot
    png_fname = main_csv_file.rsplit('.', 1)[0] + '_PB_vs_Coulomb.png'  # save as the same name as the main CSV file
    plt.savefig(png_fname)
    logging.info(f"Saved the plot of Linear Regression predictions vs actual values to {png_fname}.")

    # Train a Linear Regression model
    lr_model = LinearRegression()
    lr_model.fit(X_train_scaled, y_train)
    logging.info("Trained Linear Regression model.")
    # Save the trained model
    joblib.dump(lr_model, 'lr_model.pkl')
    logging.info("Saved the trained Linear Regression model to lr_model.pkl.")
    # Make predictions on the test set
    y_pred_lr = lr_model.predict(X_test_scaled)
    # Calculate R^2 score and RMSE for Linear Regression
    r2_lr = r2_score(y_test, y_pred_lr)
    rmse_lr = np.sqrt(mean_squared_error(y_test, y_pred_lr))
    y_range = np.ptp(y_test)  # Range of actual values
    normalized_rmse_lr = rmse_lr / y_range if y_range != 0 else rmse_lr  # Normalize RMSE by range of actual values
    logging.info(f"Linear Regression R^2: {r2_lr:.3f}, RMSE: {rmse_lr:.3f}, Normalized RMSE: {normalized_rmse_lr:.3f}")
    # Plot the predictions vs actual values
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_test, y=y_pred_lr, alpha=0.6)
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'r--', lw=2)
    plt.xlabel("Actual PB Potential")
    plt.ylabel("Predicted PB Potential")
    plt.title("Linear Regression: Actual vs Predicted")
    # print the R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2_lr:.3f} (Good if > 0.8)\nRMSE: {normalized_rmse_lr:.3f} (Good if < 0.05)", transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    plt.grid(True)
#    plt.show()  # Uncomment this line to display the plot interactively
    # save the plot
    png_fname = main_csv_file.rsplit('.', 1)[0] + '_lr_predictions_vs_actual.png'  # save as the same name as the main CSV file
    plt.savefig(png_fname)
    logging.info(f"Saved the plot of Linear Regression predictions vs actual values to {png_fname}.")


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
    plt.grid(True)
#    plt.show()  # Uncomment this line to display the plot interactively
    # save the plot
    png_fname = main_csv_file.rsplit('.', 1)[0] + '_rf_predictions_vs_actual.png'  # save as the same name as the main CSV file
    plt.savefig(png_fname)
    logging.info(f"Saved the plot of Random Forest predictions vs actual values to {png_fname}.")


    # Train a Neural Network Regressor
    nn_model = MLPRegressor(hidden_layer_sizes=(100, 50), max_iter=1000, random_state=int(time.time()))
    nn_model.fit(X_train_scaled, y_train)
    logging.info("Trained Neural Network Regressor.")
    # Save the trained model
    joblib.dump(nn_model, 'nn_model.pkl')
    logging.info("Saved the trained Neural Network model to nn_model.pkl.")
    # Make predictions on the test set
    y_pred_nn = nn_model.predict(X_test_scaled)
    # Calculate R^2 score and RMSE for Neural Network
    r2_nn = r2_score(y_test, y_pred_nn)
    rmse_nn = np.sqrt(mean_squared_error(y_test, y_pred_nn))
    y_range = np.ptp(y_test)  # Range of actual values
    normalized_rmse_nn = rmse_nn / y_range if y_range != 0 else rmse_nn  # Normalize RMSE by range of actual values
    logging.info(f"Neural Network Regressor R^2: {r2_nn:.3f}, RMSE: {rmse_nn:.3f}, Normalized RMSE: {normalized_rmse_nn:.3f}")
    # Plot the predictions vs actual values
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_test, y=y_pred_nn, alpha=0.6)
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'r--', lw=2)
    plt.xlabel("Actual PB Potential")
    plt.ylabel("Predicted PB Potential")
    plt.title("Neural Network Regressor: Actual vs Predicted")
    # print the R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2_nn:.3f} (Good if > 0.8)\nRMSE: {normalized_rmse_nn:.3f} (Good if < 0.05)", transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    plt.grid(True)
#    plt.show()  # Uncomment this line to display the plot interactively
    # save the plot
    png_fname = main_csv_file.rsplit('.', 1)[0] + '_nn_predictions_vs_actual.png'  # save as the same name as the main CSV file
    plt.savefig(png_fname)
    logging.info(f"Saved the plot of Neural Network predictions vs actual values to {png_fname}.")

    # Train a Polynomial Regression model
    poly_features = PolynomialFeatures(degree=4, include_bias=True)
    X_poly_train = poly_features.fit_transform(X_train_scaled)
    X_poly_test = poly_features.transform(X_test_scaled)
    poly_model = LinearRegression()
    
    poly_model.fit(X_poly_train, y_train)
    logging.info("Trained Polynomial Regression model.")
    # Save the trained model
    joblib.dump(poly_model, 'poly_model.pkl')
    logging.info("Saved the trained Polynomial Regression model to poly_model.pkl.")
    # Make predictions on the test set
    y_pred_poly = poly_model.predict(X_poly_test)
    # Calculate R^2 score and RMSE for Polynomial Regression
    r2_poly = r2_score(y_test, y_pred_poly)
    rmse_poly = np.sqrt(mean_squared_error(y_test, y_pred_poly))
    y_range = np.ptp(y_test)  # Range of actual values
    normalized_rmse_poly = rmse_poly / y_range if y_range != 0 else rmse_poly  # Normalize RMSE by range of actual values
    logging.info(f"Polynomial Regression R^2: {r2_poly:.3f}, RMSE: {rmse_poly:.3f}, Normalized RMSE: {normalized_rmse_poly:.3f}")
    # Plot the predictions vs actual values
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_test, y=y_pred_poly, alpha=0.6)
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'r--', lw=2)
    plt.xlabel("Actual PB Potential")
    plt.ylabel("Predicted PB Potential")
    plt.title("Polynomial Regression: Actual vs Predicted")
    # print the R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2_poly:.3f} (Good if > 0.8)\nRMSE: {normalized_rmse_poly:.3f} (Good if < 0.05)", transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    plt.grid(True)
    # plt.show()  # Uncomment this line to display the plot interactively
    # save the plot
    png_fname = main_csv_file.rsplit('.', 1)[0] + '_poly_predictions_vs_actual.png'  # save as the same name as the main CSV file
    plt.savefig(png_fname)
    logging.info(f"Saved the plot of Polynomial Regression predictions vs actual values to {png_fname}.")

    # If multiple files are provided, use the model trained on the first file to predict the electrostatic energy of the other files
    if len(args.csv_file) > 1:
        for csv_file in args.csv_file[1:]:
            logging.info(f"Processing {csv_file}...")
            df_other = pd.read_csv(csv_file)
            # Check if the required columns are present
            if not all(col in df_other.columns for col in required_columns):
                logging.error(f"CSV file {csv_file} is missing required columns: {required_columns}")
                continue
            # Split the data into features and target variable
            X_other = df_other[['Distance', 'Embedding1', 'Embedding2', 'CoulombPotential']]
            y_other = df_other['PBPotential']
            # Calculate the average embedding score between the two atoms
            X_other = X_other.copy()
            X_other["EmbeddingAverage"] = (X_other['Embedding1'] + X_other['Embedding2']) / 2.0
            # Adjust the Coulomb potential based on the average embedding and dielectric constants
            effective_dielectric = D_out - (D_out - D_in) * X_other['EmbeddingAverage']
            # Adjust the Coulomb potential
            X_other['CoulombPotentialAdjusted'] = X_other['CoulombPotential'] * (D_in / effective_dielectric)
            # Since CoulombPotentialAdjusted is derived from CoulombPotential, EmbeddingAverage, Embedding1, and Embedding2, we will drop other columns
            X_other_dropped = X_other.drop(columns=['Embedding1', 'Embedding2', 'EmbeddingAverage', 'CoulombPotential'])
            # Standardize the features
            X_other_scaled = scaler.transform(X_other_dropped)
            # Make predictions using the trained models using the best model random forest
            y_pred_rf_other = rf_model.predict(X_other_scaled)
            # Calculate R^2 score and RMSE for each model
            y_range = np.ptp(y_other)
            r2_rf_other = r2_score(y_other, y_pred_rf_other)
            rmse_rf_other = np.sqrt(mean_squared_error(y_other, y_pred_rf_other))
            normalized_rmse_rf_other = rmse_rf_other / y_range if y_range != 0 else rmse_rf_other
            logging.info(f"Random Forest Regressor R^2: {r2_rf_other:.3f}, RMSE: {rmse_rf_other:.3f}, Normalized RMSE: {normalized_rmse_rf_other:.3f}")

            # Plot the predictions vs actual values
            plt.figure(figsize=(10, 6))
            sns.scatterplot(x=y_other, y=y_pred_rf_other, alpha=0.6)
            plt.plot([y_other.min(), y_other.max()], [y_other.min(), y_other.max()], 'r--', lw=2)
            plt.xlabel("Actual PB Potential")
            plt.ylabel("Predicted PB Potential (Random Forest)")
            plt.title(f"Random Forest Regressor: Actual vs Predicted ({csv_file})")
            plt.text(0.05, 0.95, f"R^2: {r2_rf_other:.3f} (Good if > 0.8)\nRMSE: {normalized_rmse_rf_other:.3f} (Good if < 0.05)", 
                     transform=plt.gca().transAxes, fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
            plt.grid(True)
            # save the plot
            png_fname = csv_file.rsplit('.', 1)[0] + '_rf_predictions_vs_actual.png'
            plt.savefig(png_fname)
            logging.info(f"Saved the plot of Random Forest predictions vs actual values to {png_fname}.")

    plt.show()
