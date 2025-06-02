#!/usr/bin/env python
"""
Use compiled electrostatic energy to fit a model
to predict electrostatic energy based on embedding scores and distances.

The input is a CSV file with the following columns:
- conf1: Conformer ID for atom 1
- conf2: Conformer ID for atom 2
- distance: Distance between two atoms in Angstroms
- embedding1: Embedding score for atom 1
- embedding2: Embedding score for atom 2
- CoulombPotential: Coulomb potential between two atoms
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


if __name__ == "__main__":
    # Set up logging    
    logging_format = "%(asctime)s %(levelname)s: %(message)s"
    logging_datefmt='%Y-%m-%d %H:%M:%S'
    logging.basicConfig(format=logging_format, datefmt=logging_datefmt, level=logging.INFO)    

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Fit a model to predict electrostatic energy based on embedding scores and distances.")
    parser.add_argument("input_csv", help="Input CSV file with columns: conf1, conf2, distance, embedding1, embedding2, CoulombPotential, PBPotential")
    parser.add_argument("--output_model", default="electrostatic_model.pkl", help="Output model file name")
    args = parser.parse_args()

    # Load the data
    logging.info(f"Loading data from {args.input_csv} ...")
    data = pd.read_csv(args.input_csv)

    # Check if required columns are present
    required_columns = ['Conf1', 'Conf2', 'Distance', 'Embedding1', 'Embedding2', 'CoulombPotential', 'PBPotential']
    for col in required_columns:
        if col not in data.columns:
            logging.error(f"Missing required column: {col}")
            exit(1)

    # Prepare features and target variable
    X = data[['Distance', 'Embedding1', 'Embedding2', 'CoulombPotential']]
    y = data['PBPotential']

    logging.info("Features and target variable prepared.")
    # print(X.head())  # Print the first few rows of the features for debuggiSng
    # print the target name
    # print("\n  PBPotential")
    # print(y.head())  # Print the first few rows of the target variable for debugging

    # Add columns that are likely to contribute to the model
    X = X.copy()
    X['EmbeddingAverage'] = (X['Embedding1'] + X['Embedding2']) / 2
    # PBPotential = CoulombPotential * EmbeddingAverage + CoulombPotential * (1 - EmbeddingAverage) * D_in/D_out 
    X['CoulombPotentialAdjusted'] = X['CoulombPotential'] * (D_in / D_out) * (1 - X['EmbeddingAverage']) + X['CoulombPotential'] * X['EmbeddingAverage']


    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))

    # Standardize the features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Fit a linear regression model
    logging.info("Fitting linear regression model ...")
    linear_model = LinearRegression()
    linear_model.fit(X_train_scaled, y_train)

    # Predict on the test set
    y_pred_linear = linear_model.predict(X_test_scaled)

    # Evaluate the linear model
    mse_linear = mean_squared_error(y_test, y_pred_linear)
    r2_linear = r2_score(y_test, y_pred_linear)
    
    logging.info(f"Linear Regression Model - MSE: {mse_linear}, R^2: {r2_linear}")

    # Print MSA and R^2 with reference values (good or bad)
    print(f"Linear Regression Model - MSE: {mse_linear:.4f} (Good if < 0.1)")
    print(f"Linear Regression Model - R^2: {r2_linear:.4f} (Good if > 0.8)")
    print(f"Linear Regression Model - Coefficients: {linear_model.coef_}")
    print(f"Linear Regression Model - Intercept: {linear_model.intercept_}")          
    

    # plot the results
    plt.figure(figsize=(10, 6))
    plt.scatter(y_test, y_pred_linear, alpha=0.5)
    plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=2)
    plt.xlabel('True PBPotential')
    plt.ylabel('Predicted PBPotential')
    plt.title('Linear Regression Model: True vs Predicted PBPotential')
    plt.grid()
    # plt.savefig('linear_regression_results.png')
    plt.show()


    # Fit a Random Forest Regressor
    logging.info("Fitting Random Forest Regressor ...")
    rf_model = RandomForestRegressor(n_estimators=100, random_state=int(time.time()))
    rf_model.fit(X_train_scaled, y_train)
    # Predict on the test set
    y_pred_rf = rf_model.predict(X_test_scaled)
    # Evaluate the Random Forest model
    mse_rf = mean_squared_error(y_test, y_pred_rf)
    r2_rf = r2_score(y_test, y_pred_rf)
    logging.info(f"Random Forest Model - MSE: {mse_rf}, R^2: {r2_rf}")
    # Print MSA and R^2 with reference values (good or bad)
    print(f"Random Forest Model - MSE: {mse_rf:.4f} (Good if < 0.1)")
    print(f"Random Forest Model - R^2: {r2_rf:.4f} (Good if > 0.8)")
    print(f"Random Forest Model - Feature Importances: {rf_model.feature_importances_}")
    # plot the results
    plt.figure(figsize=(10, 6))
    plt.scatter(y_test, y_pred_rf, alpha=0.5)
    plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=2)
    plt.xlabel('True PBPotential')
    plt.ylabel('Predicted PBPotential')
    plt.title('Random Forest Model: True vs Predicted PBPotential')
    plt.grid()
    # plt.savefig('random_forest_results.png')
    plt.show()

    # Fit an XGBoost Regressor
    logging.info("Fitting XGBoost Regressor ...")
    xgb_model = XGBRegressor(n_estimators=100, random_state=int(time.time()), verbosity=0)
    xgb_model.fit(X_train_scaled, y_train)
    # Predict on the test set
    y_pred_xgb = xgb_model.predict(X_test_scaled)
    # Evaluate the XGBoost model
    mse_xgb = mean_squared_error(y_test, y_pred_xgb)
    r2_xgb = r2_score(y_test, y_pred_xgb)
    logging.info(f"XGBoost Model - MSE: {mse_xgb}, R^2: {r2_xgb}")
    # Print MSA and R^2 with reference values (good or bad)
    print(f"XGBoost Model - MSE: {mse_xgb:.4f} (Good if < 0.1)")
    print(f"XGBoost Model - R^2: {r2_xgb:.4f} (Good if > 0.8)")
    print(f"XGBoost Model - Feature Importances: {xgb_model.feature_importances_}")
    # plot the results
    plt.figure(figsize=(10, 6))
    plt.scatter(y_test, y_pred_xgb, alpha=0.5)
    plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=2)
    plt.xlabel('True PBPotential')
    plt.ylabel('Predicted PBPotential')
    plt.title('XGBoost Model: True vs Predicted PBPotential')
    plt.grid()
    # plt.savefig('xgboost_results.png')
    plt.show()
    # print the feature importances
    importances = xgb_model.feature_importances_
    feature_names = X.columns
    feature_importances = pd.DataFrame({'Feature': feature_names, 'Importance': importances})
    feature_importances = feature_importances.sort_values(by='Importance', ascending=False)
    logging.info("Feature Importances:")
    logging.info(feature_importances)


    # Plot PBPotential vs CoulombPotentialAdjusted
    plt.figure(figsize=(10, 6))
    plt.scatter(X['CoulombPotentialAdjusted'], y, alpha=0.5)
    plt.xlabel('CoulombPotentialAdjusted')
    plt.ylabel('PBPotential')
    plt.title('PBPotential vs CoulombPotentialAdjusted')
    plt.grid()
    # plt.savefig('pbpotential_vs_coulombpotentialadjusted.png')
    plt.show()