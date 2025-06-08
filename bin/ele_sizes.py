#!/usr/bin/env python
"""
Use Random Forest to fit a model. Compare these strategies:
1. Use medium protein as training, predict small, medium and large proteins
2. Use combined training set (small + medium + large) to predict all proteins
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

training_csv_file_medium = "medium/state_0001_compiled.csv"  # Main CSV file
training_csv_file_mix = "mix/state_0001_compiled.csv"  # Main CSV file
prediction_csv_files = [
    "small/state_0401_compiled.csv",  # Small protein dataset
    "medium/state_0401_compiled.csv",  # Medium protein dataset
    "large/state_0401_compiled.csv"    # Large protein dataset
]

def parse_arguments():
    helpmsg = "Use Random Forest to fit a model, then predict on other datasets.\n"
    parser = argparse.ArgumentParser(description=helpmsg, formatter_class=argparse.RawTextHelpFormatter)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s", datefmt='%Y-%m-%d %H:%M:%S')

    # Strategy 1: Train on medium protein, predict on small, medium, and large proteins
    logging.info("Loading training data from medium protein dataset...")
    df_medium = pd.read_csv(training_csv_file_medium)
    X_medium = df_medium[['Distance', 'AdjustedCoulombPotential']]
    y_medium = df_medium['PBPotential']
    # Split the medium dataset into training and validation sets
    X_train_medium, X_val_medium, y_train_medium, y_val_medium = train_test_split(X_medium, y_medium, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train_medium = scaler.fit_transform(X_train_medium)
    X_val_medium = scaler.transform(X_val_medium)

    logging.info("Training Random Forest model on medium protein dataset...")
    rf_medium = RandomForestRegressor(n_estimators=100, random_state=int(time.time()))
    rf_medium.fit(X_train_medium, y_train_medium)
    logging.info("Evaluating model on medium protein validation set...")
    y_pred_medium = rf_medium.predict(X_val_medium)
    rmse_medium = np.sqrt(mean_squared_error(y_val_medium, y_pred_medium))
    y_range = np.ptp(y_val_medium)  # Range of true values
    normalized_rmse_medium = rmse_medium / y_range if y_range != 0 else 0
    r2_medium = r2_score(y_val_medium, y_pred_medium)
    logging.info(f"Medium Protein - RMSE: {normalized_rmse_medium:.3f}, R2: {r2_medium:.3f}")
    # Plot the training results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val_medium, y=y_pred_medium)
    plt.xlabel("True PBPotential")
    plt.ylabel("Predicted PBPotential")
    plt.title("Training by Medium Protein: True vs Predicted PBPotential on training set")
    plt.grid(True)
    plt.xlim(y_val_medium.min(), y_val_medium.max())
    plt.ylim(y_val_medium.min(), y_val_medium.max())
    # Add a diagonal line for perfect predictions
    plt.plot([y_val_medium.min(), y_val_medium.max()], [y_val_medium.min(), y_val_medium.max()], 'k--', lw=2)
    # print the R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2_medium:.3f} Good if > 0.9\nRMSE: {normalized_rmse_medium:.3f} Good if < 0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    plt.savefig("medium_training_results.png")

    # Use the trained model to predict on small, medium, and large protein datasets
    for prediction_csv_file in prediction_csv_files:
        logging.info(f"Loading prediction data from {prediction_csv_file}...")
        df_pred = pd.read_csv(prediction_csv_file)
        X_pred = df_pred[['Distance', 'AdjustedCoulombPotential']]
        X_pred = scaler.transform(X_pred)
        y_pred = rf_medium.predict(X_pred)
        # Calculate RMSE and R^2 for predictions
        rmse_pred = np.sqrt(mean_squared_error(df_pred['PBPotential'], y_pred))
        y_range_pred = np.ptp(df_pred['PBPotential'])  # Range of true values
        normalized_rmse_pred = rmse_pred / y_range_pred if y_range_pred != 0 else 0
        r2_pred = r2_score(df_pred['PBPotential'], y_pred)
        logging.info(f"Predictions on {prediction_csv_file} - RMSE: {normalized_rmse_pred:.3f}, R2: {r2_pred:.3f}, Normalized RMSE: {normalized_rmse_pred:.3f}")
        # plot the predictions
        plt.figure(figsize=(10, 6))
        sns.scatterplot(x=df_pred['PBPotential'], y=y_pred, alpha=0.6)
        plt.xlim(df_pred['PBPotential'].min(), df_pred['PBPotential'].max())
        plt.ylim(df_pred['PBPotential'].min(), df_pred['PBPotential'].max())
        plt.grid(True)
        plt.xlabel("True PBPotential")
        plt.ylabel("Predicted PBPotential")
        plt.title(f"Prediction on {prediction_csv_file}: True vs Predicted PBPotential")
        plt.plot([df_pred['PBPotential'].min(), df_pred['PBPotential'].max()], [df_pred['PBPotential'].min(), df_pred['PBPotential'].max()], 'k--', lw=2)
        # print the R^2 and RMSE on the plot
        plt.text(0.05, 0.95, f"R^2: {r2_pred:.3f} Good if > 0.9\nRMSE: {normalized_rmse_pred:.3f} Good if < 0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
        plt.savefig(f"{prediction_csv_file.replace('_compiled.csv', '_predicted_medium.png')}")

    # show the predictions
    plt.show()
    # close the plots
    plt.close('all')

    # Strategy 2: Train on combined dataset (small + medium + large), predict on all datasets
    logging.info("Loading training data from combined small, medium, and large protein datasets...")
    df_mix = pd.read_csv(training_csv_file_mix)
    X_mix = df_mix[['Distance', 'AdjustedCoulombPotential']]
    y_mix = df_mix['PBPotential']
    # Split the combined dataset into training and validation sets
    X_train_mix, X_val_mix, y_train_mix, y_val_mix = train_test_split(X_mix, y_mix, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler_mix = StandardScaler()
    X_train_mix = scaler_mix.fit_transform(X_train_mix)
    X_val_mix = scaler_mix.transform(X_val_mix)
    logging.info("Training Random Forest model on combined dataset...")
    rf_mix = RandomForestRegressor(n_estimators=100, random_state=int(time.time()))
    rf_mix.fit(X_train_mix, y_train_mix)
    logging.info("Evaluating model on combined validation set...")
    y_pred_mix = rf_mix.predict(X_val_mix)
    rmse_mix = np.sqrt(mean_squared_error(y_val_mix, y_pred_mix))
    y_range_mix = np.ptp(y_val_mix)  # Range of true values
    normalized_rmse_mix = rmse_mix / y_range_mix if y_range_mix != 0 else 0
    r2_mix = r2_score(y_val_mix, y_pred_mix)
    logging.info(f"Combined Dataset - RMSE: {normalized_rmse_mix:.3f}, R2: {r2_mix:.3f}")
    # Plot the training results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val_mix, y=y_pred_mix)
    plt.xlim(y_val_mix.min(), y_val_mix.max())
    plt.ylim(y_val_mix.min(), y_val_mix.max())
    plt.grid(True)
    plt.xlabel("True PBPotential")
    plt.ylabel("Predicted PBPotential")
    plt.title("Training by Combined Dataset: True vs Predicted PBPotential on training set")
    plt.plot([y_val_mix.min(), y_val_mix.max()], [y_val_mix.min(), y_val_mix.max()], 'k--', lw=2)
    # print the R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2_mix:.3f} Good if > 0.9\nRMSE: {normalized_rmse_mix:.3f} Good if < 0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    plt.savefig("combined_training_results.png")
    # Use the trained model to predict on small, medium, and large protein datasets
    for prediction_csv_file in prediction_csv_files:
        logging.info(f"Loading prediction data from {prediction_csv_file}...")
        df_pred = pd.read_csv(prediction_csv_file)
        X_pred = df_pred[['Distance', 'AdjustedCoulombPotential']]
        X_pred = scaler_mix.transform(X_pred)
        y_pred = rf_mix.predict(X_pred)
        # Calculate RMSE and R^2 for predictions
        rmse_pred = np.sqrt(mean_squared_error(df_pred['PBPotential'], y_pred))
        y_range_pred = np.ptp(df_pred['PBPotential'])  # Range of true values
        normalized_rmse_pred = rmse_pred / y_range_pred if y_range_pred != 0 else 0
        r2_pred = r2_score(df_pred['PBPotential'], y_pred)
        logging.info(f"Predictions on {prediction_csv_file} - RMSE: {normalized_rmse_pred:.3f}, R2: {r2_pred:.3f}")
        # plot the predictions
        plt.figure(figsize=(10, 6))
        sns.scatterplot(x=df_pred['PBPotential'], y=y_pred, alpha=0.6)
        plt.xlim(df_pred['PBPotential'].min(), df_pred['PBPotential'].max())
        plt.ylim(df_pred['PBPotential'].min(), df_pred['PBPotential'].max())
        plt.grid(True)
        plt.xlabel("True PBPotential")
        plt.ylabel("Predicted PBPotential")
        plt.title(f"Prediction on {prediction_csv_file}: True vs Predicted PBPotential")
        plt.plot([df_pred['PBPotential'].min(), df_pred['PBPotential'].max()], [df_pred['PBPotential'].min(), df_pred['PBPotential'].max()], 'k--', lw=2)
        # print the R^2 and RMSE on the plot
        plt.text(0.05, 0.95, f"R^2: {r2_pred:.3f} Good if > 0.9\nRMSE: {normalized_rmse_pred:.3f} Good if < 0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
        plt.savefig(f"{prediction_csv_file.replace('_compiled.csv', '_predicted_mix.png')}")

    plt.show()
