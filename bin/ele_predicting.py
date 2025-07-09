#!/usr/bin/env python
"""
Use  a selected ML model to predict on another set of data in CSV format.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np
import joblib
import argparse
import logging
import time

def predict_rf(data, model, title):
    X = data[model['features']]
    y = data['PBPotential']
    rf = model['model']
    scaler = model['scaler']
    # Standardize the features
    X_scaled = scaler.transform(X)
    # Predict using the model
    logging.info(f"Predicting with {title}...")
    start_time = time.time()
    y_pred = rf.predict(X_scaled)
    logging.info(f"Prediction completed in {time.time() - start_time:.2f} seconds.")

    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    y_range = np.ptp(y)  # Range of true values
    normalized_rmse = rmse / y_range if y_range != 0 else 0
    r2 = r2_score(y, y_pred)
    logging.info(f"{title} - R2: {r2:.3f}, RMSE: {normalized_rmse:.3f}")
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y, y=y_pred)

    plt.xlabel("True PBPotential")
    plt.ylabel("Predicted PBPotential")
    plt.title(f"Prediction Results: {title}")
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)
    plt.grid(True)
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    # print the RMSE and R2 score on the plot
    plt.text(0.05, 0.95, f"R2: {r2:.3f}\nRMSE: {normalized_rmse:.3f}", transform=plt.gca().transAxes, fontsize=12,
             verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5, edgecolor='black'))
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
    parser.add_argument("-m", "--model", default="random_forest_with_scaler.pkl", help="Path to the saved model and scaler in pkl format")
    args = parser.parse_args()

    # Load the data
    logging.info(f"Loading data from {args.input_csv} ...")
    data = pd.read_csv(args.input_csv)
    data['iDistance'] = 1 / data['Distance'] # Inverse distance


    # predict using the saved model
    fname = args.model   # Load model and scaler from args.model
    title = "Pre-trained Model"
    logging.info("Loading model %s" % fname)
    try:
        model = joblib.load(fname)
    except FileNotFoundError as e:
        logging.error(f"Model not found: {e}")
        exit(1)
    predict_rf(data, model, title)


    plt.show()  # Show the plot interactively