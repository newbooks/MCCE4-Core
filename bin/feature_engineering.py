#!/usr/bin/env python3
"""
Feature Engineering
Examine the data (CSV file) and create combined features to reduce the number of features.
"""

import pandas as pd
import numpy as np
import logging
import argparse
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.neural_network import MLPRegressor
from matplotlib import pyplot as plt
import seaborn as sns
import time
import joblib


def fit_ann(X_scaled, y, title):
    """
    Fit an ANN model to the scaled features and target variable.
    This function is a placeholder for the actual ANN fitting logic.
    """
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=int(time.time()) % 1000)

    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    model = MLPRegressor(hidden_layer_sizes=(50, 20), alpha=0.01, learning_rate_init=0.001, learning_rate='adaptive', max_iter=500, random_state=int(time.time()) % 1000, early_stopping=True)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    r2 = r2_score(y_test, y_pred)

    logger.info(f"Results for {title}:")
    logger.info(f"Root Mean Squared Error: {rmse}")
    logger.info(f"R^2 Score: {r2}")

    # Plot the predictions vs actual values
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_test, y=y_pred)
    plt.xlabel("Actual Values")
    plt.ylabel("Predicted Values")
    plt.title(f"Actual vs Predicted - {title}")
    # Plot the diagonal line
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)
    plt.grid(True)
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    plt.text(0.05, 0.95, f"R^2: {r2:.3f}\nRMSE: {rmse:.3f}", transform=plt.gca().transAxes)
    plt.tight_layout()
    plt.savefig(f"{title.replace(' ', '_')}_predictions.png")
    joblib.dump(model, f"{title.replace(' ', '_')}.pkl")


if __name__ == "__main__":
    logging_format = "%(asctime)s %(levelname)s: %(message)s"
    logging_datefmt='%Y-%m-%d %H:%M:%S'
    logging.basicConfig(
        format=logging_format,
        datefmt=logging_datefmt,
        level=logging.INFO
    )

    parser = argparse.ArgumentParser(description="Test reduced feature set by combining features.")
    parser.add_argument("input_file", type=str, help="Path to the input CSV file")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # Load data
    df = pd.read_csv(args.input_file)
    logger.info(f"Loaded data from {args.input_file} with shape {df.shape}")

    df['iDistance'] = 1 / df['Distance']

    # Original features
    features = ['Distance', 'Density1_Near', 'Density1_Mid', 'Density1_Far', 'D2surface1', 'Density2_Near', 'Density2_Mid', 'Density2_Far', 'D2surface2']
    target = 'PBPotential'
    fit_ann(df[features], df[target], "Original Features")

    # Use Inversed Distance
    features = ['iDistance', 'Density1_Near', 'Density1_Mid', 'Density1_Far', 'D2surface1', 'Density2_Near', 'Density2_Mid', 'Density2_Far', 'D2surface2']
    target = 'PBPotential'
    fit_ann(df[features], df[target], "Inversed Distance")

    # Use Lower Denity and D2surface
    # df['LowerDensity_Near'] = df[['Density1_Near', 'Density2_Near']].min(axis=1)
    # df['LowerDensity_Mid'] = df[['Density1_Mid', 'Density2_Mid']].min(axis=1)
    # df['LowerDensity_Far'] = df[['Density1_Far', 'Density2_Far']].min(axis=1)
    # df['LowerD2surface'] = df[['D2surface1', 'D2surface2']].min(axis=1)
    # features = ['iDistance', 'LowerDensity_Near', 'LowerDensity_Mid', 'LowerDensity_Far', 'LowerD2surface']
    # target = 'PBPotential'
    # fit_ann(df[features], df[target], "Lower Density and D2surface") 

    # Use Average Density and D2surface. Since these numbers are averaged, the opposite side should be removed.
    df['AverageDensity_Near'] = df[['Density1_Near', 'Density2_Near']].mean(axis=1)
    df['AverageDensity_Mid'] = df[['Density1_Mid', 'Density2_Mid']].mean(axis=1)
    df['AverageDensity_Far'] = df[['Density1_Far', 'Density2_Far']].mean(axis=1)
    df['AverageD2surface'] = df[['D2surface1', 'D2surface2']].mean(axis=1)
    # show number of rows before and after removing the opposite side
    logger.info(f"Number of rows before removing opposite side: {df.shape[0]}")

    # Make a copy of DF, for the following selected features, average the target if features are within the torlerance.
    tolerance = 0.1
    features = ['iDistance', 'AverageDensity_Near', 'AverageDensity_Mid', 'AverageDensity_Far', 'AverageD2surface']
 
    # Make a copy to avoid changing original
    df_copy = df.copy()

    # Bin features by tolerance (group values that are within ±tolerance/2)
    for f in features:
        df_copy[f + '_bin'] = (df_copy[f] / tolerance).round().astype(int)

    # Group by the binned features and average both features and PBPotential
    grouped = df_copy.groupby([f + '_bin' for f in features]).agg({
        **{f: 'mean' for f in features},    # average original feature values
        'PBPotential': 'mean'              # average target
    }).reset_index(drop=True)
 
    logger.info(f"Number of rows after removing opposite side: {grouped.shape[0]}")

    target = 'PBPotential'
    fit_ann(grouped[features], grouped[target], "Average Density and D2surface")

    # Use Upper Density and D2surface
    # df['UpperDensity_Near'] = df[['Density1_Near', 'Density2_Near']].max(axis=1)
    # df['UpperDensity_Mid'] = df[['Density1_Mid', 'Density2_Mid']].max(axis=1)
    # df['UpperDensity_Far'] = df[['Density1_Far', 'Density2_Far']].max(axis=1)
    # df['UpperD2surface'] = df[['D2surface1', 'D2surface2']].max(axis=1)
    # features = ['iDistance', 'UpperDensity_Near', 'UpperDensity_Mid', 'UpperDensity_Far', 'UpperD2surface']
    # target = 'PBPotential'
    # fit_ann(df[features], df[target], "Upper Density and D2surface")

    # Use Inversed D2surface.
    grouped['iD2surface'] = 1 / grouped['AverageD2surface']
    grouped['iD2surface'] = grouped['iD2surface'].replace([np.inf, -np.inf], np.nan)  # Replace inf with NaN
    grouped['iD2surface'] = grouped['iD2surface'].fillna(0)
    features = ['iDistance', 'AverageDensity_Near', 'AverageDensity_Mid', 'AverageDensity_Far', 'iD2surface']
    target = 'PBPotential'
    fit_ann(grouped[features], grouped[target], "Average Density and Inversed D2surface")


    plt.show()
