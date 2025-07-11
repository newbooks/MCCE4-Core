#!/usr/bin/env python
"""
Tune the Neural Network model for regression tasks.
Use medium protein as test data set.
ele: medium_protein_ele.csv
rxn: medium_protein_rxn.csv
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor
import logging
import time
import matplotlib.pyplot as plt
import seaborn as sns

def tune_ann(X_scaled, y):
    param_grid = {
        'hidden_layer_sizes': [(50, 20)],
        # 'activation': ['relu', 'tanh'],
        'alpha': [0.01],
        'learning_rate_init': [0.001],
        'learning_rate': ['adaptive'],
        # 'learning_rate': ['constant', 'adaptive'],
        # 'solver': ['adam', 'sgd'],
        'early_stopping': [True]
    }
    # {'alpha': 0.01, 'early_stopping': True, 'hidden_layer_sizes': (50, 20), 'learning_rate': 'adaptive', 'learning_rate_init': 0.001}
    grid_search = GridSearchCV(MLPRegressor(max_iter=500), param_grid, cv=3, scoring='neg_mean_squared_error')
    grid_search.fit(X_scaled, y)
    logging.info(f"Best parameters: {grid_search.best_params_}")
    logging.info(f"Best cross-validation score: {-grid_search.best_score_}")
    return grid_search.best_estimator_

if __name__ == "__main__":
    print("Run this above folder ele and rxn to get the data.")

    # Set up logging
    logging_format = "%(asctime)s %(levelname)s: %(message)s"
    logging_datefmt='%Y-%m-%d %H:%M:%S'
    logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)

    # Load the data
    logging.info("Loading data...")
    ele_data = pd.read_csv('ele/medium_protein_compiled.csv')
    # rxn_data = pd.read_csv('rxn/medium_protein_rxn.csv')

    # Fit ele_data
    ele_data['iDistance'] = 1 / ele_data['Distance']
    features = ['iDistance', 'Density1_Near', 'Density1_Mid', 'Density1_Far', 'D2surface1', 'Density2_Near', 'Density2_Mid', 'Density2_Far', 'D2surface2']
    target = 'PBPotential'

    X = ele_data[features]
    y = ele_data[target]

    # standardize the features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    ele_model = tune_ann(X_scaled, y)
    logging.info("Model trained on ele_data.")

    # use the model to predict on the whole ele_data
    y_pred = ele_model.predict(X_scaled)
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    r2 = r2_score(y, y_pred)
    logging.info(f"Ele Data - RMSE: {rmse:.3f}, R2: {r2:.3f}")

    # plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y, y=y_pred)
    plt.xlabel("True PBPotential")
    plt.ylabel("Predicted PBPotential")
    plt.title("Prediction Results on Ele Data")
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)
    plt.grid(True)
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    plt.text(0.05, 0.95, f"R2: {r2:.3f}\nRMSE: {rmse:.3f}", transform=plt.gca().transAxes)
    plt.tight_layout()
    plt.savefig("ele_data_predictions.png")
    
    # Fit rxn_data
    # features = ['Density_Near', 'Density_Mid', 'Density_Far', 'D2surface']
    # target = 'PBRXN'
    # X = rxn_data[features]
    # y = rxn_data[target]
    # rxn_model = tune_xgboost(X, y)
    # logging.info("Model trained on rxn_data.")
    # # use the model to predict on the whole rxn_data
    # y_pred = rxn_model.predict(X)
    # rmse = np.sqrt(mean_squared_error(y, y_pred))
    # r2 = r2_score(y, y_pred)
    # logging.info(f"Rxn Data - RMSE: {rmse:.3f}, R2: {r2:.3f}") 
    # # plot the results
    # plt.figure(figsize=(10, 6))
    # sns.scatterplot(x=y, y=y_pred)
    # plt.xlabel("True PBRXN")
    # plt.ylabel("Predicted PBRXN")
    # plt.title("Prediction Results on Rxn Data")
    # plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)
    # plt.grid(True)
    # plt.xlim(y.min(), y.max())
    # plt.ylim(y.min(), y.max())
    # plt.text(0.05, 0.95, f"R2: {r2:.3f}\nRMSE: {rmse:.3f}", transform=plt.gca().transAxes)
    # plt.tight_layout()

    plt.show()