#!/usr/bin/env python
"""
Tune the XGBoost model for regression tasks.
Use medium protein as test data set.
ele: medium_protein_ele.csv
rxn: medium_protein_rxn.csv
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score
from xgboost import XGBRegressor
import logging
import time
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO)

def tune_xgboost(X, y):

    param_grid = {
        'max_depth': [3, 5],
        'min_child_weight': [7, 9, 11],
        'gamma': [0.2, 0.4],
        'subsample': [1.0],
        'colsample_bytree': [1.0],
        'learning_rate': [0.05, 0.1]    
    }
# ele
#  {'colsample_bytree': 1.0, 'gamma': 0.1, 'learning_rate': 0.07, 'max_depth': 5, 'min_child_weight': 3, 'subsample': 1.0}
#  {'colsample_bytree': 1.0, 'gamma': 0.05, 'learning_rate': 0.07, 'max_depth': 5, 'min_child_weight': 3, 'subsample': 1.0}
#  {'colsample_bytree': 1.0, 'gamma': 0.05, 'learning_rate': 0.07, 'max_depth': 5, 'min_child_weight': 3, 'subsample': 1.0}
# rxn
#  {'colsample_bytree': 1.0, 'gamma': 0.5, 'learning_rate': 0.1, 'max_depth': 3, 'min_child_weight': 5, 'subsample': 1.0}
#  {'colsample_bytree': 1.0, 'gamma': 0.2, 'learning_rate': 0.1, 'max_depth': 3, 'min_child_weight': 7, 'subsample': 1.0}
#  {'colsample_bytree': 1.0, 'gamma': 0.2, 'learning_rate': 0.1, 'max_depth': 3, 'min_child_weight': 9, 'subsample': 1.0}
#  {'colsample_bytree': 1.0, 'gamma': 0.2, 'learning_rate': 0.1, 'max_depth': 3, 'min_child_weight': 9, 'subsample': 1.0}
    grid_search = GridSearchCV(XGBRegressor(objective='reg:squarederror'), param_grid, cv=3, scoring='neg_mean_squared_error')
    grid_search.fit(X, y)
    logging.info(f"Best parameters: {grid_search.best_params_}")
    logging.info(f"Best cross-validation score: {-grid_search.best_score_}")
    return grid_search.best_estimator_

if __name__ == "__main__":
    print("Run this above folder ele and rxn to get the data.")

    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s - %(message)s')

    # Load the data
    logging.info("Loading data...")
    # ele_data = pd.read_csv('ele/medium_protein_compiled.csv')
    rxn_data = pd.read_csv('rxn/medium_protein_rxn.csv')

    # Fit ele_data
    # ele_data['iDistance'] = 1 / ele_data['Distance']
    # features = ['iDistance', 'Density1_Near', 'Density1_Mid', 'Density1_Far', 'D2surface1', 'Density2_Near', 'Density2_Mid', 'Density2_Far', 'D2surface2']
    # target = 'PBPotential'

    # X = ele_data[features]
    # y = ele_data[target]

    # ele_model = tune_xgboost(X, y)
    # logging.info("Model trained on ele_data.")

    # # use the model to predict on the whole ele_data
    # y_pred = ele_model.predict(X)
    # rmse = np.sqrt(mean_squared_error(y, y_pred))
    # r2 = r2_score(y, y_pred)
    # logging.info(f"Ele Data - RMSE: {rmse:.3f}, R2: {r2:.3f}")

    # # plot the results
    # plt.figure(figsize=(10, 6))
    # sns.scatterplot(x=y, y=y_pred)
    # plt.xlabel("True PBPotential")
    # plt.ylabel("Predicted PBPotential")
    # plt.title("Prediction Results on Ele Data")
    # plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)
    # plt.grid(True)
    # plt.xlim(y.min(), y.max())
    # plt.ylim(y.min(), y.max())
    # plt.text(0.05, 0.95, f"R2: {r2:.3f}\nRMSE: {rmse:.3f}", transform=plt.gca().transAxes)
    # plt.tight_layout()
    # plt.savefig("ele_data_predictions.png")
    
    # Fit rxn_data
    features = ['Density_Near', 'Density_Mid', 'Density_Far', 'D2surface']
    target = 'PBRXN'
    X = rxn_data[features]
    y = rxn_data[target]
    rxn_model = tune_xgboost(X, y)
    logging.info("Model trained on rxn_data.")
    # use the model to predict on the whole rxn_data
    y_pred = rxn_model.predict(X)
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    r2 = r2_score(y, y_pred)
    logging.info(f"Rxn Data - RMSE: {rmse:.3f}, R2: {r2:.3f}") 
    # plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y, y=y_pred)
    plt.xlabel("True PBRXN")
    plt.ylabel("Predicted PBRXN")
    plt.title("Prediction Results on Rxn Data")
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)
    plt.grid(True)
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    plt.text(0.05, 0.95, f"R2: {r2:.3f}\nRMSE: {rmse:.3f}", transform=plt.gca().transAxes)
    plt.tight_layout()

    plt.show()