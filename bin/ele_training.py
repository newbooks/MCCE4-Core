#!/usr/bin/env python
"""
Train a model to predict electrostatic energy based on distance and Local Density

The input is a CSV file with the following columns:
- Distance: Distance between the two atoms
- DensityAverage_Near: Average local density near the atom
- DensityAverage_Mid: Average local density at mid-range
- DensityAverage_Far: Average local density far from the atom
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
from xgboost import XGBRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

import joblib
import argparse
import logging
import time

def fit_rf(features, target, data, title):
    X = data[features]
    y = data[target]
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
    logging.info(f"Adjusted Coulomb Potential - R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
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
    plt.xlabel("True PB Potential")
    plt.ylabel("Predicted PB Potential")
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


def fit_xgb(features, target, data, title):
    X = data[features]
    y = data[target]
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train an XGBoost model
    logging.info(f"Training with {title}...")
    xgb_model = XGBRegressor(n_estimators=100, random_state=int(time.time()))
    xgb_model.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = xgb_model.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2 = r2_score(y_val, y_pred_adjusted)
    logging.info(f"Adjusted Coulomb Potential - R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
    # get feature importances
    feature_importances_adjusted = xgb_model.feature_importances_
    feature_names_adjusted = X.columns
    # Log the feature importances
    for name, importance in zip(feature_names_adjusted, feature_importances_adjusted):
        logging.info(f"Feature: {name}, Importance: {importance:.4f}")
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val, y=y_pred_adjusted, alpha=0.5)
    plt.grid(True)  # add grid lines
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)  # Diagonal line
    plt.xlabel("True PB Potential")
    plt.ylabel("Predicted PB Potential")
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
    joblib.dump({'model': xgb_model, 'scaler': scaler, 'features': feature_names}, model_filename)
    logging.info(f"Saved the trained model and scaler to {model_filename}.")

def fit_ann(features, target, data, title):
    X = data[features]
    y = data[target]
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train an ANN model
    logging.info(f"Training with {title}...")
    ann_model = MLPRegressor(hidden_layer_sizes=(100,), max_iter=500, random_state=int(time.time()))
    ann_model.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = ann_model.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2 = r2_score(y_val, y_pred_adjusted)
    logging.info(f"Adjusted Coulomb Potential - R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val, y=y_pred_adjusted, alpha=0.5)
    plt.grid(True)  # add grid lines
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)  # Diagonal line
    plt.xlabel("True PB Potential")
    plt.ylabel("Predicted PB Potential")
    plt.title(f"{title}")
    # print R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2:.3f} Good if > 0.9\nRMSE: {normalized_rmse_adjusted:.3E} Good if <0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    plt.savefig(f"{title}.png")
    # Save the trained model
    logging.info(f"Saving the trained model and scaler for {title} ...")
    model_filename = f"{title.replace(' ', '_').lower()}_with_scaler.pkl"
    feature_names = [f.replace(' ', '_') for f in features]  # Replace spaces with underscores in feature names
    joblib.dump({'model': ann_model, 'scaler': scaler, 'features': feature_names}, model_filename)
    logging.info(f"Saved the trained model and scaler to {model_filename}.")

def fit_linear(features, target, data, title):
    X = data[features]
    y = data[target]
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train a Linear Regression model
    logging.info(f"Training with {title}...")
    linear_model = LinearRegression()
    linear_model.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = linear_model.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2 = r2_score(y_val, y_pred_adjusted)
    logging.info(f"Adjusted Coulomb Potential - R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
    # get feature importances (coefficients for linear regression)
    feature_importances_adjusted = np.abs(linear_model.coef_)
    feature_names_adjusted = X.columns
    # Log the intercept 
    logging.info(f"Intercept: {linear_model.intercept_:.4f}")
    # Log the coefficients
    logging.info("Linear Regression Coefficients:")
    for name, coef in zip(feature_names_adjusted, linear_model.coef_):
        logging.info(f"Feature: {name}, Coefficient: {coef:.4f}")
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val, y=y_pred_adjusted, alpha=0.5)
    plt.grid(True)  # add grid lines
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)  # Diagonal line
    plt.xlabel("True PB Potential")
    plt.ylabel("Predicted PB Potential")
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
    joblib.dump({'model': linear_model, 'scaler': scaler, 'features': feature_names}, model_filename)
    logging.info(f"Saved the trained model and scaler to {model_filename}.")

def fit_ridge(features, target, data, title):
    X = data[features]
    y = data[target]
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train a Ridge Regression model
    logging.info(f"Training with {title}...")
    ridge_model = Ridge(alpha=0.1, random_state=int(time.time()))
    ridge_model.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = ridge_model.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2 = r2_score(y_val, y_pred_adjusted)
    logging.info(f"Adjusted Coulomb Potential - R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
    # get feature importances (coefficients for ridge regression)
    feature_importances_adjusted = np.abs(ridge_model.coef_)
    feature_names_adjusted = X.columns
    # Log the feature importances
    for name, importance in zip(feature_names_adjusted, feature_importances_adjusted):
        logging.info(f"Feature: {name}, Importance: {importance:.4f}")
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val, y=y_pred_adjusted, alpha=0.5)
    plt.grid(True)  # add grid lines
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)  # Diagonal line
    plt.xlabel("True PB Potential")
    plt.ylabel("Predicted PB Potential")
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
    joblib.dump({'model': ridge_model, 'scaler': scaler, 'features': feature_names}, model_filename)
    logging.info(f"Saved the trained model and scaler to {model_filename}.")

def fit_lasso(features, target, data, title):
    X = data[features]
    y = data[target]
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train a Lasso Regression model
    logging.info(f"Training with {title}...")
    lasso_model = Lasso(alpha=0.01, random_state=int(time.time()))
    lasso_model.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = lasso_model.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2 = r2_score(y_val, y_pred_adjusted)
    logging.info(f"Adjusted Coulomb Potential - R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
    # get feature importances (coefficients for lasso regression)
    feature_importances_adjusted = np.abs(lasso_model.coef_)
    feature_names_adjusted = X.columns
    # Log the feature importances
    for name, importance in zip(feature_names_adjusted, feature_importances_adjusted):
        logging.info(f"Feature: {name}, Importance: {importance:.4f}")
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val, y=y_pred_adjusted, alpha=0.5)
    plt.grid(True)  # add grid lines
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)  # Diagonal line
    plt.xlabel("True PB Potential")
    plt.ylabel("Predicted PB Potential")
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
    joblib.dump({'model': lasso_model, 'scaler': scaler, 'features': feature_names}, model_filename)
    logging.info(f"Saved the trained model and scaler to {model_filename}.")


def fit_svr(features, target, data, title):
    X = data[features]
    y = data[target]
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train a Support Vector Regression model
    logging.info(f"Training with {title}...")
    svr_model = SVR(kernel='rbf', C=1.0, epsilon=0.1)
    svr_model.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = svr_model.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2 = r2_score(y_val, y_pred_adjusted)
    logging.info(f"Adjusted Coulomb Potential - R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val, y=y_pred_adjusted, alpha=0.5)
    plt.grid(True)  # add grid lines
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)  # Diagonal line
    plt.xlabel("True PB Potential")
    plt.ylabel("Predicted PB Potential")
    plt.title(f"{title}")
    # print R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2:.3f} Good if > 0.9\nRMSE: {normalized_rmse_adjusted:.3E} Good if <0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    plt.savefig(f"{title}.png")
    # Save the trained model
    logging.info(f"Saving the trained model and scaler for {title} ...")
    model_filename = f"{title.replace(' ', '_').lower()}_with_scaler.pkl"
    feature_names = [f.replace(' ', '_') for f in features]  # Replace spaces with underscores in feature names
    joblib.dump({'model': svr_model, 'scaler': scaler, 'features': feature_names}, model_filename)
    logging.info(f"Saved the trained model and scaler to {model_filename}.")

def fit_knn(features, target, data, title):
    X = data[features]
    y = data[target]
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train a K-Nearest Neighbors model
    logging.info(f"Training with {title}...")
    knn_model = KNeighborsRegressor(n_neighbors=3)
    knn_model.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = knn_model.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2 = r2_score(y_val, y_pred_adjusted)
    logging.info(f"Adjusted Coulomb Potential - R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val, y=y_pred_adjusted, alpha=0.5)
    plt.grid(True)  # add grid lines
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)  # Diagonal line
    plt.xlabel("True PB Potential")
    plt.ylabel("Predicted PB Potential")
    plt.title(f"{title}")
    # print R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2:.3f} Good if > 0.9\nRMSE: {normalized_rmse_adjusted:.3E} Good if <0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    plt.savefig(f"{title}.png")
    # Save the trained model
    logging.info(f"Saving the trained model and scaler for {title} ...")
    model_filename = f"{title.replace(' ', '_').lower()}_with_scaler.pkl"
    feature_names = [f.replace(' ', '_') for f in features]  # Replace spaces with underscores in feature names
    joblib.dump({'model': knn_model, 'scaler': scaler, 'features': feature_names}, model_filename)
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

    # Feature Engineering: Add a new column AdjustedPotential that is PBPotential / Distance
    # data['AdjustedPotential'] = data['PBPotential'] / data['Distance']
    # features = ['DensityAverage_Near', 'DensityAverage_Mid', 'DensityAverage_Far']
    # target = 'AdjustedPotential'
    data['iDistance'] = 1 / data['Distance'] 
    # Inverse distance, by the nature of the PB potential, the closer the atoms are, the stronger the interaction scalable by the inverse of distance. 
    # So using iDistance make KNN model better but won't affect decision tree based models much.

    features = ['iDistance', 'DensityAverage_Near', 'DensityAverage_Mid', 'DensityAverage_Far']
    target = 'PBPotential'

    title = "Random Forest"
    fit_rf(features, target, data, title)
    logging.info(f"Model trained by {title}.")

    # title = "XGBoost"
    # fit_xgb(features, target, data, title)
    # logging.info(f"Model trained by {title}.")

    title = "ANN"
    fit_ann(features, target, data, title)
    logging.info(f"Model trained by {title}.")

    title = "KNN"
    fit_knn(features, target, data, title)
    logging.info(f"Model trained by {title}.")

    # title = "SVR"
    # fit_svr(features, target, data, title)
    # logging.info(f"Model trained by {title}.")

    # title = "Linear Regression"
    # fit_linear(features, target, data, title)
    # logging.info(f"Model trained by {title}.")

    # title = "Ridge Regression"
    # fit_ridge(features, target, data, title)
    # logging.info(f"Model trained by {title}.")

    # title = "Lasso Regression"
    # fit_lasso(features, target, data, title)
    # logging.info(f"Model trained by {title}.")

    plt.show()
    plt.close('all')
