#!/usr/bin/env python
"""
Train a model to predict electrostatic energy based on embedding scores and distances.

The input is a CSV file with the following columns:
- Conf1: Conformer ID for atom 1
- Conf2: Conformer ID for atom 2
- Distance: Distance between two atoms in Angstroms
- Radius1: Radius for atom 1
- Radius2: Radius for atom 2
- Embedding1: Embedding score for atom 1
- Embedding2: Embedding score for atom 2
- Density1: Density score for atom 1
- Density2: Density score for atom 2
- CoulombPotential: Coulomb potential between two atoms
- AdjustedCoulombPotential: Adjusted Coulomb potential based on embedding scores
- PBPotential: Electrostatic energy from Poisson-Boltzmann calculation
"""


# import the modules
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
from xgboost import XGBRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.linear_model import LinearRegression, Lasso, Ridge
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
import joblib
import argparse
import logging
import time

def fit_rf(features, data, title):
    X = data[features]
    y = data['PBRXN']
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train a Random Forest model
    logging.info(f"Training with {title}...")
    rf = RandomForestRegressor(n_estimators=1000, random_state=int(time.time()))
    rf.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = rf.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2 = r2_score(y_val, y_pred_adjusted)
    logging.info(f"R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
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
    plt.xlabel("True PB RXN")
    plt.ylabel("Predicted PB RXN")
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

def fit_xgb(features, data, title):
    X = data[features]
    y = data['PBRXN']
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train an XGBoost model
    logging.info(f"Training with {title}...")
    xgb_model = XGBRegressor(
        n_estimators=1000,         # Number of trees (will use early stopping)
        learning_rate=0.01,        # Low LR for smoother learning
        max_depth=3,               # Shallow trees are enough for small feature sets
        subsample=0.8,             # Random sampling of rows
        colsample_bytree=1.0,      # Use all features (since you only have 3)
        random_state=int(time.time())
    )
    xgb_model.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = xgb_model.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2 = r2_score(y_val, y_pred_adjusted)
    logging.info(f"R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
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
    plt.xlabel("True PB RXN")
    plt.ylabel("Predicted PB RXN")
    plt.title(f"{title}")
    # print R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2:.3f} Good if > 0.9\nRMSE: {normalized_rmse_adjusted:.3E} Good if <0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    # Sort and add feature importances to the plot
    feature_importance_adjusted_df = pd.DataFrame({'Feature': feature_names_adjusted, 'Importance': feature_importances_adjusted})
    feature_importance_adjusted_df = feature_importance_adjusted_df.sort_values(by='Importance', ascending=False)
    plt.text(0.05, 0.85, "\n".join([f"{row['Feature']}: {row['Importance']:.4f}" for _, row in feature_importance_adjusted_df.iterrows()]), transform=plt.gca().transAxes, fontsize=10, verticalalignment='top')
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    plt.savefig(f"{title.replace(' ', '_')}.png")
    # Save the trained model
    logging.info(f"Saving the trained model and scaler for {title} ...")
    model_filename = f"{title.replace(' ', '_').lower()}_with_scaler.pkl"
    feature_names = [f.replace(' ', '_') for f in features]  # Replace spaces with underscores in feature names
    joblib.dump({'model': xgb_model, 'scaler': scaler, 'features': feature_names}, model_filename)
    logging.info(f"Saved the trained model and scaler to {model_filename}.")


def fit_ann(features, data, title):
    X = data[features]
    y = data['PBRXN']
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train an ANN model
    logging.info(f"Training with {title}...")
    ann_model = MLPRegressor(hidden_layer_sizes=(10, 5), max_iter=1000, random_state=int(time.time()))
    ann_model.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = ann_model.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2 = r2_score(y_val, y_pred_adjusted)
    logging.info(f"R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
    # get feature importances (not directly available for ANN, so we skip this part)
    
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val, y=y_pred_adjusted, alpha=0.5)
    plt.grid(True)  # add grid lines
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)  # Diagonal line
    plt.xlabel("True PB RXN")
    plt.ylabel("Predicted PB RXN")
    plt.title(f"{title}")
    # print R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2:.3f} Good if > 0.9\nRMSE: {normalized_rmse_adjusted:.3E} Good if <0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    plt.savefig(f"{title.replace(' ', '_')}.png")
    # Save the trained model
    logging.info(f"Saving the trained model and scaler for {title} ...")
    model_filename = f"{title.replace(' ', '_').lower()}_with_scaler.pkl"
    feature_names = [f.replace(' ', '_') for f in features]  # Replace spaces with underscores in feature names
    joblib.dump({'model': ann_model, 'scaler': scaler, 'features': feature_names}, model_filename)
    logging.info(f"Saved the trained model and scaler to {model_filename}.")

def fit_linear(features, data, title):
    X = data[features]
    y = data['PBRXN']
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
    logging.info(f"R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
    
    # get feature importances (not directly available for Linear Regression, so we skip this part)
    
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val, y=y_pred_adjusted, alpha=0.5)
    plt.grid(True)  # add grid lines
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)  # Diagonal line
    plt.xlabel("True PB RXN")
    plt.ylabel("Predicted PB RXN")
    plt.title(f"{title}")
    # print R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2:.3f} Good if > 0.9\nRMSE: {normalized_rmse_adjusted:.3E} Good if <0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    plt.savefig(f"{title.replace(' ', '_')}.png")
    # Save the trained model
    logging.info(f"Saving the trained model and scaler for {title} ...")
    model_filename = f"{title.replace(' ', '_').lower()}_with_scaler.pkl"
    feature_names = [f.replace(' ', '_') for f in features]  # Replace spaces with underscores in feature names
    joblib.dump({'model': linear_model, 'scaler': scaler, 'features': feature_names}, model_filename)
    logging.info(f"Saved the trained model and scaler to {model_filename}.")

def fit_lasso(features, data, title):
    X = data[features]
    y = data['PBRXN']
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train a Lasso Regression model
    logging.info(f"Training with {title}...")
    lasso_model = Lasso(alpha=0.1)  # Adjust alpha as needed
    lasso_model.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = lasso_model.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2 = r2_score(y_val, y_pred_adjusted)
    logging.info(f"R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
    
    # get feature importances (not directly available for Lasso Regression, so we skip this part)
    
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val, y=y_pred_adjusted, alpha=0.5)
    plt.grid(True)  # add grid lines
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)  # Diagonal line
    plt.xlabel("True PB RXN")
    plt.ylabel("Predicted PB RXN")
    plt.title(f"{title}")
    # print R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2:.3f} Good if > 0.9\nRMSE: {normalized_rmse_adjusted:.3E} Good if <0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    plt.savefig(f"{title.replace(' ', '_')}.png")
    # Save the trained model
    logging.info(f"Saving the trained model and scaler for {title} ...")
    model_filename = f"{title.replace(' ', '_').lower()}_with_scaler.pkl"
    feature_names = [f.replace(' ', '_') for f in features]  # Replace spaces with underscores in feature names
    joblib.dump({'model': lasso_model, 'scaler': scaler, 'features': feature_names}, model_filename)
    logging.info(f"Saved the trained model and scaler to {model_filename}.")

def fit_ridge(features, data, title):
    X = data[features]
    y = data['PBRXN']
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train a Ridge Regression model
    logging.info(f"Training with {title}...")
    ridge_model = Ridge(alpha=1.0)  # Adjust alpha as needed
    ridge_model.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = ridge_model.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2 = r2_score(y_val, y_pred_adjusted)
    logging.info(f"R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
    
    # get feature importances (not directly available for Ridge Regression, so we skip this part)
    
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val, y=y_pred_adjusted, alpha=0.5)
    plt.grid(True)  # add grid lines
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)  # Diagonal line
    plt.xlabel("True PB RXN")
    plt.ylabel("Predicted PB RXN")
    plt.title(f"{title}")
    # print R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2:.3f} Good if > 0.9\nRMSE: {normalized_rmse_adjusted:.3E} Good if <0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    plt.savefig(f"{title.replace(' ', '_')}.png")
    # Save the trained model
    logging.info(f"Saving the trained model and scaler for {title} ...")
    model_filename = f"{title.replace(' ', '_').lower()}_with_scaler.pkl"
    feature_names = [f.replace(' ', '_') for f in features]  # Replace spaces with underscores in feature names
    joblib.dump({'model': ridge_model, 'scaler': scaler, 'features': feature_names}, model_filename)
    logging.info(f"Saved the trained model and scaler to {model_filename}.")
        
def fit_knn(features, data, title):
    X = data[features]
    y = data['PBRXN']
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=int(time.time()))
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Train a KNN model
    logging.info(f"Training with {title}...")
    knn_model = KNeighborsRegressor(n_neighbors=5)  # Adjust n_neighbors as needed
    knn_model.fit(X_train, y_train)
    # Evaluate the model
    logging.info(f"Evaluating with {title} on validation set...")
    y_pred_adjusted = knn_model.predict(X_val)
    rmse_adjusted = np.sqrt(mean_squared_error(y_val, y_pred_adjusted))
    y_range_adjusted = np.ptp(y_val)  # Range of true values
    normalized_rmse_adjusted = rmse_adjusted / y_range_adjusted if y_range_adjusted != 0 else 0
    r2 = r2_score(y_val, y_pred_adjusted)
    logging.info(f"R2: {r2:.3f}, RMSE: {normalized_rmse_adjusted:.3f}")
    
    # get feature importances (not directly available for KNN, so we skip this part)
    
    # Plot the results
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=y_val, y=y_pred_adjusted, alpha=0.5)
    plt.grid(True)  # add grid lines
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'g--', lw=2)  # Diagonal line
    plt.xlabel("True PB RXN")
    plt.ylabel("Predicted PB RXN")
    plt.title(f"{title}")
    # print R^2 and RMSE on the plot
    plt.text(0.05, 0.95, f"R^2: {r2:.3f} Good if > 0.9\nRMSE: {normalized_rmse_adjusted:.3E} Good if <0.05", transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    
    plt.xlim(y.min(), y.max())
    plt.ylim(y.min(), y.max())
    plt.savefig(f"{title.replace(' ', '_')}.png")
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
    parser.add_argument("-m", "--model", default="rf", help="Model type to use (default: rf), options: xgb, ann, linear, lasso, ridge, knn.")
    parser.add_argument("input_csv", help="Input CSV file")
    args = parser.parse_args()

    # Load the data
    logging.info(f"Loading data from {args.input_csv} ...")
    data = pd.read_csv(args.input_csv)


    # Train the model
    features = ['Density_Near', 'Density_Mid', 'Density_Far', 'Density_Variance']
    title = "rxn"
    if args.model == "rf":
        title = f"Random Forest {title}"
        fit_rf(features, data, title)
    elif args.model == "xgb":
        title = f"XGBoost {title}"
        fit_xgb(features, data, title)
    elif args.model == "ann":
        title = f"ANN {title}"
        fit_ann(features, data, title)
    elif args.model == "linear":
        title = f"Linear Regression {title}"
        fit_linear(features, data, title)
    elif args.model == "ridge":
        title = f"Ridge Regression {title}"
        fit_ridge(features, data, title)
    elif args.model == "lasso":
        title = f"Lasso Regression {title}"
        fit_lasso(features, data, title)
    elif args.model == "knn":
        title = f"KNN {title}"
        fit_knn(features, data, title)
    else:
        logging.error(f"Unknown model type: {args.model}. Supported models are: rf, xgb, ann.")
        exit(1)
    logging.info(f"Model trained by {title}.")


    plt.show()
    plt.close('all')
