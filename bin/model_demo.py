#!/usr/bin/env python
"""
Demonstration of the ANN model to:
1. capture non-linear relationships in the data
2. catpture the interaction between features
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.neural_network import MLPRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor

# Create synthetic 1D nonlinear data
np.random.seed(0)
X = np.sort(np.random.rand(200, 1) * 2 - 1, axis=0)  # X in [-1, 1]
y = X**3

# Train/test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
# Standardize the features
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# Linear regression model
linear_model = LinearRegression()
linear_model.fit(X_train, y_train)
# Predict with linear model
y_pred_linear = linear_model.predict(X_test)
# Evaluate linear model
linear_rmse = np.sqrt(mean_squared_error(y_test, y_pred_linear))
linear_r2 = r2_score(y_test, y_pred_linear)
print(f"Linear Model - RMSE: {linear_rmse:.3f}, R2: {linear_r2:.3f}")
# Plot linear model predictions
plt.figure(figsize=(24, 6))
plt.subplot(1, 4, 1)
plt.scatter(X_test, y_test, color='blue', label='True Values')
plt.scatter(X_test, y_pred_linear, color='red', label='Predicted Values')
plt.title('Linear Regression Predictions')
plt.xlabel('X')
plt.ylabel('y')
plt.legend()
plt.grid()

# Neural network model
nn_model = MLPRegressor(hidden_layer_sizes=(20,10), activation='relu', solver='adam', max_iter=2000, random_state=42)
nn_model.fit(X_train, y_train.ravel())
# Predict with neural network model
y_pred_nn = nn_model.predict(X_test)
# Evaluate neural network model
nn_rmse = np.sqrt(mean_squared_error(y_test, y_pred_nn))
nn_r2 = r2_score(y_test, y_pred_nn)
print(f"Neural Network Model - RMSE: {nn_rmse:.3f}, R2: {nn_r2:.3f}")
# Plot neural network model predictions
plt.subplot(1, 4, 2)
plt.scatter(X_test, y_test, color='blue', label='True Values')
plt.scatter(X_test, y_pred_nn, color='red', label='Predicted Values')
plt.title('Neural Network Predictions')
plt.xlabel('X')
plt.ylabel('y')
plt.legend()
plt.grid()
plt.tight_layout()


# K-Nearest Neighbors model
knn_model = KNeighborsRegressor(n_neighbors=5)
knn_model.fit(X_train, y_train.ravel())
# Predict with KNN model
y_pred_knn = knn_model.predict(X_test)
# Evaluate KNN model
knn_rmse = np.sqrt(mean_squared_error(y_test, y_pred_knn))
knn_r2 = r2_score(y_test, y_pred_knn)
print(f"KNN Model - RMSE: {knn_rmse:.3f}, R2: {knn_r2:.3f}")
# Plot KNN model predictions
plt.subplot(1, 4, 3)
plt.scatter(X_test, y_test, color='blue', label='True Values')
plt.scatter(X_test, y_pred_knn, color='red', label='Predicted Values')
plt.title('KNN Predictions')
plt.xlabel('X')
plt.ylabel('y')
plt.legend()
plt.grid()


# XGBoost model
xgb_model = XGBRegressor(objective='reg:squarederror', n_estimators=1000, learning_rate=0.1, max_depth=3, random_state=42)
xgb_model.fit(X_train, y_train.ravel())
# Predict with XGBoost model
y_pred_xgb = xgb_model.predict(X_test)
# Evaluate XGBoost model
xgb_rmse = np.sqrt(mean_squared_error(y_test, y_pred_xgb))
xgb_r2 = r2_score(y_test, y_pred_xgb)
print(f"XGBoost Model - RMSE: {xgb_rmse:.3f}, R2: {xgb_r2:.3f}")
# Plot XGBoost model predictions
plt.subplot(1, 4, 4)
plt.scatter(X_test, y_test, color='blue', label='True Values')
plt.scatter(X_test, y_pred_xgb, color='red', label='Predicted Values')
plt.title('XGBoost Predictions')
plt.xlabel('X')
plt.ylabel('y')
plt.legend()
plt.grid()
plt.tight_layout()



plt.show()