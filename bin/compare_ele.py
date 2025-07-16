#!/usr/bin/env python
"""
Compare ele energy from ML and Delphi
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging


PH2KCAL = 1.364

if __name__ == "__main__":
    logging_format = "%(asctime)s %(levelname)s: %(message)s"
    logging_datefmt='%Y-%m-%d %H:%M:%S'
    logging.basicConfig(format=logging_format, datefmt=logging_datefmt, level=logging.INFO)    


    # Load the data
    source_ml = "premade/energies"  # use raw files
    source_delphi = "4lzt/energies"  # use opp files

    # Collect energies from ML to dictionary
    ml_energies = {}
    for filename in os.listdir(source_ml):
        if filename.endswith(".raw"):
            conf_id = filename[:12]  # Extract the conformer ID from the filename
            lines = open(os.path.join(source_ml, filename), 'r').readlines()
            for line in lines:
                fields = line.split()
                conf2_id = fields[0]
                ele = float(fields[1])
                ml_energies[(conf_id, conf2_id)] = ele

    # Collect energies from Delphi to dictionary
    delphi_energies = {}
    for filename in os.listdir(source_delphi):
        if filename.endswith(".opp"):
            conf_id = filename[:3]+filename[5:14]
            lines = open(os.path.join(source_delphi, filename), 'r').readlines()
            for line in lines:
                fields = line.split()
                conf2_id = fields[1][:3]+fields[1][5:14]
                ele_average = float(fields[2])
                ele_raw = float(fields[5])
                delphi_energies[(conf_id, conf2_id)] = (ele_average, ele_raw)

    logging.info(f"Loaded ML energies: {len(ml_energies)} data points")
    logging.info(f"Loaded Delphi energies: {len(delphi_energies)} data points")

    # Prepare data for both plots
    x_raw = []
    y_raw = []
    x_avg = []
    y_avg = []
    for key in ml_energies:
        if key in delphi_energies:
            x_raw.append(delphi_energies[key][1])
            y_raw.append(ml_energies[key])
            x_avg.append(delphi_energies[key][0])
            y_avg.append(ml_energies[key])
        else:
            x_raw.append(0.0)
            y_raw.append(ml_energies[key])
            x_avg.append(0.0)
            y_avg.append(ml_energies[key])
    x_raw = np.array(x_raw)
    y_raw = np.array(y_raw)
    x_avg = np.array(x_avg)
    y_avg = np.array(y_avg)

    # Create side-by-side plots
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Plot 1: ML predicted value vs Delphi raw
    sns.scatterplot(x=x_raw, y=y_raw, alpha=0.5, ax=axes[0])
    axes[0].set_xlabel("Delphi Raw Energy (kcal/mol)")
    axes[0].set_ylabel("ML Predicted Energy (kcal/mol)")
    axes[0].set_title("ML Predicted vs Delphi Raw Energy")
    axes[0].plot([x_raw.min(), x_raw.max()], [x_raw.min(), x_raw.max()], 'g--', lw=2)
    axes[0].set_xlim(x_raw.min(), x_raw.max())
    axes[0].set_ylim(y_raw.min(), y_raw.max())
    delta = 1.0 * PH2KCAL
    axes[0].plot([x_raw.min(), x_raw.max()], [x_raw.min() - delta, x_raw.max() - delta], 'r--', lw=1)
    axes[0].plot([x_raw.min(), x_raw.max()], [x_raw.min() + delta, x_raw.max() + delta], 'r--', lw=1)
    within_delta_raw = np.sum((y_raw >= x_raw - delta) & (y_raw <= x_raw + delta))
    total_points_raw = len(x_raw)
    percentage_within_delta_raw = (within_delta_raw / total_points_raw) * 100 if total_points_raw > 0 else 0
    within_half_delta_raw = np.sum((y_raw >= x_raw - 0.5 * PH2KCAL) & (y_raw <= x_raw + 0.5 * PH2KCAL))
    percentage_within_half_delta_raw = (within_half_delta_raw / total_points_raw) * 100 if total_points_raw > 0 else 0
    axes[0].text(0.05, 0.90, f"Within ±0.5 pH unit: {percentage_within_half_delta_raw:.2f}%", transform=axes[0].transAxes,
                 fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))    
    within_quarter_delta_raw = np.sum((y_raw >= x_raw - 0.25 * PH2KCAL) & (y_raw <= x_raw + 0.25 * PH2KCAL))
    percentage_within_quarter_delta_raw = (within_quarter_delta_raw / total_points_raw) * 100 if total_points_raw > 0 else 0
    axes[0].text(0.05, 0.85, f"Within ±0.25 pH unit: {percentage_within_quarter_delta_raw:.2f}%", transform=axes[0].transAxes,
                 fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))    
    axes[0].text(0.05, 0.95, f"Within ±1 pH unit: {percentage_within_delta_raw:.2f}%", transform=axes[0].transAxes,
                 fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    axes[0].grid(True)

    # Plot 2: ML predicted value vs Delphi average
    sns.scatterplot(x=x_avg, y=y_avg, alpha=0.5, ax=axes[1])
    axes[1].set_xlabel("Delphi Average Energy (kcal/mol)")
    axes[1].set_ylabel("ML Predicted Energy (kcal/mol)")
    axes[1].set_title("ML Predicted vs Delphi Average Energy")
    axes[1].plot([x_avg.min(), x_avg.max()], [x_avg.min(), x_avg.max()], 'g--', lw=2)
    axes[1].set_xlim(x_avg.min(), x_avg.max())
    axes[1].set_ylim(y_avg.min(), y_avg.max())
    axes[1].plot([x_avg.min(), x_avg.max()], [x_avg.min() - delta, x_avg.max() - delta], 'r--', lw=1)
    axes[1].plot([x_avg.min(), x_avg.max()], [x_avg.min() + delta, x_avg.max() + delta], 'r--', lw=1)
    within_delta_avg = np.sum((y_avg >= x_avg - delta) & (y_avg <= x_avg + delta))
    total_points_avg = len(x_avg)
    percentage_within_delta_avg = (within_delta_avg / total_points_avg) * 100 if total_points_avg > 0 else 0
    within_half_delta_avg = np.sum((y_avg >= x_avg - 0.5 * PH2KCAL) & (y_avg <= x_avg + 0.5 * PH2KCAL))
    percentage_within_half_delta_avg = (within_half_delta_avg / total_points_avg) * 100 if total_points_avg > 0 else 0
    within_quarter_delta_avg = np.sum((y_avg >= x_avg - 0.25 * PH2KCAL) & (y_avg <= x_avg + 0.25 * PH2KCAL))
    percentage_within_quarter_delta_avg = (within_quarter_delta_avg / total_points_avg) * 100 if total_points_avg > 0 else 0
    axes[1].text(0.05, 0.90, f"Within ±0.5 pH unit: {percentage_within_half_delta_avg:.2f}%", transform=axes[1].transAxes,
                 fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    axes[1].text(0.05, 0.85, f"Within ±0.25 pH unit: {percentage_within_quarter_delta_avg:.2f}%", transform=axes[1].transAxes,
                 fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    axes[1].text(0.05, 0.95, f"Within ±1 pH unit: {percentage_within_delta_avg:.2f}%", transform=axes[1].transAxes,
                 fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    axes[1].grid(True)

    plt.tight_layout()
    plt.show()
