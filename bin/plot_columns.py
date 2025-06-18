#!/usr/bin/env python
"""
MCCE4 Tool: Plot Columns from space separated or CSV files
"""
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
import numpy as np

def parse_arguments():
    helpmsg = "Plot two columns in one CSV file.\n"
    parser = argparse.ArgumentParser(description=helpmsg, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("data_file", help="Data file (CSV or space separated) containing data for plotting")
    # specify the column name for x and y axis, default is empty string which maches nothing and will prompt the user to enter the column index
    parser.add_argument("--colx", default="", help="Column name for x-axis")
    parser.add_argument("--coly", default="", help="Column name for y-axis")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    # read the csv data file, assuming the first row is column names
    print(f"Reading data from {args.data_file} and the data looks like...")
    data = pd.read_csv(args.data_file)
    print(f"Data from {args.data_file}:")
    print(data.head())    

    # if colx or coly is not specified, prompt the user to enter the column index
    if not args.colx:
        print("Available columns for x-axis:")
        print(data.columns)
        colx = input("Enter the column name for x-axis: ")
    else:
        colx = args.colx
    if not args.coly:
        print("Available columns for y-axis:")
        print(data.columns)
        coly = input("Enter the column name for y-axis: ")
    else:
        coly = args.coly

    # check if the columns exist in the data
    if colx not in data.columns or coly not in data.columns:
        raise ValueError(f"Column '{colx}' or '{coly}' does not exist in the data file.")
    
    # create a scatter plot to compare the two columns
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=data, x=colx, y=coly, alpha=0.7)
    plt.title(f"Column Comparison: {colx} vs {coly}")
    plt.xlabel(colx)
    plt.ylabel(coly)
    plt.grid()
    plt.tight_layout()
    plt.show()
    plt.close()