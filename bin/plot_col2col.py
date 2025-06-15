#!/usr/bin/env python
"""
MCCE4 Tool: Plot Column to Column Comparison"""

import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import numpy as np

def parse_arguments():
    helpmsg = "Plot column comparison from two files CSV or space separated values file.\n"
    parser = argparse.ArgumentParser(description=helpmsg, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("data_file1", help="First data file (CSV or space separated) containing data for comparison")
    parser.add_argument("data_file2", help="Second data file (CSV or space separated) containing data for comparison")
    parser.add_argument("--col", default=0, type=int, help="Column number to compare, positive number counts from left, negative counts from right, default is 0 (first column)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    # read the data files and use the white space or comma as separator. The files have no column header
    # allow multiple separators by using sep=None and engine='python'
    # not all fields are numbers, so we read them as strings and convert to numeric later
    lines = open(args.data_file1).readlines()
    data1 = [line.strip().split() for line in lines]
    lines = open(args.data_file2).readlines()
    data2 = [line.strip().split() for line in lines]
    
    # determine the column index to compare
    col_index = args.col
    if col_index < 0:
        col_index = max(len(data1[0]), len(data2[0])) + col_index  # convert negative index to positive index
    elif col_index >= len(data1[0]) or col_index >= len(data2[0]):
        raise ValueError(f"Column index {args.col} is out of range for the data files.")
    # convert the selected column to numpy arrays
    col_data1 = np.array([float(row[col_index]) for row in data1 if len(row) > col_index])
    col_data2 = np.array([float(row[col_index]) for row in data2 if len(row) > col_index])

    # display the data
    print(f"Data from {args.data_file1} (Column {col_index}):")
    print(col_data1)
    print(f"Data from {args.data_file2} (Column {col_index}):")
    print(col_data2)
    
    # create a scatter plot to compare the two columns
    plt.figure(figsize=(8, 6))
    plt.scatter(col_data1, col_data2, alpha=0.7)
    plt.title(f"Column {col_index} Comparison")
    plt.xlabel(f"{args.data_file1} (Column {col_index})")
    plt.ylabel(f"{args.data_file2} (Column {col_index})")
    plt.grid()
    plt.show()
    plt.close()
