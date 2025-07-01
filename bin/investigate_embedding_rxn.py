#!/usr/bin/env python
"""
Debugging script to investigate the embedding reaction field energy prediction
"""

import argparse
import logging
import pandas as pd


def parse_arguments():
    helpmsg = "Investigate embedding reaction field energy prediction"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("rxn_csvfile", default="", help="Reaction CSV file to be processed")
    parser.add_argument("-t", default=1, type=int, help="Tolerance for filtering of local density")
    parser.add_argument("-d", default=5.0, type=float, help="Difference of PBRXN with similar local density")
    return parser.parse_args()

def investigate_embedding_rxn(rxn_csvfile, tolerance, difference):
    """
    Investigate the embedding reaction field energy prediction.
    """
    logging.info(f"Reading reaction CSV file: {rxn_csvfile}")
    df = pd.read_csv(rxn_csvfile)
  
    if 'PBRXN' not in df.columns:
        logging.error("PBRXN column not found in the CSV file.")
        return
     
    # Display first few rows
    logging.info("First few rows of the DataFrame:")
    logging.info(df.head())


    # Use the first 3 columns as local density features
    local_density_features = df.columns[:3]

    # Find local density values that are similar, but local density are different
    similar_density = df.groupby(list(local_density_features)).filter(
        lambda x: len(x) > 1 and x['PBRXN'].std() > difference
    )

    logging.info(f"Similar local density groups with different PBRXN:")
    print(similar_density)



if __name__ == "__main__":
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Parse command line arguments
    args = parse_arguments()
    
    if not args.rxn_csvfile:
        logging.error("No reaction CSV file provided.")
        exit(1)
    
    # Investigate the embedding reaction field energy prediction
    investigate_embedding_rxn(args.rxn_csvfile, args.t, args.d)
