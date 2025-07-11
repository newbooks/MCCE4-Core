#!/usr/bin/env python3
"""
All in one training script for fast force field by Machine Learning.
"""

import logging
import argparse
import os

# create a unique temporary directory to run the training
import tempfile

logging_format = "%(asctime)s %(levelname)s: %(message)s"
logging_datefmt='%Y-%m-%d %H:%M:%S'

if __name__ == "__main__":
    helpmsg = "Train a fast force field model from a pdb folder. This script attempts to reproduce PB solver delphi potential using machine learning techniques."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument('pdb_folder', type=str, help='Path to the pdb folder. The folder should contain proteins at different size and shape.')
    args = parser.parse_args()

    logging.basicConfig(format=logging_format, datefmt=logging_datefmt, level=logging.INFO)

    with tempfile.TemporaryDirectory(prefix="aiofff_training_") as temp_dir:
        logging.info(f"Created temporary directory at {temp_dir}")
        # go to the temporary directory
        os.chdir(temp_dir)
        # create a file named "a"
        with open("a", "w") as f:
            f.write("This is a temporary file.")
        
