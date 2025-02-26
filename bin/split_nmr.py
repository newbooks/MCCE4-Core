#!/usr/bin/env python

"""
Split NMR models into individual PDB files
"""

import os
import logging
import argparse

logging_format = "%(asctime)s %(levelname)s: %(message)s"
logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
logging_datefmt='%Y-%m-%d %H:%M:%S'


if __name__ == "__main__":
    helpmsg = "Split NMR models into individual PDB files"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("pdb_files", metavar="pdbs", nargs="+", default=[], help="PDB files to split")
    parser.add_argument("--debug", default=False, action="store_true", help="Print debug information")
    args = parser.parse_args()


    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format=logging_format_debug, datefmt=logging_datefmt)
    else:
        logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)

    logging.info("Split NMR models into individual PDB files")

    # Read the PDB files
    for pdb_file in args.pdb_files:
        if os.path.isfile(pdb_file):
            pdblines = open(pdb_file).readlines()

            # Split the models
            model = []
            model_number = None
            pdb_file_base = os.path.splitext(pdb_file)[0]
            fnames = []
            for line in pdblines:
                if line.startswith("MODEL"):
                    model = []
                    model_number = int(line.split()[1])
                elif line.startswith("ENDMDL"):
                    with open(f"{pdb_file_base}.model{model_number}.pdb", "w") as f:
                        f.writelines(model)
                        fnames.append(f"{pdb_file_base}.model{model_number}.pdb")
                model.append(line)
            if model_number is not None:
                with open(f"{pdb_file_base}.model{model_number}.pdb", "w") as f:
                    f.writelines(model)
                    fnames.append(f"{pdb_file_base}.model{model_number}.pdb")
                logging.info(f"Split {pdb_file} into {model_number} models: {', '.join(fnames)}")
            else:
                logging.info(f"No models found in {pdb_file}")
        else:
            logging.error(f"File \"{pdb_file}\" not found")
            exit(1)