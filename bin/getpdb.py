#!/usr/bin/env python

"""
Tool file: getpdb (uses mcce4.protinfo)
Possibly download the bioassembly of one or more proteins, else download the standard pdb file.
"""

import argparse
import requests


if __name__ == "__main__":
    helpmsg = "MCCE4 Step 1: Read PDB file and generate a mcce protein object"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("pdbid", nargs="+", default=[], help="Specify the pdb ID(s), e.g.: 1ots 4lzt 1FAT")
    args = parser.parse_args()

    pdbids = [id.lower() for id in args.pdbid]

    for pdbid in pdbids:
        url_rscb = "https://files.rcsb.org/download/" + pdbid + ".pdb"
        r = requests.get(url_rscb)
        if r.status_code == 200:
            with open(pdbid + ".pdb", "w") as f:
                f.write(r.text)
            print("Downloaded " + pdbid + ".pdb")
        else:
            print("Failed to download " + pdbid + ".pdb")

