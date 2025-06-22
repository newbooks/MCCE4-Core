#!/usr/bin/env python
"""
MCCE4 Tool: Validate Scaling Factors
Coulomb potential is propotional to the charge, we need to valid PB potential would be propotional to the charge as well.
If this is the case, we only need to scale the Coulomb potential by a factor to get the PB potential.
"""

import logging
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Constants
PDB_Template = """
ATOM    571  N   PHE A0038_000   6.547  -0.872  -7.124   1.550       0.000      BKO000_000
ATOM    572  CA  PHE A0038_000   5.411  -1.341  -7.906   1.700       0.000      BKO000_000
ATOM    573  C   PHE A0038_000   5.366  -2.842  -8.068   1.700       0.000      BKO000_000
ATOM    574  O   PHE A0038_000   4.648  -3.356  -8.936   1.520       0.000      BKO000_000
ATOM    575  H   PHE A0038_000   6.343  -0.435  -6.146   1.200       0.000      BKO000_000
ATOM    576  HA  PHE A0038_000   4.545  -1.047  -7.312   1.200       0.000      BKO000_000
ATOM    577  CB  PHE A0038_001   5.299  -0.663  -9.277   1.700       0.000      01G000_000
ATOM    578  CG  PHE A0038_001   5.205   0.831  -9.280   1.700       0.000      01G000_000
ATOM    579  CD1 PHE A0038_001   4.478   1.562  -8.355   1.700       0.000      01G000_000
ATOM    580  CD2 PHE A0038_001   5.759   1.546 -10.308   1.700       0.000      01G000_000
ATOM    581  CE1 PHE A0038_001   4.289   2.911  -8.457   1.700       0.000      01G000_000
ATOM    582  CE2 PHE A0038_001   5.580   2.900 -10.389   1.700       0.000      01G000_000
ATOM    583  CZ  PHE A0038_001   4.841   3.597  -9.503   1.700       0.000      01G000_000
ATOM    584  HB2 PHE A0038_001   4.406  -1.051  -9.766   1.200       0.000      01G000_000
ATOM    585  HB3 PHE A0038_001   6.180  -0.940  -9.856   1.200       0.000      01G000_000
ATOM    586  HD1 PHE A0038_001   4.037   1.037  -7.508   1.200       0.000      01G000_000
ATOM    587  HD2 PHE A0038_001   6.348   1.032 -11.068   1.200       0.000      01G000_000
ATOM    588  HE1 PHE A0038_001   3.700   3.439  -7.707   1.200       0.000      01G000_000
ATOM    589  HE2 PHE A0038_001   6.057   3.440 -11.207   1.200       0.000      01G000_000
ATOM    590  HZ  PHE A0038_001   4.689   4.670  -9.616   1.200       0.000      01G000_000
ATOM    591  N   ASN A0039_000   6.041  -3.574  -7.193   1.550       0.000      BKO000_000
ATOM    592  CA  ASN A0039_000   6.135  -5.036  -7.296   1.700       0.000      BKO000_000
ATOM    593  C   ASN A0039_000   5.133  -5.683  -6.322   1.700       0.000      BKO000_000
ATOM    594  O   ASN A0039_000   5.258  -5.563  -5.113   1.520       0.000      BKO000_000
ATOM    595  H   ASN A0039_000   6.545  -3.072  -6.367   1.200       0.000      BKO000_000
ATOM    596  HA  ASN A0039_000   5.901  -5.366  -8.308   1.200       0.000      BKO000_000
ATOM    597  CB  ASN A0039_001   7.554  -5.443  -6.928   1.700       0.000      01G000_000
ATOM    598  CG  ASN A0039_001   7.796  -6.898  -7.132   1.700       0.000      01G000_000
ATOM    599  OD1 ASN A0039_001   6.930  -7.719  -7.350   1.520       0.000      01G000_000
ATOM    600  ND2 ASN A0039_001   9.002  -7.313  -7.010   1.550       0.000      01G000_000
ATOM    601  HB2 ASN A0039_001   7.726  -5.203  -5.879   1.200       0.000      01G000_000
ATOM    602  HB3 ASN A0039_001   8.252  -4.881  -7.549   1.200       0.000      01G000_000
ATOM    603 HD21 ASN A0039_001   9.225  -8.376  -7.102   1.200       0.000      01G000_000
ATOM    604 HD22 ASN A0039_001   9.805  -6.603  -6.814   1.200       0.000      01G000_000
"""

CHARGE_SOURCE = "ATOM    577  CB  PHE A0038_001   5.299  -0.663  -9.277   1.700       0.000      01G000_000"      # The source of charge for the potential calculation, will change in range [-1, 1]
POTENTIAL_SITE = "ATOM    598  CG  ASN A0039_001   7.796  -6.898  -7.132   1.700       0.000      01G000_000"     # The site of potential to be reported, maintain constant charge 1
STEP3_EXE = "/home/jmao/projects/MCCE4/bin/step3.py"

def parse_arguments():
    helpmsg = "Validate scaling electrostatic potential from the source charge."
    parser = argparse.ArgumentParser(description=helpmsg, formatter_class=argparse.RawTextHelpFormatter)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    logging_format = "%(asctime)s %(levelname)s: %(message)s"
    datefmt = '%Y-%m-%d %H:%M:%S'

    # Set up logging
    logging.basicConfig(format=logging_format, datefmt=datefmt, level=logging.INFO)

    # Create a PDB file with the template
    pdb_content = PDB_Template.strip()
    pdb_lines = pdb_content.split('\n')

    # Get the source charge created potential opp file
    confname = CHARGE_SOURCE[17:20] + CHARGE_SOURCE[80:82] + CHARGE_SOURCE[21:30]
    potential_oppfile = f"energies/{confname}.opp"
    # Get the conformer name to extract the potential from potential_oppfile
    potential_confname = POTENTIAL_SITE[17:20] + POTENTIAL_SITE[80:82] + POTENTIAL_SITE[21:30]
    pdb_filename = "step2_out.pdb"

    # Set the portential site charge to 1.0
    for i, line in enumerate(pdb_lines):
        if POTENTIAL_SITE[:30] in line: # only use the first 30 characters to match the potential site
            pdb_lines[i] = line[:62] + "%12.3f" % (1.0) + line[74:]  # Set charge to 1.0
            break

    # Modify the charge source to a range of values from -1.0 to 1.0
    charge_potential_list = []  # List to store (charge, potential) tuples
    n = 10  # steps in one unit of charge, e.g. -1.0 to 1.0 with 10 steps will generate 21 charges
    for crg in range(-n, n+1, 1):
        # Set the charge in pdb_content and write to a file
        charge_value = crg / n
        for i, line in enumerate(pdb_lines):
            if CHARGE_SOURCE[:30] in line: # only use the first 30 characters to match the charge source
                pdb_lines[i] = line[:62] + "%12.3f" % (charge_value) + line[74:]
                break

        # Save PDB content to "step2_out.pdb"
        pdb_content = "\n".join(pdb_lines) + "\n"
        open(pdb_filename, 'w').write(pdb_content)
        logging.info(f"Generated PDB file: {pdb_filename} with charge {charge_value:.2f}")

        # Call mcce step3 to calculate the potential
        logging.info(f"Running MCCE step3 with charge {charge_value:.2f}...")
        result = subprocess.run([STEP3_EXE, "-s", "delphi"], capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f"step3.py failed with exit code {result.returncode}")
            logging.error(f"stderr: {result.stderr.strip()}")

        # Read the output potential from the generated file
        potential_lines = open(potential_oppfile, 'r').readlines()
        potential_value = 0.00
        for line in potential_lines:
            if potential_confname in line:
                fields = line.strip().split()
                if len(fields) >= 6 and fields[1] == potential_confname:
                    potential_value = float(fields[5])

        charge_potential_list.append((charge_value, potential_value))

        print(f"Charge: {charge_value:.2f}, Potential: {potential_value:.4f}")

    # Write the charge and potential values to a csv file
    df = pd.DataFrame(charge_potential_list, columns=['Charge', 'Potential'])
    df.to_csv('charge_potential.csv', index=False)
    logging.info("Charge and potential values saved to charge_potential.csv")
    
    # Plot the charge vs potential
    plt.figure(figsize=(10, 6))
    plt.title('Charge vs PB_Potential')
    plt.xlabel('Charge')
    plt.ylabel('PB_Potential')
    plt.grid(True)
    # add a linear fit line
    sns.scatterplot(x='Charge', y='Potential', data=df, marker='X', color='black', s=100)
    sns.regplot(x='Charge', y='Potential', data=df, scatter=False, color='green', line_kws={'lw': 2})
    plt.savefig('charge_vs_potential.png')
    plt.show()