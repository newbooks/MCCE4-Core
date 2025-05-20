#!/usr/bin/env python
"""
MCCE4 Tool: Compile Electrostatic Energy
Description: Compile electrostatic energy for a given protein structure (microstate)
This script does the Step 3 of the ele fitting

Step 0. Calculate the embedding score

Step 1. Prepare step2_out.pdb by ele_setup.py
delphi_ele.py will prepare a file from step2_out.pdb in which:
1. the atom radius will be set to the values used by embedding depth score calculation.
2. the atom charge will be set to +1 per conformer on one atom, so that the opp reports atom to atom ele

Step 2. Run delphi
After step2_out.pdb is ready, we will switch to MCCE4 to use step3.py run delphi
```
step3.py -s delphi --debug
```

Step 3. Compile ele from delphi and embedding score to a csv file
ele_compile.py will grab information from /tmp/pbe files and embedding score to make a csv file

Columns:
- distance between atom 1 and atom 2
- embedding score atom 1
- embedding score atom 2
- internal epsilon
- external epsilon 
- Columbs potential
- electrostatic energy calculated by delphi
"""