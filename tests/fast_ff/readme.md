# Modeling Electrostatic Potential Using Coulomb's Law and Embedding Depth

1. Step 1: Microstate pdb files

Before modeling, you need to have microstate pdb files. A microstate file has the same format as MCCE step2_out.pdb, except only one side chain conformer instead of multiple conformers is given to each residue.

You can create such microstate pdb file using
```
step2.py --writepdb
```
The microstate pdb files are in folder "ga_output".


2. Step 2: Local density score
Use this command to calculate the local density score on a selected microstate pdb file
```
local_density.py microstate_pdb
```

3. Step 3: step2_out.pdb
This step prepares step2_out.pdb by ele_setup.py.
```
ele_setup.py microstate_pdb
```

This script alters the atom radius and charge so 
- the redius matches the values in embedding score calculation
- the charge is +1 for one random atom in one side chain

This way the atom to atom pairwise interaction is the same as conformer to conformer interaction reported by MCCE step3.

4. Step 4: Run delphi using step3.py.

When step2_out.pdb is ready, we will switch to MCCE4 to use step3.py run delphi
```
step3.py -s delphi
```

5. Step 5: Compile ele from delphi and embedding score to a csv file
At this point, you should have an energies directory with opp files. The ele_compile.py will grab information from opp files and embedding score to make a csv file.

```
ele_compile microstate
```
Columns in CSV file:
- Conf1
- Conf2
- Distance
- Radius1
- Radius2
- Density1_Near
- Density2_Near
- Density1_Mid
- Density2_Mid
- Density1_Near
- Density2_Near
- CoulombPotential
- PBPotential

6. Modeling, compare models
Use script:
- ele_fitting.py

7. Training and Prediction
We use small(4pti), medium(4lzt) and large(2rcr) proteins to train a combined model and evaluate the performance of this combined model.
```
mkdir small
```