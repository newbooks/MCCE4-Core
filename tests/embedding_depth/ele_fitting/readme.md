# Modeling Electrostatic Potential Using Coulomb's Law and Embedding Depth

Step 0. Calculate the embedding score
```
embedding_score.py pdb
```

Step 1. Prepare step2_out.pdb by ele_setup.py
ele_setup.py will prepare a file from step2_out.pdb in which:
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
- distance
- embedding score
- internal epsilon
- external epsilon 
- Columbs potential
- delphi calculated potential