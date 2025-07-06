# Modeling Electrostatic Potential, RXN Terms, and Pairwise Interactions Using Local Density

## Local Density

### What is Local Density?
Local density measures how deeply an atom is buried by counting the number of surrounding atoms.
Local density also includes "density variance" which is defined as the standard deviation of atom counts in 8 quadrants.

### How is Local Density Calculated?
A bounding sphere is defined with Near, Mid, and Far radii. Atoms within each sphere are counted, resulting in a tuple:

LocalDensity = (near_score, mid_score, far_score)

This calculation is efficiently performed using the KDTree algorithm.

### Limitations of Local Density
The calculation does not account for atom radii. If artificial atoms are used to represent low dielectric constant regions, local density may not accurately reflect the environment. To address this, use reasonably sized atoms (radius 1.2 to 2.0) to fill these regions.

## Modeling Electrostatic Potential

In a medium with two dielectric constants, electrostatic potential is typically determined by solving the Poisson-Boltzmann (PB) equation, which is computationally intensive.

This project aims to use machine learning to fit PB-calculated potentials using features derived from local density, providing a faster method with acceptable accuracy.

### Preparation
Local density calculations require single conformer structures. To generate these, run:
```
step2.py --writepdb
```
The resulting structures are saved in the `ga_output` folder as `state_????.pdb`.

### Modeling Pairwise Interactions

To predict PB_Potential:
- Assume a point charge of +1; the model predicts a scaling factor. For other charges, multiply by q1 and q2.
- Random Forest Regressor is used.
- Input features:
    - LocalDensity with Near, Mid and Far scores
- Target:
    - PB_Potential (from PB solver) scaled down by distance

**Training Steps:**
- Select a microstate from `ga_output`
- Prepare microstate pdb files small_protein.pdb, medium_protein.pdb, large_protein.pdb
- Run `ele_setup.py microstate.pdb` to assign charge, calculate local density, run delphi to get PBPotential, and compile to a csv file
- Combine CSV files from multiple proteins for a comprehensive training set
- Train the model: `ele_training.py pairwise_data.csv`
- Some feature engineering might be necessary. For example, by the nature of the PB potential, the closer the atoms are, the stronger the interaction scalable by the inverse of distance. So using inversed Distance make KNN model better but won't affect decision tree based models much. KNN model is based on the the average of the nearest neighbors, and serves as a good reference of feature selection.

**Validation:**
- Does PB_Potential scale with the charge? This model relies on the linear replationship between Coulomb Potential and PB Potential.
- Internal validation should show good performance.
- Apply the trained model (pkl file) to new proteins using compiled CSV files from these proteins.

### Modeling RXN

**Training:**
- Amino acid set: Prepare side chains, place +1 charge on each atom, and calculate reaction field energy with Delphi.
- Small, medium, and large protein sets: For each, place a charge on one atom at a time and calculate reaction field energy with PB solver.
- These datasets are compiled into one file and used to train the RXN modifier. 
- The training and prediction takes out charge, and the actual RXN value should be quadratically scaled by the charge.

Approximation (excluding mutual polarization) for multi-charge systems:

For two charges (q1, q2) with no mutual polarization:
- In a uniform low dielectric, energy is Coulomb_Potential(q1, q2).
- In aqueous solution, energy is RXN(q1) + RXN(q2) + PB_Potential(q1, q2).
- The reaction field energy should be the difference between the above situcations:

RXN(q1, q2) = [RXN(q1) + RXN(q2) + PB_Potential(q1, q2)] - [Coulomb_potential(q1, q2)]

**Validation:**
- Validate charge scaling (should be quadratic).
- Test single charge predictions on proteins outside the training set.
- Validate multi-charge reaction field energy at residue and protein levels.
- Estimate the maximum mutual polarization effect on RXN. Theoretical maximum is 0 and 4 RXN when two charges are superimposed, while summing individual RXN gives 2 RXN.

### Limitations of This Modeling Approach
- Reaction field induced by a point charge only applies to that charge; mutual polarization is ignored.
