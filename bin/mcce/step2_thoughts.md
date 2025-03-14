# Step 2 Conformer Making Strategy

## MCCE program flow
1. Step 1: preprocess the input structure
2. Step 2: make conformers
3. Step 3: pre-calculate energy look-up table
4. Step 4: sample the microstates and link the microstates to protein's biophysical properties

## Current Strategy for Step 2: Conformer Making
Evaluating microstates is computationally expensive, so selecting the right conformer candidates in step 2 is crucial for the entire MCCE workflow.

The aim of step 2 is to quickly and inexpensively generate discrete conformers, rotamers, and ionization conformers that are likely to be accessed in step 4 which uses an accurate and costly force field. To achieve this, several strategies are employed:
- Utilize simplified force fields, such as Coulomb's law and a simplified van der Waals (VDW) force.
- Apply a divide-and-conquer approach: create conformers by reducing heavy atom clashes, searching for hydrogen bond pairs, generating the most exposed conformers, and then optimizing hydrogen atoms.
- Use clustering to select only representative conformers, further reducing the number of conformers.

While this method is relatively efficient, it still produces some inaccessible conformers and may miss low-energy conformers. One possible reason is the lack of "microstate awareness" during the process. Conformers for each residue are created considering only self-clashes, solvent exposure, and potential hydrogen bonds with neighbors. In reality, residues choose their sidechain ionization and motion in a concerted manner, which can only be captured by a method that considers the entire picture of microstates.

## Alternative Strategy