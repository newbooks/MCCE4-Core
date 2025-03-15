# Step 2 Conformer Making Strategy

## MCCE Program Flow
1. Step 1: Preprocess the input structure.
2. Step 2: Generate conformers.
3. Step 3: Pre-calculate the energy look-up table.
4. Step 4: Sample microstates and link them to the protein's biophysical properties.

## Current Strategy for Step 2: Conformer Generation
Evaluating microstates is computationally intensive, making the selection of appropriate conformer candidates in step 2 critical for the MCCE workflow.

The goal of step 2 is to efficiently generate discrete conformers, rotamers, and ionization states that are likely to be relevant in step 4, which uses a precise and computationally expensive force field. To achieve this, several strategies are employed:
- Utilize simplified force fields, such as Coulomb's law and a simplified van der Waals (VDW) force.
- Apply a divide-and-conquer approach: create conformers by reducing heavy atom clashes, searching for hydrogen bond pairs, generating the most exposed conformers, and then optimizing hydrogen atoms.
- Use clustering to select representative conformers, reducing the overall number of conformers.

While this method is relatively efficient, it can still produce some inaccessible conformers and may miss low-energy conformers.

One possible reason is the lack of "microstate awareness" during the process. Conformers for each residue are created considering only self-clashes, solvent exposure, and potential hydrogen bonds with neighbors. In reality, residues choose their sidechain ionization and motion in a concerted manner, which can only be captured by a method that considers the entire picture of microstates.

Another possible cause of missing good conformers is the discrete nature of the rotamer generation in step 2. Although the search includes stepped refinement of rotamers, the initial search uses a predefined 60-degree step, which may not cover all potential rotation angles.

## Alternative Strategy
### Force Field
Consider developing a set of quick empirical functions to account for all the forces/energies used in step 4 and apply them in step 2. These force functions don't have to be perfect.

vdw: empirical function based on distance^2 
ele: Coulomb's law
desolvation: function of SAS and charge
torsion: 1-4 vdw
hydrogen bond: angle and distance
pH and Eh energy: environment pH/Eh and pK0, Eh0


### Sampling Method
The goal of this step is to identify potential conformers, so we can use qualitative searching algorithms, such as Genetic Algorithms.

- Place H and make ionization conformers
- GA flow
   - Initialization: generate N random microstates.
   - Fitness evaluation: individual fitness, and population fitness average (PFA score)
   - Selection
   - Crossover: one-point crossover : two-point crossover = 1:1
   - Mutation: dynamically adjusted from 1% to 5% based on rolling PFA score average
   - Replacement
   - Termination condition check: convergence check (rolling average PFA score of the last 5 generations does not imporve) and max generation
   - If termination condition is met, return the top 50% population
- Do GA for pH = 4, 7, and 10

Note: PFA score variations
- PFA: whole population
- PFA90: top 90% population to exclude the buttom performers
- PFA50: top 50% population which is what we will use as final conformers