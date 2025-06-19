# Modeling Electrostatic Potential, terms RXN and pairwise interaction using Local Density

## Local Density
### What is Local Density
Local density describes how deeply buried is an atom bu counting other atoms around it.

### How to calculate Local Density
Define a bounding sphere with Near, Mid and Far radii, and search the atoms inside the spheres. The 3 counts make a tuple to represent Local Density.

LocalDensity = (near_score, mid_score, far_score)

The calculation is efficient when using KDTree algorithm.

### Shortcoming of Local Density
The atom radius is not considered in Local Desity calculation. This means if I define an artificial atom to describe the low dielectric constant bulk of material around a point charge, it won't be able to describe the local environment of this point charge accurately. To remedy this, one has to use reasonable sized (R 1.2 to 2.0) atoms to fill up the low dielectric constant material space.


## Modeling Electrostatic Potential
In a two-dielectric constant medium, the electrostatic potential can be determined by solving Poisson Boltzmann (PB) equation. This is a computationally expensive.

Instead of solving PB equation, the task of this project is to use machine learning to train a model to fit the PB solved potential by using cheaply calculated atom local environment (local density).

The goal is to develop a much faster way with acceptatble accurancy to calculate the electrostatic potential than PB solver.

### Preparation
Local Denstity only works in sigle conformer strcutures. To create single conformer structures, run 
```
step2.py --writepdb
``` 
and the single conformer structures are stored in folder `ga_output` in the name of `state_????.pdb`

### Modeling pairwise
To predict the PB_Potential, we will assume:
- Point charge is +1, and what we calculate is the potential scaling factor. To get the actual potential of a non +1 charge pair, we will multiply this factor with q1 and q2.
- Random Forest Regressor works reasonably well.
- The input features are 
    - LocalDensity_Near
    - LocalDensity_Mid
    - LocalDensity_Far
- The fitting target is
    - PB_Potential (Potential calculated by PB solver)

**Training:**
- take a microstate from `ga_output`
- run `ele_setup.py microstate.pdb` to generate a step2_out.pdb
- run `step3.py`, this creates energies directory
- run `ele_compile.py state_name` to put pairwise interaction to a csv file.
- combine multiple csv files for a training set that covers multiple protein
- train the model with a csv file: `ele_training.py pairwise_data.csv`

**Validation:**
- The internal validation of the training should look pretty good.
- Apply the trained model (pkl file) to arbitrary protein (csv file compiled in the same way as training set)


### Modeling RXN
**Training:**
- Amino acid set: Prepare amino acid side chains and place +1 charge on each atom, as a conformer, then calculate the reaction field energy by delphi
- Small protein set: Use a small protein, place charge on one atom at a time, then calculate the reaction field energy by delphi
- Medium protein set: Use a medium protein, place charge on one atom at a time, then calculate the reaction field energy by delphi
- Large protein set: Use a large protein, place charge on one atom at a time, then calculate the reaction field energy by delphi
The above will be the training set to train RXM modifier based on local density. The actual RXN should be scled by the charge.

Approximation (by excluding mutual polarization) of multi-charge reaction field energy:

Let's consider a two charge system (q1, q2), assuming no mutual polarization:

In a uniform low dielectric constant material, the system energy is Coulomb_Potential(q1, q2) which can be calculated analytically. Now we move this two charge system to aqueous solution, the new system will have energy RXN(q1) + RXN(q2) + PB_Potential(q1, q2). The energy difference is the reaction field energy of the system treated as a whole:

RXN(q1, q2) = [system in aquesou solution] - [system in uniform dielectric constant]

=>

RXN(q1, q2) = [RXN(q1) + RXN(q2) + PB_Potential(q1, q2)] - [Coulomb_potential(q1, q2)]

=>

RXN(q1, q2) = RXN(q1) + RXN(q2) + PB_Potential(q1, q2) - Coulomb_potential(q1, q2)


**Validation:**
- Validate the charge scaling. It should be quadratic with charge.
- Validation of single charge, choose a protein outside the training set.
- Validation of multi charge reaction field energy calculation at residue and protein level.
- Estimate the maximum mutual polarization effect on RXN. In theory, two superimposed charges give the maximum number: 0 and 4 RXN, while summing up the individual RXN gives 2 RXN. 

### Shortcomings of this modeling
- The point charge induced reaction field only applies to the point charge. The mutual polarization is ommited

