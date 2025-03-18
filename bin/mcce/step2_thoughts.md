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
    - Termination condition check: convergence check (rolling average PFA score of the last 5 generations does not improve) and max generation
    - If termination condition is met, return the top 50% population
- Do GA for pH = 4, 7, and 10

Note: PFA score variations
- PFA: whole population
- PFA90: top 90% population to exclude the bottom performers
- PFA50: top 50% population which is what we will use as final conformers

The whole population average fitness score is commonly used to measure the convergence, but which score is appropriate in our application is a question.

### Conformer Making Levels
In the case of GA conformer making, the rotation angles are not limited by the rotation step size, so this mechanism can potentially create a large number of ending conformers. The control of the ending number of conformers is currently dependent on the conformer generating method, in the sense of how extensively to explore the rotational space. With GA, the control of bond rotations does not exist. However, the number of ending conformers in GA can be achieved by:
1. Population size: increasing the pool size allows more diverse conformers to be created and preseved, at the cost of linearly scaled sampling expenses. 
2. Clustering criteria: after making all the conformers, clustering selects the representative conformers based on the geometry similarity. The threshold of similarity may affect the final conformers selected and used by the later steps.


### Find Rotatable Bond for H atoms
It comes natually to mind that we can use the ROTATE records to include the rotations affecting hydrogen atoms. However, in current MCCE, we already have the knowledge of finding what hydrogen atoms have rotational freedom in the place_h module, without ROTATE records. Can we use the same mechanism to figure out these rotation bonds without touching the RECORDS in ftpl files?

The mutation in GA is done by "shake" the side chain by rotating all rotatble bonds. So we need uniform rules to make rotamers, regardless heavy atom rotamer or hydrogen atom rotamers. One solution is to
- ammend the ROTATE rules with new H atom rotate rules obtained from analysis. This means the tpl object is not pure like it currently is.
- create a new set of rotate rules from two sources: ROTATE records and derived H rotate bonds

I will use combined rotate rule solution, as it requires no modification to the ftpl files, and provides easy access to affected atoms.

Proposed structure:

```
rotate_rules = {conftype: rotate_rule
                ...
}

class RotateRule:
    def __init__(self):
        self.bond = (atom1, atom2)
        self.affected_atoms = [affected_atom1, affected_atom2, ...]
```

The atoms in rules are 4-char atom names rather than the Atom() object so the rules are universal.

In the Conformer class, the following attributes will be added.
```
class Conformer:
    def __init__(self):
        ...
        self.rotatbles = [bond1, bond2, ...]  # bond is the atom pairs pointing to the actual atom instances
        self.rotated_atoms = {  # atoms point to actual atom instances
            bond1: [atom1, atom2, ...],
            bond2: [atom3, atom4, ...],
            ...
        }

    def init_rotate(self, rotaterules)
        """ Initilize the following attributes for quick access to rules and atom instances from the generic name based rotaterules. """
        self.rotatbles = ...
        self.rotated_atoms = ...
```