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
The goal of this step is to identify potential conformers using qualitative searching algorithms, such as Genetic Algorithms (GA).

The following is a modified GA approach:

- Place hydrogen atoms and generate ionization conformers based on CONFLIST conftypes.
- GA flow (Pool size: 1000, GA_preserver = 0.2):
    - Initialization: generate N random microstates.
        - Divide the protein into fixed residues and flipper residues.
        - Flipper residues (conf[0] + conf[1]) will form microstates.
        - Flipper residue side chains can only rotate during mutation.
    - Fitness evaluation: assess individual fitness and population fitness average (PFA80 score). Individuals are placed in an ordered list.
    - Selection: to avoid remaking connectivity, atoms are changed but never created during the GA cycle.
        - Use tournament selection to randomly pick a group of 5 individuals from the pool, and the top two individuals will be used for crossover.
        - Perform double crossover and single crossover alternately for a total of 20 times in a period.
        - Perform one mutation from a randomly selected individual in one period.
    - Replacement:
        - After each crossover or mutation, calculate the fitness score and insert the new individual in the correct position in the ordered list.
    - Termination condition check:
        - After each period, calculate PFA80.
        - Check for convergence (rolling average PFA80 score of the last 5 periods does not improve) and maximum periods.
    - If the termination condition is met, return the top 50% of the population.
- in GA, pH and Eh are not part of energy calculation so that all conftypes have a equal chance to be sampled.

Note:
- PFA80: top 80% of the population to exclude the bottom performers.
- The atoms in a residue will be mapped to an array, and connect12 will use the index number of the atom in the array. This way, the connectivity will be preserved when propagating new individuals without remaking the connect12 records.

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