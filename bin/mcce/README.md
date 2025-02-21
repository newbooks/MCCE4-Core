# MCCE Modules

This folder contains MCCE modules. Executables and tools can import the MCCE modules using:
```
from mcce import *
```
or 
```
from mcce.pdbio import *
```
to import functions and classes from the pdbio module.

## Importable Modules
These are the modules that can be imported by executables and tools. The variable `__all__` in file `__init__.py` defines which files can be imported using the wildcard `*`.

Importable modules should not have a leading underscore ("_"), as this is reserved for private modules.

The importable modules are categorized by their functions:
- pdbio.py: Protein and parameter input/output
- geom.py: Geometry operations
- constants.py: Constants and default values

## Private Modules
File names starting with an underscore ("_") are private modules. These modules are imported and used by importable modules and are hidden from regular developers and users.