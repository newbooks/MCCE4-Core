# Parameter Directory
This directory contains structure parameters for MCCE.
- `ftpl_[name]`: Subdirectories with structure parameters in free topology (ftpl) format.
- `ideal_structures`: Subdirectory with molecule (residue and cofactor) structure templates.
- `extra.ftpl`: File with extra energy terms and scaling factors for energy adjustment.

## Special Handling of runprm Entries
In the `run.prm.default` file, these entries specify important file/folder locations relative to the distribution:
- `FTPL_FOLDER`
- `EXTRA`
- `RENAME_RULES`
- `DELPHI_EXE`
- `APBS_EXE`
Absolute paths can also be used in the default and additional runprm files.

The `Runprm` class will have absolute paths for these attributes, prefixed with an underscore.