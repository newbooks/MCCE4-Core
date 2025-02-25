"""
Module for constants used in MCCE.
"""

# File name for run.prm.default, path from the distribution root
RUNPRM_DEFAULT = "config/run.prm.default"   # path from the distribution root
RUNPRM_DUMP = "run.prm.record"              # path from the working directory
FTPL_DUMP = "ftpl.record"                   # path from the working directory

# Special entries in run.prm that should be converted to absolute paths if seen as relative paths
RUNPRM_SPECIAL_ENTRIES = ["FTPL_FOLDER", "EXTRA", "RENAME_RULES", "DELPHI_EXE", "APBS_EXE"]

# Unit conversion factors
ROOMT = 298.15
PH2KCAL = 1.364
KCAL2KT = 1.688
