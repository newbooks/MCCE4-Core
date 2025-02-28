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

# Output file names
STEP1_OUT = "step1_out.pdb"
STEP1_HEAD = "head2.lst"
STEP2_OUT = "step2_out.pdb"


# This is for split_altloc()
# A backbone atom should have residue name AND name match the following to be considered as backbone
BACKBONE_ATOMS = {" N  ", " CA ", " C  ", " O  "}
RESIDUE_NAMES = {"ALA", "ARG", "ASN", "ASP", "CYS", "CYD", "GLN", "GLU", "GLY",
                 "HIL", "HIS", "ILE", "LEU", "LYS", "MET", "MEL", "PHE", "PRO",
                 "SER", "THR", "TRP", "TYR", "VAL"}

# This list what names qualify terminal residues and what names qualify amino acids
# They are used in make_ter_residues() in Protein class
TERMINAL_RESIDUES = {"NTR", "NTG", "CTR"}
AMINO_ACIDS = {"ALA", "ARG", "ASN", "ASP", "CYS", "CYD", "GLN", "GLU", "GLY",
               "HIL", "HIS", "ILE", "LEU", "LYS", "MET", "MEL", "PHE", "PRO",
               "SER", "THR", "TRP", "TYR", "VAL"}
NTR_ATOMS = {" N  ", " CA "}  # heavy atoms only
CTR_ATOMS = {" C  ", " O  ", " OXT"}  # heavy atoms only
