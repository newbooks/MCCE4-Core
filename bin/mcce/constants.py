"""
Module for constants used in MCCE.
"""

# File name for run.prm.default, path from the distribution root
RUNPRM_DEFAULT = "config/run.prm.default"   # path from the distribution root
RUNPRM_DUMP = "run.prm.record"              # path from the working directory
FTPL_DUMP = "ftpl.record"                   # path from the working directory

# Special entries in run.prm that should be converted to absolute paths if seen as relative paths
# The absolute path will be added with a leading underscore.
RUNPRM_SPECIAL_ENTRIES = ["FTPL_FOLDER", "EXTRA", "RENAME_RULES", "DELPHI_EXE", "APBS_EXE"]

# Unit conversion factors
ROOMT = 298.15
PH2KCAL = 1.364
KCAL2KT = 1.688


# Input and Output file names
USER_PARAM = "user_param"       # User defined ftpl files in the working directory. MCCE will use these files if found.
STEP1_OUT = "step1_out.pdb"
STEP1_HEAD = "head2.lst"
NEW_FTPL = "new.ftpl"
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
AMINO_ACIDS = {"ALA", "ARG", "ASN", "ASP", "CYS", "CYD", "CYL", "GLN", "GLU", "GLY",
               "HIL", "HIS", "ILE", "LEU", "LYS", "MET", "MEL", "PHE", "PRO",
               "SER", "THR", "TRP", "TYR", "VAL"}
NTR_ATOMS = {" N  ", " CA "}  # heavy atoms only
CTR_ATOMS = {" C  ", " O  ", " OXT"}  # heavy atoms only

# Default RADIUS values for atoms in unknow cofactors 
R_BOUNDARY = {
    " H": 1.0,
    " C": 1.7,
    " N": 1.55,
    " O": 1.52,
    " S": 1.8,
    " P": 1.8,
    " F": 1.47,
    "CL": 1.75,
    "BR": 1.85,
    " I": 1.98,
    "FE": 1.8,
    "ZN": 1.39,
    "CU": 1.4,
    "MG": 1.73,
    "CA": 1.97,
    "MN": 1.61,
    "NA": 1.73,
    "CA": 2.23,
    " X": 1.8,
}
R_VDW = {   # radius and energy well depth
    " H": (1.000, 0.020),
    " C": (2.000, 0.150),
    " N": (1.750, 0.150),
    " O": (1.600, 0.150),
    " S": (2.000, 0.200),
    " P": (2.100, 0.200),
    " F": (1.500, 0.150),
    "CL": (1.750, 0.150),
    "BR": (1.950, 0.200),
    "NA": (1.750, 0.160),
    "CA": (2.000, 0.173),
    " X": (2.000, 0.173),
}
