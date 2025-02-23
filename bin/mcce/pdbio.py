"""
Module Name: pdbio

Description:
This module provides MCCE protein and parameter data structure and input/output functions.
"""

import os
import logging
import datetime
from .constants import *
from .geom import Vector


class Atom:
    """
    Atom class
    """
    def __init__(self):
        # defined by pdb file
        self.atomname = ""          # atom name
        self.altloc = ""            # alternate location indicator
        self.resname = ""           # residue name
        self.chain = ""             # chain ID
        self.sequence = 0           # residue sequence number
        self.insertion = ""         # insertion code
        self.xyz = Vector()         # coordinates
        # extended attributes
        self.r_boundary = 0.0       # boundary radius
        self.charge = 0.0           # charge
        self.r_vdw = 0.0            # van der Waals radius
        self.e_vdw = 0.0            # van der Waals energy well depth
        self.element = ""           # element name
        self.conn12 = []            # list of atoms that are 1-2 bonded
        self.conn13 = []            # list of atoms that are 1-3 bonded
        self.conn14 = []            # list of atoms that are 1-4 bonded
        self.parent_conf = None     # parent conformer


class Conformer:
    """
    Conformer class
    """
    def __init__(self):
        self.confid = ""            # conformer name, unique ID in the protein. resname+chain+sequence+insertion+confnum
        self.altloc = ""            # alternate location indicator
        self.resname = ""           # residue name
        self.chain = ""             # chain ID
        self.sequence = 0           # residue sequence number
        self.insertion = ""         # insertion code
        self.conftype = ""          # conformer type, as defined in the ftpl file
        self.confnum = 0            # conformer number
        self.atoms = []             # list of atoms in the conformer
        self.parent_residue = None  # parent residue
        self.occ = 0.0              # occupancy
        self.history = ""           # history string
        self.charge = 0.0           # net charge
        self.calculated = False     # flag for calculated conformer


class Residue:
    """
    Residue class
    """
    def __init__(self):
        self.resname = ""           # residue name
        self.chain = ""             # chain ID
        self.sequence = 0           # residue sequence number
        self.insertion = ""         # insertion code
        self.conformers = []        # list of conformers in the residue
        self.resid = ()             # residue ID, (resname, chain, sequence, insertion)


class Protein:
    """
    Protein class
    """
    def __init__(self):
        self.residues = []          # list of residues in the protein



class Runprm:
    """
    Runprm class
    This class stores the parameters from a runprm file.
    ---------------------------------------
    Example runprm file:
    ---------------------------------------
    step 1:
    t      Step1. Label terminal residues as separate NTER and CTR?     (TERMINALS)
    2.0    Step1. Distance limit for reporting clashes                  (CLASH_DISTANCE)
    0.05   Step1. Remove water with %SAS exceeding this cutoff value    (H2O_SASCUTOFF)
    ---------------------------------------
    
    The last field in parentheses is the parameter name, which will be used as a class attribute.
    The first field is the default value, stored as a string, if the last field qualifies as a parameter name.
    Anything in between is the description.

    Runprm class attributes come from multiple sources:
    1. Default values from run.prm.default in the distribution.
    2. Optionally, the values can be updated by command line options using update_by_cmd().
    3. Optionally, the values can be updated by additional runprm files using update_by_files().
    4. Optionally, the values can be updated by direct assignment.

    Runprm class has a method dump() to print the parameters in the format of a runprm file.
    """
    class Record:
        """
        Record class
        """
        def __init__(self, value, description, set_by):
            self.value = value
            self.description = description
            self.set_by = set_by

    def __init__(self):
        dist_path = os.path.abspath(os.path.join(__file__, "../../.."))
        runprm_default = os.path.join(dist_path, RUNPRM_DEFAULT)
        if os.path.exists(runprm_default):
            logging.info(f"   Loading default runprm file {runprm_default}")
            self.update_by_files([runprm_default])
        else:
            logging.warning(f"   Default runprm file {runprm_default} not found.")
        

    def update_by_files(self, files):
        """
        Update runprm object by runprm files.
        """
        for file in files:
            if not os.path.exists(file):
                logging.warning(f"   Runprm file {file} not found.")
                continue

            with open(file) as f:
                for line in f:
                    entry_str = line.split("#")[0].strip()
                    fields = entry_str.split()
                    if len(fields) > 1 and fields[-1].startswith("(") and fields[-1].endswith(")"):
                        key = fields[-1][1:-1].strip()
                        value = fields[0]
                        description = " ".join(fields[1:-1])
                        filename = os.path.basename(file)
                        set_by = f"#Set by {filename}"
                        setattr(self, key, self.Record(value, description, set_by))
            

    def dump(self, comment=""):
        """
        Dump the parameters in the format of a runprm file.
        """
        attributes = vars(self)
        date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(RUNPRM_DUMP, "w") as f:
            f.write(f"# This runprm file is recorded on {date}\n{comment}\n")
            for key, value in attributes.items():
                if key.startswith("_"):
                    continue
                value_len = len(value.value)
                if value_len < 10:
                    line = "%-10s %-60s %-16s %s\n" % (value.value, value.description, f"({key})", value.set_by)
                else:
                    description = value.description.ljust(60-value_len+10)
                    line = "%-10s %s %-16s %s\n" % (value.value, description, f"({key})", value.set_by)
                f.write(line)
        logging.info(f"   MCCE runprm parameters are recorded in file {RUNPRM_DUMP}")
                

        