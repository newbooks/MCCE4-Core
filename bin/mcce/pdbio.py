"""
Module Name: pdbio

Description:
This module provides MCCE protein and parameter data structure and input/output functions.
"""

import sys
import os
import logging
import datetime
import copy
import glob
from collections import defaultdict

from .constants import *
from .geom import Vector


class Atom:
    """
    Atom class
    """
    def __init__(self):
        # defined by pdb file
        self.record = ""            # record name "ATOM" or "HETATM"
        self.serial = 0             # atom serial number
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

    def residue_id(self):
        """
        Return the residue ID.
        """
        return (self.resname, self.chain, self.sequence, self
                .insertion)
    
    def clone(self):
        """
        Clone the atom.
        """
        new_atom = Atom()
        new_atom.record = self.record
        new_atom.serial = self.serial
        new_atom.atomname = self.atomname
        new_atom.altloc = self.altloc
        new_atom.resname = self.resname
        new_atom.chain = self.chain
        new_atom.sequence = self.sequence
        new_atom.insertion = self.insertion
        new_atom.xyz = self.xyz.copy()
        new_atom.r_boundary = self.r_boundary
        new_atom.charge = self.charge
        new_atom.r_vdw = self.r_vdw
        new_atom.e_vdw = self.e_vdw
        new_atom.element = self.element
        new_atom.conn12 = []    # do not clone connections
        new_atom.conn13 = []    # do not clone connections
        new_atom.conn14 = []    # do not clone connections   
        new_atom.parent_conf = self.parent_conf
        return new_atom

    def load_pdbline(self, line):
        """
        Load the atom from a pdb line.
        """
        self.record = line[0:6]
        self.atomname = line[12:16]
        self.altloc = line[16]
        self.resname = line[17:20]
        self.chain = line[21]
        self.sequence = int(line[22:26])
        if line[26] == " ":
            self.insertion = "_"
        else:
            self.insertion = line[26]
        self.xyz = Vector([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        #self.element = line[76:78] # Use atom name to detemin element so it covers non standard PDB files
        if len(self.atomname.strip()) == 4 and self.atomname[0] == "H":
            self.element = " H"
        else:
            self.element = self.atomname[:2]

    def as_mccepdb_line(self):
        """
        Return the atom as a pdb line.
        """
        serial = self.serial % 100000
        return "ATOM  %5d %4s %3s %c%04d%c%03d%8.3f%8.3f%8.3f%8.3f%12.3f      %s\n" % (
            serial,
            self.atomname,
            self.resname,
            self.chain,
            self.sequence,
            self.insertion,
            self.parent_conf.confnum,
            self.xyz.x,
            self.xyz.y,
            self.xyz.z,
            self.r_boundary,
            self.charge,
            self.parent_conf.history
        )


class Conformer:
    """
    Conformer class
    """
    def __init__(self):
        self.confid = ""            # conformer name, unique ID in the protein. resname+chain+sequence+insertion+confnum
        self.conftype = ""          # conformer type, as defined in the ftpl file
        self.confnum = 0            # conformer number
        self.parent_residue = None  # parent residue
        self.occ = 0.0              # occupancy
        self.history = ""           # history string
        self.charge = 0.0           # net charge
        self.calculated = False     # flag for calculated conformer
        self.atoms = []             # list of atoms in the conformer

    def clone(self):
        """
        Clone the conformer.
        """
        new_conf = Conformer()
        new_conf.confid = self.confid
        new_conf.conftype = self.conftype
        new_conf.confnum = self.confnum
        new_conf.parent_residue = self.parent_residue
        new_conf.occ = self.occ
        new_conf.history = self.history
        new_conf.charge = self.charge
        new_conf.calculated = False # newly created conformer is not calculated
        new_conf.atoms = [atom.clone() for atom in self.atoms]
        return new_conf

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
        self.prepend_lines = []     # lines to prepend to the output
        self.append_lines = []      # lines to append to the output
        self.residues = []          # list of residues in the protein

    def make_ter_residues(self):
        """
        Split terminal residues into NTR and CTR.
        """
        # find chains
        chains = sorted(list(set(res.chain for res in self.residues)))

        aminoacids_in_chains = defaultdict(list)
        nonaminoacids_in_chains = defaultdict(list)
        for res in self.residues:
            if res.resname in AMINO_ACIDS:
                aminoacids_in_chains[res.chain].append(res)
            else:
                nonaminoacids_in_chains[res.chain].append(res)

        # within each chain group, 
        # 1. find the first amino acid and split to NTR
        # 2. find the last amino acit and split to CTR
        for chain in chains:
            res_chain = aminoacids_in_chains[chain]
            # NTR
            if res_chain[0].resname not in TERMINAL_RESIDUES:
                ntr = Residue()
                if res_chain[0].resname == "GLY":
                    ntr.resname = "NTG"
                else:
                    ntr.resname = "NTR"
                ntr.chain = chain
                ntr.sequence = res_chain[0].sequence  # Use the same sequence number as the first amino acid
                ntr.insertion = res_chain[0].insertion
                ntr.resid = (ntr.resname, ntr.chain, ntr.sequence, ntr.insertion)
                ntr.conformers = []
                for conf in res_chain[0].conformers:
                    # clone the conformer and only keep the NTR atoms
                    new_conf = conf.clone()
                    new_conf.confid = f"{ntr.resname}{ntr.chain}{ntr.sequence:04d}{ntr.insertion}000"
                    new_conf.parent_residue = ntr
                    new_conf.atoms = [atom.clone() for atom in conf.atoms if atom.atomname in NTR_ATOMS]
                    for atom in new_conf.atoms:
                        atom.parent_conf = new_conf
                        atom.resname = ntr.resname
                    # remove the NTR atoms from the original conformer
                    conf.atoms = [atom for atom in conf.atoms if atom.atomname not in NTR_ATOMS]
                    if new_conf.atoms:  # only add the conformer if it has atoms
                        ntr.conformers.append(new_conf)
                if ntr.conformers:  # only add the residue if it has conformers
                    res_chain.insert(0, ntr)
                
            # CTR
            if res_chain[-1].resname not in TERMINAL_RESIDUES:
                ctr = Residue()
                ctr.resname = "CTR"
                ctr.chain = chain
                ctr.sequence = res_chain[-1].sequence  # Use the same sequence number as the last amino acid
                ctr.insertion = res_chain[-1].insertion
                ctr.resid = (ctr.resname, ctr.chain, ctr.sequence, ctr.insertion)
                ctr.conformers = []
                for conf in res_chain[-1].conformers:
                    # clone the conformer and only keep the CTR atoms
                    new_conf = conf.clone()
                    new_conf.confid = f"{ctr.resname}{ctr.chain}{ctr.sequence:04d}{ctr.insertion}000"
                    new_conf.parent_residue = ctr
                    new_conf.atoms = [atom.clone() for atom in conf.atoms if atom.atomname in CTR_ATOMS]
                    for atom in new_conf.atoms:
                        atom.parent_conf = new_conf
                        atom.resname = ctr.resname
                    # remove the CTR atoms from the original conformer
                    conf.atoms = [atom for atom in conf.atoms if atom.atomname not in CTR_ATOMS]
                    if new_conf.atoms:
                        ctr.conformers.append(new_conf)
                if ctr.conformers:
                    res_chain.append(ctr)                    

        # assemble the residues in groups back into one array self.residues
        self.residues = [res for chain in chains for res in aminoacids_in_chains[chain] + nonaminoacids_in_chains[chain]]

    def new_ftpl(self, tpl):
        """
        Scan protein residues to find unknown cofactors, if found, create new_ftpl for them and ammend tpl database.
        This should be done before split_backbone() and split_altloc(), so that all atoms in conformer[0]
        """
        new_ftpllines = []
        for res in self.residues:
            key = ("CONFLIST", res.resname)
            if key not in tpl:
                logging.warning(f"   {res.resname} not found in the ftpl database. Creating a new ftpl entry.")
                tpl[key] = [f"{res.resname}BK, {res.resname}01"]
                new_ftpllines.append(f"CONFLIST, {res.resname}: {res.resname}BK, {res.resname}01\n")
                # CONNECT, all atoms are "ion" type of atoms with no connectivity
                for conf in res.conformers:
                    connect_lines = []
                    radius_lines = []
                    for atom in conf.atoms:
                        tpl[("CONNECT", atom.atomname, f"{res.resname}01")] = tpl.CONNECT_param('ion')
                        tpl[("RADIUS", f"{res.resname}01", atom.atomname)] = tpl.RADIUS_param(f'{R_BOUNDARY[atom.element]}, {R_VDW[atom.element][0]}, {R_VDW[atom.element][1]}')
                        connect_lines.append(f"CONNECT, {atom.atomname}, {res.resname}01: ion\n")
                        radius_lines.append(f"RADIUS, {res.resname}01, {atom.atomname}: {R_BOUNDARY[atom.element]}, {R_VDW[atom.element][0]}, {R_VDW[atom.element][1]}\n")
                    new_ftpllines.extend(connect_lines)
                    new_ftpllines.extend(radius_lines)
                new_ftpllines.append("#" + "-" * 89 + "\n")

        if new_ftpllines:
            open(NEW_FTPL, "w").writelines(new_ftpllines)
            logging.info(f"   New ftpl entries are recorded in file {NEW_FTPL}")

    def split_backbone(self, tpl):
        """
        Move backbone atoms to conformer[0], the rest atoms to conformer[1].
        At this point, the backbone atoms should only appear once, the altloc should be detected and resolved by split_altloc.py.
        At this point, we only have one conformer with all atoms in it.
        """
        for res in self.residues:
            conf = res.conformers[0]
            backbone_atoms = [atom for atom in conf.atoms if ("CONNECT", atom.atomname, f"{res.resname}BK") in tpl]
            new_conf = Conformer()
            new_conf.confid = f"{res.resname}{res.chain}{res.sequence:04d}{res.insertion}000"  # confnum=000 is the backbone conformer
            new_conf.conftype = f"{res.resname}BK"
            new_conf.parent_residue = res
            new_conf.history = "BKO000_000"
            new_conf.atoms = backbone_atoms
            for atom in backbone_atoms:
                atom.parent_conf = new_conf
            conf.atoms = [atom for atom in conf.atoms if atom not in backbone_atoms]
            res.conformers = [new_conf, conf]

    def split_altloc(self):
        """
        Split the atoms with alternate locations to different conformers.
        This function is different from split_altloc.py in that this one splits to conformers, while the other one splits to pdb files.
        This functions assumes:
        1. The backbone atoms are already split to conformer[0] by split_backbone().
        2. The atoms with alternate locations are in the same conformer.
        """
        for res in self.residues:
            conf = res.conformers[1]
            # scan all the atoms in the conformer
            # if the atom does not have an alternate location, keep it in the common atom list
            # if the atom has an alternate location, create a new conformer and move the atom to the new conformer
            # if the new conformer already exists, move the atom to the existing conformer
            # combine common atoms and new conformers to the conformer list
            common_atoms = [atom for atom in conf.atoms if atom.altloc == " "]
            new_confs = defaultdict(list)
            for atom in conf.atoms:
                if atom.altloc != " ":
                    new_confs[atom.altloc].append(atom)

            # create new conformers if new_confs is not empty
            # else do nothing (conformer[1] stays the same)
            if new_confs:
                res.conformers = [res.conformers[0]]  # keep the backbone conformer
                for altloc, atoms in new_confs.items():
                    new_conf = Conformer()
                    new_conf.conftype = conf.conftype
                    new_conf.confid = f"{res.resname}{res.chain}{res.sequence:04d}{res.insertion}000"
                    new_conf.parent_residue = res
                    new_conf.history = f"__{altloc}_______"  # The 3rd character is the altloc
                    new_conf.atoms = common_atoms + atoms
                    for atom in new_conf.atoms:
                        atom.parent_conf = new_conf
                    res.conformers.append(new_conf)
                    
    def assign_conftype(self, tpl):
        """
        Assign conformer types to conformers.
        """
        for res in self.residues:
            if len(res.conformers) > 1: # only needed for side chains, backbone conformes[0] type is always "BK"
                key = ("CONFLIST", res.resname)
                possible_types = tpl[key]
                for conf in res.conformers[1:]:
                    conf.conftype = ""  # make sure the conftype is empty
                    for t in possible_types:
                        fit_to_type = True
                        for atom in conf.atoms:
                            if ("CONNECT", atom.atomname, t) not in tpl:
                                fit_to_type = False
                                break
                        if fit_to_type:  # this type fits all atoms
                            conf.conftype = t
#                            print("before", conf.history)
                            conf.history = f"{conf.conftype[-2:]}{conf.history[2]}000_000"
#                          print("after", conf.history)
                            break   # no need to check other types
                    if conf.conftype == "": # after trying all types, if no type fits, report error
                        logging.error(f"   {res.resname} {res.chain}{res.sequence:4d}{res.insertion} has no conformer type to hold all atoms.")


    def dump(self, fname):
        """
        Dump the protein to a pdb file.
        """
        
        # serialize the protein to
        # 1. atom serial number
        # 2. conformer number
        # 3. conformer ID
        # 4. conformer history (same heavy atom type counts in a serials, H atom type counts in regard to heavy and H types)
        atom_counter = 0
        for res in self.residues:
            for conf_counter, conf in enumerate(res.conformers):
                conf.confnum = conf_counter
                conf.confid = f"{res.resname}{res.chain}{res.sequence:04d}{res.insertion}{conf.confnum:03d}"
                for atom in conf.atoms:
                    atom.serial = atom_counter
                    atom_counter += 1
        # history string is tricky, as it carries information about how the conformer is inheried from the parent
        # Sample history string: "BKO000_000" or "01R000M000"
        # history[:2] is the conformer type, 
        # history[2] is rotamer type
        # history[3:6] rotamer number of that type
        # history[6] is how H atom is placed
        # history[7:10] is the conformer number of H atom placement of the type
        for res in self.residues:
            counter_heavy_type = defaultdict(int)
            for conf in res.conformers:
                heavy_confType = conf.history[:6]
                conf.history = f"{conf.history[:3]}{counter_heavy_type[heavy_confType]:03d}{conf.history[6:]}"
                counter_heavy_type[heavy_confType] += 1

            # renumber H atom history
            heavy_id = defaultdict(int)
            for conf in res.conformers:
                id = conf.history[:7]
                conf.history = f"{id}{heavy_id[id]:03d}"
                heavy_id[id] += 1      


        lines = self.prepend_lines
        for res in self.residues:
            lines.append("#" + "=" * 89 + "\n")
            lines.append(f"# Residue: {res.resname} {res.chain}{res.sequence:4d}{' ' if res.insertion == '_' else res.insertion}\n")
            lines.append("#" + "=" * 89 + "\n")
            for conf in res.conformers:
                lines.append(f"## Conformer ID={conf.confid} Type={conf.conftype} History={conf.history}\n")
                lines.extend(atom.as_mccepdb_line() for atom in conf.atoms)
                lines.append("#" + "-" * 89 + "\n")
        lines.extend(self.append_lines)
        if fname is None:
            sys.stdout.writelines(lines)
        else:
            with open(fname, "w") as f:
                f.writelines(lines)


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
        # Load default runprm file
        self._dist_path = os.path.abspath(os.path.join(__file__, "../../.."))
        runprm_default = os.path.join(self._dist_path, RUNPRM_DEFAULT)

        if os.path.exists(runprm_default):
            logging.info(f"   Loading default runprm file {runprm_default}")
            self.update_by_files([runprm_default])
        else:
            logging.warning(f"   Default runprm file {runprm_default} not found.")

        # internal attributes, direct name and value assignment, starts with underscore and all lower case
        self._additional_loose_cofactors = []   # list of additional loose cofactors, loaded from command line


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

        # convert relative paths to absolute paths for special entries
        special_entries = RUNPRM_SPECIAL_ENTRIES
        for entry in special_entries:
            if hasattr(self, entry):
                attr = getattr(self, entry)
                if not os.path.isabs(attr.value):
                    setattr(self, f"_{entry}", copy.deepcopy(attr))
                    getattr(self, f"_{entry}").value = os.path.join(self._dist_path, attr.value)
                else:
                    setattr(self, f"_{entry}", attr)

                
    def update_by_cmd(self, args):
        """
        Update the runprm object using command line arguments.

        Command line options can set runprm key:value pairs in two ways:
        1. Directly by recognizing key and value, such as "-s cutoff".
        2. Using "-r runprm_file" to load additional runprm files.

        Notes:
        While each script handles its own command line options, the following options are recommended for consistency:
        -f ftpl_folder    Load from this ftpl folder
        -r prm [prm ...]  Load additional runprm files, in order
        --debug           Print debug information
        """
        if args.pdb_file:
            self.INPDB.value = args.pdb_file
            self.INPDB.set_by = "#Set by command line"
        
        if args.f:
            self.FTPL_FOLDER.value = args.f
            self.FTPL_FOLDER.set_by = "#Set by command line"

        if args.s:
            self.SAS_CUTOFF.value = str(args.s)
            self.SAS_CUTOFF.set_by = "#Set by command line"
        
        if args.no_ter:
            self.TERMINALS.value = "f"
            self.TERMINALS.set_by = "#Set by command line"
        
        if args.no_hoh:
            self.NO_HOH.value = "t"
            self.NO_HOH.set_by = "#Set by command line"

        # loading additional runprm files should be done at the very end to overwrite previous settings
        if args.r:
            self.update_by_files(args.r)    # set by additional runprm files specified by command line

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


class Tpl:
    """
    Tpl class
    This class stores parameters from ftpl files.
    """
    # Make Tpl a dictionary-like object
    def __init__(self):
        self._db = {}

    def __getitem__(self, key):
        return self._db[key]

    def __setitem__(self, key, value):
        self._db[key] = value

    def __delitem__(self, key):
        del self._db[key]

    def __contains__(self, key):
        return key in self._db

    def keys(self):
        return self._db.keys()

    def values(self):
        return self._db.values()

    def items(self):
        return self._db.items()

    # FTPL key is a tuple of up to 3 strings
    # FTPL value is determined by the type of record, so we have to define one by one.
    # The value types are
    # 1. list of strings (example: CONFLIST)
    # 2. float point number (example: CHARGE)
    # 3. complex object (example: CONNECT)
    class CONNECT_param:
        """
        CONNECT parameter class
        """
        def __init__(self, value_str):
            fields = value_str.split(",")
            self.orbital = fields[0].strip()
            self.connected = [f.strip().strip('"') for f in fields[1:]]

        def __str__(self):
            values = [self.orbital]
            values.extend(['"%s"'%a for a in self.connected])
            value_str = ", ".join(values)
            return value_str

    class RADIUS_param:
        """
        RADII parameter class
        """
        def __init__(self, value_str):
            fields = value_str.split(",")
            self.r_bound = float(fields[0].strip())
            self.r_vdw = float(fields[1].strip())
            self.e_vdw = float(fields[2].strip())

        def __str__(self):
            return f"{self.r_bound}, {self.r_vdw}, {self.e_vdw}"

    class CONFORMER_param:
        """
        CONFORMER parameter class
        """
        def __init__(self, value_str):
            fields = value_str.split(",")
            for f in fields:
                k, v = f.split("=")
                setattr(self, k.strip().lower(), float(v.strip()))

        def __str__(self):
            return ", ".join([f"{k}={v}" for k, v in vars(self).items()])

    class ROTATE_param:
        """
        ROTATE parameter class
        Example: ROTATE, ASP: " CA " - " CB ", " CB " - " CG "
        """
        def __init__(self, value_str):
            fields = value_str.split(",")
            rotatables = [tuple([a.strip().strip('"') for a in f.split("-")]) for f in fields]
            self.rotatables = rotatables

        def __str__(self):
            return ", ".join([f'"{a}" - "{b}"' for a, b in self.rotatables])

    class ROT_SWAP_param:
        """
        ROT_SWAP parameter class
        Example: ROT_SWAP, HIS: " ND1" - " CD2",  " CE1" - " NE2"
        """
        def __init__(self, value_str):
            fields = value_str.split(",")
            swapables = [tuple([a.strip().strip('"') for a in f.split("-")]) for f in fields]
            self.swapables = swapables

        def __str__(self):
            return ", ".join([f'"{a}" - "{b}"' for a, b in self.swapables])
    

    class LIGAND_ID_param:
        """
        LIGAND_ID parameter class
        Examples:
        LIGAND_ID, CYS, CYS: " SG " - " SG "; 2.00 +- 0.20; CYL, CYL
        LIGAND_ID, HIS, HEM: " NE2" - "FE  "; 2.50 +- 0.25; HIL, HEM
        LIGAND_ID, HIS, HEA: " NE2" - "FE  "; 2.50 +- 0.25; HIL, HEA
        LIGAND_ID, HIS, HEB: " NE2" - "FE  "; 2.50 +- 0.25; HIL, HEB
        LIGAND_ID, HIS, HEC: " NE2" - "FE  "; 2.50 +- 0.25; HIL, HEC
        LIGAND_ID, MET, HEA: " SD " - "FE  "; 2.50 +- 0.25; HIL, HEA
        LIGAND_ID, MET, HEB: " SD " - "FE  "; 2.50 +- 0.25; HIL, HEB
        LIGAND_ID, MET, HEC: " SD " - "FE  "; 2.50 +- 0.25; HIL, HEC
        """
        def __init__(self, value_str):
            fields = value_str.split(";")
            # atom pair
            atom1, atom2 = fields[0].strip().split("-")
            self.atom1 = atom1.strip().strip('"')
            self.atom2 = atom2.strip().strip('"')
            # ligand bond distance and tolerance
            distance, tolerance = fields[1].strip().split("+-")
            self.distance = float(distance.strip())
            self.tolerance = float(tolerance.strip())
            # residue rename to
            name1, name2 = fields[2].strip().split(",")
            self.res1_name = name1.strip()
            self.res2_name = name2.strip()  

        def __str__(self):
            value_str = f'"{self.atom1}" - "{self.atom2}"; {self.distance} +- {self.tolerance}; {self.res1_name}, {self.res2_name}'
            return value_str



    # MCCE specific methods            
    def load_ftpl_folder(self, ftpl_folder):
        """
        Load ftpl files from the ftpl folder.
        """
        files =glob.glob(os.path.join(ftpl_folder, "*.ftpl"))
        files.sort()    # sort the files to ensure the order, once can be used to overwrite the other
        logging.info(f"   Loading ftpl files from folder {ftpl_folder}")
        for file in files:
            self.load_ftpl_file(file)
    
    def load_ftpl_file(self, file):
        """
        Load a ftpl file.
        Sample ftpl file:
        ---------------------------------------
        # Values of the same key are appended and separated by ","
        CONFLIST, ASP: ASPBK, ASP01, ASP02, ASP-1

        # Atom definition
        CONNECT, " N  ", ASPBK: sp2, " ?  ", " CA ", " H  "
        CONNECT, " H  ", ASPBK: s, " N  "
        CONNECT, " CA ", ASPBK: sp3, " N  ", " C  ", " CB ", " HA "
        CONNECT, " HA ", ASPBK: s, " CA "
        CONNECT, " C  ", ASPBK: sp2, " CA ", " O  ", " ?  "
        ---------------------------------------
        : separates key and value
        The key is the fields (up to 3) before the first : in a line.
        The value is the fields after the first : in a line.
        """

        with open(file) as f:
            for line in f:
                entry_str = line.split("#")[0].strip()
                fields = entry_str.split(":")
                if len(fields) == 2:
                    key_str = fields[0].strip()
                    # we have up to 3 keys, separated by ","
                    keys = key_str.split(",")
                    key1 = keys[0].strip().strip('"')
                    key2 = keys[1].strip().strip('"') if len(keys) > 1 else ""
                    key3 = keys[2].strip().strip('"') if len(keys) > 2 else ""
                    
                    value_str = fields[1].strip()
                    warn_duplicate_msg = "   Duplicate key {}. Overwriting its value ..."

                    # We have to handle the value case by case here, once for all.
                    if key1 == "CONFLIST":  # value stored as a list of strings
                        #print(value_str)
                        key = (key1, key2)
                        if key in self:
                            logging.warning(warn_duplicate_msg.format(key))
                        self[key] = [v.strip() for v in value_str.split(",")]
                    elif key1 == "CONNECT":  # value stored as a complex object
                        key = (key1, key2, key3)
                        if key in self:
                            logging.warning(warn_duplicate_msg.format(key))
                        self[key] = self.CONNECT_param(value_str)
                    elif key1 == "RADIUS":
                        key = (key1, key2, key3)
                        if key in self:
                            logging.warning(warn_duplicate_msg.format(key))
                        self[key] = self.RADIUS_param(value_str)
                    elif key1 == "CONFORMER":
                        key = (key1, key2)
                        if key in self:
                            logging.warning(warn_duplicate_msg.format(key))
                        self[key] = self.CONFORMER_param(value_str)
                    elif key1 == "CHARGE":  # value stored as a float point number
                        key = (key1, key2, key3)
                        if key in self:
                            logging.warning(warn_duplicate_msg.format(key))
                        self[key] = float(value_str)
                    elif key1 == "ROTATE":
                        key = (key1, key2)
                        if key in self:
                            logging.warning(warn_duplicate_msg.format(key))
                        self[key] = self.ROTATE_param(value_str)
                    elif key1 == "ROT_SWAP":
                        key = (key1, key2)
                        if key in self:
                            logging.warning(warn_duplicate_msg.format(key))
                        self[key] = self.ROT_SWAP_param(value_str)
                    elif key1 == "TORSION":
                        logging.debug(f"   TORSION parameters are not used in this version of MCCE.")
                    elif key1 == "LIGAND_ID":
                        key = (key1, key2, key3)
                        if key in self:
                            logging.warning(warn_duplicate_msg.format(key))
                        self[key] = self.LIGAND_ID_param(value_str)
                    elif key1 == "EXTRA" or key1 == "SCALING":  # captures 2-key float value type parameters
                        key = (key1, key2)
                        if key in self:
                            logging.warning(warn_duplicate_msg.format(key))
                        self[key] = float(value_str)
                    else:
                        logging.warning(f"   {key1} parameters are not defined in {__file__}, value treated as string.")
                        key = (key1, key2, key3)
                        if key in self:
                            logging.warning(warn_duplicate_msg.format(key))
                        self[key] = value_str

            logging.debug(f"   Loaded ftpl file {file}")
    


    def dump(self, comment=""):
        """
        Dump the parameters in the format of a tpl file.
        """
        date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(FTPL_DUMP, "w") as f:
            f.write(f"# This tpl file is recorded on {date}\n{comment}\n")
            for key, value in self.items():
                # wrap double quotes if a key has leading or ending spaces, and are 4 characters long
                key_str = ", ".join(f'"{k}"' if len(k) == 4 and (k[0] == " " or k[-1] == " ") else k for k in key)
                value_str = str(value).strip("[]").replace("'", "")  # remove brackets and single quotes in case it comes from a list
                line = "%s: %s\n" % (key_str, value_str)
                f.write(line)
        logging.info(f"   MCCE ftpl parameters are recorded in file {FTPL_DUMP}")


class Pdb:
    """
    Pdb class
    This class stores the protein data from a pdb file.
    """
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.atoms = []
        self.mcce_ready = False
        self.message = ""
        self.prepend_lines = []  # lines to prepend to the pdb file, linkages, etc.
        self.append_lines = []   # lines to append to the pdb
        self.load_pdb()

    def load_pdb(self):
        """
        Load the pdb file to atoms.
        """ 
        if not os.path.exists(self.pdb_file):
            self.mcce_ready = False
            self.message = f"PDB file {self.pdb_file} not found."
            return

        with open(self.pdb_file) as f:
            # detect if this is a multi-model pdb file
            n_model = 0
            for line in f:
                if line.startswith("MODEL"):
                    n_model += 1
            if n_model > 1:
                self.mcce_ready = False
                self.message = "Multi-model pdb file is not supported. Use split_nmr.py to divide the file into individual models."
                return

            # read the pdb file
            f.seek(0)
            for line in f:
                if line.startswith("ATOM  ") or line.startswith("HETATM"):
                    atom = Atom()
                    atom.load_pdbline(line)
                    if not is_H(atom.atomname):
                        self.atoms.append(atom)

            # detect is the altloc happens on backbone atoms
            altloc_on_backbone = set()
            for atom in self.atoms:
                if atom.atomname in BACKBONE_ATOMS and atom.resname in RESIDUE_NAMES and atom.altloc != " ":
                    altloc_on_backbone.add(atom.altloc)
                    
            if len(altloc_on_backbone) > 1:
                self.mcce_ready = False
                self.message = f"Backbone atoms with altLoc are not supported. Use split_altloc.py to split the pdb file."
                return

            self.mcce_ready = True

    def rename(self, rules):
        """
        Rename pdb lines according to the rules in runprm.
        """

        def match_rule2string(rule, string):
            return all(r == "*" or r == s for r, s in zip(rule, string))

        def rename_rule2string(rule, string):
            return "".join([r if r != "*" else s for r, s in zip(rule, string)])

        if os.path.exists(rules):
            with open(rules) as f:
                rename_rules = [(line[:14], line[16:30]) for line in f if len(line.split("#")[0]) >= 30]

            for atom in self.atoms:
                for rule, newname in rename_rules:
                    if (
                        match_rule2string(rule[:4], atom.atomname) and
                        (rule[4] == "*" or rule[4] == atom.altloc) and
                        match_rule2string(rule[5:8], atom.resname) and
                        (rule[9] == "*" or rule[9] == atom.chain) and
                        match_rule2string(rule[10:14], f"{atom.sequence:4d}")
                    ):
                        atom.atomname = rename_rule2string(newname[:4], atom.atomname)
                        atom.altloc = newname[4] if rule[4] != "*" else atom.altloc
                        atom.resname = rename_rule2string(newname[5:8], atom.resname)
                        atom.chain = newname[9] if rule[9] != "*" else atom.chain
                        atom.sequence = int(rename_rule2string(newname[10:14], f"{atom.sequence:4d}"))
            logging.info(f"   Atoms are renamed according to the rules in {rules}")
        else:
            logging.warning(f"   Rename rules file {rules} not found. Nothing is renamed.")
            
    def remove_hoh(self):
        """
        Remove water atoms from the pdb file.
        """
        self.atoms = [atom for atom in self.atoms if atom.resname != "HOH"]
        

    def dump_pdb(self, fname):
        """
        Dump atoms to a pdb file.
        """
        with open(fname, "w") as f:
            counter = 0
            for atom in self.atoms:
                counter += 1
                line = (
                    f"{atom.record:6}{counter:5d} {atom.atomname:4}{atom.altloc:1}{atom.resname:3} "
                    f"{atom.chain:1}{atom.sequence:4}{atom.insertion:1}   {atom.xyz.x:8.3f}{atom.xyz.y:8.3f}{atom.xyz.z:8.3f}  "
                    f"1.00  0.00          {atom.element:2}\n"
                )
                f.write(line)

    def identify_ligands(self, tpl):
        """
        Identify ligands in the pdb file. The ligands detection rules are in ligand_detect_rules.ftpl
        Sample rules:
        ---------------------------------------
        LIGAND_ID, CYS, CYS: " SG " - " SG "; 2.03 +- 0.90; CYD, CYD
        LIGAND_ID, CYS, HEC: " SG " - " CA*"; 1.90 +- 1.00; CYL, HEC
        LIGAND_ID, CYS, HEM: " SG " - " CA*"; 1.90 +- 1.00; CYL, HEM
        """
        # This is an internal function to compare two strings with wildcard "*"
        def match_strings(s1, s2):
            return all([r == "*" or s=="*" or r == s for r, s in zip(s1, s2)]+[len(s1) == len(s2)])


        logging.info("   Identifying ligands in the pdb file.")


        link_lines = []
        ssbond_serial = 1  # counter for SSBOND serial number starting from 1
        ssbond_fmt = "SSBOND %3d CYS %c %4d%c   CYS %c %4d%c %s  1555   1555 %5.2f\n"
        link_fmt = "LINK        %4s %3s %c%4d%c %s  %4s %3s %c%4d%c    1555   1555 %5.2f\n"

        # group atoms by residue, the residues should be in the order of appearance in the pdb file after Python 3.7
        residue_atoms_dict = defaultdict(list)
        for atom in self.atoms:
            residue_atoms_dict[(atom.resname, atom.chain, atom.sequence, atom.insertion)].append(atom)
        residue_ids = list(residue_atoms_dict.keys())
        
        # loop over residues from 1 to the second last residue
        for i, res1_id in enumerate(residue_ids[:-1]):
            for res2_id in residue_ids[i+1:]:
                # check if the two residues can form a ligand
                key1 = ("LIGAND_ID", res1_id[0], res2_id[0])
                key2 = ("LIGAND_ID", res2_id[0], res1_id[0])
                if key1 in tpl or key2 in tpl:      # key is found, potential ligand pair
                    res1, res2 = (res1_id, res2_id) if key1 in tpl else (res2_id, res1_id)
                    ligand_param = tpl[key1 if key1 in tpl else key2]
                    # check if the two residues are close enough to form a ligand
                    atom1_name = ligand_param.atom1
                    atom2_name = ligand_param.atom2
                    distance = ligand_param.distance
                    tolerance = ligand_param.tolerance
                    # find atom1 and atom2 in the residues
                    atom1 = [atom for atom in residue_atoms_dict[res1] if match_strings(atom.atomname, atom1_name)]
                    atom2 = [atom for atom in residue_atoms_dict[res2] if match_strings(atom.atomname, atom2_name)]
                    if atom1 and atom2:
                        # calculate the distance between the two atoms
                        for a1 in atom1:
                            for a2 in atom2:
                                d = a1.xyz.distance(a2.xyz)
                                if distance - tolerance <= d <= distance + tolerance:
                                    # create a link line
                                    logging.debug(f"   Ligand detected between {a1.atomname} {res1} and {a2.atomname} {res2} with distance {d:.2f}")
                                    for atom in residue_atoms_dict[res1]:
                                        atom.resname = ligand_param.res1_name
                                    for atom in residue_atoms_dict[res2]:
                                        atom.resname = ligand_param.res2_name

                                    if res1[0] == "CYS" and res2[0] == "CYS":     # SSBOND
                                        line = ssbond_fmt % (ssbond_serial, res1[1], res1[2], res1[3], res2[1], res2[2], res2[3], " "*22, d)
                                        link_lines.append(line)
                                        ssbond_serial += 1
                                    else:  # LINK
                                        line = link_fmt % (a1.atomname, res1[0], res1[1], res1[2], res1[3], " "*12,
                                            a2.atomname, res2[0], res2[1], res2[2], res2[3], d)
                                        link_lines.append(line)
        self.atoms = [atom for res in residue_atoms_dict.values() for atom in res]
        self.prepend_lines.extend(link_lines)
        logging.info(f"   {len(link_lines)} ligands are identified in the pdb file.")


    def convert_to_protein(self, tpl):
        """
        Convert pdb atoms to Protein object.
        Protein is a hirarchy of Residue -> Conformer -> Atom.
        """

        protein = Protein()
        for atom in self.atoms:
            # create a new conformer if the residue does not exist
            residue = None
            for res in protein.residues:
                if res.resid == (atom.resname, atom.chain, atom.sequence, atom.insertion):
                    residue = res
                    break
            if residue is None:
                residue = Residue()
                residue.resname = atom.resname
                residue.chain = atom.chain
                residue.sequence = atom.sequence
                residue.insertion = atom.insertion
                residue.resid = (atom.resname, atom.chain, atom.sequence, atom.insertion)
                protein.residues.append(residue)

            # create a new conformer if the conformer does not exist
            conformer = None
            for conf in residue.conformers:
                if conf.confid == f"{atom.resname}{atom.chain}{atom.sequence:04d}{atom.insertion}000": # conformer number assumed to be 0
                    conformer = conf
                    break
            if conformer is None:
                conformer = Conformer()
                conformer.confnum = 0
                conformer.confid = f"{atom.resname}{atom.chain}{atom.sequence:04d}{atom.insertion}000" # use 000 as conformer number
                conformer.conftype = "NA"
                conformer.history = "__O000_000"
                conformer.altloc = atom.altloc
                conformer.resname = atom.resname
                conformer.chain = atom.chain
                conformer.sequence = atom.sequence
                conformer.insertion = atom.insertion
                conformer.parent_residue = residue
                residue.conformers.append(conformer)

            # add this atom to conformer
            atom.parent_conf = conformer
            conformer.atoms.append(atom)

        # update prepend and append lines
        protein.prepend_lines = self.prepend_lines
        protein.append_lines = self.append_lines

        return protein

def is_H(atom_name):
    # It is an H atom when full 4 chars starts with H, or less than 4-char but the 2nd char is H
    return (len(atom_name.strip()) == 4 and atom_name[0] == "H") or (len(atom_name.strip()) < 4 and atom_name[1] == "H" )
