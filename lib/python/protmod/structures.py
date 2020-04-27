# Protein structure functions

import re
import math
import numpy as np
from Bio import PDB

################################################################################
def calculate_dihedral(A, B, C, D, units="degrees"):
    """
    Calculate the dihedral/torsion angle between four points.
    
    --PARAMETERS--
    A,B,C,D: the four points used to calculate the angle. These can be in the
        form of coordinates (a numpy array) or Biopython atoms.
    
    --RETURNS--
    The dihedral angle A-B-C-D.
    """

    # Get coordinates if user provided Biopython atoms
    if str(type(A)) == "<class 'Bio.PDB.Atom.Atom'>":
        A = A.coord

    if str(type(B)) == "<class 'Bio.PDB.Atom.Atom'>":
        B = B.coord

    if str(type(C)) == "<class 'Bio.PDB.Atom.Atom'>":
        C = C.coord

    if str(type(D)) == "<class 'Bio.PDB.Atom.Atom'>":
        D = D.coord

    b1 = B-A
    b2 = C-B
    b3 = D-C

    n1 = np.cross(b1,b2)
    n2 = np.cross(b2,b3)

    n1_unit = n1/np.linalg.norm(n1)
    n2_unit = n2/np.linalg.norm(n2)

    b2_unit = b2/np.linalg.norm(b2)
    m1 = np.cross(n1,b2_unit)

    x = np.dot(n1,n2)
    y = np.dot(m1,n2)

    if units == "degrees":
        angle = -math.atan2(y,x) * 180 / math.pi
    elif units == "radians":
        angle = -math.atan2(y,x)

    return angle

################################################################################
def biopython_resid_to_str(resid, remove_whitespace=True):
    """
    Convert a Biopython-style residue ID to a simple string (as would be found
    in a PDB file).
    
    --PARAMETERS--
    resid: a residue ID in Biopython format - for example (' ', 100, 'A')
    remove_whitespace: if True, blank insertion codes are removed (i.e. for the
        residue (' ', 100, ' '), instead of returning ' 100 ', returns '100'.
    
    --RETURNS--
    A string version of the residue ID - e.g. (' ', 100, 'A') becomes '100A'.
    """
    
    resid_str = resid[0] + str(resid[1]) + resid[2]
    
    if remove_whitespace:
        resid_str = resid_str.strip()
    
    return resid_str
    
################################################################################
def to_biopython_resid(resid, het=" "):
    """
    Convert a residue ID (int or str) to an ID in Biopython format.
    
    --PARAMETERS--
    resid: a string or integer to be converted to a Biopython residue ID.
    het: specify a 'HETATM' type; e.g. 'W' for water molecules.
    
    --RETURNS--
    A residue ID in Biopython format. For example, the string '100A' becomes 
        (' ', 100, 'A').
    
    """
    
    inscode = " "
    if type(resid) != int:
        try:
            # If the input can be converted directly to an integer, then there
            # is no insertion code
            resno = int(resid)
        except:
            # Split string into numbers + letters
            split_string = re.split("(\d+)", resid)
            split_string = [i for i in split_string if i != ""]
            resno = int(split_string[0])
            
            # Check whether there is an insertion code present
            if len(split_string) > 1:
                inscode = split_string[1]
    else:
        resno = resid
                
    biopy_resid = (het, resno, inscode)

    return biopy_resid
    
################################################################################
def assign_secondary_structure(pdbfile, modelno=0):
    """
    Uses DSSP via Biopython to assign secondary structure.
    Requires DSSP to be installed.
    
    --PARAMETERS--
    pdbfile: the path to the structure (in PDB format) for which you wish to
        get secondary structure assignments.
    modelno: specify which model in the PDB file should be analysed. Default is 
        the first.
        
    --RETURNS--
    A dictionary of secondary structure assignments, using DSSP codes (i.e. E =
        beta strand, H = alpha helix, etc.). 
        Format is dict[chain_id][res_id] = DSSP assignment
    """
    
    # Load structure using Biopython and select specified model
    structure = PDB.PDBParser(QUIET=True).get_structure("structure", pdbfile)
    model = structure[modelno]
    
    # Run DSSP
    dssp_result = PDB.DSSP(model, pdbfile)
    
    # Extract SS assignments
    result_dict = {}
    for chain, res in dssp_result.keys():
        resid = biopython_resid_to_str(res)
        if chain in result_dict:
            result_dict[chain][resid] = dssp_result[(chain, res)][2]
        else:
            result_dict[chain] = {resid: dssp_result[(chain, res)][2]}

    return result_dict
    
################################################################################
def get_fragment(structure, model, chain, start, end, anchorlength=0, 
                 returntype="atom", atoms=["N", "CA", "C", "O"]):
    """
    Extracts a section of a protein.
    
    --PARAMETERS--
    structure: either a path to a PDB file, or a BioPython PDB Structure object.
    model: the number of the model from which to extract the fragment.
    chain: the chain of the protein from which to extract the fragment.
    start,end: the IDs of the first and last residues of your desired fragment 
        (can be ints, or strings if an insertion code is used).
    anchorlength: the number of residues to be added to each end of the 
        fragment.
    returntype: either 'atom' or 'residue'; specifies whether to return a list
        of Biopython 'Atoms' or Biopython 'Residues'.
    atoms: a list of required atom types (default is backbone atoms).

    --RETURNS--
    A list of Biopython Atoms/Residues for the requested fragment, and if 
    returntype="atoms", a list of atoms that are missing from the structure.    
    """

    # Convert start/end inputs to Biopython format
    start = to_biopython_resid(start)
    end   = to_biopython_resid(end)

    # Load structure, get a list of residues
    if type(structure) == str:
        structure = PDB.PDBParser(QUIET=True).get_structure("protein", structure)
    residues = structure[model][chain].child_list

    # Look in residue list for the start residue
    if start in structure[model][chain]:
        start_index = residues.index(structure[model][chain][start])
    else:
        # Couldn't find it - maybe it has a HETATM code
        try:
            res = [i for i in residues if i.id[1] == start[1] and i.id[2] == start[2]][0]
            start_index = residues.index(res)
        except:
            print "Error: start residue not found. May be a modified residue"
            return

    # Look in residue list for the end residue
    if end in structure[model][chain]:
        end_index = residues.index(structure[model][chain][end])
    else:
        # Couldn't find it - maybe it has a HETATM code
        try:
            res = [i for i in residues if i.id[1] == end[1] and i.id[2] == end[2]][0]
            end_index = residues.index(res)
        except:
            print "Error: end residue not found. May be a modified residue"
            return

    # Extract fragment
    fragment_res = residues[start_index-anchorlength:end_index+anchorlength+1]

    # Convert to requested format (residue or atom list)
    if returntype.lower() == "residue":
        return fragment_res
    elif returntype.lower() == "atom":
        fragment_atom = []
        missing = []
        for res in fragment_res:
            for atom in atoms:
                try:
                    fragment_atoms.append(res[atom])
                except:
                    fragment_atoms.append(None)
                    missing.append(res.resname+"_"+atom)
        return fragment_atoms, missing
    else:
        print "Argument 'returntype' not recognised, must be either 'atom' or 'residue'"
        return
