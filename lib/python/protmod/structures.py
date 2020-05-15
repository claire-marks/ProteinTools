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
    A dictionary containing secondary structure assignments, using DSSP codes 
        (i.e. E = beta strand, H = alpha helix, etc.). Dictionary keys are the
        protein chain IDs, and values are a list of tuples giving residue 
        ss assignments: (resid, ss_assignment)
    """
    
    # Load structure using Biopython and select specified model
    structure = PDB.PDBParser(QUIET=True).get_structure("structure", pdbfile)
    model = structure[modelno]
    
    # Run DSSP
    dssp_result = PDB.DSSP(model, pdbfile)
    
    # Extract SS assignments
    # List of tuples format maintains the correct residue order
    result = {}
    for chain, res in dssp_result.keys():
        if chain not in result:
            result[chain] = []
        resid = biopython_resid_to_str(res)
        k = (chain, res)
        result[chain].append((resid, dssp_result[k][1], dssp_result[k][2]))

    return result
    
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
            raise Exception("Error: start residue not found. May be a modified residue")

    # Look in residue list for the end residue
    if end in structure[model][chain]:
        end_index = residues.index(structure[model][chain][end])
    else:
        # Couldn't find it - maybe it has a HETATM code
        try:
            res = [i for i in residues if i.id[1] == end[1] and i.id[2] == end[2]][0]
            end_index = residues.index(res)
        except:
            raise Exception("Error: end residue not found. May be a modified residue")

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
        raise Exception("Argument 'returntype' not recognised, must be either 'atom' or 'residue'")

################################################################################
def check_consecutive(fragment):
    """
    Takes a list of BioPython residues and checks the residues are consecutive.
    
    In the first instance, based on the distance between the C of one residue 
    and the N of the next; if less than 2A then the residues are assumed to be 
    consecutive. If either of those two atoms are not present, we try the CA-CA
    distance (non-consecutive if over 4.5A).

    --PARAMETERS--
    fragment: a list of Biopython residues.
    
    --RETURNS--
    True if the residues in the input are consecutive, False if not.
    """

    consecutive = True
    for resno in range(len(fragment)-1):
        # Check the C-N distance between the two residues in the first instance
        if fragment[resno].has_id("C") and fragment[resno+1].has_id("N"):
            CN_distance = fragment[resno]["C"] - fragment[resno+1]["N"]
            if CN_distance > 2:
                consecutive = False
                return consecutive
        # If one or other of the atoms is not present, check CA-CA distance
        elif fragment[resno].has_id("CA") and fragment[resno+1].has_id("CA"):
            CA_distance = fragment[resno]["CA"] - fragment[resno+1]["CA"]
            if CA_distance > 4.5:
                consecutive = False
                return consecutive
        else:
            print("Warning: missing atoms; cannot check consecutiveness of some residues")
            continue

    return consecutive
    
################################################################################
def calculate_rmsd(fragment1, fragment2, atom_types=["N", "CA", "C", "O"]):
    """
    Calculates the RMSD (root mean square difference) between two protein 
    fragments. It does NOT superimpose the two, so do that first if required!
    
    --PARAMETERS--
    fragment1, fragment2: two lists of either Biopython Residues, Biopython
        Atoms, or simple coordinates.
    atom_types: if provided with lists of residues, which atom types should be
        included in the RMSD calculation. This has no effect if the fragments
        given are lists of atoms or coordinates (in that case everything is 
        included).
    
    --RETURNS--
    The RMSD between the two structures.
    """
    
    # Check that the two lists are the same length
    assert len(fragment1) == len(fragment2), "Number of atoms/residues doesn't match!"
    
    # If the lists contain Biopython residues:
    if str(type(fragment1[0])) == "<class 'Bio.PDB.Residue.Residue'>":

        try:
            sum_sqdiff = 0
            atom_count = 0
            for i in range(len(fragment1)):
                for a in atom_types:
                    if fragment1[i].has_id(a) and fragment2[i].has_id(a):
                        sum_sqdiff += ((fragment1[i][a].coord[0] - fragment2[i][a].coord[0])**2 +
                                       (fragment1[i][a].coord[1] - fragment2[i][a].coord[1])**2 +
                                       (fragment1[i][a].coord[2] - fragment2[i][a].coord[2])**2)
                        atom_count += 1

            RMSD = (sum_sqdiff/atom_count)**0.5

            return RMSD
        except:
            "Error: Something went wrong! Expected Biopython Residues?"
            return None

    # If the lists contain Biopython atoms:
    if str(type(fragment1[0])) == "<class 'Bio.PDB.Atom.Atom'>" or str(type(fragment1[0]) == "<class 'Bio.PDB.Atom.DisorderedAtom'>"):

        try:
            sum_sqdiff = 0
            for i in range(len(fragment1)):
                sum_sqdiff += ((fragment1[i].coord[0] - fragment2[i].coord[0])**2 +
                               (fragment1[i].coord[1] - fragment2[i].coord[1])**2 +
                               (fragment1[i].coord[2] - fragment2[i].coord[2])**2)

            RMSD = (sum_sqdiff/len(fragment1))**0.5

            return RMSD
        except:
            "Error: Something went wrong! Expected Biopython atoms?"
            return None

    # If not Biopython atoms, try just list of coordinates
    if str(type(fragment1[0])) == "<class 'Bio.PDB.Atom.Atom'>":

        try:
            sum_sqdiff = 0
            for i in range(len(fragment1)):
                sum_sqdiff += ((fragment1[i][0] - fragment2[i][0])**2 +
                               (fragment1[i][1] - fragment2[i][1])**2 +
                               (fragment1[i][2] - fragment2[i][2])**2)

            RMSD = (sum_sqdiff/len(fragment1))**0.5

            return RMSD
        except:
            "Error: Something went wrong! Expected a list of coordinates?"
            return None
