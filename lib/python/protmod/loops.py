# Functions for tasks involving protein loops

from numpy import split, diff, where
from .structures import assign_secondary_structure, get_fragment, check_consecutive

################################################################################
def find_loops(pdbfile, ss_assignments=["H", "G", "I", "E"], max_ss_len=3, 
    min_length=3, max_length=30, remove_termini=True, check_missing=True):
    """
    Find which regions in a protein structure correspond to loops.
    
    --PARAMETERS--
    pdbfile: the structure in which to search for loops
    ss_assignments: which DSSP types should be considered as non-loops
    max_ss_len: the maximum number of consecutive non-loop residues that can 
        appear within a loop (i.e. will not split the loop into two).
    min_length: the minimum length of loops to be returned.
    max_length: the maximum length of loops to be returned.
    remove_termini: if True, the function will not return the termini of the 
        protein.
    check_missing: if True, each loop will be checked for consecutive-ness to
        establish whether the loop is complete. If False, the function will be 
        faster, but the check will not be done. If it is indicated that a loop
        is incomplete, the sequence and length given will not reflect the true
        loop structure, only what was present in the PDB file.
    
    --RETURNS--
    A list of dictionaries containing the details of each loop.
    """

    # Run DSSP to get secondary structure assignments
    ss = assign_secondary_structure(pdbfile)
    
    result = {}
    for chain in ss:
        result[chain] = []
        # Find which residues are not part of secondary structure elements such as 
        # helices and sheets
        loop_res = [r for r in zip(range(len(ss[chain])), ss[chain]) if r[1][2] not in ss_assignments]

        # Split into loop regions, split by secondary structure elements of at least 3 residues
        # So loops can have up to 2 consecutive residues of secondary structure within them
        loop_regions = split(loop_res, where(diff([r[0] for r in loop_res])>max_ss_len)[0]+1)

        # Get rid of loops that are too short or too long        
        loop_regions = [[list(r) for r in loop] for loop in loop_regions if len(loop) >= min_length]
        loop_regions = [loop for loop in loop_regions if len(loop) <= max_length]

        # Get rid of protein termini if requested
        if remove_termini:
            loop_regions = [loop for loop in loop_regions if loop[0][0] != 0 and loop[-1][0] != len(ss[chain])]
            
        # Iterate through selected loops and check whether the residues are consecutive
        # Hence checking for missing residues (loop will be longer than reported)
        for loop in loop_regions:
            start = loop[0][1][0]
            end = loop[-1][1][0]
            sequence = "".join([r[1][1] for r in loop])
            
            result[chain].append({"start"   : start, 
                                  "end"     : end,
                                  "length"  : len(loop),
                                  "sequence": sequence})
            
            if check_missing:
                loop_structure = get_fragment(pdbfile, 0, chain, start, end, returntype="residue") 
                if check_consecutive(loop_structure):
                    result[chain][-1]["complete"] = True
                else:
                    result[chain][-1]["complete"] = False

    return result