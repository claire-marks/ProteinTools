# Contains basic functions for dealing with proteins

################################################################################
def get_AAs(codelength=1):
    """
    Get a list of the 20 standard amino acid types. 
    
    --PARAMETERS--
    codelength: specifies whether to return 1-letter or 3-letter codes. Default
    is 1. If a value other than 1 or 3 is given, defaults to 1.
        
    --RETURNS--
    A list of amino acid types (in alphabetical order).
    """

    oneletter   = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

    threeletter = ["ALA", "ARG", "ASN", "ASP", "CYS", 
                   "GLN", "GLU", "GLY", "HIS", "ILE",
                   "LEU", "LYS", "MET", "PHE", "PRO",
                   "SER", "THR", "TRP", "TYR", "VAL"]

    if codelength == 1:
        return oneletter
    elif codelength == 3:
        return threeletter
    else:
        print("Can 't return %d-letter codes! Returning 1-letter codes instead" %(codelength))
        return oneletter
      
################################################################################  
def code_convert(aa, non_aa_return=None):
    """
    Convert a 3-letter amino acid code to a 1-letter one, or vice versa.
    
    --PARAMETERS--
    aa: the amino acid code you wish to convert
    non_aa_return: what the function should return if the provided amino acid 
        code is not one of the standard 20. Default is None.
        
    --RETURNS--
    The equivalent 1- or 3-letter code for the input.
    """

    conversions =  {"ALA": "A", "GLU": "E", "GLN": "Q", "ASP": "D", "ASN": "N", 
                    "LEU": "L", "GLY": "G", "LYS": "K", "SER": "S", "VAL": "V", 
                    "ARG": "R", "THR": "T", "PRO": "P", "ILE": "I", "MET": "M", 
                    "PHE": "F", "TYR": "Y", "CYS": "C", "TRP": "W", "HIS": "H",
                    
                    "A": "ALA", "E": "GLU", "Q": "GLN", "D": "ASP", "N": "ASN", 
                    "L": "LEU", "G": "GLY", "K": "LYS", "S": "SER", "V": "VAL", 
                    "R": "ARG", "T": "THR", "P": "PRO", "I": "ILE", "M": "MET", 
                    "F": "PHE", "Y": "TYR", "C": "CYS", "W": "TRP", "H": "HIS"}

    if aa in conversions:
        return conversions[aa]
    else:
        return non_aa_return