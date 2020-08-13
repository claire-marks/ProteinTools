# Functions that deal with downloading data from the web

import requests
from .sequences import parse_fasta

################################################################################
def download(url, filename="download.txt"):
    """
    Download a file from the internet.
    
    --ARGUMENTS--
    url: the address to attempt to download from.
    filename: where the file should be saved.
    
    --RETURNS--
    The path to the downloaded file. If the download failed, this is None.
    """
        
    r = requests.get(url, allow_redirects=True)
    
    # Check file was downloaded successfully; if not return None
    if r.ok:
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename
    else:
        print("Error downloading file %s. Returning None" %(url))
        return None

################################################################################
def download_pdb_structure(pdbcode, pdbfile=None):
    """
    Download a PDB structure from the PDB website.
    
    --ARGUMENTS--
    pdbcode: the four-letter PDB ID for the entry that you wish to download.
    pdbfile: where the file should be saved. Default is to save the file in
        the current directory, filename <pdbcode>.pdb.
    
    --RETURNS--
    The path to the downloaded PDB file. If the download failed, this is None.
    """

    # Ensure the provided PDB code is only 4 letters long
    if len(pdbcode) != 4:
        raise Exception("%s is not a valid PDB code (should be four letters long)." %(pdbcode))
    
    # Establish url and output filename
    url = "http://files.rcsb.org/download/%s.pdb" %(pdbcode)
    
    if not pdbfile:
        pdbfile = "%s.pdb" %(pdbcode)
        
    outfile = download(url, pdbfile)
        
    return outfile
    
################################################################################
def fetch_sequence_from_pdb(pdbcode, chain=None):
    """
    Fetch a sequence of a PDB entry.
    
    --ARGUMENTS--
    pdbcode: the four-letter PDB ID for the entry that you wish to download.
    chain: specify the chain you are interested in (optional).
        
    --RETURNS--
    The sequence for the specified PDB entry. If the chain is not specified, the 
        function will return all sequences stored in a dictionary, with the
        chain ID as the key.
    """
    
    # Ensure the provided PDB code is only 4 letters long
    if len(pdbcode) != 4:
        raise Exception("%s is not a valid PDB code (should be four letters long)." %(pdbcode))
    
    # Establish url and output filename
    url = "https://www.rcsb.org/fasta/entry/%s/display" %(pdbcode.upper())
    
    # Load fasta file from web
    r = requests.get(url, allow_redirects=True)
    fasta = r.content.decode().split("\n")
    sequences = parse_fasta(fasta)

    # Change dictionary keys to just chain IDs
    by_chain = {}
    for seq in sequences:
        # Relies on PDB fasta headings remaining in the same format!
        chainid = seq.split("|")[1].split()[-1]
        if "," in chainid:
            chainids = chainid.split(",")
            for chainid in chainids:
                by_chain[chainid] = sequences[seq]
        else:
            by_chain[chainid] = sequences[seq]
        
    if not chain:
        return by_chain
    else:
        if chain in by_chain:
            return by_chain[chain]
        else:
            raise Exception("Chain %s not found for PDB entry %s." %(chain, pdbcode))
            return
    