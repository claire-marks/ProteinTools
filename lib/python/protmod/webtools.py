# Functions that deal with downloading data from the web

import requests

################################################################################
def download(url, filename="download.txt"):
    """
    Download a file from the internet.
    
    --PARAMETERS--
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
    
    --PARAMETERS--
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