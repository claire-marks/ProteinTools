# Functions to deal with protein sequences

import os

################################################################################
def write_fasta(seqid, sequence, outpath):
    """
    Write a sequence to a file with fasta format.
    
    --ARGUMENTS--
    seqid: the identifier for the sequence (will be stored in the line starting
        with '>' in the fasta file, above the sequence)
    sequence: the amino acid sequence to be written
    outpath: the location to which the fasta file should be written.
    """

    write_str = ">%s\n" %(seqid)

    i = 0
    while sequence[i:i+80] != "":
        write_str += sequence[i:i+80] + "\n"
        i += 80

    with open(outpath, "w") as f:
        f.write(write_str)

    return
    
################################################################################
def parse_fasta(fasta):
    """
    Parse a set of sequences stored in a fasta file.
    
    --ARGUMENTS--
    fastafile: the contents of a fasta file, or the path to a fasta file.
    
    --RETURNS--
    A dictionary of sequences, where the key is the ID for the sequence (given 
    on the line that starts with '>' before each sequence).
    """
    
    # If provided with a path to a fasta file, load data
    if isinstance(fasta, str):
        if os.path.isfile(fasta):
            fasta = open(fasta).readlines()
        else:
            raise Exception("File not found.")
      
    # Read sequences and add to dictionary line by line  
    sequences = {}
    for l in fasta:
        if l.startswith(">"):
            seqid = l[1:].strip()
            sequences[seqid] = ""
        else:
            sequences[seqid] += l.strip()
            
    return sequences
