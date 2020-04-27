# Functions to deal with protein sequences

################################################################################
def write_fasta(seqid, sequence, outpath):

    write_str = ">%s\n" %(seqid)

    i = 0
    while sequence[i:i+80] != "":
        write_str += sequence[i:i+80] + "\n"
        i += 80

    with open(outpath, "w") as f:
        f.write(write_str)

    return