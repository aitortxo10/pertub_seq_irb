import subprocess, sys, os, copy, gzip
from Bio import SeqIO
from io import StringIO
from fuzzywuzzy import process

def main(read2,gz=True):
    """Helper function that contains the overall process. This function will then be passed to multiprocessing commands"""

    bc_dict=parse_barcodes(bc_file) #Parse the barcode fasta into a barcode dict
    bc_list=list(bc_dict.values()) #Generate the list of barcodes sequence and ids
    bc_ids=list(bc_dict.keys())

    mutant_dict={i:[] for i in bc_ids} #Empty dictionary with all the bc ids as keys
#    lineage_dict=copy.deepcopy(mutant_dict) #Deep copy it to store the lineage counts

    read2_bc_dict=run_cutadapt(read2,"CGAGCTCGAATTCATCGAT","CTACGAGACCGACACCG")#Extract the mutant bc from the reads
    sys.stderr.write(str(len(read2_bc_dict.keys())))
    for r2 in read2_bc_dict:
        if len(read2_bc_dict[r2])==0:
            continue

        highest_match=process.extractOne(read2_bc_dict[r2],bc_dict.values())

        if not highest_match:
            continue
        elif not highest_match[1]>=95:
            continue

        id_location=bc_list.index(highest_match[0])
        id=bc_ids[id_location]

        mutant_dict[id].append(r2)

    mutant_dict=dict([(k,v) for k,v in mutant_dict.items() if len(v)>0])
    return mutant_dict

    lineage_dict=run_cutadapt(read2,"AACGCCGCCATCCAGTGTCGAAAACGAGCTCGAATTCATCGAT")

    mutant_lineages=dict()

    for mutant in mutant_dict:
        lineages=[lineage_dict[x] for x in mutant_dict[mutant]]
        lineages_unique=set(lineages)
        for lin in lineages_unique:
            new_id=mutant+lin
            mutant_lineages[new_id]=lineages.count(lin)

    return results

def run_cutadapt(target_file,target1,target2=""):

    if target2:
        cmd="cutadapt -j 10 --rc --discard-untrimmed -a {s1}...{s2} {file}".format(s1=target1,s2=target2,file=target_file)
    else:
        cmd="cutadapt -j 10 --rc --discard-untrimmed -a {s1} {file}".format(s1=target1,file=target_file)

    cutadapt_res=subprocess.run(cmd, shell=True, capture_output=True)
    reparse=cutadapt_res.stdout.decode("utf-8")
    return_dict=dict()

    with StringIO(reparse) as fastqIO:
        parser=SeqIO.parse(fastqIO,"fastq")

        for bc in parser:
            return_dict[bc.id]=str(bc.seq)

    return return_dict

def parse_barcodes(bc):
    """Given a fasta file with all the bcs, return a dictionary with the bc name as key and the seq as values. In this case, the bc is considered to be placed 17 bases up the end
    of the seq and 20 bp long"""

    parser=SeqIO.parse(bc,"fasta")
    bc_dict=dict()

    for sequence in parser:
        bc_dict[sequence.id]=str(sequence.seq)[-37:-17]

    return bc_dict
