import subprocess, sys, os, copy
from Bio import SeqIO
from io import StringIO
from fuzzywuzzy import process
from multiprocessing import Pool
from functools import partial

def main(read2,bc_file):
    """Helper function that contains the overall process. This function will then be passed to multiprocessing commands"""

    bc_dict=parse_barcodes(bc_file) #Parse the barcode fasta into a barcode dict
    bc_list=list(bc_dict.values()) #Generate the list of barcodes sequence and ids
    bc_ids=list(bc_dict.keys())

    mutant_dict={i:[] for i in bc_ids} #Empty dictionary with all the bc ids as keys
#    lineage_dict=copy.deepcopy(mutant_dict) #Deep copy it to store the lineage counts

    read2_bc_dict=run_cutadapt(read2,"CGAGCTCGAATTCATCGAT","CTACGAGACCGACACCG")#Extract the mutant bc from the reads
    sys.stderr.write(str(len(read2_bc_dict.keys()))+"\n")

    with Pool(10) as pool:
        assign_mutant_partial=partial(assign_mutant,mutant_dict,read2_bc_dict,bc_dict,bc_list,bc_ids)
        res=pool.map(assign_mutant_partial,read2_bc_dict.keys())
        pool.close()
        pool.join()
    return res
    # for dicti in res:
    #     for k,v in dicti.items():
    #         if len(dicti[k])>0:
    #             mutant_dict[k]=v
    # mutant_dict={k:v for x in res for k,v in x.items()}

    mutant_dict=dict([(k,v) for dicti in res for k,v in dicti.items() if len(v)>0])
    mutant_dict=dict([(k,v) for k,v in mutant_dict.items() if len(v)>0])
    return mutant_dict

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

def assign_mutant(mutant_dict,read2_bc_dict,bc_dict,bc_list,bc_ids,key):
    if len(read2_bc_dict[key])==0:
        return None
    highest_match=process.extractOne(read2_bc_dict[key],bc_dict.values())

    if not highest_match:
        return None
    elif not highest_match[1]>=95:
        return None

    id_location=bc_list.index(highest_match[0])
    id=bc_ids[id_location]

    mutant_dict[id].append(key)
    return mutant_dict
