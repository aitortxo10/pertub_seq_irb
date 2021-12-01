import subprocess, sys
from Bio import SeqIO
from io import StringIO

def detect_mutants(read2,bc):
    """Using cutadapt, extract for each barcode the sequences that belong"""

#Iniate the variables that will be needed (mutant_dict will be the dictionary containing the mutant id --> sequences id relationship)
    barcodes=SeqIO.parse(bc,"fasta")
    barcode_dict=dict()
    mutant_dict=dict()
    read2_dict=dict()

    for sequence in barcodes:
        barcode_dict[sequence.id]=str(sequence.seq)[-37:-17] #In the noempty fasta files, the barcode is between these positions (17 last bp are endogenous terminator)
    sys.stderr.write("Barcode parsing done\n")

    for bc in barcode_dict:
        cmd="cutadapt -j 10 --rc --discard-untrimmed -a " + barcode_dict[bc] +" --action=none "+ read2
        cutadapt_res=subprocess.run(cmd, shell=True, capture_output=True)
        reparse=cutadapt_res.stdout.decode("utf-8")
        
        with StringIO(reparse) as fastqIO:
            tmp=SeqIO.parse(fastqIO,"fastq")

            for sequence in tmp:
                if not bc in mutant_dict:
                    mutant_dict[bc]=[sequence.id]
                else:
                    mutant_dict[bc].append(sequence.id)
            sys.stderr.write("Barcode {0} done\n".format(bc))
    return mutant_dict
