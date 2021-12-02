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

#Parse the barcode sequences into a dictionary where the id is the key and the values are the sequences
    for sequence in barcodes:
        barcode_dict[sequence.id]=str(sequence.seq)[-37:-17] #In the noempty fasta files, the barcode is between these positions (17 last bp are endogenous terminator)
    sys.stderr.write("Barcode parsing done\n")

#Run cutadapt for each barcode
    for bc in barcode_dict:
        cmd="cutadapt -j 10 --rc --no-indels --discard-untrimmed -a " + barcode_dict[bc] +" --action=none "+ read2 #Generate the command that will run
        cutadapt_res=subprocess.run(cmd, shell=True, capture_output=True) #Run cutadapt from python

        if not cutadapt_res.stdout: #If there are no matches, skip the iteration
            continue

        reparse=cutadapt_res.stdout.decode("utf-8") #Get the stdout, which in this case is the sequences that contain the barcode, and transform them to a string

        with StringIO(reparse) as fastqIO: #Small cheat to be able to run Biopython parse over a string
            tmp=SeqIO.parse(fastqIO,"fastq") #Generate the parser

            mutant_dict[bc]=list() #Initiate the dictionary entry

            for sequence in tmp: #Store the sequences in the dictionary
                mutant_dict[bc].append(sequence.id)

            # for sequence in tmp: #Store the sequences in the dictionary
            #     if not bc in mutant_dict:
            #         mutant_dict[bc]=[sequence.id]
            #     else:
            #         mutant_dict[bc].append(sequence.id)

            sys.stderr.write("Barcode {0} done\n".format(bc))
    return mutant_dict
