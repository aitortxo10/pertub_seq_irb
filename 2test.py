from Bio import SeqIO
import regex, sys
from fuzzywuzzy import process

def detect_mutants_from_transcriptomic(reads2,bc,upstream="CATCGAT", downstream="CTACGAGA"):
    """Given a read file sequenced from transcriptomic data and a barcode mapping fasta, return a dictionary with the genotype as key and the ids of such mutant as values"""

#Iniate the variables that will be needed (mutant_dict will be the dictionary containing the mutant id --> sequences id relationship)
    barcodes=SeqIO.parse(bc,"fasta")
    read2=SeqIO.parse(reads2,"fastq")
    barcode_dict=dict()
    mutant_dict=dict()
    read2_dict=dict()

    for sequence in barcodes:
        barcode_dict[sequence.id]=str(sequence.seq)[-37:-17] #In the noempty fasta files, the barcode is between these positions (17 last bp are endogenous terminator)
    sys.stderr.write("Barcode parsing done\n")

    bc_list=list(barcode_dict.values())
    mutant_id_list=list(barcode_dict.keys())

    for sequence in read2:
        read2_seq=str(sequence.reverse_complement().seq) #In genomic DNA, the read bc is located between bp 21 and, but this is the reverse complementary
        try:
            up_bc=regex.search("("+upstream+"){s<=1}",read2_seq).span()[1] #Search for the end of the region immediately up the genotype bc (allows up to 1 substitution)
        except AttributeError:
            sys.stderr.write("No upstream detected for read {0}\n{1}\n".format(sequence.id,read2_seq))
            continue
        try:
            down_bc=regex.search("("+downstream+"){s<=1}",read2_seq).span()[0] #Look for the start of the region downstream the mutant bc (allows up to 1 substitution)
        except AttributeError:
            sys.stderr.write("No downstream detected for read {0}\n".format(sequence.id))
            continue
        read2_dict[sequence.id]=read2_seq[up_bc:down_bc] #Extract the genotype bc
    sys.stderr.write("Read parsing done\n")

    for read_bc in read2_dict:
        highest_match=process.extractOne(read2_dict[read_bc], bc_list)

        if not highest_match:
            continue
        elif not highest_match[1]>=95:
            continue

        id_location=bc_list.index(highest_match[0])
        id=mutant_id_list[id_location]

        if not id in mutant_dict:
            mutant_dict[id]=[read_bc]
        else:
            mutant_dict[id].append(read_bc)

    return mutant_dict
