from Bio import SeqIO
import regex, concurrent.futures, sys
#
# parser=argparse.ArgumentParser(description="Extract R1 paired to clean R2 reads\n")
#
# parser.add_argument('-r2','--read2',
#                     type=str,
#                     required=True,
#                     action="store",
#                     dest="read2",
#                     help= 'path to the R2 file')
#
# parser.add_argument('-r1','--read1',
#                     type=str,
#                     required=True,
#                     action="store",
#                     dest="read1",
#                     help="path to the R1 file")
#
# parser.add_argument('-o','--outfile',
#                     type=str,
#                     action="store",
#                     dest="outname",
#                     default="lineage_clone_count",
#                     help='output file name prefix')
#
# parser.add_argument('-b','--barcodes',
#                     type=str,
#                     required=True,
#                     action="store",
#                     dest="barcodes",
#                     help="path to the barcode file")
#
# parser.add_argument('-x','--bc_reference',
#                     type=str,
#                     default=" AACGCCGCCATCCAGTGTCGAAAACGAGCTCGAATTCATCGAT",
#                     action="store",
#                     dest="ref_point",
#                     help="DNA string that will be sought after to localise the barcodes.")
#
# parser.add_argument('-g','--genomic',
#                     default=False,
#                     action="store_true",
#                     dest="geno",
#                     help="Use this flag in case the read comes from genomic data.")
#
# args=parser.parse_args()

def detect_mutants_from_genomic(reads,bc):
    """Given a read file sequenced from genomic data and a barcode mapping fasta, return a dictionary with the genotype as key and the ids of such mutant as values"""

#Iniate the variables that will be needed (mutant_dict will be the dictionary containing the mutant id --> sequences id relationship)
    barcodes=SeqIO.parse(bc,"fasta")
    read2=SeqIO.parse(reads,"fastq")
    barcode_dict=dict()
    mutant_dict=dict()
    read2_dict=dict()

#Save the barcode sequences in a dictionary
    for sequence in barcodes:
        barcode_dict[sequence.id]=str(sequence.seq)[-37:-17] #In the noempty fasta files, the barcode is between these positions (17 last bp are endogenous terminator)

    sys.stderr.write("Barcode parsing done")

#Process read2 and generate the mutant_dict
    for sequence in read2:
        read2_seq=str(sequence.reverse_complement().seq)[11:31] #In genomic DNA, the read bc is located between bp 21 and, but this is the reverse complementary.
        read2_dict[sequence.id]=read2_seq

    sys.stderr.write("Read parsing done")

    bc_seq_list=''.join(barcode_dict.values()) #This is a dirty way to have all the barcodes together and not have to iterate over them. CAUTION: it may generate artifactual bc
    bc_values=list(barcode_dict.values()) #A list that will be used to later extract the index of the mutant id based on the sequence
    bc_keys=list(barcode_dict.keys())#A list of the mutant ids (same order as the list above)

    for key in read2_dict.keys():
        tmp=regex.findall(r"("+read2_dict[key]+"){s<=1}", bc_seq_list) #Search in the barcode list for all possible coincidences of the read bc with up to 1 substitution
        if tmp: #If the bc is found anywhere

            for hit in tmp:
                if not hit in barcode_dict.values(): #If there is any match at all, check it is not an artifact (just in case, this bit can perhaps be removed)
                    continue
                bc_ind=bc_values.index(hit) #Search for the mutant id
                if not bc_keys[bc_ind] in mutant_dict: #Save the read id as a new key-list variable if it does not exist. Otherwise, append it where it belongs
                    mutant_dict[bc_keys[bc_ind]]=[key]
                else:
                    mutant_dict[bc_keys[bc_ind]].append(key)
                sys.stderr.write("Key {0} done".format(key))

    return mutant_dict

        # for mutant_bc in barcode_dict.keys():
        #     print("Analysis of mutant bc: "+mutant_bc)
        #     if regex.findall(r"("+str(barcode_dict[mutant_bc])+"){s<=1}",read2_seq[11:31]):
        #         if not mutant_bc in mutant_dict:
        #             mutant_dict[mutant_bc]=[sequence.id]
        #         else:
        #             mutant_dict[mutant_bc].append(sequence.id)

    # return mutant_dict

def count_lineages(mutant_dict,read1,downstream=""):
    """
    Given a dictionary with the relationship between read2 mutant id and read id, extract and count the lineages from read1 file. Return a dictionaty with
    mutant id_n(lineage) as keys and sequence ids as values.
    """
    #Create the iterator that will allow us to read the read1 file and initiate the variables that will be needed
    read1_iter=SeqIO.parse(read1,"fastq")
    lineage_count=dict()
    read2_ids=list(mutant_dict.values())
    mutant_bcs=list(mutant_dict.keys())

    for sequence in read1_iter: #Loop through every sequence in the file
        if sequence.id in mutant_dict.values(): # If the paired read2 was associated to a mutant
            read_seq=str(sequence.seq) #Save the sequence
            tmp=regex.search("("+downstream+"){s<=1}",read_seq) #Search for the sequence downstream the lineage barcode with up to 1 substitution

            if tmp: #If there is a match
                lineage_start=tmp.span()[0]-5 #Variable that defines the starting position of the lineage bc in the read
                lineage_bc=read_seq[lineage_start:tmp.span()[0]] #Extract the lineage bc
                mutant=read2_ids.index(str(sequence.id)) #This and the next line are to identify the mutant the read corresponds to
                mutant=mutant_bcs[mutant]
                if not mutant in lineage_count: #If it is the first time the mutant is seen
                    lineage_count[mutant]=[lineage_bc] #Initiate the dictionary entry
                else:
                    if not lineage_bc in lineage_count[mutant]: #If there is no record of the lineage barcode in this mutant, save the lineage bc
                        lineage_count[mutant].append(lineage_bc)

                n_mutant=lineage_count[mutant].index(lineage_bc) #Get the lineage number for this mutant and generate a key
                new_id=mutant+str(n_mutant)

                if not new_id in lineage_dict: #Check if this lineage is already present. If it is, add the sequence id. Otherwise, create an entry
                    lineage_dict[new_id]=[sequence.id]
                else:
                    lineage_dict.append(sequence.id)

    return lineage_dict

def detect_mutants_from_transcriptomic(reads2,bc,upstream="AATTCATCGAT", downstream="CTACGAGA"):
    """Given a read file sequenced from transcriptomic data and a barcode mapping fasta, return a dictionary with the genotype as key and the ids of such mutant as values"""

#Iniate the variables that will be needed (mutant_dict will be the dictionary containing the mutant id --> sequences id relationship)
    barcodes=SeqIO.parse(bc,"fasta")
    read2=SeqIO.parse(reads2,"fastq")
    barcode_dict=dict()
    mutant_dict=dict()
    read2_dict=dict()

#Save the barcode sequences in a dictionary
    for sequence in barcodes:
        barcode_dict[sequence.id]=str(sequence.seq)[-37:-17] #In the noempty fasta files, the barcode is between these positions (17 last bp are endogenous terminator)
    sys.stderr.write("Barcode parsing done\n")

#Process read2 and
    for sequence in read2:
        read2_seq=str(sequence.reverse_complement().seq) #In genomic DNA, the read bc is located between bp 21 and, but this is the reverse complementary.
        up_bc=regex.search("("+upstream+"){s<=1}",read2_seq).span()[1] #Search for the region immediately up the genotype bc
        down_bc=regex.search("("+downstream+"){s<=1}",read2_seq).span()[0] #Search for the downstream region start
        read2_dict[sequence.id]=read2_seq[up_bc:down_bc] #Extract the genotype bc
    sys.stderr.write("Read parsing done\n")

    bc_seq_list=''.join(barcode_dict.values()) #This is a dirty way to have all the barcodes together and not have to iterate over them. CAUTION: it may generate artifactual bc
    bc_values=list(barcode_dict.values()) #A list that will be used to later extract the index of the mutant id based on the sequence
    bc_keys=list(barcode_dict.keys())#A list of the mutant ids (same order as the list above)

    for key in read2_dict.keys():
        tmp=regex.findall("("+read2_dict[key]+"){s<=1,d<=1}", bc_seq_list) #Search in the barcode list for all possible coincidences of the read bc with up to 1 substitution
        # print(tmp)
        if tmp:
            for hit in tmp:
                if not hit in barcode_dict.values(): #If there is any match at all, check it is not an artifact (just in case, this bit can perhaps be removed)
                    continue
                bc_ind=bc_values.index(hit) #Search for the mutant id
                if not bc_keys[bc_ind] in mutant_dict: #Save the read id as a new key-list variable if it does not exist. Otherwise, append it where it belongs
                    mutant_dict[bc_keys[bc_ind]]=[key]
                else:
                    mutant_dict[bc_keys[bc_ind]].append(key)
                sys.stderr.write("Key {0} done\n".format(key))

    return mutant_dict

# if __name__=="__main__":
#
#     if geno == False:
#         step1=detect_mutants_from_transcriptomic(args.read2, args.barcodes,args.ref_point)
#
#     elif geno==True:
#         step1=detect_mutants_from_genomic(args.read2,args.barcodes)
#
#     else:
#         raise ValueError("Assigment of transcriptomic/genomic does not work \n")
#
#     step2=count_lineages(step1,args.read1,args.ref_point)
    # with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
    #     mutants=executor.map(detect_mutants(read2,barcodes))
