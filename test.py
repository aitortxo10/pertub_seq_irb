from Bio import SeqIO
import regex, sys

def detect_mutants_from_transcriptomic(reads2,bc,upstream="CATCGAT", downstream="CTACGAGA"):
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
        try:
            up_bc=regex.search("("+upstream+"){s<=1}",read2_seq).span()[1] #Search for the end of the region immediately up the genotype bc (allows up to 1 substitution)
        except AttributeError:
            sys.stderr.write("No upstream detected for read {0}\n{1}".format(sequence.id,read2_seq))
            continue
        try:
            down_bc=regex.search("("+downstream+"){s<=1}",read2_seq).span()[0] #Look for the start of the region downstream the mutant bc (allows up to 1 substitution)
        except AttributeError:
            sys.stderr.write("No downstream detected for read {0}\n".format(sequence.id))
            continue
        read2_dict[sequence.id]=read2_seq[up_bc:down_bc] #Extract the genotype bc
    sys.stderr.write("Read parsing done\n")

    bc_seq_list=''.join(barcode_dict.values()) #This is a dirty way to have all the barcodes together and not have to iterate over them. CAUTION: it may generate artifactual bc
    bc_values=list(barcode_dict.values()) #A list that will be used to later extract the index of the mutant id based on the sequence
    bc_keys=list(barcode_dict.keys())#A list of the mutant ids (same order as the list above)

    for key in read2_dict:
        tmp=regex.findall("("+read2_dict[key]+"){s<=1,d<=1}", bc_seq_list) #Search in the bc list for possible mutant bc coincidences with up to 1 substitution and 1 insertion
        # print(key,tmp)
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

def count_lineages(mutant_dict,read1,upstream="GCCAGCAAAACTAA",downstream="AACGCCGCCATCC"):
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
            tmp=regex.search("("+upstream+"){s<=1}",read_seq) #Search for the sequence downstream the lineage barcode with up to 1 substitution

            if tmp: #If there is a match
                lineage_start=tmp.span()[1] #Variable that defines the starting position of the lineage bc in the read
                lineage_end=regex.search("("+downstream+"){s<=1}",read_seq).span()[0] #define the start of the region downstream the lineage barcode
                lineage_bc=read_seq[lineage_start:lineage_end] #Extract the lineage bc
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
