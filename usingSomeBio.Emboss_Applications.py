from Bio import SeqIO #Allows me to use SeqIO.parse to easily read a FASTA file, as well as record.seq and record.ID to extract the sequence and ID portion of a FASTA file, respectively.
from Bio import Entrez #Allows me to use e-utilities like elink, esearch, and efetch in order to get sequences (nucleotide and amino acid) from GenBank.
from Bio.Blast import NCBIWWW #Allows me to access BLAST.
from Bio.Emboss.Applications import FDNADistCommandline #Allows me to create a distance matrix, which is required for FNeighborCommandline.
from Bio.Emboss.Applications import FNeighborCommandline #Allows me to create a phylogentic tree based on the distance matrix created with FDNADistCommandline.
from Bio.Emboss.Applications import NeedleallCommandline #Allows me to align nucleotide and amino acid sequences.
from Bio.Emboss.Applications import IepCommandline #Allows me to use Bio.Emboss to calculate the isoelectric point of proteins.

##The first part of the code uses Entrez to extract nucleotide sequences from GenBank. 
Entrez.email = "desnoyej@uoguelph.ca" #This line is required in order to use Entrez functions. 
query = input("What would you like to search for in Entrez? ") #This line allows the user to search GenBank for whatever they would like at the time. 
result = Entrez.esearch(db="nuccore",term=query) #Entrez.esearch uses the above query to search GenBank and pull the ID of sequences matching the query. 
record = Entrez.read(result) #In order to actually aqcuire the IDs for sequences matching the query, the program must use Entrez.read.
IDs = record["IdList"] #This line stores the IDs that were acquired above in variable called "IDs" so that they can be used later. 

query_sequence_dict = {} #Creating an empty dictionary to store IDs as the keys, and their respective sequence as the values.

for i in range(6): #The next block of code iterates through six of the IDs acquired above and retrieves their sequences using Entrez.efetch().
    query_seqs = Entrez.efetch(db="nuccore",id=IDs[i],rettype="fasta",retmode="text") #This line actually pulls the sequences based on the IDs provided. It uses a nucleotide database
    #in order to retrieve nucleotide sequences, and returns these sequences as fasta. 
    query_record = SeqIO.parse(query_seqs,"fasta") #This line parses through the sequences 'fetched' above and stores them in a variable called "query_record" so that they can be accessed later.
    query_key = IDs[i] #This line stores the IDs in the 'IDs' variable to another variable named 'query_key'. The 'query_key' will later be used to define the keys in the dictionary
    #called 'query_sequence_dict'. 
    for j in query_record: #This line iterates through the sequences stored in the 'query_record' variable. 
        query_value = j.seq #This line stores the sequences in the 'query_record' variable to another variable named 'query_value'. The 'query_value' will later be used to define the values
        #in the dictionary called 'query_sequence_dict'. 
        query_sequence_dict[query_key] = query_value #This line of code actually assigns the keys and values defined above to the 'query_sequence_dict' so that the dictionary will contain IDs as the keys
        #and their associated sequences as the values.

with open("query_searchSequences_dna.fasta",'w') as outfile_sequences: #This line opens the file called "query_searchSequences_dna.fasta" and assigns it to a variable called 'outfile_sequences'.
    #The 'with' portion will make sure that the "query_searchSequences_dna.fasta" file is closed after the code below this line has been executed. 
    for key, value in query_sequence_dict.items(): #This line iterates through each key:value pair from the 'query_sequence_dict' dictionary one at a time.
        outfile_sequences.write('>' + str(key) + '\n' + str(value) + '\n' + '\n') #This line writes a '>' character followed by a key from the 'query_sequence_dict' dictionary on one line and then
        #the corresponding value on the next line and then a blank line to separate key:value pairs to the file opened above. 

##The next part of the code creates the input file for FDNADistCommandline.
speciesList = [] #This is an empty list that will eventually store species names (or FASTA ID's aqcuired using record.id).
with open("query_searchSequences_dna.fasta",'r') as infile_sequences:
    for record in SeqIO.parse("query_searchSequences_dna.fasta", "fasta"): 
        speciesList.append(record.id) #This line uses the SeqIO module's record.id to exract the FASTA ID for the current record and then appends this ID to the list called 'speciesList'.

numberOfSpecies = len(speciesList) #This line finds the number of species using the length function on the 'speciesList' list. Knowing the number of species is required by FDNADistCommandline as the
#first part of the first line of the input file. 

        
#The next few lines get the number of positions to be analyzed. This will be set to the shortest sequence that was pulled from Entrez. All other sequences will be trimmed to this length
#because most Bio.Emboss applications require that all of the sequences be the same length. This part of the code will also restrict the length of the sequences written to the file so that the
#pipeline is less computationally intensive.
sequenceList = [] #An empty list that will eventually hold sequences. 
with open("query_searchSequences_dna.fasta",'r') as infile_sequences:
    for record in SeqIO.parse("query_searchSequences_dna.fasta", "fasta"):
        sequenceList.append(record.seq) #This line uses the SeqIO module's record.seq to exract the FASTA sequence for the current record, and then appends the sequence to the 'sequenceList'.

#The next few lines use a for loop to iterate through each sequence in 'sequenceList', and for each sequence, uses the length function to get the length and then appends to 'sequenceLengths' list.
sequenceLengths = [] 
for i in sequenceList:  
    sequenceLengths.append(len(i))

smallestSeqLength = min(sequenceLengths) #Uses the min function to find the smallest length in 'sequenceLengths'. This will become the maximum length of all sequences put into the
#'FormattedInputFile.txt' file.

#The next few lines use an if/else statement to make sure that if the smallest sequence length is less than 1000, then that will be the maximum length used. And if the smallest sequence is greater
#than 1000 then all sequences will be cut off at 1000 nucleotides. Again, this is to make the pipeline less computationally intensive. 
if int(smallestSeqLength) < 1000:
    numberOfPositions = smallestSeqLength
else:
    numberOfPositions = 1000

#The next few lines create the input file required by FDNADistCommandline. The first line is four spaces, number of species (stored in the variable 'numberOfSpecies'), another four spaces and then
#number of positions to be analysed (stored in the variable 'numberOfPositions'). Each line after that contains the species' name (which can only be a maximum of 10 characters, and if the name is less
#than 10, then the remaining characters must be filled with spaces), followed by the species sequence. 
with open("FormattedInputFile.txt",'w') as infile_Formatted:
    infile_Formatted.write(' ' + ' ' + ' ' + ' ')
    infile_Formatted.write(str(numberOfSpecies)) #The str() function is required because the elements in 'numberOfSpecies' are Seq objects, which Python will not write to a file. 
    infile_Formatted.write(' ' + ' ' + ' ' + ' ')
    infile_Formatted.write(str(numberOfPositions))
    infile_Formatted.write("\n")
    with open("query_searchSequences_dna.fasta",'r') as infile_sequences:
        for record in SeqIO.parse("query_searchSequences_dna.fasta", "fasta"):
            speciesName = '{:<10}'.format(record.id[:10]) #This line ensures that any species' name that is greater than 10 characters, gets cut to 10 characters. And any species' name that is less than
            #10 characters, gets increased to 10 characters using spaces. 
            speciesSequence = record.seq[:numberOfPositions] #This line cuts each sequence to the length stored in 'numberOfPositions' so that each sequence is the same length.
            infile_Formatted.write(str(speciesName))
            infile_Formatted.write(str(speciesSequence))
            infile_Formatted.write("\n")


##The next part of the code uses FDNADistCommandline and the "FormattedInputFile.txt" file created above to create a distance matrix.
#The first part is "/usr/local/bin/fdnadist", which tells Python where it can find this application. The sequence parameter tells FDNADistCommandline which file you would like to use that contains
#the sequences you want to compute a distance matrix for. The method parameter tells FDNADistCommandline which distance matrix algorithm you would like to use; in this case, f means FDNADistCommandline
#will use the F84 distance model. The outfile parameter specifies the name of the outfile in which FDNADistCommandline will write the results to.
FDNADist_matrix = FDNADistCommandline("/usr/local/bin/fdnadist",sequence="FormattedInputFile.txt",method="f",outfile="distanceMatrix.fdnadist")

stdout, stderr = FDNADist_matrix() #This line is required in order for FDNADistCommandline to actually write anything to the outfile. 


##The next part of the code uses FNeighborCommandline and the "distanceMatrix.fdnadist" file created above to create a phylogenetic tree.
#The first part is "/usr/local/bin/fneighbor", which tells Python where it can find this application. The datafile parameter tells FNeighborCommandline which file you would like to use that contains
#the distance matrix that will be used to create the phylogenetic tree. The outfile parameter specifies the name of the outfile in which FNeighborCommandline will write the results to.
FNeighbor_tree = FNeighborCommandline("/usr/local/bin/fneighbor",datafile="distanceMatrix.fdnadist",outfile="treeFile.fneighbor")

stdout, stderr = FNeighbor_tree() #This line is required in order for FNeighborCommandline to actually write anything to the outfile. 


##The next part of the code creates alignments between each and every one of the sequences in the "query_searchSequences_dna.fasta" using the NeedleallCommandline.
#The next few lines create an input file for NeedleallCommandline. The FASTA file called 'query_searchSequences_dna.fasta' cannot be used because it contains spaces between each record and for
#some reason the spaces affect the Needleall alignemnt. So the next few lines create a FASTA without spaces between each record called 'needlallAlignmentInput_Nucleotide.fasta'. The sequences are also all
#made the same length using the slie function -> [:numberOfPositions].
with open("needlallAlignmentInput_Nucleotide.fasta",'w') as outfile_alignment:
    with open("query_searchSequences_dna.fasta",'r') as infile_sequences:
        for record in SeqIO.parse("query_searchSequences_dna.fasta", "fasta"):
            sequence = record.seq[:numberOfPositions]
            outfile_alignment.write('>' + str(record.id) + '\n' + str(sequence) + '\n') 


#The next few lines uses the NeedleallCommandline application. The first parameter is "/usr/local/bin/needleall", which tells Python where it can find this application. The asequence parameter
#tells NeedleallCommandline which file you would like to use that contains the sequences to be aligned with the sequences in the bsequence parameter. The bsequence parameter tells NeedleallCommandline
#which file you would like to use that contains the sequences to be aligned with the sequences in the asequence parameter. The outfile parameter specifies the name of the outfile in which
#NeedleallCommandline will write the results to.
needleall_align = NeedleallCommandline("/usr/local/bin/needleall",asequence="needlallAlignmentInput_Nucleotide.fasta",bsequence="needlallAlignmentInput_Nucleotide.fasta",gapopen=10,gapextend=0.5,outfile="nucleotide_Alignment.needleall")

stdout,stderr = needleall_align() #This line is required in order for NeedleallCommandline to actually write anything to the outfile. 


##The next part of the code gets the linked protein sequences for each of the above queried nucleotide sequences.
link_protein = Entrez.elink(db='protein',dbfrom='nuccore',id=IDs) #This line uses Entrez.elink() to find the protein IDs that are linked to the sequence ID's from part one. 
record_link = Entrez.read(link_protein) #This line extracts the results that ELink has found.

query_protein_IDs = [] #An empty list to store protein IDs acquired below.

#The next few lines assign the IDs extracted above to variables, and then appends these variables one by one to the 'query_protein_IDs' list. 
proteinID1 = record_link[1]['LinkSetDb'][0]['Link'][0]['Id']
query_protein_IDs.append(int(proteinID1))
proteinID2 = record_link[2]['LinkSetDb'][0]['Link'][0]['Id']
query_protein_IDs.append(int(proteinID2))
proteinID3 = record_link[3]['LinkSetDb'][0]['Link'][0]['Id']
query_protein_IDs.append(int(proteinID3))
proteinID4 = record_link[4]['LinkSetDb'][0]['Link'][0]['Id']
query_protein_IDs.append(int(proteinID4))
proteinID5 = record_link[5]['LinkSetDb'][0]['Link'][0]['Id']
query_protein_IDs.append(int(proteinID5))
proteinID6 = record_link[6]['LinkSetDb'][0]['Link'][0]['Id']
query_protein_IDs.append(int(proteinID6))

query_protein_dict = {} #Creates an empty dictionary that will later store the linked protein IDs as the keys and the corresponding amino acid sequences as the values
#so that they can later be written to a file.


#For the next two blocks of code, refer to the nucleotide section for what each code line is doing. This part is essentially the same except the program is using linked protein IDs to
#acquire the corresponding amino acid sequences, rather than using ID's to acquire the corresponding nucleotide sequences of the CDS.
for i in range(6):
    query_proteins = Entrez.efetch(db="protein",id=query_protein_IDs[i],rettype="fasta",retmode="text") 
    query_protein_record = SeqIO.parse(query_proteins,"fasta")
    query_protein_key = query_protein_IDs[i]
    for j in query_protein_record:
        query_protein_value = j.seq
        query_protein_dict[query_protein_key] = query_protein_value

with open("query_searchSequences_linked_aa.fasta",'w') as outfile_proteins:
    for key, value in query_protein_dict.items():
        outfile_proteins.write('>' + str(key) + '\n' + str(value) + '\n' + '\n')

#The next few lines get the number of positions to be analyzed. This will be set to the shortest sequence that was pulled from Entrez. All other sequences will be trimmed to this length
#because most Bio.Emboss applications require that all of the sequences be the same length. This part of the code will also restrict the length of the sequences written to the file so that the
#pipeline is less computationally intensive.
aa_sequenceList = [] #An empty list that will eventually hold sequences. 
with open("query_searchSequences_linked_aa.fasta",'r') as infile_sequences:
    for record in SeqIO.parse("query_searchSequences_linked_aa.fasta", "fasta"):
        aa_sequenceList.append(record.seq) #This line uses the SeqIO module's record.seq to exract the FASTA sequence for the current record, and then appends the sequence to 'aa_sequenceList'.

#The next few lines use a for loop to iterate through each sequence in 'aa_sequenceList', and for each sequence, uses the length function to get the length
#and then appends this length to the 'aa_sequenceLengths' list.
aa_sequenceLengths = [] 
for i in aa_sequenceList:  
    aa_sequenceLengths.append(len(i))

aa_smallestSeqLength = min(aa_sequenceLengths) #Uses the min function to find the smallest length in 'aa_sequenceLengths'. This will become the maximum length of all sequences put into the
#'needlallAlignmentInput_AminoAcid.fasta' file

#The next few lines use an if/else statement to make sure that if the smallest sequence length is less than 1000, then that will be the maximum length used. And if the smallest sequence is greater
#than 1000 then all sequences will be cut off at 1000 amino acids. Again, this is to make the pipeline less computationally intensive. 
if int(aa_smallestSeqLength) < 1000:
    aa_numberOfPositions = aa_smallestSeqLength
else:
    aa_numberOfPositions = 1000

##The next part of the code creates alignments between each and every one of the amino acid sequences in the "query_searchSequences_linked_aa.fasta" file using the NeedleallCommandline.
#Everything from creating the input file to running NeedleallComandline is the same as what was done with the nucleotide sequences above. 
with open("needlallAlignmentInput_AminoAcid.fasta",'w') as outfile_alignment:
    with open("query_searchSequences_linked_aa.fasta",'r') as infile_sequences:
        for record in SeqIO.parse("query_searchSequences_linked_aa.fasta", "fasta"):
            sequence = record.seq[:aa_numberOfPositions]
            outfile_alignment.write('>' + str(record.id) + '\n' + str(sequence) + '\n') 

needleall_align = NeedleallCommandline("/usr/local/bin/needleall",asequence="needlallAlignmentInput_AminoAcid.fasta",bsequence="needlallAlignmentInput_AminoAcid.fasta",gapopen=10,gapextend=0.5,outfile="aminoAcid_Alignment.needleall")

stdout,stderr = needleall_align()


##The next part of the code calculates the isoelectric point of each protein in the 'query_searchSequences_linked_aa.fasta' file using the IepCommandline.
#The first parameter is "/usr/local/bin/iep", which tells Python where to look for the Iep application. The sequence parameter tells IepCommandline which file you would like to use that contains the
#proteins that you want to calculate isoelectric points for. The outfile parameter tells IepCommandline the name of the file in which you would like it to write the results to.

iep_isoelectric = IepCommandline("/usr/local/bin/iep",sequence='query_searchSequences_linked_aa.fasta',outfile='isoelectricPoints.iep')

stdout,stderr = iep_isoelectric() #This line is required in order for IepCommandline to actually write anything to the outfile. 










        
