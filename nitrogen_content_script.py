
"""
Created on Wed Oct  21 18:50:23 2020

author: Thibaut Goldsborough

Please email: tg76@st-andrews.ac.uk for any queries 

This code calculates the nitrogen content of the 
genome, CDS, non-CDS, mRNA, protein, introns, exons and tRNA of Genlisea aurea, Arabidopsis thaliana and
Trifolium pratense. 

This script takes about 20 mins to run on a 1,4 GHz Intel Core i5. (The reader might prefer to run the code in sections, rather than all at once.)

In this script, the word 'domain' refers to a genome contig. In the 
complete genome of A. thaliana, this refers to individual chromosomes.

Details of how information is stored in gff files: 
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

"""


import numpy as np         #Numpy is used for calculating mean and variances of lists efficiently
from Bio.Seq import Seq    #Bio.Seq is used to help handle fasta files and converting sequences
import scipy.stats         #scipy.stats is used to calculate confidence intervals



#This function was taken from stackoverflow user gcamargo, https://stackoverflow.com/questions/15033511/compute-a-confidence-interval-from-sample-data
#The function calculates the mean and a 95% confidence interval of a list  
def mean_confidence_interval(data, confidence=0.95):  
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, h


#This function calculates the nitrogen content of DNA, RNA or protein.
#sequence must be 'string' when sequence_stype is 'DNA' or 'RNA'
#sequence must be 'list' of 'strings' when sequence_stype is 'Protein' where each element of the list is one gene to be translated

def Calculate_Nitrogen_Content(sequence,sequence_type): 
     

    if sequence_type=="DNA": #Dealing with DNA sequences
        
         nitrogen_sequence=sequence.replace("A",'7').replace("T",'7').replace("C",'8').replace("G",'8') 
         #Make a new string of numbers corresponding to the number of nitrogen atoms for each DNA base pair. 
          
     
    if sequence_type=="RNA": #Dealing with RNA sequences
         
         gene_seq=Seq(sequence) #Convert string sequence into a Bio.Seq sequence

         sequence=str(gene_seq.transcribe()) #Convert to RNA using Bio.Seq transcribe() function
        
         nitrogen_sequence=sequence.replace("A",'5').replace("U",'2').replace("C",'3').replace("G",'5')
         #Make a new string of numbers corresponding to the number of nitrogen atoms for each RNA nucleotide.
    
    if sequence_type=="Protein": #Dealing with RNA sequences
        
        #Dictionary that assigns the number of nitrogen atoms to each amino acid
        nitrogen_table={'R': 4,'H': 3,'K': 2,'D': 1,'E': 1,'S': 1,'T': 1,'N': 2,'Q': 2,'C': 1,'U': 1,'G': 1,'P': 1,'A': 1,'V': 1,'I': 1,'L': 1,'M': 1,'F': 1,'Y': 1,'W': 2}

        nitrogen_sequence=str() #Create a new string that will contain the number of nitrogen atoms for each amino acid
        for gene in sequence: #sequence is a list and gene is a string
              
            gene_seq=Seq(gene) #Convert string sequence into a Bio.Seq sequence
    
            mRNA_sequence=gene_seq.transcribe()  #Convert to RNA using Bio.Seq transcribe() function
            
            protein_sequence=mRNA_sequence.translate() #Convert to protein using Bio.Seq translate() function
            
            for amino_acid in protein_sequence: #iterate over each amino acid in the protein sequence
            #Unfortunately for an unclear reason some CDS sequences contain stop codons within the sequence.
            #This prevents the use of the faster .replace() method used for DNA and RNA
            
                if amino_acid=='*': #Stop at stop codon
                    break
    
                if amino_acid in nitrogen_table.keys(): #Make sure amino acid corresponds to one of the keys of nitrogen_table, this gets rid of unknown amino acids 
                    nitrogen_sequence+=str(nitrogen_table[amino_acid]) 
                    #update the nitrogen_sequence string with the number of nitrogen atoms
                    #corresponding to the amino_acid


            

         
    cleaned_sequence = [ int(base) for base in nitrogen_sequence if str(base).isdigit() ]
    #This removes all elements that were not converted to nitrogen counts
    #This is useful because fasta sequences are not only composed of A,G,C,T (e.g. 'Y', 'W' etc...)
    #This line also converts string to a list of integers
    
    #Example given for nitrogen_sequence = "5587Y533"
    #The corresponding cleaned_sequence =[5,5,8,7,5,3,3]
    
    
    relative_nitrogen,interval=mean_confidence_interval(cleaned_sequence)
    #Calculate the mean number of nitrogen atoms and its 95% confidence interval
    #relative_nitrogen is the average number of nitrogen atoms per element
    #relative_nitrogen + or - interval is the 95% confidence interval of the mean

    total_nitrogen=np.sum(cleaned_sequence) #Calculate the sum of all nitrogen atoms in the sequence
    
    length=len(cleaned_sequence) #Get the number of elements in the sequence 

    return(total_nitrogen,relative_nitrogen,interval,length) 


        


Species_Names=("Genlisea","Arabidopsis","Trifolium") #List of the three species names to study

for species in Species_Names: #Iterate over each species 
    
    if species=="Genlisea": #Load the appropriate gff, fna and cds files

        gff_file="./Sequences/GCA_000441915.1_GenAur_1.0_genomic.gff"
        fna_file="./Sequences/GCA_000441915.1_GenAur_1.0_genomic.fna"
        cds_file="./Sequences/GCA_000441915.1_GenAur_1.0_cds_from_genomic.fna"
        trnas_file="./Sequences/genlisea_trnas"
       
            
    if species=="Arabidopsis":
        
        
       gff_file="./Sequences/GCA_904420315.1_AT9943.Cdm-0.scaffold_genomic.gff"
       fna_file="./Sequences/GCA_904420315.1_AT9943.Cdm-0.scaffold_genomic.fna"
       cds_file="./Sequences/GCA_904420315.1_AT9943.Cdm-0.scaffold_cds_from_genomic.fna"
       trnas_file="./Sequences/arabidopsis_trnas"
        
    if species=="Trifolium":

        gff_file="./Sequences/GCA_000583005.2_Tp1.0_genomic.gff"
        fna_file="./Sequences/GCA_000583005.2_Tp1.0_genomic.fna"
        cds_file="./Sequences/GCA_000583005.2_Tp1.0_cds_from_genomic.fna"
        trnas_file="./Sequences/trifolium_trnas"

    print("") #New line
    print("")
    print("#########################",species,"#########################")
  
        
    
    #PART I
    
    ########## Calculate total genomic nitrogen content ############
    
    #In this part, the number of nitrogen atoms in the genomic fasta file will be determined
     
    with open(fna_file,'r') as f: 
        raw_genome_str = f.read() #Open file and store contents in genome_str 
    
    
    raw_genome=raw_genome_str.split('>') #Split the genome at every occurence of '>' (corresponding to a new sequence domain)
    
    genomic_str=str()  #This string will contain the 'cleaned' fasta sequences
    genomic_dict={}    #This dict will assign a string DNA sequence to each domain name
    for line in raw_genome: #Iterate over each genomic domain
        if line!="": #Some lines were found to be empty
        
            split_line=list(line) #convert line (a string) to a list
            
            nucleotide_sequence=line[split_line.index("\n"):] 

            #Find where the first "\n" element occurs, anything before corresponds to the 
            #header and anything after corresponds to the fasta sequence
            
            nucleotide_sequence=nucleotide_sequence.replace('\n','')  #Remove all "\n" elements
            
            genomic_str+=nucleotide_sequence.upper() #Capitalize the sequence (e.g. "aaAAtA" converts to the homogeneous "AAAATA")
            
            domain_name=''.join(map(str, split_line[:10])) #The ''.join(map(str,list) method converts a list to a string
            #Unless the domain name starts with an "A", all domain names are 10 letters long (this line may fail for other genomes)
            
            if domain_name[0]=='A': 
                domain_name=''.join(map(str, split_line[:14])) #When the domain name starts with an "A" the name is 14 letters long
            
            if  domain_name not in genomic_dict:
                genomic_dict[domain_name]=[] #If domain name is a new domain name, then define a new empty list 
            genomic_dict[domain_name].append(nucleotide_sequence.upper()) #Fill the list with the corresponding string sequences. 
            
            #We now have genomic_dict which contains a string sequence for every domain in the scaffold genome
            #We also have genomic_str which is a string containing all sequences in the genome
    
        
    
    #In the next line, we calculate the nitrogen content of genomic_str
    genomic_total_nitrogen,relative_nitrogen,interval,genomic_length=Calculate_Nitrogen_Content(genomic_str,"DNA") 
    
    
    print("")
    print("Length genomic DNA:",genomic_length)
    print("The genomic nitrogen was:",genomic_total_nitrogen,"nitrogen atoms")
    print("The relative genomic nitrogen was:",relative_nitrogen,"+ or -",interval,"nitrogen atoms per basepair")

    genomic_str=[] #Free up some RAM (this variable is no longer used)

    
    #Part II
    
    ########## Calculate total CDS nitrogen content ############
    
    
    #In this part, the nitrogen content of the CDS is calculated. 
    
    with open(cds_file,'r') as f:
        cds_raw_str = f.read().replace("\n","") #Open file and store contents in cds_raw_str, remove "\n" elements 
    
        
        cds_seq=cds_raw_str.split(">lcl") #Split the raw string at each ">lcl" element corresponding to a new gene, obtained is a new list
    
    cds_str=str() #This variable will contain the sequences of all the CDS
    
    cds_list=[] #This variable will be a list of strings of each gene (used for protein nitrogen calculation later on)
    
    for cds_sequence in cds_seq:   #Iterate of ever every gene in the list   
        pos = cds_sequence.find('[gbkey=CDS]') + 11 #Find where the actual protein sequence starts 
        cds_str+=(cds_sequence[pos:]) #Only add the coding sequence part to the variable cds_str (i.e. remove the header)
        cds_list.append((cds_sequence[pos:])) #Do the same for the list version of the CDS (used for protein nitrogen content later)
    
    
     #In the next line, we calculate the nitrogen content of cds_str
    cds_total_nitrogen,relative_nitrogen,interval,cds_length=Calculate_Nitrogen_Content(cds_str,"DNA") 
    print("")
    print("Length of CDS:",cds_length)
    print("The CDS nitrogen was:",cds_total_nitrogen,"nitrogen atoms")
    print("The relative CDS nitrogen was:",relative_nitrogen,"+ or -",interval,"nitrogen atoms per basepair")
  
    
    #Part III

    ########## Calculate total protein nitrogen content ############
    
    #In the next line, we calculate the nitrogen content of protein corresponding to each gene in cds_list (list of all the sequences of each gene), these sequences are converted 
    #to protein in the function Calculate_Nitrogen_Content() 

    protein_total_nitrogen,relative_nitrogen,interval,protein_length=Calculate_Nitrogen_Content(cds_list,"Protein") 
    print("")
    print("Length of all protein:",protein_length)
    print("The protein nitrogen was:",protein_total_nitrogen,"nitrogen atoms")
    print("The relative protein nitrogen was:",relative_nitrogen,"+ or -",interval,"nitrogen atoms per basepair")
   

    cds_str=[] #Free up some RAM (these variables are no longer used and are rather big)
    
    
    #Part IV
    
    ########## Calculate total nitrogen content of non coding DNA ############
    
    nc_length=genomic_length-cds_length #Define the non coding DNA length as the genomic DNA length minus the CDS length

    nc_nitrogen=genomic_total_nitrogen-cds_total_nitrogen #Define nc_nitrogen as the total number of nitrogen atoms in the non coding region 
    
    print("")
    print("Length of the non-coding genomic DNA:",nc_length)
    print("The Non-Coding DNA nitrogen was:",nc_nitrogen,"nitrogen atoms")
    print("The relative Non-Coding DNA nitrogen was:",nc_nitrogen/nc_length,"nitrogen atoms per basepair") 
    #Unfortunately this method is very fast but prevents the calculation of a 95% confidence interval
    
    
    
    #Part V
    
    ################## Introns and exons ########################
    
    
    gff_dict={} 
    #This dictionary will assign all the information of a domain in the gff file to the domain name
    #gff_dict[domain_name]=gff information of domain_name
    

    

    
    opened_gff_file=open(gff_file,'r') #read the gff_file
    
    #Generate the gff_dictionary:
        
    for line in opened_gff_file: #Iterate over every domain in the gff_file
        if len(line)>100: #This is a clumsy way of removing the multiple lines at the start and end of the file that have nothing to do with DNA sequences
            gff_str=line.rstrip().split("\t")  #Remove all the "\t" and convert string to list
            
            #gff_str[0] contains domain name
            if gff_str[0] not in gff_dict: 
                gff_dict[gff_str[0]]=[] #Initiate a new list for each domain name
         
                        
            element=gff_str
            #element[0] contains domain name
            #element[2] contains feature type (e.g. 'gene','exon')
            #element[3] and element[4] contain start and end positions of the feature 
            element[3],element[4]=int(element[3]),int(element[4])
    
            gff_dict[element[0]].append(element) #Assign these elements to the appropriate list in gff_dict
     
        
    all_exons=str() #This string will contain all the exonic sequences of the entire CDS
    all_introns=str() #This string will contain all the intronic sequences of the entire CDS
           
        
    introns_str=str() #This string will contain all the exonic sequences of a gene
    exons_str=str()   #This string will contain all the intronic sequences of a gene
          
    
    
    for domain in gff_dict.keys(): #Iterate over each domain contained in the gff dictionary 


        for element in gff_dict[domain]: #Iterate over every element contained in the domain
            if element[0] in genomic_dict:      #Check whether the domain name was found in genomic_dict defined in part I 
                                                #(recall that genomic_dict from part I assigns a domain sequence to a domain name)
                                                
                if len(gff_dict[domain])>1: #This is equal to 1 when no genes or features were annotated in this domain
                
                    if element[2]=='gene': #If a gene is found

                        gene_length=element[4]-element[3]            
                        
                        #element[3] and element[4] contain start and end positions of the gene

                        gene_start,gene_end=element[3],element[4]
                        
                        #In the following two lines, all_exons and all_introns are updated with the sequences of the PREVIOUS gene, if it is the first gene then an empty string is added
                        all_exons+=exons_str                       
                        all_introns+=introns_str
                        
                        
                        exons_str=str() #Reset exons_str
                        
                        #Define the introns as the entire gene sequence then remove the exons one by one in the following steps    
                        introns_str=genomic_dict[element[0]][0][gene_start:gene_end] 
                       
                    if element[2]=='exon':
                        #Update exons_str with the exon sequence ranging from element[3] to element[4]
                        exons_str+=genomic_dict[element[0]][0][element[3]:element[4]]
                        
                       
                        #Update introns_str by removing the exon sequence ranging from element[3] to element[4]                    
                        introns_str=introns_str.replace(genomic_dict[element[0]][0][element[3]:element[4]],"")
        
            else:
                #This line shouldn't come up
                print("WARNING:",element[0],"from the gff file was not found in the genomic fasta file") # if element[0] not in genomic_dict
                    
    #For completeness, the exonic and intronic sequences of the last gene are added. 
    all_exons+=exons_str                
    all_introns+=introns_str
            
    #mRNA is defined as the sum of the intronic sequence and the exonic sequences
    
    mRNA_total_nitrogen,relative_nitrogen,interval,mRNA_length=Calculate_Nitrogen_Content(all_introns+all_exons,"RNA") 

    print("")
    print("Length of all mRNA:",mRNA_length)
    print("The total mRNA nitrogen was:",mRNA_total_nitrogen,"nitrogen atoms")
    print("The relative mRNA nitrogen was:",relative_nitrogen,"+ or -",interval,"nitrogen atoms per basepair")
 
    
    #Calculate the nitrogen content of all the introns
    intron_total_nitrogen,relative_nitrogen,interval,intron_length=Calculate_Nitrogen_Content(all_introns,"RNA") 
  
    print("")
    
    print("Length of all introns:",intron_length)
    print("The intron nitrogen was:",intron_total_nitrogen,"nitrogen atoms")
    print("The relative intron nitrogen was:",relative_nitrogen,"+ or -",interval,"nitrogen atoms per basepair")

    
    
    #Calculate the nitrogen content of all the exons
    
    exon_total_nitrogen,relative_nitrogen,interval,exon_length=Calculate_Nitrogen_Content(all_exons,"RNA") 
 
    print("")
    
    print("Length of all exons:",exon_length)
    print("The exon nitrogen was:",exon_total_nitrogen,"nitrogen atoms")
    print("The relative exon nitrogen was:",relative_nitrogen,"+ or -",interval,"nitrogen atoms per basepair")
 
    
 
    # #Part VI
    
    # ################## tRNA sequences ########################
 
    # In this part, the nitrogen content of tRNAs identified using the tRNAscan-SE software (not available for python, 
    # the following unix command was run in the terminal:  tRNAscan-SE -o genlisea_trnas -m genlisea_stats ../GCA_000441915.1_GenAur_1.0_genomic.fna)
    
    
    with open(trnas_file,'r') as f: #The obtained trnas results file was opened
        trnas = f.read()
        
    trnas=trnas.split("\n") #Split the data at every new line
    
    trnas_list=[] #This list will contain all the information obtained from the results file such as tRNA position and anticodon
    
    trna_str=str() #This string will contain all the tRNA sequences as one single string
    
    trna_dict={} #This dictionary will assign trna sequence(s) to each corresponding codon
    
    for line in trnas:  
        trnas_list.append(line.split("\t"))  #Remove all the "\t" elements 
        
    #An example of an element in trnas_list:
        
    # Domain name          Start     End         Anticodon 
    #['KE526704.1 ', '1', '97097', '97170', 'His', 'GTG', '0', '0', '49.4', '']
    
    
    
    for trna in trnas_list[3:]: #The first 2 elements of the list is a header
        if len(trna)==10: #Some of the lines were found to be empty
            
            
            if min(int(trna[2]),int(trna[3]))==int(trna[2]): 
                # If the start pos of the gene is smaller than the end pos, then the gene is in the 5' to 3' direction.  
                # Recall that genomic_dict[domain_name][0]=fasta sequence of the domain
                # Here the domain name is trna[0]
                # trna[0]='KE526704.1 ' so trna[0][:-1]]='KE526704.1' (note the removal of the space at the end)
                # so genomic_dict[trna[0][:-1]][0][int(trna[2]):int(trna[3])] is the DNA sequence between the start and end pos of the tRNA gene. 
                
                    trna_seq=genomic_dict[trna[0][:-1]][0][int(trna[2]):int(trna[3])]
                    
                   
            else:
                # If the start pos of the gene is larger than the end pos, then the gene is in the 3' to 5' direction (reverse direction). 
                # Note the extra [::-1] term at the end that reverses the string sequence to obtain the gene sequence in the 5' to 3'        
                trna_seq=genomic_dict[trna[0][:-1]][0][int(trna[3]):int(trna[2])][::-1]
            
            
            trna_str+=trna_seq #Add the tRNA gene sequence to trna_str
              
            codon=str(Seq(trna[5]).complement().transcribe()) #trna[5] corresponds to the anticodon of the trna, in this line, a codon is obtained from the anticodon
                   
            if codon not in trna_dict: #update the dictionary to assign a trna sequence to its corresponding codon
                trna_dict[codon]=[]
            trna_dict[codon].append(trna_seq)
          
    #The following loop calculates the nitrogen content of the tRNAs and assigning them to their complimentary codon 
    trna_nitrogen_dict={}
    for codon in trna_dict.keys():
        trna_nitrogen=[]
        for trna in trna_dict[codon]:
            trna_nitrogen.append(Calculate_Nitrogen_Content(trna,"RNA")[0])
        trna_nitrogen_dict[codon]=np.mean(trna_nitrogen) #When there are multiple tRNAs for a single codon, the mean value of the nitrogen content is taken. 
             
        
    #Calculate the nitrogen content of all the trna sequences
        
    trna_total_nitrogen,relative_nitrogen,interval,trna_length=Calculate_Nitrogen_Content(trna_str,"RNA") 
 
    print("")
    
    print("Length of all tRNAs:",trna_length)
    print("The exon nitrogen was:",trna_total_nitrogen,"nitrogen atoms")
    print("The relative exon nitrogen was:",relative_nitrogen,"+ or -",interval,"nitrogen atoms per basepair")
 
    
 
    # #Part VII
    
    # ######################## Codon Usage Bias ########################
 
    
    #The following dictionary was taken from https://pythonforbiologists.com/dictionaries
    
    genetic_code={ 
    'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
    'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
    'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
    'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
    'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
    'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
    'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
    'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
    'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
    'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W'}
    
    
    codon_usage={}
    
    # Create a new dictionary that has amino acids as keys, and for every amino acid, 
    # the dictionary contains a subdictionary containing the codons for that amino acid and their counts
    # e.g. codon_usage['C']= {'UGC': O, 'UGU': O}
    
    
    for codon in genetic_code:  
        if genetic_code[codon] in codon_usage: 
            codon_usage[genetic_code[codon]][codon]=0
    
        else:
            codon_usage[genetic_code[codon]]={}
            codon_usage[genetic_code[codon]][codon]=0
    
    
    
    #This dictionary is somewhat more practical as it will have 'codon' as a key and will contain the number of occurences of each codon. 
    #e.g. dict_codons['UGC']=0  
    dict_codons = {}
    
    for gene in cds_list: #cds_list was generated in part II
        
        gene_seq=Seq(gene) #Load the gene sequence as a biopython sequence
        
        mRNA=gene_seq.transcribe() #Transcribe the gene sequence to RNA
        
        codons = list(str(mRNA)[n:n+3] for n in range(0,len(mRNA),3)) #Obtain a list of codons in the RNA sequence 
        
        
        for codon in codons:
            if len(codon)==3: #Make sure the codon length is 3 to avaid any errors when calling the dictionary 
                if codon.count("A")+codon.count("U")+codon.count("G")+codon.count("C")==len(codon): #Make sure the codon doesn't contain undefined nucleotides such as 'N' for example
                    if codon in dict_codons:
                        dict_codons[codon] += 1 #Count the occurence of each codon
                    else:
                        dict_codons[codon] = 1
             

     

   
    #update the codon_usage dictionary based on the counts stored in the dict_codons dictionary  
    for codon in dict_codons:
        codon_usage[genetic_code[codon]][codon]+=dict_codons[codon]
     #This translates to:
     #  codon_usage[amino acid][codon]+=count of codon
        
    
    
    #This new dictionary will contain the percentage occurence of each codon coding for each amino acid:
    # e.g. codon_usage_percentage['C']= {'UGC': 40, 'UGU': 60}
    codon_usage_percentage={}
    
    
    
    for amino_acid in codon_usage: #Iterate through every amino acid
        total_amino_acid_usage=0 #Number of times an amino acid is encoded
        for codon in codon_usage[amino_acid]:
            total_amino_acid_usage+=codon_usage[amino_acid][codon] #Count the number of times an amino acid is encoded by each codon
        codon_usage_percentage[amino_acid]={} #Intitiate dictionary
        for codon in codon_usage[amino_acid]:
            try:
                codon_usage_percentage[amino_acid][codon]=round(codon_usage[amino_acid][codon]/total_amino_acid_usage*100,2) #Calculate the percentages of codon usage for each amino acid.
            except:
                pass
                
    
    #Here the user has the option to display the obtained dictionary, for aesthetic purposes, the following line is commented out but the dictionary was used for analyses featuring in the report. 
    #print(codon_usage_percentage)
 
    
 

         
 

     
     
"""




The script may take some time to run, here is the output of the script:



######################### Genlisea #########################

Length genomic DNA: 43356312
The genomic nitrogen was: 321101118 nitrogen atoms
The relative genomic nitrogen was: 7.4060985168664715 + or - 0.00014618231939446562 nitrogen atoms per basepair

Length of CDS: 17092254
The CDS nitrogen was: 128226799 nitrogen atoms
The relative CDS nitrogen was: 7.5020415095633375 + or - 0.000237036294680822 nitrogen atoms per basepair

Length of all protein: 5689403
The protein nitrogen was: 7786394 nitrogen atoms
The relative protein nitrogen was: 1.368578390386478 + or - 0.0006572932084183225 nitrogen atoms per basepair

Length of the non-coding genomic DNA: 26264058
The Non-Coding DNA nitrogen was: 192874319 nitrogen atoms
The relative Non-Coding DNA nitrogen was: 7.3436602599643965 nitrogen atoms per basepair

Length of all mRNA: 25415294
The total mRNA nitrogen was: 94718483 nitrogen atoms
The relative mRNA nitrogen was: 3.726830112608573 + or - 0.0005138304855715966 nitrogen atoms per basepair

Length of all introns: 8403199
The intron nitrogen was: 30904986 nitrogen atoms
The relative intron nitrogen was: 3.6777643847301484 + or - 0.0009232759112659433 nitrogen atoms per basepair

Length of all exons: 17012095
The exon nitrogen was: 63813497 nitrogen atoms
The relative exon nitrogen was: 3.7510663442685925 + or - 0.0006171562057590841 nitrogen atoms per basepair

Length of all tRNAs: 17747
The exon nitrogen was: 66691 nitrogen atoms
The relative exon nitrogen was: 3.757874570349918 + or - 0.0019021942776001596 nitrogen atoms per basepair

######################### Arabidopsis #########################

Length genomic DNA: 125938000
The genomic nitrogen was: 927282241 nitrogen atoms
The relative genomic nitrogen was: 7.363005931490099 + or - 8.398355529394707e-05 nitrogen atoms per basepair

Length of CDS: 33598839
The CDS nitrogen was: 250052111 nitrogen atoms
The relative CDS nitrogen was: 7.442284270596374 + or - 0.00016793575517931323 nitrogen atoms per basepair

Length of all protein: 10924993
The protein nitrogen was: 14893472 nitrogen atoms
The relative protein nitrogen was: 1.3632477384653703 + or - 0.00046012864754889284 nitrogen atoms per basepair

Length of the non-coding genomic DNA: 92339161
The Non-Coding DNA nitrogen was: 677230130 nitrogen atoms
The relative Non-Coding DNA nitrogen was: 7.3341594472577025 nitrogen atoms per basepair

Length of all mRNA: 108645137
The total mRNA nitrogen was: 402034069 nitrogen atoms
The relative mRNA nitrogen was: 3.7004331726324757 + or - 0.0002529736947185023 nitrogen atoms per basepair

Length of all introns: 24704888
The intron nitrogen was: 90592683 nitrogen atoms
The relative intron nitrogen was: 3.666994280646 + or - 0.0005420664866365375 nitrogen atoms per basepair

Length of all exons: 83940249
The exon nitrogen was: 311441386 nitrogen atoms
The relative exon nitrogen was: 3.710274745551446 + or - 0.0002858966685049005 nitrogen atoms per basepair

Length of all tRNAs: 58337
The exon nitrogen was: 219666 nitrogen atoms
The relative exon nitrogen was: 3.765466170697842 + or - 0.0010397512571922063 nitrogen atoms per basepair

######################### Trifolium #########################

Length genomic DNA: 304880418
The genomic nitrogen was: 2234200013 nitrogen atoms
The relative genomic nitrogen was: 7.328119095533384 + or - 5.270420051936954e-05 nitrogen atoms per basepair

Length of CDS: 48632115
The CDS nitrogen was: 360612689 nitrogen atoms
The relative CDS nitrogen was: 7.415114251148649 + or - 0.00013848599925682422 nitrogen atoms per basepair

Length of all protein: 16156565
The protein nitrogen was: 22009618 nitrogen atoms
The relative protein nitrogen was: 1.3622708787418614 + or - 0.00037060953283520704 nitrogen atoms per basepair

Length of the non-coding genomic DNA: 256248303
The Non-Coding DNA nitrogen was: 1873587324 nitrogen atoms
The relative Non-Coding DNA nitrogen was: 7.311608709463337 nitrogen atoms per basepair

Length of all mRNA: 93703227
The total mRNA nitrogen was: 345056310 nitrogen atoms
The relative mRNA nitrogen was: 3.6824378524338335 + or - 0.0002753800125201837 nitrogen atoms per basepair

Length of all introns: 45268711
The intron nitrogen was: 165448012 nitrogen atoms
The relative intron nitrogen was: 3.654798388228903 + or - 0.0004024233759553797 nitrogen atoms per basepair

Length of all exons: 48434516
The exon nitrogen was: 179608298 nitrogen atoms
The relative exon nitrogen was: 3.7082707299067468 + or - 0.0003771710019850335 nitrogen atoms per basepair

Length of all tRNAs: 20185
The exon nitrogen was: 76226 nitrogen atoms
The relative exon nitrogen was: 3.776368590537528 + or - 0.0017541559006970124 nitrogen atoms per basepair

"""

