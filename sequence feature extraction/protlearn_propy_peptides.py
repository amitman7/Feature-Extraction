import re
import protlearn
import numpy as np
import pandas as pd
from protlearn.features import length
from protlearn.features import aac
from protlearn.features import aaindex1
from protlearn.features import ngram
from protlearn.features import entropy
from protlearn.features import atc
from protlearn.features import binary
from protlearn.features import cksaap
from protlearn.features import ctd
from protlearn.features import ctdc
from protlearn.features import ctdt
from protlearn.features import ctdd
from protlearn.features import moreau_broto
from protlearn.features import moran
from protlearn.features import geary
from protlearn.features import paac
from protlearn.features import apaac
from propy import AAComposition
from propy import Autocorrelation
from propy import CTD
from propy import PyPro
import peptides
import matlab.engine
import csv
from protlearn.utils.validation import check_natural



def read_fasta(file_path):
    with open(file_path, 'r') as file:
        return file.read()

def extract_sequences(fasta_data):
    sequences = []
    current_sequence = ''
    current_id = None

    for line in fasta_data.split('\n'):
        line = line.strip() #removing \n  
        if line.startswith('>'): #sequence start with >
            if current_id is not None and current_sequence:
                sequences.append((current_id, current_sequence))
                current_sequence = ''
            current_id = extract_uniprotID(line)
        else:
            current_sequence += line

    if current_id is not None and current_sequence:
        sequences.append((current_id, current_sequence))
    
    return sequences

def extract_uniprotID(header_line):
    match = re.search(r'\|(\w+)\|', header_line)
    if match:
        return match.group(1)
    else:
        return None  

def min_len_sequence(sequences_list,sequences_num):
     min_len = len(sequences_list[0])
     for i in range (sequences_num):
         if (len(sequences_list[i]) < min_len):
             min_len = len(sequences_list[i])
     return min_len

def max_len_sequence(sequences_list,sequences_num):
     max_len = len(sequences_list[0])
     for i in range (sequences_num):
         if (len(sequences_list[i]) > max_len):
             max_len = len(sequences_list[i])
     return max_len 


##### protlearn #####   

def vector_to_dataframe_protlearn(uniprotID_list,vectors_list):

    columns =  [str(i) for i in range(1, 10298)]

    df = pd.DataFrame(vectors_list, columns=columns)
    df['uniprotID'] = uniprotID_list

    return df

def protlearn_to_vector(protein_sequence,min):
    vector = []

    #length - return 1 value
    vector.append(length(protein_sequence))

    #amino acid composition - return 20 values
    vector.append((aac(protein_sequence))[0][0])
    
    #aaIndex1 - return 553 values
    vector.append((aaindex1(protein_sequence))[0][0])

    #ngram - for dipeptide composition return 400 values
    vector.append((ngram(protein_sequence,n=2))[0][0])

    #ngram - for dipeptide composition return 8000 values
    vector.append((ngram(protein_sequence,n=3))[0][0])

    #entropy - return 1 value
    x = np.array(entropy(protein_sequence))
    vector.append(x)

    #atc - Array containing atomic compositions. 5 values
    vector.append(atc(protein_sequence)[0][0])

    #Array containing bond compositions. 3 values
    vector.append(atc(protein_sequence)[1][0])

    #cksaap - 400 values
    vector.append(cksaap(protein_sequence)[0][0])

    vector.append(ctd(protein_sequence)[0][0])
    
    #ctdc - 39 values
    vector.append(ctdc(protein_sequence)[0][0])

    #ctdt - 39 values
    vector.append(ctdt(protein_sequence)[0][0])

    #ctdd - 39 values
    vector.append(ctdd(protein_sequence)[0][0])

    #moreau_broto - 8 values
    vector.append(moreau_broto(protein_sequence)[0])

    #moran - 8 values
    vector.append(moran(protein_sequence)[0])

    #geary - 8 values
    vector.append(geary(protein_sequence)[0])

    #paac 
    vector.append(paac(protein_sequence, lambda_ = min -1)[0][0])

    #apaac 
    vector.append(apaac(protein_sequence,lambda_ = min -1)[0][0])


    
   

    return vector
 

##### propy #####   

def vector_to_dataframe_propy(uniprotID_list,vectors_list):

    columns =  [str(i) for i in range(1, 25204)]

    df = pd.DataFrame(vectors_list, columns=columns)
    df['uniprotID'] = uniprotID_list

    return df


def propy_to_vector(protein_sequence,max):
    vector = []

    #AAComposition - return 20 values
    res = AAComposition.CalculateAAComposition(protein_sequence)
    vector.append(list(res.values()))


    #DipeptideComposition - return 400 values
    res = AAComposition.CalculateDipeptideComposition(protein_sequence)
    vector.append(list(res.values()))

    #SpectrumDict - return 8000 values
    res = AAComposition.GetSpectrumDict(protein_sequence)
    vector.append(list(res.values()))

    #GearyAutoTotal - return 240 values
    res = Autocorrelation.CalculateGearyAutoTotal(protein_sequence)
    vector.append(list(res.values()))


    #NormalizedMoreauBrotoAutoTotal - return 240 values
    res = Autocorrelation.CalculateNormalizedMoreauBrotoAutoTotal(protein_sequence)
    vector.append(list(res.values()))


    #MoranAutoTotal - return 240 values
    res = Autocorrelation.CalculateMoranAutoTotal(protein_sequence)
    vector.append(list(res.values()))

    #CalculateC - return 21 values
    res = CTD.CalculateC(protein_sequence)
    vector.append(list(res.values()))

    #CalculateT - return 21 values
    res = CTD.CalculateT(protein_sequence)
    vector.append(list(res.values()))

    #CalculateD - return 105 values
    res = CTD.CalculateD(protein_sequence)
    vector.append(list(res.values()))


    #GetSOCN - return 7938 values. in the end there is padding
    res = PyPro.GetProDes(protein_sequence).GetSOCN(maxlag = max)
    vector.append(list(res.values()))

    #GetQSO - return 7978 values. in the end there is padding
    res = PyPro.GetProDes(protein_sequence).GetQSO(maxlag = max)
    vector.append(list(res.values()))

   
    return vector



##### peptides #####

def peptides_init(sequences_list,sequences_num,peptide_list):
    for i in range(sequences_num):
        peptide = peptides.Peptide(sequences_list[i])
        peptide_list.append(peptide)      

def vector_to_dataframe_peptides(uniprotID_list,vectors_list):

    columns =  [str(i) for i in range(1, 293)]

    df = pd.DataFrame(vectors_list, columns=columns)
    df['uniprotID'] = uniprotID_list

    return df

def peptides_to_vector(peptide):
    vector = []

    #aliphatic_index - 1 value
    result = peptide.aliphatic_index()
    vector.append([result])

    #auto_correlation - 45 values
    tables_types_len = len(peptides.tables.HYDROPHOBICITY.keys())
    result = []
    for j in range(tables_types_len):
        table =  peptides.tables.HYDROPHOBICITY[list(peptides.tables.HYDROPHOBICITY.keys())[j]]
        result.append(peptide.auto_correlation(table=table))
    vector.append(result)

    #auto_covariance - 45 values
    tables_types_len = len(peptides.tables.HYDROPHOBICITY.keys())
    result = []
    for j in range(tables_types_len):
        table =  peptides.tables.HYDROPHOBICITY[list(peptides.tables.HYDROPHOBICITY.keys())[j]]
        result.append(peptide.auto_covariance(table=table))
    vector.append(result)
    

    #blosum_indices - 10 value
    result = list(peptide.blosum_indices())
    vector.append(result)

    #boman - 1 value
    result = peptide.boman()
    vector.append([result])

    #charge - 9 values
    tables_types_len = len(peptides.tables.PK.keys())
    result = []
    for j in range(tables_types_len):
        result.append(peptide.charge(pKscale = list(peptides.tables.PK.keys())[j]))
    vector.append(result)

    
    #counts - 26 value
    result = list(peptide.counts().values())
    vector.append(result)

    #cruciani_properties - 3 value
    result = list(peptide.cruciani_properties())
    vector.append(result)

    #descriptors - 88 value
    result = list(peptide.descriptors().values())
    vector.append(result)


    #hydrophobic_moment - 2 value
    result =[]
    result.append(peptide.hydrophobic_moment(angle =100))
    result.append(peptide.hydrophobic_moment(angle =160))
    vector.append(result)


    #hydrophobicity - 45 values
    scale_len = len(peptides.tables.HYDROPHOBICITY.keys())
    result = []
    for j in range(scale_len):
        result.append(peptide.hydrophobicity(scale = list(peptides.tables.HYDROPHOBICITY.keys())[j]))
    vector.append(result)


    #instability_index - 1 value
    result = peptide.instability_index()
    vector.append([result])


    #isoelectric_point - 9 values
    pKscale_len = len(peptides.tables.PK.keys())
    result = []
    for j in range(pKscale_len):
        result.append(peptide.isoelectric_point(pKscale = list(peptides.tables.PK.keys())[j]))
    vector.append(result)

    #mass_shift - 3 values
    aa_shift_len = len(peptides.tables.MASS_SHIFT.keys())
    result = []
    for j in range(aa_shift_len):
       result.append(peptide.mass_shift(aa_shift = list(peptides.tables.MASS_SHIFT.keys())[j]))
    vector.append(result)


    #molecular_weight - 3 values
    average_len = len(peptides.tables.MOLECULAR_WEIGHT.keys())
    result = []
    for j in range(average_len):
       result.append(peptide.molecular_weight(average = list(peptides.tables.MOLECULAR_WEIGHT.keys())[j]))
    vector.append(result)


    #mz - 1 value
    result = peptide.mz()
    vector.append([result])


    return vector


def main():
    file_path = r'your path for the fasta files'
    fasta_data = read_fasta(file_path)
    sequences_and_id = extract_sequences(fasta_data)
    sequences_num = len(sequences_and_id) 
    sequences_list = [i[1] for i in sequences_and_id]
    uniprotID_list = [i[0] for i in sequences_and_id]
    min = min_len_sequence(sequences_list,sequences_num)
    max = max_len_sequence(sequences_list,sequences_num)


    
    #### protlearn vector ####
    vectors_list = []
    for i in range (sequences_num):
        x = protlearn_to_vector(sequences_list[i],min)
        flattened_array = np.concatenate([arr.flatten() for arr in x])
        flattened_list = flattened_array.tolist()
        vectors_list.append(flattened_list)

    df_protlearn = vector_to_dataframe_protlearn(uniprotID_list,vectors_list)
    df_protlearn.to_csv('protlearn.csv', index=False)

    
   

#### propy vector ####
    vectors_list = []
    for i in range (sequences_num):
        x = propy_to_vector(sequences_list[i],max)
        one_dim_list = [item for sublist in x for item in sublist]
        vectors_list.append(one_dim_list)
       
    df_propy = vector_to_dataframe_propy(uniprotID_list,vectors_list)
    df_propy.to_csv('propy.csv', index=False)
    



###### peptides vector ##### 
    peptide_list = []
    peptides_init(sequences_list,sequences_num,peptide_list)


    vectors_list = []
    for i in range (sequences_num):
        x = peptides_to_vector(peptide_list[i])
        one_dim_list = [item for sublist in x for item in sublist]
        vectors_list.append(one_dim_list)
        
    df_peptides = vector_to_dataframe_peptides(uniprotID_list,vectors_list)
    df_peptides.to_csv('peptides.csv', index=False)




    #combine all the  csv files to one file
    """
    df1 = pd.read_csv("protlearn_without_C5J4T0.csv")
    df1 = df1.iloc[:, :-1] 

    df2 = pd.read_csv("propy_without_C5J4T0.csv")
    df2 = df2.iloc[:, :-1] 
    

    df3 = pd.read_csv("peptides_without_C5J4T0.csv")
    df2 = df2.iloc[:, :-1] 

    combined_df = pd.concat([df1,df2,df3], axis=1,ignore_index=True)
    combined_df.to_csv("combined_file.csv", index=False)
    """

  

    
   

    



   
    
    
    
    
    

    

    return None




if __name__ == "__main__":
    main()
