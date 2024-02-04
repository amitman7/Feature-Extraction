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


#create the data frame
def feature_to_dataframe(sequences,sequences_num):
    data = []
    for i in range(sequences_num):
        data.append(sequences[i])
    df = pd.DataFrame(data, columns=['Uniprot_ID','sequence'])
    return df


def main():
    file_path = # write your fasta file path
    fasta_data = read_fasta(file_path)
    sequences = extract_sequences(fasta_data)
    sequences_dict= dict(sequences)
    sequences_num = len(sequences) 
    sequences_list = [i[1] for i in sequences]
    my_dataframe = feature_to_dataframe(sequences,sequences_num)


    ### protlearn ###
    length_feature_add(my_dataframe,sequences_list)
    aac_feature_add(my_dataframe,sequences_num)
    aaindex1_feature_add(my_dataframe,sequences_num)
    ngram_feature_add(my_dataframe,sequences_num)
    entropy_feature_add(my_dataframe,sequences_list)
    atc_feature_add(my_dataframe,sequences_num)
    binary_feature_add(my_dataframe,sequences_num,sequences_list)
    cksaap_feature_add(my_dataframe,sequences_num)
    ctd_feature_add(my_dataframe,sequences_num)
    ctdc_feature_add(my_dataframe,sequences_num)
    ctdt_feature_add(my_dataframe,sequences_num)
    ctdd_feature_add(my_dataframe,sequences_num)
    moreau_broto_feature_add(my_dataframe,sequences_num)
    moran_feature_add(my_dataframe,sequences_num)
    geary_feature_add(my_dataframe,sequences_num)
    paac_feature_add(my_dataframe,sequences_num,sequences_list)
    apaac_feature_add(my_dataframe,sequences_num,sequences_list)

    ### propy ###
    AAComposition_feature_add(my_dataframe,sequences_num)
    DipeptideComposition_feature_add(my_dataframe,sequences_num)
    SpectrumDict_feature_add(my_dataframe,sequences_num)
    GearyAutoTotal_feature_add(my_dataframe,sequences_num)
    Moreau_BrotoAutoTotal_feature_add(my_dataframe,sequences_num)
    MoranAutoTotal_feature_add(my_dataframe,sequences_num)
    CalculateC_feature_add(my_dataframe,sequences_num)
    CalculateT_feature_add(my_dataframe,sequences_num)
    CalculateD_feature_add(my_dataframe,sequences_num)
    GetSOCN_feature_add(my_dataframe,sequences_num,sequences_list)
    GetQSO_feature_add(my_dataframe,sequences_num,sequences_list)


    ### peptides ###
    peptide_list = []
    peptides_init(my_dataframe,sequences_num,peptide_list)

    aliphatic_index_feature_add(my_dataframe,sequences_num,peptide_list)
    auto_correlation_feature_add(my_dataframe,sequences_num,peptide_list)
    auto_covariance_feature_add(my_dataframe,sequences_num,peptide_list)
    blosum_indices_feature_add(my_dataframe,sequences_num,peptide_list)
    boman_feature_add(my_dataframe,sequences_num,peptide_list)
    charge_feature_add(my_dataframe,sequences_num,peptide_list)
    counts_feature_add(my_dataframe,sequences_num,peptide_list)
    cruciani_properties_feature_add(my_dataframe,sequences_num,peptide_list)
    descriptors_feature_add(my_dataframe,sequences_num,peptide_list)
    hydrophobic_moment_feature_add(my_dataframe,sequences_num,peptide_list)
    hydrophobicity_feature_add(my_dataframe,sequences_num,peptide_list)
    instability_index_feature_add(my_dataframe,sequences_num,peptide_list)
    isoelectric_point_feature_add(my_dataframe,sequences_num,peptide_list)
    mass_shift_point_feature_add(my_dataframe,sequences_num,peptide_list)
    molecular_weight_feature_add(my_dataframe,sequences_num,peptide_list)
    mz_feature_add(my_dataframe,sequences_num,peptide_list)


    ### FEGS  ###
    fegs_path = r'C:\Users\97252\Desktop\עמית\אוניברסיטה\שנה ג\התמחות\FEGS'  
    FEGS_add(my_dataframe,sequences_num,fasta_data,fegs_path)


 
    #dictionary with info for each feature
    feature_info = {}
    feature_info_add(my_dataframe,feature_info)


    return None
   

""" 
Here is the implementation for the functions above, organized in the order of APIs.
"""


### protlearn  ###

def length_feature_add(my_dataframe,sequences_list):
    my_dataframe['length'] = length(sequences_list)

def aac_feature_add(my_dataframe,sequences_num):
    my_dataframe['aac'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = aac(my_dataframe['sequence'][i])
        list_toAdd.append(result[0]) 
    my_dataframe['aac'] = list_toAdd 

def aaindex1_feature_add(my_dataframe,sequences_num):
    my_dataframe['aaindex1'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = aaindex1(my_dataframe['sequence'][i],standardize='zscore')
        list_toAdd.append(result[0][0]) 
    my_dataframe['aaindex1'] = list_toAdd  

def ngram_feature_add(my_dataframe,sequences_num):
    my_dataframe['ngram'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = ngram(my_dataframe['sequence'][i])
        list_toAdd.append(result[0][0]) 
    my_dataframe['ngram'] = list_toAdd  

def entropy_feature_add(my_dataframe,sequences_list):
    my_dataframe['entropy'] = entropy(sequences_list)

def atc_feature_add(my_dataframe,sequences_num):
    my_dataframe['atc'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = atc(my_dataframe['sequence'][i])
        list_toAdd.append(result) 
    my_dataframe['atc'] = list_toAdd 

def binary_feature_add(my_dataframe,sequences_num,sequences_list):
   binary_data = binary(sequences_list)
   my_dataframe['binary'] = None
   list_toAdd = []
   for i in range(sequences_num): 
        list_toAdd.append(binary_data[i])
   my_dataframe['binary'] = list_toAdd    

def cksaap_feature_add(my_dataframe,sequences_num):
    my_dataframe['cksaap'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = cksaap(my_dataframe['sequence'][i])
        list_toAdd.append(result[0][0]) 
    my_dataframe['cksaap'] = list_toAdd 

def ctd_feature_add(my_dataframe,sequences_num):
    my_dataframe['ctd'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = ctd(my_dataframe['sequence'][i])
        list_toAdd.append(result[0]) 
    my_dataframe['ctd'] = list_toAdd 

def ctdc_feature_add(my_dataframe,sequences_num):
    my_dataframe['ctdc'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = ctdc(my_dataframe['sequence'][i])
        list_toAdd.append(result[0][0]) 
    my_dataframe['ctdc'] = list_toAdd     

def ctdt_feature_add(my_dataframe,sequences_num):
    my_dataframe['ctdt'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = ctdt(my_dataframe['sequence'][i])
        list_toAdd.append(result[0][0]) 
    my_dataframe['ctdt'] = list_toAdd     

def ctdd_feature_add(my_dataframe,sequences_num):
    my_dataframe['ctdd'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = ctdd(my_dataframe['sequence'][i])
        list_toAdd.append(result[0][0]) 
    my_dataframe['ctdd'] = list_toAdd     

def moreau_broto_feature_add(my_dataframe,sequences_num):
    my_dataframe['moreau_broto'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = moreau_broto(my_dataframe['sequence'][i])
        list_toAdd.append(result[0]) 
    my_dataframe['moreau_broto'] = list_toAdd 

def moran_feature_add(my_dataframe,sequences_num):
    my_dataframe['moran'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = moran(my_dataframe['sequence'][i])
        list_toAdd.append(result[0]) 
    my_dataframe['moran'] = list_toAdd         

def geary_feature_add(my_dataframe,sequences_num):
    my_dataframe['geary'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = geary(my_dataframe['sequence'][i])
        list_toAdd.append(result[0]) 
    my_dataframe['geary'] = list_toAdd

def paac_feature_add(my_dataframe,sequences_num,sequences_list):
    my_dataframe['paac'] = None
    list_toAdd = []
    lanbda = min_len_sequence(sequences_list,sequences_num)
    for i in range(sequences_num):
        result = paac(my_dataframe['sequence'][i],lambda_= lanbda -1)
        list_toAdd.append(result[0][0]) 
    my_dataframe['paac'] = list_toAdd         

def apaac_feature_add(my_dataframe,sequences_num,sequences_list):
    my_dataframe['apaac'] = None
    list_toAdd = []
    lanbda = min_len_sequence(sequences_list,sequences_num)
    for i in range(sequences_num):
        result = apaac(my_dataframe['sequence'][i],lambda_= lanbda -1)
        list_toAdd.append(result[0][0]) 
    my_dataframe['apaac'] = list_toAdd 



### propy ###

def AAComposition_feature_add(my_dataframe,sequences_num):
    my_dataframe['AAComposition'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = AAComposition.CalculateAAComposition(my_dataframe['sequence'][i])
        list_result = list(result.values())
        list_toAdd.append(list_result) 
    my_dataframe['AAComposition'] = list_toAdd     

def DipeptideComposition_feature_add(my_dataframe,sequences_num):
    my_dataframe['DipeptideComposition'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = AAComposition.CalculateDipeptideComposition(my_dataframe['sequence'][i])
        list_result = list(result.values())
        list_toAdd.append(list_result) 
    my_dataframe['DipeptideComposition'] = list_toAdd   
   
def SpectrumDict_feature_add(my_dataframe,sequences_num):
    my_dataframe['SpectrumDict'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = AAComposition.GetSpectrumDict(my_dataframe['sequence'][i])
        list_result = list(result.values())
        list_toAdd.append(list_result) 
    my_dataframe['SpectrumDict'] = list_toAdd

def GearyAutoTotal_feature_add(my_dataframe,sequences_num):
    my_dataframe['GearyAutoTotal'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = Autocorrelation.CalculateGearyAutoTotal(my_dataframe['sequence'][i])
        list_result = list(result.values())
        list_toAdd.append(list_result) 
        #list_toAdd.append(result) 
    my_dataframe['GearyAutoTotal'] = list_toAdd

def Moreau_BrotoAutoTotal_feature_add(my_dataframe,sequences_num):
    my_dataframe['Moreau_BrotoAutoTotal'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = Autocorrelation.CalculateNormalizedMoreauBrotoAutoTotal(my_dataframe['sequence'][i])
        list_result = list(result.values())
        list_toAdd.append(list_result) 
        #list_toAdd.append(result) 
    my_dataframe['Moreau_BrotoAutoTotal'] = list_toAdd

def MoranAutoTotal_feature_add(my_dataframe,sequences_num):
    my_dataframe['MoranAutoTotal'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = Autocorrelation.CalculateMoranAutoTotal(my_dataframe['sequence'][i])
        list_result = list(result.values())
        list_toAdd.append(list_result) 
        #list_toAdd.append(result) 
    my_dataframe['MoranAutoTotal'] = list_toAdd

def CalculateC_feature_add(my_dataframe,sequences_num):
    my_dataframe['CalculateC'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = CTD.CalculateC(my_dataframe['sequence'][i])
        list_result = list(result.values())
        list_toAdd.append(list_result) 
        #list_toAdd.append(result) 
    my_dataframe['CalculateC'] = list_toAdd

def CalculateT_feature_add(my_dataframe,sequences_num):
    my_dataframe['CalculateT'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = CTD.CalculateT(my_dataframe['sequence'][i])
        list_result = list(result.values())
        list_toAdd.append(list_result) 
        #list_toAdd.append(result) 
    my_dataframe['CalculateT'] = list_toAdd

def CalculateD_feature_add(my_dataframe,sequences_num):
    my_dataframe['CalculateD'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = CTD.CalculateD(my_dataframe['sequence'][i])
        list_result = list(result.values())
        list_toAdd.append(list_result) 
        #list_toAdd.append(result) 
    my_dataframe['CalculateD'] = list_toAdd

def GetSOCN_feature_add(my_dataframe,sequences_num,sequences_list):
    my_dataframe['GetSOCN'] = None
    list_toAdd = []
    lanbda = max_len_sequence(sequences_list,sequences_num)
    for i in range(sequences_num):
        seq = my_dataframe['sequence'][i]
        result = PyPro.GetProDes(seq).GetSOCN(maxlag = lanbda)
        list_result = list(result.values())
        list_toAdd.append(list_result) 
        #list_toAdd.append(result) 
    my_dataframe['GetSOCN'] = list_toAdd

def GetQSO_feature_add(my_dataframe,sequences_num,sequences_list):
    my_dataframe['GetQSO'] = None
    list_toAdd = []
    lanbda = max_len_sequence(sequences_list,sequences_num)
    for i in range(sequences_num):
        seq = my_dataframe['sequence'][i]
        result = PyPro.GetProDes(seq).GetQSO(maxlag = lanbda)
        #list_result = list(result.values())
        #list_toAdd.append(list_result) 
        list_toAdd.append(result) 
    my_dataframe['GetQSO'] = list_toAdd



###### peptides #####

#init for the sequnces in to type peptide
def peptides_init(my_dataframe,sequences_num,peptide_list):
    for i in range(sequences_num):
        peptide = peptides.Peptide(my_dataframe['sequence'][i])
        peptide_list.append(peptide) 


def aliphatic_index_feature_add(my_dataframe,sequences_num,peptide_list):
    my_dataframe['aliphatic_index'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = peptide_list[i].aliphatic_index()
        list_toAdd.append(result)
    my_dataframe['aliphatic_index'] = list_toAdd

def auto_correlation_feature_add(my_dataframe,sequences_num,peptide_list):
    peptides.tables.HYDROPHOBICITY.keys()
    my_dataframe['auto_correlation'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = []
        tables_types_len = len(peptides.tables.HYDROPHOBICITY.keys())
        for j in range(tables_types_len):
            table =  peptides.tables.HYDROPHOBICITY[list(peptides.tables.HYDROPHOBICITY.keys())[j]]
            result.append(peptide_list[i].auto_correlation(table=table))
        list_toAdd.append(result)
    my_dataframe['auto_correlation'] = list_toAdd

def auto_covariance_feature_add(my_dataframe,sequences_num,peptide_list):
    peptides.tables.HYDROPHOBICITY.keys()
    my_dataframe['auto_covariance'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = []
        tables_types_len = len(peptides.tables.HYDROPHOBICITY.keys())
        for j in range(tables_types_len):
            table =  peptides.tables.HYDROPHOBICITY[list(peptides.tables.HYDROPHOBICITY.keys())[j]]
            result.append(peptide_list[i].auto_covariance(table=table))
        list_toAdd.append(result)
    my_dataframe['auto_covariance'] = list_toAdd

def blosum_indices_feature_add(my_dataframe,sequences_num,peptide_list):
    my_dataframe['blosum_indices'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = list(peptide_list[i].blosum_indices())
        list_toAdd.append(result)
    my_dataframe['blosum_indices'] = list_toAdd

def boman_feature_add(my_dataframe,sequences_num,peptide_list):
    my_dataframe['boman'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = peptide_list[i].boman()
        list_toAdd.append(result)
    my_dataframe['boman'] = list_toAdd

def charge_feature_add(my_dataframe,sequences_num,peptide_list):
    my_dataframe['charge'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = []
        pKscale_len = len(peptides.tables.PK.keys())
        for j in range(pKscale_len):
            result.append(peptide_list[i].charge(pKscale = list(peptides.tables.PK.keys())[j]))
        list_toAdd.append(result)
    my_dataframe['charge'] = list_toAdd
    
def counts_feature_add(my_dataframe,sequences_num,peptide_list):
    my_dataframe['counts'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = list(peptide_list[i].counts().values())
        list_toAdd.append(result)
    my_dataframe['counts'] = list_toAdd

def cruciani_properties_feature_add(my_dataframe,sequences_num,peptide_list):
    my_dataframe['cruciani_properties'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result =list(peptide_list[i].cruciani_properties()) 
        list_toAdd.append(result)
    my_dataframe['cruciani_properties'] = list_toAdd

def descriptors_feature_add(my_dataframe,sequences_num,peptide_list):
    my_dataframe['descriptors'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = list(peptide_list[i].descriptors().values())
        list_toAdd.append(result)
    my_dataframe['descriptors'] = list_toAdd

def hydrophobic_moment_feature_add(my_dataframe,sequences_num,peptide_list):
    my_dataframe['hydrophobic_moment'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = []
        result.append(peptide_list[i].hydrophobic_moment(angle =100))
        result.append(peptide_list[i].hydrophobic_moment(angle =160))
        list_toAdd.append(result)
    my_dataframe['hydrophobic_moment'] = list_toAdd

def hydrophobicity_feature_add(my_dataframe,sequences_num,peptide_list):
    my_dataframe['hydrophobicity'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = []
        scale_len = len(peptides.tables.HYDROPHOBICITY.keys())
        for j in range(scale_len):
            result.append(peptide_list[i].hydrophobicity(scale = list(peptides.tables.HYDROPHOBICITY.keys())[j]))
        list_toAdd.append(result)
    my_dataframe['hydrophobicity'] = list_toAdd
    
def instability_index_feature_add(my_dataframe,sequences_num,peptide_list):
    my_dataframe['instability_index'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = peptide_list[i].instability_index()
        list_toAdd.append(result)
    my_dataframe['instability_index'] = list_toAdd

def isoelectric_point_feature_add(my_dataframe,sequences_num,peptide_list):
    my_dataframe['isoelectric_point'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = []
        pKscale_len = len(peptides.tables.PK.keys())
        for j in range(pKscale_len):
            result.append(peptide_list[i].isoelectric_point(pKscale = list(peptides.tables.PK.keys())[j]))
        list_toAdd.append(result)
    my_dataframe['isoelectric_point'] = list_toAdd

def mass_shift_point_feature_add(my_dataframe,sequences_num,peptide_list):
    my_dataframe['mass_shift'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = []
        aa_shift_len = len(peptides.tables.MASS_SHIFT.keys())
        for j in range(aa_shift_len):
            result.append(peptide_list[i].mass_shift(aa_shift = list(peptides.tables.MASS_SHIFT.keys())[j]))
        list_toAdd.append(result)
    my_dataframe['mass_shift'] = list_toAdd

def molecular_weight_feature_add(my_dataframe,sequences_num,peptide_list):
    my_dataframe['molecular_weight'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = []
        average_len = len(peptides.tables.MOLECULAR_WEIGHT.keys())
        for j in range(average_len):
            result.append(peptide_list[i].molecular_weight(average = list(peptides.tables.MOLECULAR_WEIGHT.keys())[j]))
        list_toAdd.append(result)
    my_dataframe['molecular_weight'] = list_toAdd

def mz_feature_add(my_dataframe,sequences_num,peptide_list):
    my_dataframe['mz'] = None
    list_toAdd = []
    for i in range(sequences_num):
        result = peptide_list[i].mz()
        list_toAdd.append(result)
    my_dataframe['mz'] = list_toAdd


##### FEGS ########

def FEGS_add(my_dataframe,sequences_num,fasta_data,fegs_path):
    # Start MATLAB engine
    eng = matlab.engine.start_matlab()
    
    # Set path to FEGS
    eng.addpath(fegs_path)
    
    # Execute FEGS
    result = eng.FEGS(fasta_data)

    # Stop MATLAB engine
    eng.quit()

    result_array = list(result)
    my_dataframe['FEGS'] = None
    list_toAdd = []
    for i in range(sequences_num):
        list_toAdd.append(result_array[i])   
    my_dataframe['FEGS'] = list_toAdd 





### info for the features ###
"len(my_dataframe.columns)-2 is the number of features"
def feature_info_add(my_dataframe,feature_info):
    info = []
    
    #length
    info.append("""name of the feature: length
feature role: calculate the length of each sequence
returns: int - value of the length
================================================================================= 
""")
    #aac
    info.append("""name of the feature: aac
feature role: calculate the frequencies of amino acids for each sequence in the dataset.
returns: a list of size 2.
[0] - the matrix of the frequencies.
[1] - the order of the amino acids in the matrix - ACDEFGHIKLMNPQRSTVWY  

the dataFrame contains just the values, to get the order - print(result[1])
=================================================================================""")
    #aaindex1
    info.append("""name of the feature: aaindex1
feature role: aaindex1_feature: AAIndex1-based physicochemical properties.
AAindex1 ver.9.2 (release Feb, 2017) is a set of 20 numerical values representing various physicochemical and biological properties of amino acids. Currently, it contains 566 indices, of which 553 contain no NaNs. The indices will be collected for each amino acid in the sequence, then averaged across the sequence.
returns: a list of size 2.
[0] - matrix with 553 cell with the values calculated. 
[1] - the order of the features which the values calculated by.  

the dataFrame contains just the values, to get the order - print(result[1])
=================================================================================""")
    #ngram
    info.append("""name of the feature:  ngram (N-gram composition)
feature role: This function computes the di- or tripeptide composition of amino acid sequences. Therefore, the function parameter n can only take on the arguments 2 and 3 - otherwise, it will raise a ValueError.
returns: a list of size 2. 
[0] - Depending on n, the returned array will be of size:
- (n_samples, 400) for dipeptide composition
- (n_samples, 8000) for tripeptide composition . 
[1] - the order of the features which the values calculated by.  

the dataFrame contains just the values, to get the order - print(result[1])
=================================================================================""")
    #entropy
    info.append("""
name of the feature: entropy
feature role: This function computes the Shannon entropy for each sequence in the dataset 
returns: int
=================================================================================""")
    #atc
    info.append("""name of the feature: Atomic and bond composition - atc
feature role: This function returns the sum of atomic and bond compositions for each amino acid sequence. The atomic features are comprised of five atoms (C, H, N, O, and S), and the bond features are comprised of total bonds, single bonds, and double bonds.
atc is calculated as follows:
atoms(i)=niN
where i denotes the type of atoms, ni is the total number of atoms of type i, and N is the total number of atoms in the sequence.
bonds(j)=nj
where j denotes the type of bond and nj is the total number of bonds of type j.
returns: tuple of two arrays - 
1. Array containing atomic compositions len = 5
2. Array containing bond compositions len = 3(how many first, how many seconds..)
=================================================================================""")
    #binary
    info.append("""name of the feature: binary
feature role: create a binary descrption for each amino acid in the sequence and merge them by their order to one vector 
returns : list containing binary vectors which padd with zeros to the biggest sequence in the data set 
""")
    #cksaap
    info.append("""name of the feature: cksaap
feature role: This function returns the k-spaced amino acid pair composition of each sequence in the dataset.
Since there are 20 natural amino acids, there are 400 possible amino acid pairs. 
The parameter 'k' represents the gap between the amino acid pair. An example for k=1 would be AxY, where 'x' can be any amino acid. Similary, an example for k=2 would be AxxY. If k=0, the function returns the dipeptide composition of each sequence.
returns: list of size 2
[0] - the values of thr pairs.
[1] - to the order of the pairs

to extract the order of the pairs in the matrix print(result[1])
=================================================================================""")
    #ctd
    info.append("""name of the feature: ctd
feature role: These descriptors were initially developed to model protein-protein interactions. Amino acids can be grouped into 7 different classes based on their dipoles and side chain volumes, which reflect their electrostatic and hydrophobic interactions.
1 Dipole Scale (Debye): -, Dipole < 1.0; +, 1.0 < Dipole < 2.0; ++, 2.0 < Dipole < 3.0; +++, Dipole > 3.0; +'+'+', Dipole > 3.0 with opposite orientation.
2 Volume Scale (A∘3): -, Volume < 50; +, Volume > 50.                
3 Cys is separated from class 3 because of its ability to form disulfide bonds.
                
After grouping, these class triads are computed and normalized.
returns: list of size 3               

to extract the order of the 3 digits in the matrix print(result[1])
=================================================================================""")
    #ctdc
    info.append("""name of the feature: ctdc. Composition/Transition/Distribution - Composition
feature role: Amino acids are categorized into 3 groups based on their physicochemical properties. The properties used here include hydrophobicity, normalized van der Waals volume, polarity, polarizability, charge, secondary structure, and solvent accessibility. For hydrophobicity, we use seven different groupings based on different studies, which can all be found in AAIndex1.
After grouping, the frequency of each class will be calculated for each physicochemical property per sequence in the dataset as follows:
c(i)=niNi=1,2,3
where ni denotes the frequency of group i in the sequence, and N is the total sequence length.
For instance, the sequence "ARKLY" translates to "23311" with respect to polarity groups. Thus, for polarity, the outcome will be P1 = 2/5, P2 = 1/5, and P3 = 2/5.
As there are 13 different properties and 3 groups for each, the total dimension of this descriptor is 39. 
returns: list of size 2.
[0] - the values .
[1] - the order of properties for example: 'Polarity-G1' = it means the descriptor value corresponding to the frequency of amino acids in the first group based on their polarity in the given sequence dataset.

to extract the order of the properties in the matrix print(result[1])
=================================================================================""")
    #ctdt
    info.append("""name of the feature: ctdt. Composition/Transition/Distribution - Transition
feature role: Amino acids are categorized into 3 groups based on their physicochemical properties. The properties used here include hydrophobicity, normalized van der Waals volume, polarity, polarizability, charge, secondary structure, and solvent accessibility. For hydrophobicity, we use seven different groupings based on different studies, which can all be found in AAIndex1.
This descriptor computes the frequency of transitions between groups as follows:
t(ij)=nij+njiN-1ij=12,13,23
where nij and nji denote the frequency of transitions from group i to j and j to i, respectively, in the sequence, and N is the total sequence length.
For instance, if the encoded sequence is "32132223311311222222", then there are 2 transitions from groups 1 to 2 and 2 to 1. Therefore, the descriptor for this particular group transition will be 2/19. Similarly, there are 3 transitions from groups 2 to 3 and 3 to 2, so the descriptor for this transition is 3/19.
As with CTDC, the dimensionality of this feature is 39.
return value: list of size 2.
[0] - the values . 
[1] - the order of properties for example: 'Hydrophobicity_ARGP820101-T1221' = it means the descriptor value corresponding to the frequency of transitions between hydrophobicity group 1 and group 2 based on the ARG820101 hydrophobicity scale in

to extract the order of the properties in the matrix print(result[1])
=================================================================================""")
    #ctdd
    info.append("""name of the feature: ctdd. Composition/Transition/Distribution - Distribution
feature role: Amino acids are categorized into 3 groups based on their physicochemical properties. The properties used here include hydrophobicity, normalized van der Waals volume, polarity, polarizability, charge, secondary structure, and solvent accessibility. For hydrophobicity, we use seven different groupings based on different studies, which can all be found in AAIndex1.
There are five distribution descriptors for each physicochemical property and they are the position percents in the whole sequence for the first residue, 25 residues, 50 residues, 75 residues, and 100 residues for a certain encoded class. For instance, if the encoded sequence is '32132223311311222222', then there are 10 residues encoded as 2. The positions for the first residue 2, the 2nd residue 2 (25% * 10 = 2), the 5th 2 residue (50% * 10 = 5), the 7th 2 (75% * 10 = 7) and the 10th residue 2 (100% * 10) in the encoded sequence are 2, 5, 15, 17, 20, so that the distribution descriptors for 2 are: 10.0 (2/20 * 100), 25.0 (5/20 * 100), 75.0 (15/20 * 100), 85.0 (17/20 * 100), 100.0 (20/20 * 100).
Since there are 13 physicochemical properties, 3 groups (see ctdc), and 5 distribution descriptors, the total dimensionality of this feature is 195.
return value: list of size 2.
[0] - the values . 
[1] - the order of properties (as in ctdc and ctdt)

to extract the order of the properties in the matrix print(result[1])
=================================================================================""")
    #moreau_broto
    info.append("""name of the feature: Moreau-Broto
role of the feature: Normalized Moreau-Broto autocorrelation based on AAIndex1.
Moreau-Broto autocorrelation descriptors are defined based on the distribution of AAIndex1-based amino acid properties along the sequence. All indices are standardized before computing the descriptors:
returns: list with float values
=================================================================================""")
    #moran
    info.append("""name of the feature: moran
role of the feature: Moran's I based on AAIndex1.
Moran's I autocorrelation descriptors are defined based on the distribution of AAIndex1-based amino acid properties along the sequence. All indices are standardized before computing the descriptors.
returns: list with float values
=================================================================================""")
    #geary
    info.append("""name of the feature: geary
role of the feature:Geary's C based on AAIndex1.
Geary's C autocorrelation descriptors are defined based on the distribution of AAIndex1-based amino acid properties along the sequence. All indices are standardized before computing the descriptors
is the same as in Moran's I.
returns: list with float values
=================================================================================""")
    #paac
    info.append("""name of the feature: paac
role of the feature: Pseudo amino acid composition.
Similar to the vanilla amino acid composition, this feature characterizes the protein mainly using a matrix of amino-acid frequencies, which helps with dealing with proteins without significant sequence homology to other proteins. However, additional information are also included in the matrix to represent some local features, such as correlation between residues of a certain distance.
First, the original hydrophobicity, hydrophilicity, and side chain mass values are converted 
returns: return value: 2-D list.
[0] - the values.
[1] - the order of properties.
                
to extract the order of the properties in the matrix print(result[1])
=================================================================================""")
    #apaac
    info.append("""name of the feature: apaac
role of the feature: Amphiphilic pseudo amino acid composition.
This feature has the same form as the vanilla amino acid composition, but contains much more information that is related to the sequence order of a protein and the distribution of the hydrophobic and hydrophilic amino acids along its chain.
Using H1 (i) and H2 (i) as defined above in the pseudo amino acid composition function, the correlation functions for hydrophobicity and hydrophilicity 
returns: list of size 2.
[0] - the values . [1] - the order of properties 
                
the data frame contains just the values             
to extract the order of the properties in the matrix print(result[1]) 
=================================================================================""")
    #AAComposition
    info.append("""name of the feature: AAComposition
role of the feature: return the percentage of each amino acid in the sequnce. 
the function return dicaonry as the keys are the amino acids.
return value: contains the composition of 20 amino acids.
return type: Dict[str, float]
                
the data frame contains just the percentage
the order of it is ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
==============================================================================""")
    #DipeptideComposition
    info.append("""name of the feature: DipeptideComposition
role of the feature: Calculate the composition of dipeptidefor a given protein sequence.
Parameters:	ProteinSequence (a pure protein sequence) -
return value: result - contains the composition of 400 dipeptides
Return type: Dict[str, float]

the data frame contains just the percentage
the order of it is ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'] * ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
=================================================================================""")
    #SpectrumDict
    info.append("""name of the feature: SpectrumDict
role of the feature: Calcualte the spectrum descriptors of 3-mers for a given protein.
Parameters:	proteinsequence (a pure protein sequence) -
Returns : Dict[str, int] contains the composition values of 8000 3-mers 
                
the data frame contains just the values
the order of it is ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'] times 3
to get the order = print(AAComposition.Getkmers())
=================================================================================""")
    #GearyAutoTotal
    info.append("""name of the feature: GearyAutoTotal
role of the feature: Compute all autocorrelation descriptors based on 8 properties of AADs.
returns: Dict[str, int] contains 30*8=240 Geary
autocorrelation descriptors based on the given properties(i.e.,
_AAPropert).
                
the data frame contains just the values. 
=================================================================================""")
    #Moreau_BrotoAutoTotal
    info.append("""name of the feature: Moreau_BrotoAutoTotal
role of the feature: Compute all autocorrelation descriptors based on 8 properties of AADs.
result contains 30*8=240 Moreau_Broto
autocorrelation descriptors based on the given properties(i.e.,
_AAPropert).

the data frame contains just the values. 
=================================================================================""")
    #MoranAutoTotal
    info.append("""name of the feature: MoranAutoTotal
role of the feature: Compute all autocorrelation descriptors based on 8 properties of AADs.
result contains 30*8=240 moran
autocorrelation descriptors based on the given properties(i.e.,
_AAPropert).

the data frame contains just the values. 
=================================================================================""")
    #CalculateC
    info.append("""name of the feature: CalculateC
role of the feature: Calculate all composition descriptors based seven different properties of AADs.
returns: list of values
=================================================================================""")
    #CalculateT
    info.append("""name of the feature: CalculateT
role of the feature: Calculate all Transition descriptors based seven different properties of AADs.
returns: list of values
=================================================================================""")
    #CalculateD
    info.append("""name of the feature: CalculateD
role of the feature: Calculate all Distribution  descriptors based seven different properties of AADs.
returns: list of values
=================================================================================""")
    #GetSOCN
    info.append("""name of the feature: GetSOCN
role of the feature: Sequence-order-coupling number.
This feature is derived from the distance matrix between the 20 amino acids. Here, we use two different distance matrices based on studies by Grantham and Schneider-Wrede, respectively.
returns: 
[0] - Array containing SOCN based on the Schneider-Wrede distance matrix.
[1] - Array containing SOCN based on the Grantham distance matrix.
=================================================================================""")
    #GetQSO
    info.append("""name of the feature: GetQSO
role of the feature: Quasi-sequence-order.
This feature is derived from the distance matrix between the 20 amino acids. Here, we use two different distance matrices based on studies by Grantham and Schneider-Wrede, respectively.
For each amino acid i, a quasi-sequence-order descriptor 
returns: 
[0]:Array containing QSO based on the Schneider-Wrede distance matrix.
[1]:Array containing QSO based on the Grantham distance matrix.
[2]: Order of QSO descriptors corresponding to columns in arr_sw and arr_g.

[2] : is not in the data frame
=================================================================================""")
    #aliphatic_index
    info.append("""name of the feature: aliphatic_index 
role of the feature: Compute the aliphatic index of the peptide.
The aliphatic index of a protein was proposed in Ikai (1980). It is defined as the relative volume occupied by aliphatic side chains (Alanine, Valine, Isoleucine, and Leucine). It may be regarded as a positive factor for the increase of thermostability of globular proteins.
Returns: float - The computed aliphatic index for the peptide sequence, between 0 and 100.
=================================================================================""")
    #auto_correlation
    info.append("""name of the feature: auto_correlation 
role of the feature: Compute the auto-correlation index of a peptide sequence.
returns:
i took the all the table types from the api.
=================================================================================""")
    #auto_covariance
    info.append("""name of the feature: auto_covariance 
role of the feature: Compute the auto_covariance index of a peptide sequence.
 all the table types from the api.
=================================================================================""")
    #blosum_indices
    info.append("""name of the feature: blosum_indices 
role of the feature:Compute the BLOSUM62-derived indices of the peptide.
The BLOSUM62-derived indices of a peptide.
BLOSUM indices were derived of physicochemical properties that have been subjected to a VARIMAX analysis and an alignment matrix of the 20 natural AAs using the BLOSUM62 matrix
Returns: a list. peptides.BLOSUMIndices - The computed average BLOSUM indices for all the amino acids in the peptide.
=================================================================================""")
    #boman
    info.append("""name of the feature: boman 
role of the feature: Compute the Boman (potential peptide interaction) index.
The potential interaction index proposed by Boman (2003) is an index computed by averaging the solubility values for all residues in a sequence. It can be used to give an overall estimate of the potential of a peptide to bind to membranes or other proteins.
Returns: float - The Boman index for the peptide. A value greater than 2.48 indicates that a protein has high binding potential.
=================================================================================""")
    #charge
    info.append("""name of the feature: charge 
role of the feature: Compute the theoretical net charge of a peptide sequence.
This function computes the theoretical net charge of a peptide sequence, based on the Henderson-Hasselbach equation described by Dexter S. Moore (1985). The net charge can be computed at a given pH using one of the 9 pKa scales available.
Parameters:
pH (float) - The pH value for which to compute the charge.
pKscale (str) - The name of the pKa scale to be used. A list of all the allowed values can be retrieved from the keys of the peptides.tables.PK dictionary.
Returns:
float - The net charge of the peptide.

there are 9 pKscale 
=================================================================================""")
    #counts
    info.append("""name of the feature: counts 
role of the feature: Return a table of amino-acid counts in the peptide.
Returns: dict - A dictionary mapping each amino-acid code to the number of times it occurs in the peptide sequence.

i made it a list of just the numbers. the order is the regular.
=================================================================================""")
    #cruciani_properties
    info.append("""Compute the Cruciani properties of the peptide.
The Cruciani properties of a peptide.
The Cruciani properties are a collection of scaled principal component scores that summarize a broad set of descriptors calculated based on the interaction of each amino acid residue with several chemical groups (or “probes”), such as charged ions, methyl, hydroxyl groups, and so forth.
Returns:peptides.CrucianiProperties - The computed average Cruciani properties of all the amino acids in the corresponding peptide sequence.
=================================================================================""")
    #descriptors
    info.append("""Create a dictionary containing every protein descriptor available.
i took all the values
to get all the descriptors: peptide_list[0].descriptors().keys())
=================================================================================""")
    #hydrophobic_moment
    info.append("""Compute the maximal hydrophobic moment of a protein sequence.
This function computes the hydrophobic moment based on Eisenberg et al (1984). Hydrophobic moment is a quantitative measure of the amphiphilicity perpendicular to the axis of any periodic peptide structure, such as the α-helix or β-sheet.
Parameters:angle (int) - A protein rotational angle, in degrees. Usual values are 100 for α-helix, and 160 for β-sheet.
window (int) - The size of the sliding window for which to compute the local hydrophobic moment.
Returns: float - The maximal hydrophobic moment of the peptide.

i made it a list of 2: [0] - angle =100, [1] - angel = 160
=================================================================================""")
    #hydrophobicity
    info.append("""name of the feature: hydrophobicity
role of the feature: Compute the hydrophobicity index of a protein sequence.
This function calculates the hydrophobicity index of an amino acid sequence by averaging the hydrophobicity values of each residue using one of the 39 scales from different sources.
Parameters:scale (str) - The name of the hydrophobicity scale to be used. A list of all the allowed values can be retrieved from the keys of the peptides.tables.HYDROPHOBICITY dictionary.
Returns: float - The hydrophobicity index of the peptide.

the data frame contains for each sequence a list with values for all the table of 45 types.
=================================================================================""")
    #instability_index
    info.append("""name of the feature: instability_index
role of the feature: Compute the instability index of a protein sequence.
This function calculates the instability index proposed by Guruprasad et al (1990). This index predicts the stability of a protein based on its dipeptide composition.
Returns: float - The instability index of the peptide. A protein whose instability index is smaller than 40 is predicted as stable, a value above 40 predicts that the protein may be unstable.
=================================================================================""")
    #isoelectric_point
    info.append("""
isoelectric_point(pKscale: str = 'EMBOSS') → float
Compute the isoelectric point of a protein sequence.
The isoelectric point (pI), is the pH at which a particular molecule or surface carries no net electrical charge.
Parameters:
pKscale (str) - The name of the pKa scale to be used. A list of all the allowed values can be retrieved from the keys of the peptides.tables.PK dictionary.
Returns:
float - The pH at which the peptide has a neutral net charge.
=================================================================================""")
    #mass_shift
    info.append("""name of the feature: mass_shift
role of the feature: Compute the mass difference of modified peptides.
This function calculates the mass difference of peptides introduced by chemical modifications or heavy isotope labelling.
Parameters:aa_shift (str or dict) - Either the key to a pre-defined isotope label (see peptides.tables.MASS_SHIFT), or a dictionary mapping each amino acid to it mass difference in Dalton (use nTer and cTer keys for N-terminal and C-terminal modifications).
monoisotopic (bool) - Flag whether monoisotopic weights of amino-acids should be used.
Returns: list of 3 floats - The mass difference of the modified peptide. each value is outcome of different type from aa_shift table
=================================================================================""")
    #molecular_weight
    info.append("""name of the feature: molecular_weight
role of the feature: Compute the molecular weight of a protein sequence.
This function calculates the molecular weight of a protein sequence. It is calculated as the sum of the mass of each amino acid using one of the 3 available scales. It also supports mass calculation of proteins with predefined or custom stable isotope mass labels.
Parameters: average (str) - The name of the average amino acid average weight scale. See peptides.tables.MOLECULAR_WEIGHT for a list of appropriate values.
aa_shift (str, dict or None) - Either an appropriate shift value to pass to Peptide.mass_shift, or None to get the unmodified weight.
Returns: float - The molecular weight of the peptide, in Dalton.
=================================================================================""")
    #mz
    info.append("""name of the feature: mz
role of the feature:Compute the m/z (mass over charge) ratio for a peptide.
This function calculates the (monoisotopic) mass over charge ratio (m/z) for peptides, as measured in mass spectrometry.
Parameters:charge (int) - The net charge for which the m/z should be computed.
aa_shift (str, dict or None) - Either an appropriate shift value to pass to Peptide.mass_shift, or None to get the unmodified weight.
cysteins (float) - The mass shift (in Dalton) of blocked cysteins. Default corresponds to cysteins blocked by iodoacetamide.
Returns: float - The m/z ratio of the peptide.
=================================================================================""")
    #FEGS
    info.append("""name of the feature: FEGS
role of the feature: The passage describes the FEGS model, a feature extraction model for protein sequences that aims to represent peptides or proteins using a combination of physicochemical properties and statistical information.
Graphical Feature Encoding:
Physicochemical Properties:
158 physicochemical properties of amino acids are used, selected from the AAindex database.
Amino acids are ranked based on their physicochemical indices.
The ranked amino acids are positioned on the circumference of the bottom of a right circular cone.
Amino acid pairs (400 in total) are arranged on the underside of the cone.
3D Graphical Curve Construction:
A 3D graphical curve for a peptide or protein sequence is constructed based on the right circular cone.
Each amino acid in the sequence corresponds to a point on the curve.
Coordinates of each point are determined by specific formulas involving physicochemical properties and amino acid pairs frequencies.
Graphical Curve Representation:
A nonnegative symmetric matrix M is computed for the graphical curve.
Off-diagonal entries of M are quotients of Euclidean distances between points on the curve.
Diagonal elements are set to zero.
The largest eigenvalue of M is computed and divided by the length of the sequence to characterize the graphical curve.
This process is repeated for each physicochemical property, resulting in a 158-dimensional vector representing graphical features for a peptide or protein.
Statistical Features:
Amino Acid Composition:
Reflects the frequency of 20 different amino acids in the sequence.
Generates a 20-dimensional feature vector.
Computed using the formula given.
Dipeptide Composition:
Calculates the frequency of each type of amino acid pair in the sequence.
Describes the fraction of amino acids and their local order.
Generates a 400-dimensional feature vector.
Computed using the formula provided.
Overall Feature Vector:
The combination of the graphical features (158-dimensional vector) and the statistical features (20-dimensional + 400-dimensional vectors) results in a final feature vector of 578 dimensions, which is used to represent a peptide or protein sequence in the FEGS model.
Returns: list of 578
=================================================================================""")
   

    for i in range(len(my_dataframe.columns)-2):
        feature_info[my_dataframe.columns[i+2]] = info[i]
 
   
if __name__ == "__main__":
    main()



