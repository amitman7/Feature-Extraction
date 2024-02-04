import re
import numpy as np
import pandas as pd
import os
import iFeatureOmegaCLI
import os
import requests

def download_alphafold_pdb(uniprot_id):
    """
    Downloads the PDB file for a given UniProt ID from the AlphaFold database and saves it
    in a folder named 'pdb_files' in directory.

    :param uniprot_id: UniProt ID of the protein
    """

    # Base URL for AlphaFold API
    base_url = "https://alphafold.ebi.ac.uk/api/prediction"

    # API key for requests
    api_key = "AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94"

    # Directory where the PDB files will be saved
    save_directory = r"C:\Users\97252\Desktop\עמית\אוניברסיטה\שנה ג\התמחות\structure_extarction\pdb"
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)

    # File path for the PDB file
    file_path = os.path.join(save_directory, f"{uniprot_id}.pdb")

    # Building the full API request URL
    request_url = f"{base_url}/{uniprot_id}?key={api_key}"

    # Fetching prediction details
    response = requests.get(request_url, headers={'accept': 'application/json'})
    if response.status_code != 200:
        print(f"Failed to fetch prediction details for {uniprot_id}: HTTP {response.status_code}")
        return

    # Parsing the JSON response
    try:
        prediction_data = response.json()[0]  # Taking the first element from the list
    except (IndexError, KeyError):
        print(f"No data found for {uniprot_id}")
        return

    # Extracting PDB file URL
    pdb_file_url = prediction_data.get("pdbUrl")
    if not pdb_file_url:
        print(f"No PDB file URL found for {uniprot_id}")
        return

    # Downloading the PDB file
    pdb_response = requests.get(pdb_file_url)
    if pdb_response.status_code != 200:
        print(f"Failed to download PDB file for {uniprot_id}: HTTP {pdb_response.status_code}")
        return

    # Saving the PDB file
    with open(file_path, 'w') as file:
        file.write(pdb_response.text)
    print(f"PDB file for {uniprot_id} saved to {file_path}")

def feature_to_dataframe(uniprots_id_list):
    data = []
    for i in range(len(uniprots_id_list)):
        data.append(uniprots_id_list[i])
    df = pd.DataFrame(data, columns=['uniprots_id'])
    return df

def import_pdb_to_ifeature(list_of_uniprots):
    structure_list = []
    for pdb_id  in list_of_uniprots:
        structure = iFeatureOmegaCLI.iStructure(rf"C:\Users\97252\Desktop\עמית\אוניברסיטה\שנה ג\התמחות\structure_extarction\pdb\{pdb_id}.pdb")
        structure_list.append(structure)
    return structure_list


#Amino acids content type 1
def AAC_type1_add(structure_list,df):
    df['AAC_type1'] = None
    list_toAdd = []
    for i in range(len(structure_list)):
        structure_list[i].get_descriptor("AAC_type1")
        res = structure_list[i].encodings
        list_toAdd.append(res)
    df['AAC_type1'] = list_toAdd 

#Amino acids content type 
def AAC_type2_add(structure_list,df):
    df['AAC_type2'] = None
    list_toAdd = []
    for i in range(len(structure_list)):
        structure_list[i].get_descriptor("AAC_type2")
        res = structure_list[i].encodings
        list_toAdd.append(res)
    df['AAC_type2'] = list_toAdd

#Grouped amino acids content type 1
def GAAC_type1_add(structure_list,df):
    df['GAAC_type1'] = None
    list_toAdd = []
    for i in range(len(structure_list)):
        structure_list[i].get_descriptor("GAAC_type1")
        res = structure_list[i].encodings
        list_toAdd.append(res)
    df['GAAC_type1'] = list_toAdd     
    
#Grouped amino acids content type 2
def GAAC_type2_add(structure_list,df):
    df['GAAC_type2'] = None
    list_toAdd = []
    for i in range(len(structure_list)):
        structure_list[i].get_descriptor("GAAC_type2")
        res = structure_list[i].encodings
        list_toAdd.append(res)
    df['GAAC_type2'] = list_toAdd     
    
#Secondary structure elements (3) type 1
def SS3_type1_add(structure_list,df):
    df['GAAC_type2'] = None
    list_toAdd = []
    for i in range(len(structure_list)):
        structure_list[i].get_descriptor("SS3_type1")
        res = structure_list[i].encodings
        list_toAdd.append(res)
    df['SS3_type1'] = list_toAdd  

#Secondary structure elements (3) type 2
def SS3_type2_add(structure_list,df):
    df['SS3_type2'] = None
    list_toAdd = []
    for i in range(len(structure_list)):
        structure_list[i].get_descriptor("SS3_type2")
        res = structure_list[i].encodings
        list_toAdd.append(res)
    df['SS3_type2'] = list_toAdd  

# Secondary structure elements (8) type 1
def SS8_type1_add(structure_list,df):
    df['SS8_type1'] = None
    list_toAdd = []
    for i in range(len(structure_list)):
        structure_list[i].get_descriptor("SS8_type1")
        res = structure_list[i].encodings
        list_toAdd.append(res)
    df['SS8_type1'] = list_toAdd  

# Secondary structure elements (8) type 2
def SS8_type2_add(structure_list,df):
    df['SS8_type2'] = None
    list_toAdd = []
    for i in range(len(structure_list)):
        structure_list[i].get_descriptor("SS8_type2")
        res = structure_list[i].encodings
        list_toAdd.append(res)
    df['SS8_type2'] = list_toAdd  

# Half sphere exposure alpha
def HSE_CA_add(structure_list,df):
    df['HSE_CA'] = None
    list_toAdd = []
    for i in range(len(structure_list)):
        structure_list[i].get_descriptor("HSE_CA")
        res = structure_list[i].encodings
        list_toAdd.append(res)
    df['HSE_CA'] = list_toAdd  

# Half sphere exposure β
def HSE_CB_add(structure_list,df):
    df['HSE_CB'] = None
    list_toAdd = []
    for i in range(len(structure_list)):
        structure_list[i].get_descriptor("HSE_CB")
        res = structure_list[i].encodings
        list_toAdd.append(res)
    df['HSE_CB'] = list_toAdd  


# Residue depth no "mems" in windows cans use it
"""def Residue_depth_add(structure_list,df):
    df['Residue depth'] = None
    list_toAdd = []
    for i in range(len(structure_list)):
        structure_list[i].get_descriptor("Residue depth")
        res = structure_list[i].encodings
        list_toAdd.append(res)
    df['Residue depth'] = list_toAdd  """


# Atom content type 1
def AC_type1_add(structure_list,df):
    df['AC_type1'] = None
    list_toAdd = []
    for i in range(len(structure_list)):
        structure_list[i].get_descriptor("AC_type1")
        res = structure_list[i].encodings
        list_toAdd.append(res)
    df['AC_type1'] = list_toAdd  


# Atom content type 2
def AC_type2_add(structure_list,df):
    df['AC_type2'] = None
    list_toAdd = []
    for i in range(len(structure_list)):
        structure_list[i].get_descriptor("AC_type2")
        res = structure_list[i].encodings
        list_toAdd.append(res)
    df['AC_type2'] = list_toAdd  


# Network-based index
def Network_based_index_add(structure_list,df):
    df['Network-based index'] = None
    list_toAdd = []
    for i in range(len(structure_list)):
        structure_list[i].get_descriptor("Network-based index")
        res = structure_list[i].encodings
        list_toAdd.append(res)
    df['Network-based index'] = list_toAdd  



def main():

    """list_of_uniprots = ["Q92800", "O75530", "Q6ZN18", "Q09028"]
    for uniprot in list_of_uniprots:
        download_alphafold_pdb(uniprot)"""
    
    list_of_uniprots = ["Q92800", "O75530", "Q6ZN18", "Q09028"]
    df = feature_to_dataframe(list_of_uniprots)
    structure_list = import_pdb_to_ifeature(list_of_uniprots)
    
    structure_list[0].import_parameters(#write your path to parameters file)
    
    AAC_type1_add(structure_list,df)  
    AAC_type2_add(structure_list,df)
    GAAC_type1_add(structure_list,df)
    GAAC_type2_add(structure_list,df)
    SS3_type1_add(structure_list,df)
    SS3_type2_add(structure_list,df)
    SS8_type1_add(structure_list,df)
    SS8_type2_add(structure_list,df)
    HSE_CA_add(structure_list,df)
    HSE_CB_add(structure_list,df)
    AC_type1_add(structure_list,df)
    AC_type2_add(structure_list,df)
    Network_based_index_add(structure_list,df)

           
   
    return None









if __name__ == "__main__":
    main()
