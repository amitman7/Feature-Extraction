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


def import_pdb_to_ifeature(list_of_uniprots):
    structure_list = []
    for pdb_id  in list_of_uniprots:
        structure = iFeatureOmegaCLI.iStructure(rf"C:\Users\97252\Desktop\עמית\אוניברסיטה\שנה ג\התמחות\structure_extarction\pdb\{pdb_id}.pdb")
        structure_list.append(structure)
    for i in range (len(structure_list)):
        structure_list[i].import_parameters(r'C:\Users\97252\Desktop\עמית\אוניברסיטה\שנה ג\התמחות\structure_extarction\Structure_parameters_setting.json') 

    return structure_list


def protein_features_to_vector(protein):
    vector = []

    protein.get_descriptor("AAC_type1")
    res = protein.encodings
    len_seq, len_features = res.shape
    for i in range(len_features):
        x= res.iloc[:, i].mean()
        vector.append(x)  

    protein.get_descriptor("AAC_type2")
    res = protein.encodings
    len_seq, len_features = res.shape
    for i in range(len_features):
        x= res.iloc[:, i].mean()
        vector.append(x)  

    protein.get_descriptor("GAAC_type1")
    res = protein.encodings
    len_seq, len_features = res.shape
    for i in range(len_features):
        x= res.iloc[:, i].mean()
        vector.append(x)  

    protein.get_descriptor("GAAC_type2")
    res = protein.encodings
    len_seq, len_features = res.shape
    for i in range(len_features):
        x= res.iloc[:, i].mean()
        vector.append(x)  

    protein.get_descriptor("SS3_type1")
    res = protein.encodings
    len_seq, len_features = res.shape
    for i in range(len_features):
        x= res.iloc[:, i].mean()
        vector.append(x)  

    protein.get_descriptor("SS3_type2")
    res = protein.encodings
    len_seq, len_features = res.shape
    for i in range(len_features):
        x= res.iloc[:, i].mean()
        vector.append(x)  

    protein.get_descriptor("SS8_type1")
    res = protein.encodings
    len_seq, len_features = res.shape
    for i in range(len_features):
        x= res.iloc[:, i].mean()
        vector.append(x)  

    protein.get_descriptor("SS8_type2")
    res = protein.encodings
    len_seq, len_features = res.shape
    for i in range(len_features):
        x= res.iloc[:, i].mean()
        vector.append(x)  

    protein.get_descriptor("HSE_CA")
    res = protein.encodings
    len_seq, len_features = res.shape
    for i in range(len_features):
        x= res.iloc[:, i].mean()
        vector.append(x)  

    protein.get_descriptor("HSE_CB")
    res = protein.encodings
    len_seq, len_features = res.shape
    for i in range(len_features):
        x= res.iloc[:, i].mean()
        vector.append(x)  

    protein.get_descriptor("AC_type1")
    res = protein.encodings
    len_seq, len_features = res.shape
    for i in range(len_features):
        x= res.iloc[:, i].mean()
        vector.append(x)  

    protein.get_descriptor("AC_type2")
    res = protein.encodings
    len_seq, len_features = res.shape
    for i in range(len_features):
        x= res.iloc[:, i].mean()
        vector.append(x)  

    protein.get_descriptor("Network-based index")
    res = protein.encodings
    len_seq, len_features = res.shape
    for i in range(len_features):
        x= res.iloc[:, i].mean()
        vector.append(x)  



    return vector   
        


def main():

    """list_of_uniprots = ["Q92800", "O75530", "Q6ZN18", "Q09028"]
    for uniprot in list_of_uniprots:
        download_alphafold_pdb(uniprot)"""
    
    list_of_uniprots = ["Q92800", "O75530", "Q6ZN18", "Q09028"]
    structure_list = import_pdb_to_ifeature(list_of_uniprots)
    


    vectors_list = []
    for i in range (len(structure_list)):
       vectors_list.append(protein_features_to_vector(structure_list[i]))
        
    

    return None






if __name__ == "__main__":
    main()



# first lets make each feature a vector of size length * num of features without len
   #then lets take all the vectors of all feature and combuend tehm
