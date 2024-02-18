import re
import numpy as np
import pandas as pd
import os
import requests
import matlab.engine



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

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        return file.read()


##### FEGS ########

def FEGS_add(fasta_data,fegs_path):
    # Start MATLAB engine
    eng = matlab.engine.start_matlab()
    
    # Set path to FEGS
    eng.addpath(fegs_path)
    print("im here")
    
    # Execute FEGS
    result = eng.FEGS(fasta_data)

    # Stop MATLAB engine
    eng.quit()

    result_array = list(result)
    print("hi")

    return result_array


def vector_to_dataframe(uniprotID_list,vectors_list):

    columns =  [str(i) for i in range(1, 579)]

    df = pd.DataFrame(vectors_list, columns=columns)
    df['uniprotID'] = uniprotID_list
 
    return df



def main():

    fegs_path = 'write your path for fegs'
    fasta_data = read_fasta('write your path to your fasta file') 
    vectors_list= FEGS_add(fasta_data,fegs_path)

    sequences = extract_sequences(fasta_data)
    uniprotID_list = [i[0] for i in sequences]
    uniprotID_num = len(uniprotID_list)

    df = vector_to_dataframe(uniprotID_list,vectors_list)

  


    # Save the DataFrame to a CSV file
    df.to_csv('fegs_features.csv', index=False)



    return None






if __name__ == "__main__":
    main()


