import re
import numpy as np
import pandas as pd
import iFeatureOmegaCLI


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

def vector_to_dataframe(vectors_list, list_of_uniprots):
    data = {
    'uniprotId': list_of_uniprots,
    'vector': vectors_list
}
    df = pd.DataFrame(data)
 
    return df

def proteins_instance_create(sequences,sequences_num):
    proteins = []
    for i in range(sequences_num):
        file = open(r"C:\Users\97252\Desktop\עמית\אוניברסיטה\שנה ג\התמחות\חלבונים\protein.txt", "w")
        file.write(">")
        file.write(sequences[i][0])
        file.write("\n")
        file.write(sequences[i][1])
        file.write("\n")
        file.close()
        protein = iFeatureOmegaCLI.iProtein(r"C:\Users\97252\Desktop\עמית\אוניברסיטה\שנה ג\התמחות\חלבונים\protein.txt")
        proteins.append(protein)

    return proteins    
       
def protein_features_to_vector(protein):
    vector = []

    #20 values 
    protein.get_descriptor("AAC")
    res = protein.encodings
    num_seq, len_features = res.shape
    vector.append(res.iloc[0].values)  
    
    #one value cause the length of the sequence is not the same. i made mean
    protein.get_descriptor("EAAC")
    res = protein.encodings
    arr = np.array([res.iloc[0].values.mean()])
    vector.append(arr)  

    #1600 values 
    protein.get_descriptor("CKSAAP type 1")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #1600 values 
    protein.get_descriptor("CKSAAP type 2")
    res = protein.encodings
    vector.append(res.iloc[0].values) 

    #400 values 
    protein.get_descriptor("DPC type 1")
    res = protein.encodings
    vector.append(res.iloc[0].values) 

    #400 values 
    protein.get_descriptor("DPC type 2")
    res = protein.encodings
    vector.append(res.iloc[0].values) 

    #8000 values 
    protein.get_descriptor("TPC type 1")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #8000 values 
    protein.get_descriptor("TPC type 2")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #39 values 
    protein.get_descriptor("CTDC")
    res = protein.encodings
    vector.append(res.iloc[0].values) 

    #39 values 
    protein.get_descriptor("CTDT")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #195 values 
    protein.get_descriptor("CTDD")
    res = protein.encodings
    vector.append(res.iloc[0].values)       
        
    #343 values 
    protein.get_descriptor("CTriad")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #1379 values 
    protein.get_descriptor("KSCTriad")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #400 values 
    protein.get_descriptor("ASDC")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #20 values 
    protein.get_descriptor("DistancePair")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #20 values 
    protein.get_descriptor("GAAC")
    res = protein.encodings
    vector.append(res.iloc[0].values)


    #1 value cause the length of the sequence is not the same. i made mean
    protein.get_descriptor("EGAAC")
    res = protein.encodings
    arr = np.array([res.iloc[0].values.mean()])
    vector.append(arr)

    #100 values 
    protein.get_descriptor("CKSAAGP type 1")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #100 values 
    protein.get_descriptor("CKSAAGP type 2")
    res = protein.encodings
    vector.append(res.iloc[0].values)


    #25 values 
    protein.get_descriptor("GDPC type 1")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #25 values 
    protein.get_descriptor("GDPC type 2")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #125 values 
    protein.get_descriptor("GTPC type 1")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #125 values 
    protein.get_descriptor("GTPC type 1")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #24 values 
    protein.get_descriptor("Moran")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #24 values 
    protein.get_descriptor("Geary")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #24 values 
    protein.get_descriptor("NMBroto")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #24 values 
    protein.get_descriptor("AC")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #6 values 
    protein.get_descriptor("SOCNumber")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #46 values 
    protein.get_descriptor("QSOrder")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #23 values 
    protein.get_descriptor("PAAC")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #26 values 
    protein.get_descriptor("APAAC")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #4 values 
    protein.get_descriptor("PseKRAAC type 1")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #4 values 
    protein.get_descriptor("PseKRAAC type 2")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #4 values 
    protein.get_descriptor("PseKRAAC type 3A")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #4 values 
    protein.get_descriptor("PseKRAAC type 3B")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #25 values 
    protein.get_descriptor("PseKRAAC type 4")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #9 values 
    protein.get_descriptor("PseKRAAC type 5")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #16 values 
    protein.get_descriptor("PseKRAAC type 6A")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #25 values 
    protein.get_descriptor("PseKRAAC type 6B")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #25 values 
    protein.get_descriptor("PseKRAAC type 6C")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #4 values 
    protein.get_descriptor("PseKRAAC type 7")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #4 values 
    protein.get_descriptor("PseKRAAC type 8")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #4 values 
    protein.get_descriptor("PseKRAAC type 9")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #4 values 
    protein.get_descriptor("PseKRAAC type 10")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #4 values 
    protein.get_descriptor("PseKRAAC type 11")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #4 values 
    protein.get_descriptor("PseKRAAC type 12")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #4 values 
    protein.get_descriptor("PseKRAAC type 13")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #4 values 
    protein.get_descriptor("PseKRAAC type 14")
    res = protein.encodings
    vector.append(res.iloc[0].values)

    #4 values 
    protein.get_descriptor("PseKRAAC type 15")
    res = protein.encodings
    vector.append(res.iloc[0].values)
    
    #4 values 
    protein.get_descriptor("PseKRAAC type 16")
    res = protein.encodings
    vector.append(res.iloc[0].values)




     #1 value cause the length of the sequence is not the same. i made mean
    protein.get_descriptor("binary")
    res = protein.encodings
    arr = np.array([res.iloc[0].values.mean()])
    vector.append(arr)

     #1 value cause the length of the sequence is not the same. i made mean 
    protein.get_descriptor("binary_6bit")
    res = protein.encodings
    arr = np.array([res.iloc[0].values.mean()])
    vector.append(arr)

    #1 value cause the length of the sequence is not the same. i made mean 
    protein.get_descriptor("binary_5bit type 1")
    res = protein.encodings
    arr = np.array([res.iloc[0].values.mean()])
    vector.append(arr)

    #1 value cause the length of the sequence is not the same. i made mean 
    protein.get_descriptor("binary_5bit type 2")
    res = protein.encodings
    arr = np.array([res.iloc[0].values.mean()])
    vector.append(arr)

    #1 value cause the length of the sequence is not the same. i made mean 
    protein.get_descriptor("binary_3bit type 1")
    res = protein.encodings
    arr = np.array([res.iloc[0].values.mean()])
    vector.append(arr)

    #1 value cause the length of the sequence is not the same. i made mean 
    protein.get_descriptor("binary_3bit type 2")
    res = protein.encodings
    arr = np.array([res.iloc[0].values.mean()])
    vector.append(arr)

    #1 value cause the length of the sequence is not the same. i made mean
    protein.get_descriptor("binary_3bit type 3")
    res = protein.encodings
    arr = np.array([res.iloc[0].values.mean()])
    vector.append(arr)

    #1 value cause the length of the sequence is not the same. i made mean 
    protein.get_descriptor("binary_3bit type 4")
    res = protein.encodings
    arr = np.array([res.iloc[0].values.mean()])
    vector.append(arr)

    #1 value cause the length of the sequence is not the same. i made mean 
    protein.get_descriptor("binary_3bit type 5")
    res = protein.encodings
    arr = np.array([res.iloc[0].values.mean()])
    vector.append(arr)
    #1 value cause the length of the sequence is not the same. i made mean
    protein.get_descriptor("binary_3bit type 6")
    res = protein.encodings
    arr = np.array([res.iloc[0].values.mean()])
    vector.append(arr)

    #1 value cause the length of the sequence is not the same. i made mean 
    protein.get_descriptor("binary_3bit type 7")
    res = protein.encodings
    arr = np.array([res.iloc[0].values.mean()])
    vector.append(arr)

    #3 value cause the features return 3*num of the sequence length. i made mean for each feature 
    protein.get_descriptor("AESNN3")
    res = protein.encodings
    num_seq_seq, len_features = res.shape
    sum0 =0
    sum1 =0
    sum2 = 0
    x = res.iloc[0].values
    i=0
    for j in range(len_features):
        if (i%3 == 0):
            sum0 +=  x[i]
        if (i%3 == 1):
            sum1 +=  x[i]
        if (i%3 == 2):
            sum2 += x[i]

        i += 1
    z = [sum0,sum1,sum2] 

    vector.append(np.array(z))
    



    #10 value cause the features return 3*num of the sequence length. i made mean for each feature 
    protein.get_descriptor("OPF_10bit")
    res = protein.encodings
    num_seq_seq, len_features = res.shape
    sum0=sum1=sum2=sum3=sum4=sum5=sum6=sum7=sum8=sum9 =0
  
    x = res.iloc[0].values
   
    i=0
    for j in range(len_features):
        if (i%10 == 0):
            sum0 +=  x[i]
        if (i%10 == 1):
            sum1 +=  x[i]
        if (i%10 == 2):
            sum2 += x[i]
        if (i%10 == 3):
            sum3 +=  x[i]
        if (i%10 == 4):
            sum4 +=  x[i]
        if (i%10 == 5):
            sum5 += x[i] 
        if (i%10 == 6):
            sum6 +=  x[i]
        if (i%10 == 7):
            sum7 +=  x[i]
        if (i%10 == 8):
            sum8 += x[i]  
        if (i%10 == 9):
            sum9 += x[i]        

        i += 1
    z = [sum0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9]         
    vector.append(np.array(z))



    #7 value cause the features return 3*num of the sequence length. i made mean for each feature
    protein.get_descriptor("OPF_7bit type 1") 
    res = protein.encodings
    num_seq_seq, len_features = res.shape
    sum0=sum1=sum2=sum3=sum4=sum5=sum6 =0
  
    x = res.iloc[0].values
   
    i=0
    for j in range(len_features):
        if (i%7 == 0):
            sum0 +=  x[i]
        if (i%7 == 1):
            sum1 +=  x[i]
        if (i%7 == 2):
            sum2 += x[i]
        if (i%7 == 3):
            sum3 +=  x[i]
        if (i%7 == 4):
            sum4 +=  x[i]
        if (i%7 == 5):
            sum5 += x[i]
        if (i%7 == 6):
            sum6 +=  x[i]

        i += 1     
    z= [sum0,sum1,sum2,sum3,sum4,sum5,sum6]    
    vector.append(np.array(z))
    
    
    
     #7 value cause the features return 3*num of the sequence length. i made mean for each feature
    protein.get_descriptor("OPF_7bit type 2") 
    res = protein.encodings
    num_seq_seq, len_features = res.shape
    sum0=sum1=sum2=sum3=sum4=sum5=sum6 =0
  
    x = res.iloc[0].values
   
    i=0
    for j in range(len_features):
        if (i%7 == 0):
            sum0 +=  x[i]
        if (i%7 == 1):
            sum1 +=  x[i]
        if (i%7 == 2):
            sum2 += x[i]
        if (i%7 == 3):
            sum3 +=  x[i]
        if (i%7 == 4):
            sum4 +=  x[i]
        if (i%7 == 5):
            sum5 += x[i]
        if (i%7 == 6):
            sum6 +=  x[i]

        i += 1     
    z= [sum0,sum1,sum2,sum3,sum4,sum5,sum6]    
    vector.append(np.array(z))
    
    

    #7 value cause the features return 3*num of the sequence length. i made mean for each feature
    protein.get_descriptor("OPF_7bit type 3") 
    res = protein.encodings
    num_seq_seq, len_features = res.shape
    sum0=sum1=sum2=sum3=sum4=sum5=sum6 =0
  
    x = res.iloc[0].values
   
    i=0
    for j in range(len_features):
        if (i%7 == 0):
            sum0 +=  x[i]
        if (i%7 == 1):
            sum1 +=  x[i]
        if (i%7 == 2):
            sum2 += x[i]
        if (i%7 == 3):
            sum3 +=  x[i]
        if (i%7 == 4):
            sum4 +=  x[i]
        if (i%7 == 5):
            sum5 += x[i]
        if (i%7 == 6):
            sum6 +=  x[i]

        i += 1  
    z= [sum0,sum1,sum2,sum3,sum4,sum5,sum6]    
    vector.append(np.array(z))
    


    #8 value cause the features return 3*num of the sequence length. i made mean for each feature
    protein.get_descriptor("AAIndex") 
    res = protein.encodings
    num_seq_seq, len_features = res.shape
    sum0=sum1=sum2=sum3=sum4=sum5=sum6=sum7 =0
  
    x = res.iloc[0].values
   
    i=0
    for j in range(len_features):
        if (i%8 == 0):
            sum0 +=  x[i]
        if (i%8 == 1):
            sum1 +=  x[i]
        if (i%8 == 2):
            sum2 += x[i]
        if (i%8 == 3):
            sum3 +=  x[i]
        if (i%8 == 4):
            sum4 +=  x[i]
        if (i%8 == 5):
            sum5 += x[i]
        if (i%8 == 6):
            sum6 +=  x[i]
        if (i%8 == 7):
            sum7 +=  x[i]    

        i += 1  
    z = [sum0,sum1,sum2,sum3,sum4,sum5,sum6,sum7]   
    vector.append(np.array(z))


    #20 value cause the features return 3*num of the sequence length. i made mean for each feature 
    protein.get_descriptor("BLOSUM62")
    res = protein.encodings
    num_seq_seq, len_features = res.shape
    sum0=sum1=sum2=sum3=sum4=sum5=sum6=sum7=sum8 =sum9=sum10=sum11=sum12=sum13=sum14=sum15=sum16=sum17=sum18=sum19 =0
  
    x = res.iloc[0].values
   
    i=0
    for j in range(len_features):
        if (i%20 == 0):
            sum0 +=  x[i]
        if (i%20 == 1): 
            sum1 +=  x[i]
        if (i%20 == 2):
            sum2 += x[i]    
        if (i%20 == 3):
            sum3 +=  x[i]
        if (i%20 == 4):
            sum4 +=  x[i]
        if (i%20 == 5):
            sum5 += x[i]
        if (i%20 == 6):
            sum6 +=  x[i]
        if (i%20 == 7):
            sum7 +=  x[i]
        if (i%20 == 8):
            sum8 +=  x[i]
        if (i%20 == 9):
            sum9 +=  x[i]
        if (i%20 == 10):
            sum10 +=  x[i]
        if (i%20 == 11):
            sum11 +=  x[i]
        if (i%20 == 12):
            sum12 +=  x[i]
        if (i%20 == 13):
            sum13 +=  x[i]
        if (i%20 == 14):
            sum14 +=  x[i]
        if (i%20 == 15):
            sum15 +=  x[i]
        if (i%20 == 16):
            sum16 +=  x[i]
        if (i%20 == 17):
            sum17 +=  x[i]
        if (i%20 == 18):
            sum18 +=  x[i]
        if (i%20 == 19):
            sum19 +=  x[i]



        i += 1
    z= [sum0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,sum17,sum18,sum19]     
    vector.append(np.array(z))





    #5 value cause the features return 3*num of the sequence length. i made mean for each feature
    protein.get_descriptor("ZScale") 
    res = protein.encodings
    num_seq_seq, len_features = res.shape
    sum0=sum1=sum2=sum3=sum4 =0
  
    x = res.iloc[0].values
   
    i=0
    for j in range(len_features):
        if (i%5 == 0):
            sum0 +=  x[i]
        if (i%5 == 1):
            sum1 +=  x[i]
        if (i%5 == 2):
            sum2 += x[i]
        if (i%5 == 3):
            sum3 +=  x[i]
        if (i%5 == 4):
            sum4 +=  x[i] 

        i += 1   
    z = [sum0,sum1,sum2,sum3,sum4] 
    vector.append(np.array(z))


    return vector



        
        


def main():

    file_path = r'write your path for the uniprots'
    fasta_data = read_fasta(file_path)
    sequences = extract_sequences(fasta_data)
    sequences_num = len(sequences) 
    uniprotID_list = [i[0] for i in sequences]
    proteins = proteins_instance_create(sequences,sequences_num)


  
    vectors_list = []
    for i in range (sequences_num):
       x = protein_features_to_vector(proteins[i])
       x= np.array(x)
       vectors_list.append(np.concatenate(x))
       
    
    df = vector_to_dataframe(vectors_list, uniprotID_list)
    

    df.to_csv('seq_features.csv', index=False) there are to much values - 23378 to each vector the csv cant hold it
    
  




   
    
    



   

    return None




if __name__ == "__main__":
    main()



