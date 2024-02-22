# Feature-Extraction
Python feature extraction project from protein sequences and protein structures.

During my internship, I undertook a project focused on protein sequence feature extraction, with the aim of subsequently applying these features in machine and deep learning algorithms to address protein-protein interactions. Throughout the project, I leveraged various APIs, including ProtLearn, Propy, Peptides, FEGS and iFeature.

# Sequence feature extraction

The initial phase involved downloading a FASTA file containing sequences of 100 proteins from UniProt. I developed a dedicated function to extract both the sequence and UniProt ID for each protein.

this diractory contains 3 files:

### protlearn,propy,peptides API's - feature extraction:
I gathered over 50 sequential features and amalgamated them into a vector. This was followed by the construction of a comprehensive data frame using the Pandas library. This data frame encompasses more than 10,000 columns, each representing a distinct value from one of the features extracted from the aforementioned APIs. The final column contains the UniProt IDs.

### iFeature Omega API -  feature extraction
The method similar to protlearn, propy, and peptides APIs was utilized, resulting in the extraction of over 20,000 values.

### FEGS API - feature extraction
contains 578 values.



# Structures feature extraction

I've took several UniProt IDs to extract their corresponding PDB files from the AlphaFold database. Subsequently, I generated a one-dimensional vector for each protein. Using the iFeature Omega API, I obtained 13 structural features and integrated them into the vector. Each resulting vector comprises 722 values.

