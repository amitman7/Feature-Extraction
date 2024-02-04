# Feature-extraction
Python feature extraction project from protein sequences and protein structures.

During my internship, I undertook a project focused on protein sequence feature extraction, with the aim of subsequently applying these features in machine and deep learning algorithms to address protein-protein interactions. Throughout the project, I leveraged various APIs, including ProtLearn, Propy, Peptides, FEGS, iFeature.

# Sequence feature extraction
The initial phase involved downloading a FASTA file containing sequences of 100 proteins from UniProt. I developed a dedicated function to extract both the sequence and UniProt ID for each protein, followed by the creation of a comprehensive data frame using the Pandas library. This data frame encompasses over 50 columns, each representing a distinct feature extracted from the aforementioned APIs.

Furthermore, the data frame provides detailed expansions for each feature, elucidating the underlying rules and functionalities associated with them.

# Structures feature extraction

I've took several UniProt IDs to extract their corresponding PDB files from the AlphaFold database, followed by the creation of a comprehensive data frame using the Pandas library. Utilizing the iFeature Omega API, I fetched 14 structural features, each represented as a column in this comprehensive data frame. 

