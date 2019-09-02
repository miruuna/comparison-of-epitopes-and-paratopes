# How to obtain contacting resdiues and visualise the data
1. **Get all pdbs from IMGT associated with lysozyme and their chain annotation**
- run **get_annotation.py** to store the list of pdb codes and the chain annotation in json format files
2. **Get interacting residues**
- script **get_contacts_list_lysozyme.py** will be used in **step 4** to store the interacting residues
3. **MSA of lysozyme sequences and antibodies**
- download the sequences of the lysozyme proteins in each complex in fasta format by running **get_fasta_files()** of **seq_antibody.py**
- process the fasta files by running **process_fasta()** and **edit_fasta_for_msa()** of **seq_antibody.py**. *you need to create a new folder fasta_antigen or fasta_antibody. The chains in **get_annotation_string()** have to be changed depending on whether you would like to obtain the antibody or antigen sequences.
- run **write_new_file()** to store the sequences ready to be submitted for MSA on EBI's Clustal Omega tool. The alignment is then obtained by clicking on *Get Alignment file* and the link will be used in the next step for computing the equivalent positions.
4. **Store the equivalent positions of contacting residues( names and positions) 
- run **equiv_pos_lyso_all_contacts.py** to store the contacting residues for the chain/chains of choice in json format
5. **Filter the complexes with unique Heavy Chain CDRS**
- run **get_unique_cdrs("heavy")** of **remove_similar.py** to filter out the complexes with identical CDRs
- run steps 2, 3 and 4 again for the new list of pdbs of complexes of antibodies with unique CDRs
6. **Visualise data**
- create different heatmpas using the functions in **visualise_data.py**
7. **Run statistical tests**
- run the different functions in **compute_similarity.py**
