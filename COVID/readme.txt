# Steps for preprocessing 
Script1 steps: 
1) Extract the DNA and Protein sequence as per gene from the NCBI file into a data frame. 
2) Merge the DNA and Protein sequence to the missense variants file from E.coli LTEE.

Script2 steps: 
3) Separate the negative and positive strands to be processed separately. 
4) Regenerate the mutant DNA sequence from the WT DNA sequence in order to get only bp change
5) Compare the mutant DNA sequence with the WT DNA sequence to identify differences, which should correspond to a single change and match the reference and alternative sequences.

Script1 steps: 
6) Translate the WT DNA sequence to the protein sequence and match it with NCBI-given sequences. Takes the ones that perfectly matches.
7) Create a Seq ID for each variant along with its WT protein sequence and give it as an input to the ESM1b model for prediction.
8) Now, translate the mutant DNA to the mutant protein sequence. 

Sript2 steps:
9) Find the difference between the WT protein and mutant protein sequence. It should be one difference. 
10) Now save this final data frame

Script1 steps: 
11) Take the ESM1b output which has seq_id , mut_id and score, create a universal ID with seq_id and mut_id.
12) create a universal ID with seq_id and SNP column in the final data frame.
13) Use this universal ID to merge the scores with the missense variants. 
14) You have a final sheet with all variants and their ESM1b score.

Script3 is for generating plots.