# Steps for preprocessing and calculation of LLR score 
LLR_score_generation steps:
1) Remove duplicates from all the input files

Script_R_ESM1b steps:
2) Get the wildytpe and mutated DNA sequences of each SARS-CoV-2 variant

LLR_score_generation steps:
3) Generate protein sequences from wildtype and mutated DNA sequences

Script_R_ESM1b steps:
4) Get the SNP between wildtype and mutated protein sequences
5) Compare the mutant DNA sequence with the WT DNA sequence as well as protein sequences to identify differences, which should correspond to a single change and match the reference and alternative sequences.

LLR_score_generation steps:
6) Take the protein sequences from above results
7) Create a Seq ID for each variant along with its WT protein sequence and give it as an input to the ESM1b model for prediction.
8) Get the ESM1b output which has seq_id , mut_name and esm_score.
9) Create a universal ID with seq_id and mut_name in the final data frame.
10) Use this universal ID to merge the scores with the missense variants. 
11) You have a final sheet with all variants and their LLR score.
