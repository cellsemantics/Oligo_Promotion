# Steps for preprocessing and calculation of kGain score 
Script_flank_generation steps:
1) Get the SNP positions
2) Read the reference file for SARS-CoV-2 virus
3) Get the reference DNA using reference file
4) Specify the length of the flanking region
5) Extract the flanking region from reference sequences using SNP position
6) Similarly, extract the mutated sequences from reference sequences using SNP position and extract flanking regions.

kGain_script steps:
7) Read the flanks
8) Define kmer length
9) Generate fcgr_sequence_covid.pickle file and generate density plot 
10) Add altered and reference sliding window
11) Calculate gain score and using that, calculate accumulated gain which is kGain score.