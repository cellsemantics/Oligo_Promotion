# Selective promotion of oligonucleotides in the course of evolution 
This repository offers a collection of resources, tools, and scripts created for the research paper titled "**Selective promotion of oligonucleotides in the course of evolution**" by Bernadette Mathew, Abhishek Halder, Nancy Jaiswal, Smruti Panda, Debjit Pramanik, Sreeram Chandra Murthy Peela, Abhishek Garg, Sadhana Tripathi,Prashant Gupta, Vandana Malhotra, Gaurav Ahuja, Debarka Sengupta.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

**Contents :** 

1. Environment setup: The environment1.yml file is for Python environment generation
2. Score generations : Folder named " score_generation" contains the scripts for generating two scores essential for this analysis .i.e., the kGain and the LLR scores 
3. Figures: The "figure" folder contains the scripts and data associated with the figures presented in the manuscript and supplementary materials.

______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
**Work-flow for implementing the kGain scores :**

**Oligo Promotion** 

![snip_of_kgain](https://github.com/user-attachments/assets/59b4fcda-7535-43e7-9b07-cc9f18a87141)


Here are the steps for computing the kGain score:

### 1. **Generate Sequences with SNV and Flanks:**
   - For each variant, create two sequences, each 21 nucleotides long (assuming a k-mer length of 10).
   - One sequence includes the variant allele at the center (11th position), while the other includes the reference allele at the same position.
   - The left and right flanking regions are taken from the corresponding reference genome, depending on the organism.

### 2. **Generate k-mers Using Rolling Windows:**
   - For each variant, generate a total of **k** (k = 10) rolling windows, where each window contains either the reference or alternate allele.
   - Apply the rolling window method twice: once for the reference sequence and once for the sequence with the variant.
   - This results in sets of k-mers for each variant, with each k-mer containing the position of the variant.

### 3. **Compute kGain Score:**
   - Track the occurrence of each k-mer across the reference genome.
   - For each k-mer, calculate the fold change between the genomic frequencies of the k-mers containing the alternate allele ($F^{\text{alt}}_{i(v)}$) and the reference allele (F^Ref(i(v)))  .
   - Compute the kGain score for each variant (kGain_v) by summing the natural logarithm of the fold changes for each k-mer across all windows.

The kGain score is mathematically represented as:

![kgainformula](https://github.com/user-attachments/assets/fe11a306-95f2-49c7-a2ba-d42e7ec9ceae)

Where:
- ![gain](https://github.com/user-attachments/assets/db2b8ed7-7d4d-4d84-bf54-e07046144576) is the score for variant (v).
  
- ![alt](https://github.com/user-attachments/assets/1564acf3-1c16-4d47-be39-f578865f997f) is the frequency of the k-mer containing the alternate allele in the \(i\)th window.

- ![ref](https://github.com/user-attachments/assets/c91fc951-0810-498a-ab9e-60d1d839c124) is the frequency of the k-mer containing the reference allele in the \(i\)th window.

- _**k**_ is the total number of k-mers generated using the rolling window method.


