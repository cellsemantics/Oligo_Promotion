#===============================================================================
#                      Organised Code 
#===============================================================================
# Step 1 - Load library 
#Load library 
library(Biostrings)
library(dplyr)

#step 2 - Read fasta files & preprocess 
# File path to the FASTA file
file_path <- "G:/My Drive/PhD Project/Mutational pseudotime_from local_main/Chaos Game/E.Coli_LTEE/Data/GCF_000017985.1_ASM1798v1_genomic.fna"
# Read the file as a character vector
fasta_lines <- readLines(file_path)
# Remove the first line (header line)
sequence_lines <- fasta_lines[-1]
# Concatenate the sequence lines into a single string
sequence <- paste(sequence_lines, collapse = "")
# Print the resulting sequence
print(sequence)

#Step 3 - Read VCF files & preprocess it 
#Read VCF files & filter it 
# Read VCF 
data=read_excel("G:/My Drive/PhD Project/Mutational pseudotime_from local_main/Chaos Game/E.Coli_LTEE/Data/SNP file upto 20K genome.xlsx")
View(data)
# Filter rows based on SNP.indel
#SNP_data <- data[data$SNP.Indel == "S", ]
#write.csv(SNP_data, file = "G:/My Drive/PhD Project/Mutational pseudotime_from local_main/Chaos Game/MegaPlate_paper_seq/SNP_data.csv", row.names = FALSE)

#Step 4 - Index the position of genes 
#Index the gene position as list 
snp_position <- data$Position

#Step 5 - process to obtain ref seq
# Specify the length of the flanking region
flanking_length <- 10
# Specify the length of the k-mers
k <- ((flanking_length) * 2) + 1
# For ref 
#Index the gene position as list 
snp_position <- data$`Genome position`
# Convert snp_position to numeric
snp_position <- as.numeric(snp_position)

# Create an empty vector to store the flanking regions
flanking_regions <- vector("character", length = length(snp_position))
# Convert flanking_length to numeric
flanking_length <- as.numeric(flanking_length)

# Iterate over SNP positions and extract the flanking regions
for (i in seq_along(snp_position)) {
  # Extract the flanking region from the reference sequence
  start <- snp_position[i] - flanking_length
  end <- snp_position[i] + flanking_length
  flanking_region <- substring(sequence, start, end)
  
  # Store the flanking region in the vector
  flanking_regions[i] <- flanking_region
}

# Print the flanking regions
print(flanking_regions)
ref_seq<-as.data.frame(flanking_regions)
#Save in local
write.csv(ref_seq, file = "G:/My Drive/PhD Project/Mutational pseudotime_from local_main/Chaos Game/E.Coli_LTEE/Data/ref.csv", row.names = FALSE)

#Step 5 - process to obtain alt seq
# Alt 
# Index the gene position as list 
snp_position <- data$`Genome position`
# Create an empty vector to store modified regions
alternate_regions <- vector("character", length = length(snp_position))
# Convert flanking_length to numeric
alternate_regions <- as.numeric(alternate_regions)
# Iterate over SNP positions and substitute the alternate allele
for (i in seq_along(snp_position)) {
  # Extract the flanking region from the reference sequence
  start <- snp_position[i] - flanking_length
  end <- snp_position[i] + flanking_length
  flanking_region <- substring(sequence, start, end)
  
  # Get the alternate allele sequence for the current SNP
  alternate_seq <- data$`Evolved nucleotide`[i]
  
  # Substitute the alternate allele in the flanking region
  modified_region <- paste0(substr(flanking_region, 1, flanking_length), alternate_seq, substr(flanking_region, flanking_length + 2, 2 * flanking_length + 1))
  
  # Store the modified region in the vector
  alternate_regions[i] <- modified_region
}

# Print the modified regions
print(alternate_regions)
alt_seq<-as.data.frame(alternate_regions)
#save in local
write.csv(ref_seq, file = "G:/My Drive/PhD Project/Mutational pseudotime_from local_main/Chaos Game/E.Coli_LTEE/Data/alt.csv", row.names = FALSE)

#===============================================================================
#               To verify if the loop ran well 
#===============================================================================
# Specify the SNP position and alternate allele sequence
snp_position <- 161041
# Specify the length of the flanking region
flanking_length <- 10

# Specify the length of the k-mers
k <- ((flanking_length) * 2) + 1
# Extract the flanking region from the reference sequence
start <- snp_position - flanking_length
print (start)
end <- snp_position + flanking_length
print (end)
flanking_region <- substring(sequence, start, end)
print (flanking_region)
# Extract a substring
substring <- substr(text_file, start = 161031, stop = 161051)
