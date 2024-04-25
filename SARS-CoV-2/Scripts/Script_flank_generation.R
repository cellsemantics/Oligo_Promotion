#Step 1 ----->  Take all merged tav from the fasta directories and make 1 tsv for 1 fasta file.
# Load required libraries
# library(data.table)
# 
# # Set your directory path where the folders containing merged.tsv files are located
# parent_folder <- "/home/smrutip/smruti/smruti_project/gisaid_project/strainflow_manuscript/data/Segmented/1660645735959.sequences.fasta/"
# 
# # List all subdirectories in the parent folder
# subdirs <- list.dirs(parent_folder, recursive = FALSE)
# 
# # Initialize an empty list to store dataframes
# list_of_data <- list()
# 
# # Iterate through each subdirectory
# for (subdir in subdirs) {
#   # List all files named merged.tsv in the current subdirectory
#   file_list <- list.files(path = subdir, pattern = "merged.tsv", full.names = TRUE)
#   
#   # Read and merge files into a single dataframe for the current subdirectory
#   if (length(file_list) > 0) {
#     data <- rbindlist(lapply(file_list, fread))
#     list_of_data[[subdir]] <- data
#   }
# }

# # Combine all dataframes into a single dataframe
# merged_data <- rbindlist(list_of_data)

#Step 2 ----->  Filter the data as per only SNPs 
df <- read_excel('~/smruti/smruti_project/gisaid_project/Oligo_analysis /all_genes_esm_scores_covid.xlsx')
filtered_data_split <- df
# # Assuming 'merged_data' is your dataframe, and you want to filter 
# filtered_data <- df[nchar(as.character(REF)) <= 1 & nchar(as.character(ALT)) <= 1]

# Step 3 ---> Read the WT ref file 
# step 3.1 - Load library
#Load library
library(Biostrings)
library(dplyr)
library(readxl)
library(writexl)

# step 3.2 - Read fasta files & preprocess
# Specify the path to your FASTA file
fasta_file <- "/home/smrutip/smruti/smruti_project/gisaid_project/reference_covid_sequence.fasta"

# Read the FASTA file
fasta_sequences <- readDNAStringSet(fasta_file)

# Print the first few sequences to check
head(fasta_sequences)

# Step 4 ----> Ref Flank 
# Split Col1 
# Load required libraries
library(tidyr)
# filtered_data$`POS_Gene POS`
# # Assuming 'merged_data' is your dataframe and 'col2' is the column to split
# filtered_data_split <- separate(filtered_data, col = `POS_Gene POS`, into = c("POS_Gene", "POS"), sep = " ")
# as.numeric(filtered_data_split$POS)
# write_xlsx(filtered_data_split, path = "~/smruti/smruti_project/gisaid_project/Oligo_analysis /Data_filtered_1660645735959.xlsx")
# # Specify the length of the flanking region
flanking_length <- 10
# For ref
#Index the gene position as list
snp_position <- filtered_data_split$POS
# Create an empty vector to store the flanking regions
flanking_regions <- vector("character", length = length(snp_position))
# snp_position
# Convert snp_position to numeric
snp_position <- as.numeric(snp_position)
# Convert flanking_length to numeric
flanking_length <- as.numeric(flanking_length)
snp_position
start <-filtered_data_split$Start
end <-filtered_data_split$End

# Iterate over SNP positions and extract the flanking regions
for (i in seq_along(snp_position)) {
  # Extract the flanking region from the reference sequence
  start <- snp_position[i] - flanking_length
  end <- snp_position[i] + flanking_length
  flanking_region <- substring(fasta_sequences, start, end)
  # Store the flanking region in the vector
  flanking_regions[i] <- flanking_region
}
# Print the flanking regions
print(flanking_regions)
ref_seq<-as.data.frame(flanking_regions)
View(ref_seq)

# Step 5 ---> Generate alternative seq 

# Create an empty vector to store modified regions
alternate_regions <- vector("character", length = length(snp_position))
# Convert flanking_length to numeric
alternate_regions <- as.numeric(alternate_regions)
# Iterate over SNP positions and substitute the alternate allele
for (i in seq_along(snp_position)) {
  # Extract the flanking region from the reference sequence
  start <- snp_position[i] - flanking_length
  end <- snp_position[i] + flanking_length
  flanking_region <- substring(fasta_sequences, start, end)
  # Get the alternate allele sequence for the current SNP
  alternate_seq <- filtered_data_split$ALT[i]
  # Substitute the alternate allele in the flanking region
  modified_region <- paste0(substr(flanking_region, 1, flanking_length), alternate_seq, substr(flanking_region, flanking_length + 2, 2 * flanking_length + 1))
  # Store the modified region in the vector
  alternate_regions[i] <- modified_region
}
alt_seq<-as.data.frame(alternate_regions)
View(alt_seq)

# Step 5 ---> write all in local 
# Write the data frame to an Excel file
write_xlsx(ref_seq, path = "~/smruti/smruti_project/gisaid_project/Oligo_analysis /ref_flank.xlsx")
write_xlsx(alt_seq, path = "~/smruti/smruti_project/gisaid_project/Oligo_analysis /alt_flank.xlsx")



