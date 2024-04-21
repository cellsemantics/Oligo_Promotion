#===============================================================================
# Actual Data from python analysis continued here 
#===============================================================================
# Just add the WT and mutated DNA here.

# Step 1
# Read fasta file 
# Load the Biostrings package
library(Biostrings)

# Specify the path to your FASTA file
fasta_file <- "/home/smrutip/smruti/smruti_project/gisaid_project/reference_covid_sequence.fasta"

# Read the FASTA file
fasta_sequences <- readDNAStringSet(fasta_file)

# Print the first few sequences to check
head(fasta_sequences)

# Step 2 
# To get the WT DNA seq from start and end retrieved after pre-processing from python script 
# loop this for all rows 
# Read the df 
df = read_excel("/home/smrutip/smruti/smruti_project/gisaid_project/E_gene_data.xlsx")

# Assuming 'start' and 'end' are columns in your data frame 'df'
df$Start <- df$Start
start_positions <- df$Start
df$End <- as.integer(df$End)
end_positions <- df$End

# Initialize an empty list to store the flanks
flanks <- list()

# Loop through each row in the data frame
for (i in 1:nrow(df)) {
  start_position <- start_positions[i]
  end_position <- end_positions[i]
  
  # Extract the flank
  flank <- subseq(fasta_sequences, start_position, end_position)
  flank <- as.character(flank)
  
  # Append the flank to the list
  flanks[[i]] <- flank
}

# Add the 'flanks' list as a new column in the data frame
df$flank <- unlist(flanks)

# Now, 'df' will have a new column 'flank' containing the extracted WT DNA seq 

# Step 3 
# Create Mut DNA from WT DNA (Flanks)
mut <- character()

# Iterate over rows of the data frame 'pos'
for (i in 1:nrow(df)) {
  start <- df$Start[i]
  end <- df$End[i]
  posi <- df$POS[i]
  alt <- df$ALT[i]
  seq <- df$flank[i]
  
  # Check for missing values in posi, start, or seq
  if (is.na(posi) || is.na(start) || is.na(seq)) {
    warning("Skipping row ", i, " due to missing values")
    next  # Skip to the next iteration
  }
  
  # Initialize an empty string for the mutated sequence
  Sequence_Mut <- ""
  
  # Iterate over the characters in the wild-type sequence
  for (k in 1:nchar(seq)) {
    if (k == posi - start + 1) {  # Check if we've reached the mutation position
      Sequence_Mut <- paste0(Sequence_Mut, alt)
    } else {
      Sequence_Mut <- paste0(Sequence_Mut, substr(seq, k, k))
    }
  }
  
  # Append the mutated sequence to the 'mut' vector
  mut <- c(mut, Sequence_Mut)
}

# Print the 'mut' vector
print(mut)

# Add the mutated sequences to the 'pos' data frame
df$mut_flank <- mut

# Step 4 
# check if WT and Mut has the difference as per in alt and position 
# Difference between the WT and Mut DNA seq and match it with the original cols that did they change rightly 

# Create a function that calculates character differences
character_difference <- function(str1, str2, start) {
  diff_count <- sum(str1 != str2)
  temp_pos <- 0  # Start indexing from 1
  diff <- character()
  diff_pos <- 0
  
  for (i in 1:nchar(str1)) {
    char1 <- substr(str1, i, i)
    char2 <- substr(str2, i, i)
    
    if (char1 != char2) {
      diff <- c(diff, paste0(char1, char2))
      diff_pos <- temp_pos + start
    }
    
    temp_pos <- temp_pos + 1
  }
  
  return(list(diff_count, diff, diff_pos))
}

# Apply the character_difference function to all rows of the data frame 'pos'
# Apply the character_difference function to all rows of the data frame 'pos'
results <- apply(df, 1, function(row) {
  character_difference(row['flank'], row['mut_flank'], as.numeric(row['Start']))
})


# Create new columns in 'pos' to store the results
df$Diff_Count_DNA <- sapply(results, function(result) result[[1]])
df$Differences_DNA <- sapply(results, function(result) paste(result[[2]], collapse = ', '))
df$Diff_Pos_DNA <- sapply(results, function(result) result[[3]])

# Write in local
library(openxlsx)

# Specify the file path and name for the Excel file
file_path <- "/home/smrutip/smruti/smruti_project/gisaid_project/df_orf1ab_10_2.xlsx"

# Write the data frame to Excel
write.xlsx(df, file_path, rowNames = FALSE)

# Take this sheet to python to translate into amino acid (protein sequence)
#==============================================================================
# Now got the WT prt and Mut prt seq from Python
#===============================================================================
# Now we have to just know which position changed to which

library(readxl)

# converted the DNA of WT and mut to prt now find the difference 
df = read_excel("/home/smrutip/smruti/smruti_project/gisaid_project/orf1ab.xlsx")

# Difference between the WT and Mut DNA seq and match it with the original cols that did they change rightly 
# Create a function that calculates character differences
character_difference <- function(str1, str2, start) {
  diff_count <- sum(str1 != str2)
  temp_pos <- 0  # Start indexing from 1
  diff <- character()
  diff_pos <- 0
  
  for (i in 1:nchar(str1)) {
    char1 <- substr(str1, i, i)
    char2 <- substr(str2, i, i)
    
    if (char1 != char2) {
      diff <- c(diff, paste0(char1, char2))
      diff_pos <- temp_pos + start
    }
    
    temp_pos <- temp_pos + 1
  }
  
  return(list(diff_count, diff, diff_pos))
}

# Apply the character_difference function to all rows of the data frame 'pos'
# Apply the character_difference function to all rows of the data frame 'pos'
results <- apply(df, 1, function(row) {
  character_difference(row['WT_aa'], row['Mut_aa'], as.numeric(row['Start']))
})

# Create new columns in 'pos' to store the results
df$Diff_Count_aa <- sapply(results, function(result) result[[1]])
df$Differences_aa <- sapply(results, function(result) paste(result[[2]], collapse = ', '))
df$Diff_Pos_aa <- sapply(results, function(result) result[[3]])

# Get the SNP from wildtype and mutated protein 
# Filter out rows with missing values in either column
df <- df[complete.cases(df$WT_aa, df$Mut_aa), ]

# Function to compare amino acid sequences and generate SNP column
compare_amino_acid_sequences <- function(WT_aa, Mut_aa) {
  
  # Initialize an empty vector to store SNP changes
  snp_changes <- character()
  
  # Iterate through each position and compare amino acids
  for (i in 1:nchar(WT_aa)) {
    if (substr(WT_aa, i, i) != substr(Mut_aa, i, i)) {
      # Amino acid change detected, add to SNP changes
      snp_changes <- c(snp_changes, paste0(substr(WT_aa, i, i), i, substr(Mut_aa, i, i)))
    }
  }
  
  # Return "no_change" if there are no differences, otherwise return the SNP changes
  if (length(snp_changes) == 0) {
    return("no_change")
  } else {
    return(paste(snp_changes, collapse = ","))
  }
}

# Apply the function to create the SNP column
df$SNP <- mapply(compare_amino_acid_sequences, df$WT_aa, df$Mut_aa)

# Specify the file path and name for the Excel file
file_path <- "/home/smrutip/smruti/smruti_project/gisaid_project/orf1ab_4.xlsx"

# Write the data frame to Excel
write.xlsx(df, file_path, rowNames = FALSE)

# Take this sheet to python for input preparation for ESM1b model
#==============================================================================
