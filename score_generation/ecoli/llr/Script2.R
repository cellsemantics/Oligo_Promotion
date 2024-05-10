#===============================================================================
# Actual Data from python analysis continued here 
#===============================================================================
# Just change the WT to mutant DNA here.
library(readxl)
library(Biostrings)
#Read the sheet 
df <- read_excel('/home/bernadettem/bernadettenotebook/E.Coli_LTEE/Prt_final/esm1b_redo/Correspondance/Results/DNA_Prt_Strand_variant.xlsx')
# Separate the - strand and + strand for change to WT to mutant 
pos = df[df$Strand == "+",]
table(pos$Strand)
neg = df[df$Strand == "-",]
table(neg$Strand)

# For pos strand 
# Convert the WT DNA seq to mutant DNA seq 
mut <- character()

# Iterate over rows of the data frame 'pos'
for (i in 1:nrow(pos)) {
  start <- pos$Start[i]
  end <- pos$End[i]
  posi <- pos$Position[i]
  alt <- pos$Alt_allele[i]
  seq <- pos$DNA_Seq[i]
  
  # Initialize an empty string for the mutated sequence
  Sequence_Mut <- ""
  
  # Iterate over the characters in the wild-type sequence
  for (k in 1:nchar(seq)) {
    if (k == posi - start + 1 ) {  # Check if we've reached the mutation position
      Sequence_Mut <- paste0(Sequence_Mut, alt)
    } else {
      Sequence_Mut <- paste0(Sequence_Mut, substr(seq, k, k))
    }
  }
  
  # Append the mutated sequence to the 'mut' vector
  mut <- c(mut, Sequence_Mut)
}

# Add the mutated sequences to the 'pos' data frame
pos$Sequence_Mut <- mut
# After conversion to mut DNA , we should know the difference between WT and Mut DNA and match. 

### To know the difference 
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
results <- apply(pos, 1, function(row) {
  character_difference(row['DNA_Seq'], row['Sequence_Mut'], as.numeric(row['Start']))
})


# Create new columns in 'pos' to store the results
pos$Diff_Count <- sapply(results, function(result) result[[1]])
pos$Differences <- sapply(results, function(result) paste(result[[2]], collapse = ', '))
pos$Diff_Pos <- sapply(results, function(result) result[[3]])

# Print the updated data frame 'pos'
print(pos)

# Check if there is a space in each element of 'col2'
pos$has_space <- grepl(" ", pos$DNA_Seq) # check if space in WT seq

# Check 
table(pos$Annotation)

table(pos$Diff_Count)
#==============================================================================
# For the negative strand the given sequence is reverse complementary therefore complement the alternate allele and replace the position in mut seq with the complementary alternative allele
# Convert the alt seq 
# Function to complement a single character
complement_char <- function(char) {
  if (char == "A") {
    return("T")
  } else if (char == "T") {
    return("A")
  } else if (char == "G") {
    return("C")
  } else if (char == "C") {
    return("G")
  } else {
    return(char)
  }
}

# Complement alt allele using the function
neg$Alt_allele <- sapply(strsplit(neg$Alt_allele, NULL), function(chars) {
  paste0(sapply(chars, complement_char), collapse = "")
})

#For negative strand 
# Convert the WT DNA seq to mutant DNA seq 
mut <- character()

# Iterate over rows of the data frame 'neg'
for (i in 1:nrow(neg)) {
  start <- neg$End[i]
  end <- neg$Start[i]
  pos <- neg$Position[i]
  alt <- neg$Alt_allele[i]
  seq <- neg$DNA_Seq[i]
  
  # Initialize an empty string for the mutated sequence
  Sequence_Mut <- ""
  
  # Iterate over the characters in the wild-type sequence
  for (k in 1:nchar(seq)) {
    if (k == start - pos + 1 ) {  # Check if we've reached the mutation position
      Sequence_Mut <- paste0(Sequence_Mut, alt)
    } else {
      Sequence_Mut <- paste0(Sequence_Mut, substr(seq, k, k))
    }
  }
  
  # Append the mutated sequence to the 'mut' vector
  mut <- c(mut, Sequence_Mut)
}

# Add the mutated sequences to the 'neg' data frame
neg$Sequence_Mut <- mut



### To know the difference 
# Difference between the WT and Mut DNA seq and match it with the oriinal cols that did they change rightly 

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

# Apply the character_difference function to all rows of the data frame 'neg'
# Apply the character_difference function to all rows of the data frame 'neg'
results <- apply(neg, 1, function(row) {
  character_difference(row['DNA_Seq'], row['Sequence_Mut'], as.numeric(row['Start']))
})


# Create new columns in 'neg' to store the results
neg$Diff_Count <- sapply(results, function(result) result[[1]])
neg$Differences <- sapply(results, function(result) paste(result[[2]], collapse = ', '))
neg$Diff_Pos <- sapply(results, function(result) result[[3]])

# Print the updated data frame 'neg'
print(neg)

# Check if there is a space in each element of 'col2'
neg$has_space <- grepl(" ", neg$DNA_Seq)

# Check 
table(neg$Annotation)

table(neg$Diff_Count)

# Concatenate df2 below df1
concatenated_df <- rbind(pos, neg)
# Specify the file path and name for the Excel file
file_path <- "~/bernadettenotebook/E.Coli_LTEE/Prt_final/esm1b_redo/Correspondance/Results/DNA_Prt_Strand_Variant_MutDNA.xlsx"
# Write the data frame to Excel
write.xlsx(concatenated_df, file_path, rowNames = FALSE)
# Take this sheet to python to translate into amino acid 
#==============================================================================
# Now got the WT prt and Mut prt seq from Python
#===============================================================================
# Now we have to just know which position changed to which 
# read the df 
# data = read_excel("~/bernadettenotebook/E.Coli_LTEE/Prt_final/esm1b_redo/DNA_PRT_Ecolimeta_true_fresh.xlsx")
data = read_excel("~/bernadettenotebook/E.Coli_LTEE/Prt_final/esm1b_redo/Correspondance/Results/DNA_Prt_Strand_Variant_MutDNA_True.xlsx")
# Code given by nancy for this 
# Filter out rows with missing values in either column
data <- data[complete.cases(data$Prt_Seq, data$Mut_aa_made), ]
# Function to compare amino acid sequences and generate SNP column
compare_amino_acid_sequences <- function(WT_aa_made, Mut_aa_made) {
  # Initialize an empty vector to store SNP changes
  snp_changes <- character()
  
  # Iterate through each position and compare amino acids
  for (i in 1:nchar(WT_aa_made)) {
    if (substr(WT_aa_made, i, i) != substr(Mut_aa_made, i, i)) {
      # Amino acid change detected, add to SNP changes
      snp_changes <- c(snp_changes, paste0(substr(WT_aa_made, i, i), i, substr(Mut_aa_made, i, i)))
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
data$SNP <- mapply(compare_amino_acid_sequences, data$Prt_Seq, data$Mut_aa_made)
table(data$SNP)
# Write in local 
# Specify the file path and name for the Excel file
file_path <- "~/bernadettenotebook/E.Coli_LTEE/Prt_final/esm1b_redo/Correspondance/Results/DNA_Prt_Strand_Variant_MutDNA_True_SNP.xlsx"

# Write the data frame to Excel
write.xlsx(data, file_path, rowNames = FALSE)
