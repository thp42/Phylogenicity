# Load the required namespaces and packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DECIPHER")
BiocManager::install("BiocParallel")

if (!requireNamespace("Biostrings", quietly = TRUE))
    BiocManager::install("Biostrings")

library(DECIPHER)
library(Biostrings)
library(BiocParallel)
library(parallel)

# Register a parallel backend
register(MulticoreParam(workers = detectCores()))

# Function to modify headers
modify_headers <- function(seqs) {
  new_names <- sapply(names(seqs), function(name) {
    # Replace spaces and special characters with underscores
    new_name <- gsub("[^A-Za-z0-9]", "_", name)
    return(new_name)
  })
  names(seqs) <- new_names
  return(seqs)
}

# Load your amino acid sequences
fasta_file <- "/home/pythagoras/Documents/PhD/Evolution/FL/CGN/CGN.fasta"
seqs <- readAAStringSet(fasta_file)

# Modify headers
seqs <- modify_headers(seqs)

# Align the sequences using MAFFT with parallel processing
aligned_seqs <- AlignSeqs(seqs, iterations = 100, refinements = 200, processors = detectCores())

# Save the aligned sequences to the same directory
output_file <- file.path(dirname(fasta_file), "CGN_aligned_protein_sequences.fasta")
writeXStringSet(aligned_seqs, filepath=output_file)

# Visualize the alignment in the browser
BrowseSeqs(aligned_seqs)
