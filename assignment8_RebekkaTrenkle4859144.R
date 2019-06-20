# import package
library(Biostrings) 

# create strings for the sequences to be aligned
s1 <- "GAATC" 
s2 <- "CATAC" 

# create the nucleotide substitution matrix called sigma. This can only be symmetric
sigma <- nucleotideSubstitutionMatrix(match = 10, mismatch = -5, baseOnly = TRUE) 
sigma 
# change entries in sigma to zero, that it matches the one from the lecture
sigma[1,3] <- 0
sigma[2,4] <- 0
sigma[3,1] <- 0
sigma[4,2] <- 0
sigma

# perform the alignment using the function from the package Biostrings
# https://www.rdocumentation.org/packages/Biostrings/versions/2.40.2/topics/pairwiseAlignment
# change the values for gapOpening and gapExtension to -4 to match the lecture
alignment <- pairwiseAlignment(s1, s2, substitutionMatrix = sigma, gapOpening = -4, gapExtension = -4, scoreOnly = FALSE) 
alignment # shows the aligned sequences


##########
# install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")
library(DECIPHER)

# find example sequence file from package
# https://www.rdocumentation.org/packages/base/versions/3.6.0/topics/system.file
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
# load dna sequences (not aligned), store in the variable "dna"
dna <- readDNAStringSet(fas)
dna


# as the example sequences are protein coding genes, it makes sense to align the 
# translated sequences, meaning the amino acid sequence, instead of aligning the 
# nucleotide sequences. The aligned sequences then get reverse translated
DNA <- AlignTranslation(dna)
BrowseSeqs(DNA, highlight=1) # view the alignment
# write aligned sequence in fasta file
writeXStringSet(DNA, file="multipleSequenceAlignment.fasta")








