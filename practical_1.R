# ---------- SET THE WORKING DIRECTORY ----------
# i.e. the place where you put the files for the practical

setwd("~/Documents/2019/biol0030/BIOL3010-ZY-resources/Practical")


# ---------- MULTIPLE SEQUENCE ALIGNMENT ----------
# Align sequences using the 'msa' package

# Load the required package
library(msa)

# Load the unaligned sequences
sequences <- readDNAStringSet('NADH6.fas', format='fasta')
# sequences <- readAAStringSet('~/Documents/2018/tree_annotator/etc/examples/flu.pb2.aminoacids/PB2.fasta', format='fasta')

# The msa package provides several alignment tools. By default it uses the 
# program "ClustalW". To align the sequences using the default options
aligned_sequences <- msa(sequences)

# Take a look at the alignment
aligned_sequences

# Look at the entire alignment, not truncated
print(aligned_sequences, show='complete')

# More:
# 
# * Try two other algorithms: ClustalOmega and MUSCLE
#  `nadh6_alignment_muscle <- msa(nadh6_sequences, "Muscle")`
# * Options for the Muscle aligner: gapOpen, gapExtend


# ---------- VISUALISE TREES ----------
# https://4va.github.io/biodatasci/r-ggtree.html
  
  
# ---------- CALCULATE PAIRWISE DISTANCES ----------

# The 'ape' library provides many tools for phylogenetic analysis. 
  
# Load the package
library(ape)

# Convert our multiple sequence alignment into the format required by ape:

# TODO: this doesn't work
aligned_sequences <- msaConvert(aligned_sequences, type='ape::DNAbin')

# Save the 'msa' file and read it back using 'ape'
writeXStringSet(unmasked(aligned_sequences), file='test.fasta')
aligned_sequences <- read.FASTA('test.fasta', type='DNA')

# Calculate pairwise distances between sequences
distances <- dist.dna(aligned_sequences)

# Look at the distance matrix
distances

# Visualise the distances (requires loading the 'ade4' library)
library(ade4)


# ---------- CREATE A NEIGHBOUR-JOINING TREE USING DISTANCES ----------

# Use the distances we calculated above
nj_tree <- nj(distances)

# View the tree
plot(nj_tree)

# By default the tree is unrooted. Even if it looks OK, _always_ explicitly root the tree
rooted_nj_tree <- root(nj_tree, outgroup=c('pig', 'cow'))

# Plot the rooted tree. how does it look different?

# ---------- ESTIMATE TREE BY MAXIMUM LIKELIHOOD ----------
library(phangorn)
dna2 <- as.phyDat(aligned_sequences)
tre.ini <- nj(dist.dna(aligned_sequences, model = "TN93"))
fit.ini <- pml(tre.ini, dna2, k = 4)
fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE)
ml.tree <- fit$tree
ml.tree <- root(ml.tree, outgroup=c('pig', 'cow'))
ml.tree <- ladderize(ml.tree)
plot(ml.tree)

install.packages(c('ape', 'adegenet', 'phangorn', 
