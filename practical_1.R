# ---------- DOWNLOAD THE PRACTICAL ----------
#
# 1. Download the .zip file for the practical from
#
#    https://github.com/tamuri/BIOL0033/archive/master.zip
#
# 2. Extract the files

# ---------- SET THE WORKING DIRECTORY ----------
# i.e. the place where you put the files for the practical

setwd('~/Documents/2019/BIOL0033')

# ---------- LOAD THE REQUIRED PACKAGES ----------

# 'ape' and 'phangorn' provide vast array of phylogenetic analyses
library(ape)
library(phangorn)

# 'msa' provides multiple sequence alignment
library(msa)

# ---------- BASIC PHYLOGENY & PLOTTING -----------

# Trees are commonly defined using the Newick format
my_newick <- '(a, (b, c), d);'

# Plot the tree defined by the string above
my_tree <- ape::read.tree(text=my_newick)
print(my_tree)
ape::plot.phylo(x=my_tree, type='unrooted')

# Before we can make any evolutionary inference, we need to root the tree using an outgroup
rerooted_tree <- ape::root(phy=my_tree, outgroup = 'a', resolve.root = TRUE)
ape::plot.phylo(x=rerooted_tree, type='phylogram')

# EXERCISE 1:
#
# 1. Plot the following trees:
#
#   i) (frog, human, bird);
#  ii) (gorilla:1, (human:2, chimp:3):4);
#
# 2. What is the Newick string for the following tree. (Plot it to make sure)
#
#  /-----+ sus_scrofa                                                   
#  |                                                                     
# =+      /-----+ tursiops_truncatus                             
#  |      |                                                      
#  \------+           /-----+ capra_hircus     
#         |     /-----+                        
#         \-----+     \-----+ ovis_aries       
#               |                                       
#               \-----+ bos_taurus                      
#
# 3. Write the Newick strings for each of the 3 unrooted trees for 4 species (plot to check)
#
# 4. In the 'my_newick' example tree above, what is taxon A most closely related to 
#    if the unrooted tree is rooted by
#    i) taxon a
#   ii) taxon b
#  iii) taxon c
#   iv) taxon d

# Read in a tree from a file (take a look at the file too)
mammals_72sp <- ape::read.tree(file='72sp.tree')

# Get summary of tree
mammals_72sp

# Get tips
mammals_72sp$tip.label

# Get the tree length (sum of branch lengths)
sum(mammals_72sp$edge.length)

# Explore different ways of plotting the tree (one at a time)
ape::plot.phylo(x=mammals_72sp, type='cladogram', use.edge.length=FALSE, cex=0.3)

ape::plot.phylo(x=mammals_72sp, type='phylogram')
add.scale.bar()

# Q: What's the difference between a phylogram and a cladogram?

ape::plot.phylo(x=mammals_72sp, type='unrooted')

ape::plot.phylo(x=mammals_72sp, type='fan')

ape::plot.phylo(x=mammals_72sp, type='radial')

# you can pass most R plot arguments e.g. adjust size of text using 'cex'
ape::plot.phylo(x=mammals_72sp, type='phylogram', align.tip.label=TRUE, cex=0.7)
add.scale.bar()

# Add the node numbers for internal nodes
nodelabels()

# Plot with specific node highlighted
ape::plot.phylo(x=mammals_72sp, type='phylogram', use.edge.length=FALSE)
nodelabels('Carnivora', 132)
nodelabels('Glires', 79)


# ---------- MULTIPLE SEQUENCE ALIGNMENT ----------
# Align sequences using the 'msa' package

# Load the unaligned sequences
seqs <- Biostrings::readDNAStringSet(filepath='nadh6.8apes.fasta', format='fasta')

# The msa package provides several alignment tools. By default it uses the 
# program 'ClustalW'. To align the sequences using the default options
aln <- msa::msa(inputSeqs=seqs)

# Take a look at the alignment
aln

# Look at the entire alignment, not truncated
print(aln, show='complete')

# Save the multiple sequence alignment
Biostrings::writeXStringSet(Biostrings::unmasked(x=aln), file='nadh6.8apes.aln.fasta')

# EXERCISE 2:
#
# Try different options. For each variation, save the MSA to a different file and compare in AliView
#
# 1. Try a different algorithm: ClustalOmega or MUSCLE
#    e.g. aln <- msa(seqs, 'Muscle')
#
# 2. Try varying the gapOpening and gapExtension values
#    e.g. msa(seqs, gapOpening=0, gapExtension=0)
#    e.g. msa(seqs, gapOpening=1000, gapExtension=10000)
#    What's happening here? How do the resulting MSAs differ?

  
# ---------- CALCULATE PAIRWISE DISTANCES ----------

# The 'ape' library provides many tools for phylogenetic analysis. 

# Read the aligned sequences we had saved to file above
aln <- ape::read.dna(file='nadh6.8apes.aln.fasta', format='fasta')

# Lets see what the alignment looks like
ape::checkAlignment(x=aln)

# What kind of things should we look for to check the quality of our alignment?

# Calculate pairwise distances between sequences
D <- ape::dist.dna(x=aln, model='JC')

# Look at the distance matrix
D

# EXERCISE 3:
#
# 1. Try a few other substitutions models to calculate the pairwise distances
#   e.g. D2 <- dist.dna(aln, model='F84')
#   how do the distances compare under different models?
#
# 2. The Gamma model of among site rate variation (ASRV)

  # Plot a Gamma distribution for high-degree of site rate variation,
  # e.g. alpha parameter is 0.5 
  alpha <- 0.5
  curve(dgamma(x, shape=alpha, rate=alpha), from=0, to=10)
  
  # Plot a Gamma distribution for very little site rate variation
  alpha <- 100
  curve(dgamma(x, shape=alpha, rate=alpha), from=0, to=2)

# Calculate the pairwise distances using a model with ASRV (supply the alpha parameter)
# e.g. dist.dna(aln, model='F84', gamma=0.5)
# How do the distances compare for a model with and without ASRV?


# ---------- CREATE A UPGMA TREE USING DISTANCE MATRIX ----------

upgma_tree <- phangorn::upgma(D=D)
upgma_tree <- ape::ladderize(phy=upgma_tree)
plot(upgma_tree)

# Q: What's the distinctive feature of a UPGMA tree? (look at the plot). Is it suitable for our data set?

# look at the branch lengths
upgma_tree$edge.length

# look at the newick representation
ape::write.tree(phy=upgma_tree)

# save the tree to file
ape::write.tree(phy=upgma_tree, file='nadh6_upgma.tree')

# Q: Open the tree in Figtree (if you installed it)? Does it render the tree the same way?

# ---------- CREATE A NEIGHBOUR-JOINING TREE USING DISTANCE MATRIX ----------

# state the outgroup for the sequences
outg <- c('pig', 'cow')

# Use the distances we calculated above
nj_tree <- ape::nj(X=D)

# View the tree
ape::plot.phylo(x=nj_tree)
add.scale.bar()

# By default the tree is unrooted. Even if it looks OK, _always_ root the tree
rooted_nj <- ape::ladderize(phy=ape::root(phy=nj_tree, outgroup=outg, resolve.root=TRUE))
ape::plot.phylo(x=rooted_nj)
add.scale.bar()

# To save a tree
ape::write.tree(phy=nj_tree, file = 'nadh6_nj.tree')

# EXERCISE 4:
#
# 1. Use a distance matrix from a different model (from exercise 3) to construct
#    a NJ tree. How does the tree differ?



