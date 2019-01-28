# ---------- SET THE WORKING DIRECTORY ----------
# i.e. the place where you put the files for the practical

setwd("~/Documents/2019/BIOL0033")

# ---------- LOAD THE REQUIRED PACKAGES ----------

# 'ape' and 'phangorn' provide vast array of phylogenetic analyses
library(ape)
library(phangorn)

# 'msa' provides multiple sequence alignment
library(msa)

# ---------- BASIC PHYLOGENY & PLOTTING -----------

# Trees are commonly using the Newick format
my_newick <- '(a, (b, c), d);'

# Plot the tree defined by the string above
my_tree <- read.tree(text=my_newick)
print(my_tree)
plot(my_tree, type='unrooted')

# Before we can make any evolutionary inference, we need to root the tree using an outgroup
rerooted_tree <- root(my_tree, outgroup = 'a', resolve.root = TRUE)
plot(rerooted_tree, type='phylogram')

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
# 4. In the "my_newick" example tree above, what is taxon A most closely related to 
#    if the unrooted tree is rooted by
#    i) taxon a
#   ii) taxon b
#  iii) taxon c
#   iv) taxon d

# Read in a tree from a file (take a look at the file too)
mammals_72sp <- read.tree(file='72sp.tree')

# Get summary of tree
mammals_72sp

# Get tips
mammals_72sp$tip.label

# Get the tree length (sum of branch lengths)
sum(mammals_72sp$edges)

# Explore different ways of plotting the tree (one at a time)
plot(mammals_72sp, type='cladogram')

plot(mammals_72sp, type='cladogram', use.edge.length = FALSE)

plot(mammals_72sp, type='unrooted')

plot(mammals_72sp, type='fan')

plot(mammals_72sp, type='radial')

plot(mammals_72sp, type='phylogram')
add.scale.bar()

plot(mammals_72sp, type='phylogram', align.tip.label=TRUE)
add.scale.bar()

# Add the node numbers for internal nodes
nodelabels()

# Plot with specific node highlighted
plot(mammals_72sp, type='phylogram', use.edge.length=FALSE)
nodelabels("Carnivora", 132)
nodelabels("Glires", 79)


# ---------- MULTIPLE SEQUENCE ALIGNMENT ----------
# Align sequences using the 'msa' package

# Load the unaligned sequences
seqs <- readDNAStringSet('nadh6.8apes.fasta', format='fasta')

# The msa package provides several alignment tools. By default it uses the 
# program "ClustalW". To align the sequences using the default options
aln <- msa(seqs)

# Take a look at the alignment
aln

# Look at the entire alignment, not truncated
print(aln, show='complete')

# Save the multiple sequence alignment
writeXStringSet(unmasked(aln), file='nadh6.8apes.aln.fasta')

# EXERCISE 2:
#
# Try different options. For each variation, save the MSA to a different file and compare in AliView
#
# 1. Try a different algorithm: ClustalOmega or MUSCLE
#    e.g. aln <- msa(seqs, "Muscle")
#
# 2. Try varying the gapOpening and gapExtension values
#    e.g. msa(seqs, gapOpening=0, gapExtension=0)
#    e.g. msa(seqs, gapOpening=1000, gapExtension=10000)
#    What's happening here? How do the resulting MSAs differ?

  
# ---------- CALCULATE PAIRWISE DISTANCES ----------

# The 'ape' library provides many tools for phylogenetic analysis. 

# Read the aligned sequences we had saved to file above
aln <- read.dna('nadh6.8apes.aln.fasta', format='fasta')

# Lets see what the alignment looks like
checkAlignment(aln)

# What kind of things should we look for to check the quality of our alignment?

# Calculate pairwise distances between sequences
D <- dist.dna(aln, model='JC')

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

# state the outgroup for the sequences
outg <- c('pig', 'cow')

upgma_tree <- upgma(D)
upgma_tree <- ladderize(root(upgma_tree, outgroup=outg, resolve.root = TRUE))
plot(upgma_tree)

# Q: What's the distinctive feature of a UPGMA tree? (look at the plot)

# look at the branch lengths
upgma_tree$edge.length

# look at the newick representation
write.tree(upgma_tree)

# save the tree to file
write.tree(upgma_tree, file="nadh6_upgma.tree")

# Open the tree in Figtree (if you installed it)? Does it render the tree the same way?


# ---------- CREATE A NEIGHBOUR-JOINING TREE USING DISTANCE MATRIX ----------

# Use the distances we calculated above
nj_tree <- nj(D)

# View the tree
plot(nj_tree)
add.scale.bar()

# By default the tree is unrooted. Even if it looks OK, _always_ explicitly root the tree
rooted_nj <- ladderize(root(nj_tree, outgroup=outg, resolve.root=TRUE))
plot(rooted_nj)
add.scale.bar()

# To save a tree
write.tree(nj_tree, file = 'nadh6_nj.tree')


# ---------- ESTIMATE TREE BY MAXIMUM PARSIMONY ----------

# Use the phangorn library
aln <- read.phyDat('nadh6.8apes.aln.fasta', format = "fasta")

# Use a random starting tree
rnd_tree <- rtree(length(aln), tip.label = names(aln))
plot(rnd_tree)

# the parsimony score for the random starting tree
parsimony(rnd_tree, aln)

# optimise the tree using maximum parsimony
pars_tree <- optim.parsimony(rnd_tree, aln)
pars_tree <- ladderize(root(pars_tree, outgroup=outg, resolve.root=TRUE))
plot(pars_tree)

# what's missing? look at the newick representation
write.tree(pars_tree)

# we can look at the changes on the tree for any given site. look at the 12th site
anc.pars <- ancestral.pars(pars_tree, aln)
plotAnc(pars_tree, anc.pars, attr(anc.pars, "index")[12])

# ---------- ESTIMATE TREE BY MAXIMUM LIKELIHOOD ----------

# Package up the start tree (nj tree), alignment and model information into a likelihhod
fit_ini <- pml(nj_tree, aln, k = 4)

# Optimise the tree+alignment+model using maximum likelihood (ML)
fit <- optim.pml(fit_ini, optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE)

# Get the ML tree
ml_tree <- fit$tree

# Root the ML tree using the outgroup species
ml_tree <- ladderize(root(ml_tree, outgroup=outg, resolve.root=T))

# Plot the ML tree
plot(ml_tree)
add.scale.bar()

# Maybe a bit clearer without the outgroup
plot(drop.tip(ml_tree, outg))

# ---------- COMPARING TREES ----------

# Let's compare the NJ tree with the ML tree
# i.e. nj_tree & ml.tree

# We want a plot with plots on 2 rows, 1 column (mfrow). Use small margins (mai)
par(mfrow=c(2,1), mai=c(0.2,0.2,0.2,0.2))

# Plot the NJ tree, then the ML tree and add a scale bar
plot(rooted_nj)
edgelabels(round(rooted_nj$edge.length, digits=3))
add.scale.bar()

plot(ml_tree)
edgelabels(round(ml_tree$edge.length, digits=3))
add.scale.bar()

# Calculate the Robinson-Foulds distance between the two trees
# If they are identical the RF distance is 0
RF.dist(pars_tree, ml_tree)
RF.dist(nj_tree, ml_tree)

# Compare the phylogenies visually
comparePhylo(nj_tree, ml_tree, plot=TRUE)

# ----------- LARGER EXAMPLE (NOT COMPLETE!) ----------

# Compare trees for two mitochondrial genes

outg <- 'sus_scrofa'

# cyb msa
seqs1 <- readDNAStringSet('cyb.30prim.fasta', format='fasta')
aln1 <- msa(seqs1)
writeXStringSet(unmasked(aln1), file='cyb.30prim.aln.fasta')

# estimate a neighbour joining tree for cyb
aln1 <- read.dna('cyb.30prim.aln.fasta', format='fasta')
dist1 <- dist.dna(aln1)
tree1 <- nj(dist1)
tree1 <- ladderize(root(tree1, outgroup=outg, resolve.root=TRUE))

# nd6 msa
seqs2 <- readDNAStringSet('nd6.30prim.fasta', format='fasta')
aln2 <- msa(seqs2)
writeXStringSet(unmasked(aln2), file='nd6.30prim.aln.fasta')

# estimate a nj tree for nd6
aln2 <- read.dna('nd6.30prim.aln.fasta', format='fasta')
dist2 <- dist.dna(aln2)
tree2 <- nj(dist2)
tree2 <- ladderize(root(tree2, outgroup=outg, resolve.root=TRUE))

# compare the two trees
comparePhylo(tree1, tree2, plot=TRUE)

association <- cbind(tree1$tip.label, tree2$tip.label)
cophyloplot(tree1, tree2, assoc = association)

# estimate cyb tree using maximum likelihood
aln <- read.phyDat('cyb.30prim.aln.fasta', format='fasta')
fit_ini <- pml(tree1, aln, k = 4)
fit <- optim.pml(fit_ini, optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE)
ml_tree <- ladderize(root(fit$tree, outgroup=outg, resolve.root=TRUE))

# compare with cyb nj tree
comparePhylo(tree1, ml_tree, plot=TRUE)

# plot one tree on top of the other
par(mfrow=c(2,1), mai=c(0.2,0.2,0.2,0.2))
plot(tree1)
add.scale.bar()
plot(ml_tree)
add.scale.bar()

association <- cbind(tree1$tip.label, ml_tree$tip.label)
cophyloplot(tree1, ml_tree, assoc = association)

# is nj tree a good fit to the observed pairwise distances?
dist2 <- as.dist(cophenetic(tree1))
plot(dist1, dist2)
abline(0,1)


# compare the distances of nj and ml trees
d_nj <- cophenetic(tree1)
d_ml <- cophenetic(ml_tree)
plot(d_nj, d_ml)
abline(0,1)

d_upgma <- as.dist(cophenetic(upgma_tree))
plot(D, d_upgma)
abline(0,1)

