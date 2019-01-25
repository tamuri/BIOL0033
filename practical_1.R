# ---------- SET THE WORKING DIRECTORY ----------
# i.e. the place where you put the files for the practical

setwd("~/Documents/2019/biol0030/BIOL0033")

# ---------- LOAD ALL THE REQUIRED PACKAGES ----------

# 'ape' and 'phangorn' provide vast array of phylogenetic analyses
library(ape)
library(phangorn)

# 'msa' provides multiple sequence alignment
library(msa)

# ---------- BASIC PHYLOGENY PLOTTING -----------

# Trees are typically stored sing the Newick format
my_newick <- '(a, (b, c));'

# Plot the tree
# Learn how the Newick format represents the plotted tree
my_tree <- read.tree(text=my_newick)
print(my_tree)
plot(my_tree, type='phylogram')

# Try the following Newick strings and plot each
# 1. (frog, human, bird);
# 2. (gorilla:1, (human:2, chimp:3):4);

# What is the Newick string for the following tree. Test it to make sure
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

# Read in a tree from a file (take a look at the file too)
mammals_72sp <- read.tree(file='72sp.tree')

# Get summary of tree
mammals_72sp

# Get tips
mammals_72sp$tip.label

# Make the margins a bit smaller for the plots
par(mai=c(0.1, 0.1, 0.1, 0.1))

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

# Get the node numbers for internal nodes
plot(mammals_72sp, type='phylogram', use.edge.length=FALSE)
nodelabels()

# Plot with specific node highlighted
plot(mammals_72sp, type='phylogram', use.edge.length=FALSE)
nodelabels("Carnivora", 132)
nodelabels("Glires", 79)


# ---------- MULTIPLE SEQUENCE ALIGNMENT ----------
# Align sequences using the 'msa' package

# Load the unaligned sequences
sequences <- readDNAStringSet('nadh6.8apes.fasta', format='fasta')

# The msa package provides several alignment tools. By default it uses the 
# program "ClustalW". To align the sequences using the default options
aligned_sequences <- msa(sequences)

# Take a look at the alignment
aligned_sequences

# Look at the entire alignment, not truncated
print(aligned_sequences, show='complete')

# Save the multiple sequence alignment
writeXStringSet(unmasked(aligned_sequences), file='nadh6.8apes.aln.fasta')

# Try some variations below. For each, save the MSA to a different file and compare in AliView

# Try a different algorithm: ClustalOmega or MUSCLE
# e.g. aligned_sequences <- msa(aligned_sequences, "Muscle")`

# Try varying the gapOpening and gapExtension values
# e.g. aligned_sequences <- msa(aligned_sequences, "Muscle", gapOpening=0, gapExtension=0)
# e.g. aligned_sequences <- msa(aligned_sequences, "Muscle", gapOpening=1000, gapExtension=1000)


# ---------- VISUALISE TREES ----------
# https://4va.github.io/biodatasci/r-ggtree.html
  
  
# ---------- CALCULATE PAIRWISE DISTANCES ----------

# The 'ape' library provides many tools for phylogenetic analysis. 
  
# Load the package
library(ape)

# Read the aligned sequences we had saved to file above
aligned_sequences <- read.dna('nadh6.8apes.aln.fasta', format='fasta')

# Lets see what the alignment looks like
checkAlignment(aligned_sequences)

# What kind of things should we look for?

# Calculate pairwise distances between sequences
distances <- dist.dna(aligned_sequences, model='JC')

# Look at the distance matrix
distances


# ---------- CREATE A UPGMA TREE USING DISTANCES ----------

upgma_tree <- upgma(distances)
upgma_tree <- ladderize(root(upgma_tree, outgroup=c('pig', 'cow'), resolve.root = TRUE))
plot(upgma_tree)

# look at the branch lengths
upgma_tree$edge.length

# look at the newick representation
write.tree(upgma_tree)

# save the tree to file
write.tree(upgma_tree, file="nadh6_upgma.tree")

# ---------- CREATE A NEIGHBOUR-JOINING TREE USING DISTANCES ----------

# Use the distances we calculated above
nj_tree <- nj(distances)

# View the tree
plot(nj_tree)
add.scale.bar()

# By default the tree is unrooted. Even if it looks OK, _always_ explicitly root the tree
rooted_nj <- ladderize(root(nj_tree, outgroup=c('pig', 'cow'), resolve.root=TRUE))
plot(rooted_nj)
add.scale.bar()

# To save a tree
write.tree(nj_tree, file = 'nadh6_nj.tree')


# ---------- ESTIMATE TREE BY MAXIMUM PARSIMONY ----------

# Convert the aligned sequences loaded using 'ape' into the format required by phangorn
aligned_sequences_pd <- as.phyDat(aligned_sequences)

random_starting_tree <- random.addition(aligned_sequences_pd)

# calculate the parsimony score for the NJ tree
parsimony(random_starting_tree, aligned_sequences_pd)

# optimise the tree using maximm parsimony
pars_tree <- optim.parsimony(random_starting_tree, aligned_sequences_pd)
pars_tree <- ladderize(root(pars_tree, outgroup=c('cow', 'pig'), resolve.root=TRUE))
plot(pars_tree)

# what's missing? look at the newick representation
write.tree(pars_tree)

# we can look at the changes on the tree for any given site. look at the 12th site
anc.pars <- ancestral.pars(pars_tree, aligned_sequences_pd)
plotAnc(pars_tree, anc.pars, attr(anc.pars, "index")[12])

# ---------- ESTIMATE TREE BY MAXIMUM LIKELIHOOD ----------

# Get an initial tree to optimise
tre.ini <- nj_tree

# Package up the tree, alignment and model information into a likelihhod
fit.ini <- pml(tre.ini, aligned_sequences_pd, k = 4)

# Optimise the tree+alignment+model using maximum likelihood (ML)
fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE)

# Get the ML tree
ml.tree <- fit$tree

# Root the ML tree using the outgroup species
ml.tree <- ladderize(root(ml.tree, outgroup=c('pig', 'cow'), resolve.root=T))

# Plot the ML tree
plot(ml.tree)
add.scale.bar()



# ---------- COMPARING TREES ----------

# Let's compare the NJ tree with the ML tree
# i.e. nj_tree & ml.tree

# We want a plot with plots on 2 rows, 1 column (mfrow). Use small margins (mai)
par(mfrow=c(2,1), mai=c(0.2,0.2,0.2,0.2))

# Plot the NJ tree, then the ML tree and add a scale bar
plot(rooted_nj)
edgelabels(round(rooted_nj$edge.length, digits=3))
add.scale.bar()

plot(ml.tree)
edgelabels(round(ml.tree$edge.length, digits=3))
add.scale.bar()

# Calculate the Robinson-Foulds distance between the two trees
# If they are identical the RF distance is 0
rf_dist <- RF.dist(nj_tree, ml.tree)

# Compare the phylogenies visually
comparePhylo(nj_tree, ml.tree, plot=TRUE)

# ----------- LARGER EXAMPLE ----------

# Compare trees for two mitochondrial genes

# cyb msa
cyb_seqs <- readDNAStringSet('cyb.30prim.fasta', format='fasta')
cyb_aln <- msa(cyb_seqs)
writeXStringSet(unmasked(cyb_aln), file='cyb.30prim.aln.fasta')

# estimate a neighbour joining tree for cyb
cyb_aln <- read.dna('cyb.30prim.aln.fasta', format='fasta')
cyb_dist <- dist.dna(cyb_aln)
cyb_tree <- nj(cyb_dist)
cyb_tree <- ladderize(root(cyb_tree, outgroup='sus_scrofa', resolve.root=TRUE))

# nd6 msa
nd6_seqs <- readDNAStringSet('nd6.30prim.fasta', format='fasta')
nd6_seqs <- msa(nd6_seqs)
writeXStringSet(unmasked(nd6_seqs), file='nd6.30prim.aln.fasta')

# estimate a nj tree for nd6
nd6_aln <- read.dna('nd6.30prim.aln.fasta', format='fasta')
nd6_dist <- dist.dna(nd6_aln)
nd6_tree <- nj(nd6_dist)
nd6_tree <- ladderize(root(nd6_tree, outgroup='sus_scrofa', resolve.root=TRUE))

# compare the two trees
comparePhylo(cyb_tree, nd6_tree, plot=TRUE)

association <- cbind(cyb_tree$tip.label, nd6_tree$tip.label)
cophyloplot(cyb_tree, nd6_tree, assoc = association)

# estimate cyb tree using maximum likelihood
aligned_sequences_pd <- as.phyDat(cyb_aln)
tre.ini <- nj(dist.dna(cyb_aln, model = "TN93"))
fit.ini <- pml(tre.ini, aligned_sequences_pd, k = 4)
fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE)
ml.tree <- fit$tree
ml.tree <- root(ml.tree, outgroup=c('sus_scrofa'), resolve.root=TRUE)
ml.tree <- ladderize(ml.tree)

# compare with cyb nj tree
comparePhylo(cyb_tree, ml.tree, plot=TRUE)

# plot one tree on top of the other
par(mfrow=c(2,1), mai=c(0.2,0.2,0.2,0.2))
plot(cyb_tree)
add.scale.bar()
plot(ml.tree)
add.scale.bar()

association <- cbind(cyb_tree$tip.label, ml.tree$tip.label)
cophyloplot(cyb_tree, ml.tree, assoc = association)



