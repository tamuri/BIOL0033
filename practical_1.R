# ---------- SET THE WORKING DIRECTORY ----------
# i.e. the place where you put the files for the practical

setwd("~/Documents/2019/BIOL0033")

# ---------- BASIC PHYLOGENY PLOTTING -----------

# Load the 'ape' phylogenetics package
library(ape)

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
nodelabels("Carnivora", PUT-NODE-NUMBER)
nodelabels("Glires", PUT-NODE-NUMBER)

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
# aligned_sequences <- msaConvert(aligned_sequences, type='ape::DNAbin')

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
plot(ladderize(nj_tree))
add.scale.bar()

# By default the tree is unrooted. Even if it looks OK, _always_ explicitly root the tree
nj_tree <- root(nj_tree, outgroup=c('pig', 'cow'), resolve.root=TRUE)
nj_tree <- ladderize(nj_tree)
plot(nj_tree)
add.scale.bar()

# Plot the rooted tree. how does it look different?

# ---------- ESTIMATE TREE BY MAXIMUM LIKELIHOOD ----------

# The 'phangorn' library provides several methods of pylogeny estimation
library(phangorn)

# Convert the aligned sequences loaded using 'ape' into the format required by phangorn
aligned_sequences_pd <- as.phyDat(aligned_sequences)

# Get an initial tree to optimise
tre.ini <- nj(dist.dna(aligned_sequences, model = "JC"))

# Package up the tree, alignment and model information into a likelihhod
fit.ini <- pml(tre.ini, aligned_sequences_pd, k = 4)

# Optimise the tree+alignment+model using maximum likelihood (ML)
fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE)

# Get the ML tree
ml.tree <- fit$tree

# Root the ML tree using the outgroup species
ml.tree <- root(ml.tree, outgroup=c('pig', 'cow'), resolve.root=TRUE)
ml.tree <- ladderize(ml.tree)

# Plot the ML tree
plot(ml.tree)
add.scale.bar()



# ---------- COMPARING TREES ----------

# Let's compare the NJ tree with the ML tree
# i.e. nj_tree & ml.tree

# We want a plot with plots on 2 rows, 1 column (mfrow). Use small margins (mai)
par(mfrow=c(2,1), mai=c(0.2,0.2,0.2,0.2))

# Plot the NJ tree, then the ML tree and add a scale bar
plot(nj_tree, x.lim=c(0,1))
plot(ml.tree, x.lim=c(0,1))
add.scale.bar()

# Calculate the Robinson-Foulds distance between the two trees
# If they are identical the RF distance is 0
rf_dist <- RF.dist(nj_tree, ml.tree)

# Compare the phylogenies visually
comparePhylo(nj_tree, ml.tree, plot=TRUE)

