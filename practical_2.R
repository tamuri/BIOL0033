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

# set the outgroup for our sequences
outg <- c('cow', 'pig')

# ---------- ESTIMATE TREE BY MAXIMUM PARSIMONY ----------

# Use the phangorn library
aln <- read.phyDat('nadh6.8apes.aln.fasta', format = 'fasta')

# Use a random starting tree
rnd_tree <- rtree(length(aln), tip.label = names(aln))
plot(rnd_tree)

# Q: What's the parsimony score for your starting tree?
parsimony(rnd_tree, aln)

# optimise the tree using maximum parsimony
pars_tree <- optim.parsimony(rnd_tree, aln)
pars_tree <- ladderize(root(pars_tree, outgroup=outg, resolve.root=TRUE))
plot(pars_tree)

# Q: What's the parsimony score of the maximum parsimony tree?
# Q: How many operations did it take to find the MP tree?

# what's missing? look at the newick representation
write.tree(pars_tree)

# we can look at the changes on the tree for any given site. look at the 12th site
anc.pars <- ancestral.pars(pars_tree, aln)
plotAnc(pars_tree, anc.pars, attr(anc.pars, 'index')[12])

# EXERCISE 1:
#
# 1. What's the parsimony score for site 12 (assume uniform cost matrix)?
#
# 2. Are the following sites parsimony informative: 5, 12, 57 and 492?

# ---------- ESTIMATE TREE BY MAXIMUM LIKELIHOOD ----------

# Package up a start tree (NJ), alignment and model information into a likelihood object
fit_ini <- pml(nj(dist.hamming(aln)), aln, model='JC')

# Optimise the tree+alignment+model using maximum likelihood (ML)
fit <- optim.pml(fit_ini, optNni = TRUE, model='JC')

# Get the ML tree
ml_tree <- fit$tree

# Root the ML tree using the outgroup species
ml_tree <- ladderize(root(ml_tree, outgroup=outg, resolve.root=T))

# Plot the ML tree
plot(ml_tree)
add.scale.bar()

# Plotting without the outgroup
plot(drop.tip(ml_tree, outg))
add.scale.bar()

# EXERCISE 2:
#
# 1. What's the ML when using the K80 model? What is the ML estimate of the kappa parameter?
#
# 2. What's the likelihood ratio test statistic? Does the K80 model offer a significantly 
#    better fit to the data?
#    (Chi-sq with 1 degree of freedom has critical values of 5%: 3.84; 1%: 6.63)
#
# 3. Try a more complex model and answer the following.
#
#    e.g. try HKY+G (transition/transversion bias, unequal base frequencies, site rate variation)
#    
#    create a new initial likelihood object that includes the gamma rate categories:
#
#        fit_ini <- pml(nj(dist.hamming(aln)), aln, model='HKY', k=4)
#
#    then
#
#        optim.pml(fit_ini, optNni=TRUE, model='HKY', optGamma=TRUE)
#
#    a) What is the MLE of the transition/transversion bias?
#    b) What is the MLE gamma shape parameter? Is this high or low site rate variation?
#    c) How does the total tree length compare with the simple JC model above?


# ---------- GETTING BOOTSTRAP SUPPORT VALUES OF TREE ----------

# Bootstrap allows us to assign confidence scores to splits in our phylogenetic tree

# Load our real alignment (ape-package requires use of the `read.dna` function)
aln <- read.dna('nadh6.8apes.aln.fasta', format='fasta')

# Calculate the NJ tree for this alignment
nj_tree <- nj(dist.dna(aln))
nj_tree <- root(nj_tree, outgroup=outg)

# Look at the tree
plot(nj_tree)

# Calculate the bootstrap support for each split in the ML tree.
# The `boot.phylo` function is provided by the ape package.
# The function pass to `FUN` should be exactly what you did on the real data to estimate the tree
boots <- boot.phylo(phy=nj_tree,
                    x=aln,
                    FUN=function(bs_aln) nj(dist.dna(bs_aln)),
                    trees=TRUE)

# Plot the tree with bootstrap support
plot(nj_tree, main='nj tree w/ bootstrap support')
add.scale.bar()
nodelabels(boots$BP)

# HINT: Save the tree and look at it in Figtree. It might be clearer.
nj_tree$node.label <- boots$BP
write.tree(nj_tree, file='nadh6.8apes.bs.tree')

# Get the majority-rule consensus tree
con1 <- consensus(boots$trees, p=0.5, check.labels=TRUE)
plot(root(con1, outgroup=outg, resolve.root=TRUE), main='majority-rule consensus')

# Get the strict consensus tree
con2 <- consensus(boots$trees, p=1, check.labels=TRUE)
plot(root(con2, outgroup=outg, resolve.root=TRUE), main='strict consensus')

# EXERCISE 3:
#
# 1. Why are the majority-rule and strict consensus trees different?
#
# 2. Which split has the lowest support according to bootstrap analysis?

# ---------- LRT OF MOLECULAR CLOCK ----------

# Test the molecular clock hypothesis
aln <- read.phyDat('nadh6.8apes.aln.fasta', format = 'fasta')

# Optimise the branch lengths under the molecular clock
H0_ini <- pml(upgma(dist.hamming(aln)), aln, model='JC')
H0_opt <- optim.pml(H0_ini, optNni = TRUE, optRooted=TRUE, model='JC')

# Optimise the branch lengths without the molecular clock
H1_ini <- pml(upgma(dist.hamming(aln)), aln, model='JC')
H1_opt <- optim.pml(H1_ini, optNni = TRUE, optRooted=FALSE, model='JC')

# The likelihood:
H1_opt$logLik

# The number of parameters in the model:
H1_opt$df

# Q: Plot the clock tree and unconstrained tree

# EXERCISE 4:
#
# Perform the likelihood ratio test
# 
# 1. How do we calculate the test statistic? What is the value?
#
# 2. Calculate the degrees of freedom between the simpler and more complex model.
#
# 3. What is the P value given the value of test statistic? (use pchisq)
#
# 4. Do we accept or reject the hypothesis of the molecular clock?


# --------- EXTRA: USE PHYML ONLINE TO GET THE ML TREE WITH BOOTSTRAP ----------

# EXERCISE 4:
#
# Visit http://www.atgc-montpellier.fr/phyml/ in your browser
#
# First, convert the FASTA-format alignment to PHYLIP-format

aln <- read.dna('nadh6.8apes.aln.fasta', format='fasta')
write.dna(aln, file='nadh6.8apes.aln.phylip', format='interleaved')

# Have a look at the PHYLIP file and compare with FASTA format

# Then, upload the PHYLIP file as input data
#
# Select the options:
#    a) Let PhyML automatically select the substitution model
#    b) Get bootstrap support for the tree
#    c) Leave the rest as defaults
#
# Finally, submit the job and wait for results in your email

# 1. What was the best model according to the PhyML selection criteria?
#
# 2. What was the maximum log-likelihood of the data given the model?
#
# 3. Is there any nucleotide bias?
#
# 4. How does the PhyML tree & support values compare with the tree you estimated above? 
#    (You can use FigTree to open the tree file - add the bootstrap support 'label' to the nodes)
#    Which split has lowest bootstrap support?


