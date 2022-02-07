# ---------- DOWNLOAD THE PRACTICAL ----------
#
# 1. Download the .zip file for the practical from
#
#    https://github.com/tamuri/BIOL0033/archive/master.zip
#
# 2. Extract the files

# ---------- SET THE WORKING DIRECTORY ----------

# This is the place where you put the files for the practical,
# so you should modify the path within brackets below with 
# the path to the directory where you have these files on 
# your PC
setwd( '~/Documents/2019/BIOL0033' )

# Alternatively, you can also automatically find the path
# to the R script you are using if you load a specific 
# function from the `rstudioapi` library. Follow the 
# instructions below:
#
# 1. Load the package `rstudioapi`. If you do not have 
#    it installed, then uncomment and run the
#    command below
# install.packages( "rstudioapi" )
library( rstudioapi ) 
# 2. Get the path to current open R script
path_to_file <- rstudioapi::getActiveDocumentContext()$path
# 3. Now just keep the path to the directory where the 
#    R script is
wd           <- paste( dirname( path_to_file ), "/", sep = "" )
# 4. Set the working directory
setwd( wd )

# ---------- LOAD THE REQUIRED PACKAGES ----------

# `ape` and `phangorn` provide vast array of phylogenetic analyses
library( ape )
library( phangorn )

# Set the outgroup for our sequences
outg <- c( 'cow', 'pig' )

# ---------- ESTIMATE TREE BY MAXIMUM PARSIMONY ----------

# Use the `phangorn` library to read the alignment you generated 
# in practical 1
aln <- phangorn::read.phyDat( file = 'nadh6.8apes.aln.fasta', format = 'fasta' )

# Use a random starting tree
# HINT: You might want to set a seed number before generating the random 
# tree in case you want to reproduce the same results!
set.seed( 12345 )
rnd_tree <- ape::rtree( n = length( aln ), tip.label = names( aln ) )
ape::plot.phylo( rnd_tree )

# Q: What's the parsimony score for your starting tree?
phangorn::parsimony( tree = rnd_tree, data = aln )

# Optimise the tree using maximum parsimony
pars_tree_unr <- phangorn::optim.parsimony( tree = rnd_tree, data = aln )
pars_tree     <- ape::ladderize( ape::root( phy = pars_tree_unr, outgroup = outg, resolve.root = TRUE ) )
ape::plot.phylo( pars_tree )

# Q: What's the parsimony score of the maximum parsimony tree?
# Q: How many operations did it take to find the MP tree?
# Q: What's missing? Look at the Newick representation
write.tree( pars_tree )

# We can look at the changes on the tree for any given site.
# We will be looking at the 12th site!
anc_pars <- phangorn::ancestral.pars( tree = pars_tree, data = aln )
phangorn::plotAnc( tree = pars_tree, data = anc_pars, i = attr( x = anc_pars, which = 'index' )[12],
                   col = c( "green", "blue", "black", "red" ) )

# EXERCISE 1:
#
# 1. What's the parsimony score for site 12 (assume uniform cost matrix)?
#
# 2. Are the following sites parsimony informative: 5, 12, 57 and 492?

# ---------- ESTIMATE TREE BY MAXIMUM LIKELIHOOD ----------

# Package up a start tree (NJ), alignment and model information into a likelihood object
fit_ini <- phangorn::pml( tree = nj( phangorn::dist.hamming( aln ) ), data = aln, model = 'JC' )

# Optimise the tree+alignment+model using maximum likelihood (ML)
fit <- phangorn::optim.pml( object = fit_ini, optNni = TRUE, model = 'JC' )

# Get the ML tree
ml_tree <- fit$tree

# The maximum log-likelihood
fit$logLik

# The number of parameters in the model
fit$df

# Root the ML tree using the outgroup species
ml_tree <- ape::ladderize( phy = ape::root( phy = ml_tree, outgroup = outg, resolve.root = TRUE ) )

# Plot the ML tree
ape::plot.phylo( ml_tree )
ape::add.scale.bar()


# Plotting without the outgroup
ape::plot.phylo( x = ape::drop.tip( phy = ml_tree, tip = outg ) )
ape::add.scale.bar()

# EXERCISE 2:
#
# 1. What's the ML when using the K80 model? What is the ML estimate of the kappa parameter?
#
# 2. What's the likelihood ratio test statistic? Does the K80 model offer a significantly 
#    better fit to the data?
#    (Chi-sq with 1 degree of freedom has critical values of 5%: 3.84; 1%: 6.63)
#
# 3. Try a more complex model and answer the following.
#    E.g., try HKY+G (transition/transversion bias, unequal base frequencies, site rate variation)
#    
#    Now, create a new initial likelihood object that includes the gamma rate categories:
#
#        fit_ini <- pml(nj(dist.hamming(aln)), aln, model='HKY', k=4)
#
#    Then, run the following:
#
#        optim.pml(fit_ini, optNni=TRUE, model='HKY', optGamma=TRUE)
#
#    a) What is the MLE of the transition/transversion bias?
#    b) What is the MLE gamma shape parameter? Is this high or low site rate variation?
#    c) How does the total tree length compare with the simple JC model above?


# ---------- GETTING BOOTSTRAP SUPPORT VALUES OF TREE (NJ) ----------

# Bootstrap allows us to assign confidence scores to splits in our phylogenetic tree

# Load our real alignment (ape-package requires use of the `read.dna` function)
aln <- ape::read.dna( file = 'nadh6.8apes.aln.fasta', format = 'fasta' )

# Calculate the NJ tree for this alignment
nj_tree <- ape::nj( X = ape::dist.dna( x = aln ) )
nj_tree <- ape::root( phy = nj_tree, outgroup = outg )

# Look at the tree
ape::plot.phylo( x = nj_tree )

# Calculate the bootstrap support for each split in the tree.
# The `boot.phylo` function is provided by the ape package.
# The function pass to `FUN` should be exactly what you did on the real data to estimate the tree
boots <- ape::boot.phylo( phy = nj_tree,
                          x = aln,
                          FUN = function( bs_aln ) ape::nj( X = ape::dist.dna( x = bs_aln ) ),
                          trees = TRUE)

# Plot the tree with bootstrap support
ape::plot.phylo( x = nj_tree, main='nj tree w/ bootstrap support')
ape::add.scale.bar()
# ape::nodelabels( boots$BP ) # SAC-210226: Not used anymore as this prints
#                               an extra bootstrap value for the root 
#                               on the tree, which should not be printed.
#                               Use command below:
ape::drawSupportOnEdges( value = boots$BP )
# SAC-210226: It is also recommended that node labels are added
# before rooting the tree. Therefore:
nj_tree_unr            <- ape::nj( X = ape::dist.dna( x = aln ) )
nj_tree_unr$node.label <- boots$BP
nj_tree_bs             <- ape::root( phy = nj_tree, outgroup = outg,
                                     edgelabel = TRUE )
ape::plot.phylo( x = nj_tree_bs, main='nj tree w/ bootstrap support' )
ape::drawSupportOnEdges( value = nj_tree_bs$node.label )

# HINT: Save the tree and look at it in `Figtree` - it might be clearer!
nj_tree$node.label <- boots$BP
ape::write.tree( phy = nj_tree, file='nadh6.8apes.bs.tree')

# Get the majority-rule consensus tree
con1 <- ape::consensus( boots$trees, p = 0.5, check.labels = TRUE)
ape::plot.phylo( x = ape::root( phy = con1, outgroup = outg, resolve.root = TRUE ),
                 main = 'majority-rule consensus' )

# Get the strict consensus tree
con2 <- ape::consensus( boots$trees, p = 1, check.labels = TRUE )
ape::plot.phylo( x = ape::root( phy = con2, outgroup = outg, resolve.root = TRUE ),
                 main = 'strict consensus' )

# EXERCISE 3:
#
# 1. Why are the majority-rule and strict consensus trees different?
#
# 2. Which split has the lowest support according to bootstrap analysis?

# ---------- GETTING BOOTSTRAP SUPPORT VALUES OF TREE (ML) ----------

# Load the alignment (ape-package requires use of the `read.dna` function)
aln <- ape::read.dna( file = 'nadh6.8apes.aln.fasta', format = 'fasta' )

# A function to estimate the ML tree for this alignment
# See section 'ESTIMATE TREE BY MAXIMUM LIKELIHOOD' above for details
# (model is JC in this example, but you can change)
get_ml_tree <- function(x) {
  x_pd    <- phangorn::as.phyDat( x = x )
  fit_ini <- phangorn::pml( tree = ape::nj( phangorn::dist.hamming( x = x_pd ) ),
                            data = x_pd, model = 'JC' )
  fit     <- phangorn::optim.pml( object = fit_ini, optNni = TRUE, model = 'JC' )
  return ( fit$tree )
}

# Run the function to get the tree for the real data
ml_tree <- get_ml_tree( x = aln )

# Look at the tree
ape::plot.phylo( x = ml_tree )

# Calculate the bootstrap support for each split in the tree.
# The `boot.phylo` function is provided by the ape package.
# The function passed to `FUN` should be exactly what you did
# on the real data to estimate the tree.
# Here, we pass the `get_ml_tree` function we defined above:
boots <- ape::boot.phylo( phy = ml_tree,
                          x = aln,
                          FUN = function( bs_aln ) get_ml_tree( x = bs_aln ),
                          trees = TRUE )

# Plot the tree with bootstrap support
ape::plot.phylo( x = ml_tree, main = 'ml tree w/ bootstrap support' )
ape::add.scale.bar()
# SAC-220207: Not used anymore as this prints an extra bootstrap value for the root 
#             on the tree, which should not be printed.
#ape::nodelabels( boots$BP )
ape::drawSupportOnEdges( boots$BP )

# Save the tree with BS values and look at it in `Figtree`.
ml_tree$node.label <- boots$BP
ape::write.tree( phy = ml_tree, file = 'nadh6.8apes.ml.bs.tree' )

# ---------- LRT OF MOLECULAR CLOCK ----------

# Test the molecular clock hypothesis
aln <- phangorn::read.phyDat( file = 'nadh6.8apes.aln.fasta', format = 'fasta')

# Optimise the branch lengths under the molecular clock
H0_ini <- phangorn::pml( tree = phangorn::upgma( phangorn::dist.hamming(x = aln ) ),
                         data = aln, model = 'JC' )
H0_opt <- phangorn::optim.pml( object = H0_ini, optNni = TRUE,
                               optRooted = TRUE, model = 'JC' )

# Optimise the branch lengths without the molecular clock
H1_ini <- phangorn::pml( tree = phangorn::upgma( phangorn::dist.hamming( x = aln ) ),
                         data = aln, model = 'JC' )
H1_opt <- phangorn::optim.pml( object = H1_ini, optNni = TRUE,
                               optRooted=FALSE, model = 'JC' )

# Q: Plot the clock tree and unconstrained tree

# EXERCISE 4:
#
# Perform the likelihood ratio test and answer the questions below:
# 
# 1. How do we calculate the test statistic? What is the value?
#
# 2. Calculate the degrees of freedom between the simpler and more complex model.
#
# 3. What is the p-value given the value of test statistic?
#    HINT: you can use `pchisq`!
#
# 4. Do we accept or reject the hypothesis of the molecular clock?


# --------- EXTRA: USE PHYML ONLINE TO GET THE ML TREE WITH BOOTSTRAP ----------

# EXERCISE 4:
#
# Visit http://www.atgc-montpellier.fr/phyml/ in your browser
#
# First, convert the FASTA-format alignment to PHYLIP-format

aln <- ape::read.dna( file = 'nadh6.8apes.aln.fasta', format = 'fasta' )
ape::write.dna( x = aln, file = 'nadh6.8apes.aln.phylip', format = 'interleaved' )

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


