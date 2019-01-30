
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

