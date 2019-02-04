# scratch

# bootstrapping
ml_opt <- function(a) {
  d <- dist.dna(a)
  t <- nj(d)
  fini <- pml(t, as.phyDat(a), k=4)
  fit <- optim.pml(fini, optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE)
  return(fit$tree)
}

boots <- boot.phylo(ml_tree, aln, function(x) ml_opt(x))
nodelabels(boots)



boots <- boot.phylo(tree1, aln1, function(e) root(nj(dist.dna(e)), outgroup='sus_scrofa'))
nodelabels(boots)



# removing similar taxa

comparePhylo(cyb_tree, ml.tree, plot=TRUE)

otuPhylo <- function(phy, cutoff=0.05) {
  
  dists <- cophenetic(phy)
  diag(dists) <- NA
  pruned <- vector()
  
  for (i in phy$tip.label) {
    if (!(i %in% pruned)) {
      d <- na.omit(dists[i,])
      pruned <- c(pruned, names(which(d <= cutoff)))
    }
  }
  
  if (length(pruned) > 0) {
    print("Dropping the following taxa:")
    print(pruned)
    return(drop.tip(phy, pruned))
  }
  else {
    return(phy)
  }
  
}



drop <- cyb_tree$tip.label[!(cyb_tree$tip.label %in% apes)]
cyb_tree <- drop.tip(cyb_tree, drop)


# get subset of sequences
cyb_aln <- cyb_aln[which(labels(cyb_aln) %in% apes),]
nd6_aln <- nd6_aln[which(labels(nd6_aln) %in% apes),]

cyb_tree <- otuPhylo(cyb_tree, 0.15)

drop <- nd6_tree$tip.label[!(nd6_tree$tip.label %in% cyb_tree$tip.label)]
nd6_tree <- drop.tip(nd6_tree, drop)

cyb_aln <- cyb_aln[which(labels(cyb_aln) %in% cyb_tree$tip.label),]
nd6_aln <- nd6_aln[which(labels(nd6_aln) %in% cyb_tree$tip.label),]

write.FASTA(cyb_aln, file='cyb.30prim.aln.fasta')
write.FASTA(nd6_aln, file='nd6.30prim.aln.fasta')



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