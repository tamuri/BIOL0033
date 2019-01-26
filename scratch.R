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