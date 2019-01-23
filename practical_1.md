Set the "working directory", i.e. the directory where you put all the files for the practical

```
setwd("~/Documents/2019/biol0030/BIOL3010-ZY-resources/Practical")
```

## First multiple sequence alignment

We will use the package "msa" for multiple sequence alignment. First load the package:

```
library(msa)
```

Then load the unaligned sequences:

```
nadh6_sequences <- readDNAStringSet('NADH6.fas', format='fasta')
nadh6_sequences <- readAAStringSet('~/Documents/2018/tree_annotator/etc/examples/flu.pb2.aminoacids/PB2.fasta', format='fasta')
```

The msa package provides several alignment tools. The default is the program "ClustalW". To align the sequences using the default options, type:

```
nadh6_alignment_clustal <- msa(nadh6_sequences)
```

Take a look at the alignment by simply typing:

```
nadh6_alignment_clustal
```

To see the entire alignment, type:

```
print(nadh6_alignment_clustal, show='complete')
```

More:

* Try two other algorithms: ClustalOmega and MUSCLE

  `nadh6_alignment_muscle <- msa(nadh6_sequences, "Muscle")`
  
* Options for the Muscle aligner: gapOpen, gapExtend


## Calculate the distances between sequences

The "ape" library provides many tools for phylogenetic analysis. Load the package in the usual way:

```
library(ape)
```

We need to convert our multiple sequence alignment into the format required by ape:

TODO: this doesn't work
```
nadh6_aligned <- msaConvert(nadh6_alignment_clustal, type='ape::DNAbin')
```

Save and read the file in instead:

```
writeXStringSet(unmasked(nadh6_alignment_muscle), file='test.fasta')
nadh6_aligned <- read.FASTA('test.fasta', type='DNA')
```

To calculate the distances between sequences:

```
nadh6_distances <- dist.dna(nadh6_aligned)
```

Have a look at the distance matrix:

```
nadh6_distances
```

Plot a visual representation of the distances. This requires the "ade4" library

```
library(ade4)


```

Generate the neighbour joining tree using the distance matrix:

```
njTree <- nj(nadh6_distances)
```

View the tree

```
plot(njTree)
```

By default the tree is unrooted. Even if it looks right, we should _always_ explicitly root the tree.

```
rooted_nj_tree <- root(njTree, outgroup=c('pig', 'cow'))
```

Plot the rooted tree. How does it look different?


## Likelihood tree

```
dna2 <- as.phyDat(nadh6_aligned)
tre.ini <- nj(dist.dna(nadh6_aligned, model = "TN93"))
fit.ini <- pml(tre.ini, dna2, k = 4)
fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE)
ml.tree <- fit$tree
ml.tree <- root(ml.tree, outgroup=c('pig', 'cow'))
ml.tree <- ladderize(ml.tree)
plot(ml.tree)
```


