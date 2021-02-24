# ---------- SET THE WORKING DIRECTORY ----------
# i.e. the place where you put the files for the practical

setwd( '~/Documents/2019/BIOL0033' )

# load the APE package
library( ape )

# ---------- EXTRACTING CODON POSITIONS FROM A GENE ALIGNMENT ----------

# read the alignment
aln <- ape::read.dna( file = 'nadh6.8apes.aln.fasta', format = 'fasta' )

# ncol() gives the number of columns in alignment
# get the site index of each of the codon positions
codon1_ix <- seq( from = 1, to = ncol( aln ), by = 3 )
codon2_ix <- seq( from = 2, to = ncol( aln ), by = 3 )
codon3_ix <- seq( from = 3, to = ncol( aln ), by = 3 )

# alignment of each codon position
aln_codon1 <- aln[, codon1_ix]
aln_codon2 <- aln[, codon2_ix]
aln_codon3 <- aln[, codon3_ix]

# alignment without 3rd codon position
aln_nocodon3 <- aln[, -codon3_ix]

# ---------- SAVE ALIGNMENT IN NEXUS FORMAT ----------

# read alignment from fasta format
aln <- ape::read.dna( file = 'nadh6.8apes.aln.fasta', format = 'fasta' )

# write alignment in nexus format
ape::write.nexus.data( x = aln, file = 'nadh6.8apes.aln.nexus', format = 'dna' )

