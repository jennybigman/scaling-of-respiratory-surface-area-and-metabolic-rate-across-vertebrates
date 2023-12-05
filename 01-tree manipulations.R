# Bigman et al. metabolic rate and respiratory surface area scaling Science Advances

# 01 - tree manipulations

library(MCMCglmm)
library(phytools)
library(phangorn)
library(here)


## read in tree ##

tree <- read.tree("final_vert_tree.tre")

# # check to make sure tree is binary and ultrametric ##
is.ultrametric(tree) #this works
is.binary(tree) #no polytomes

#plot tree
plotTree(tree, ftype = "i", fsize = 0.6, lwd = 1)

# now extract phylogeny in tabular form
inv_tree <- inverseA(tree, nodes = "TIPS", scale = TRUE) 
A <- solve(inv_tree$Ainv) # is this reversing the inverse so the matrix is in the form of a covariance matrix and not its inverse (required by brms)
rownames(A) <- rownames(inv_tree$Ainv) # assigning rownames of species ( = to the tips of the tree) to the covariance matrix

## transform into correlation matrix for stan ##
vcov_mat <- as.matrix(A)

