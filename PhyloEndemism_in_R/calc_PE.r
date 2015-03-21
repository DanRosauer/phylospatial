##################################################
##  A script to calculate phylogenetic endemism ##
##  Dan Rosauer                                 ##
##  dan.rosauer@anu.edu.au                      ##
##  revised March 2015                          ##
##################################################

## the calc_PE() function in this script calculates phylogenetic endemism from a
## tree and a sites by tips (species) matrix.
  
##   tree                 - the names on tree's tip labels need to match the column names on the matrix
##   sites x tips matrix  - the sites should represent equal areas (presumably grid cells).   There is no need to include unoccupied grid cells.  One easy way to get this, is to simply round the coordinates to the appropriate grid resolution, and then group occurrences at the same rounded location together using aggregate or dplyr
##   presence             - specifies what the values in the matrix cells mean, and how to calculate PE
##                             presence uses 0 or 1 for presence / absence.  This is PE exactly as in Rosauer et al 2009
##                             abundance is an amount (could be number of individuals, proportion of cell occupied etc).  With abundance, PE is equivalent to Caddotte & Davies BED)
##                             probability is a value from 0 to 1, for example from an SDM.  Probability then propagates to internal branches at the probability that any of the descendent branches are present.  This method is described in a paper of mine in, which is being pending minor revisions.

library(phylobase)

calc_PE <- function(tree, sites_x_tips,presence=c("presence","abundance","probability")) {
  
  # add code to check that the values are correct for the presence type:
  # 0 or 1 for presence - this calculates PE (Rosauer et al 2009)
  # from 0 to 1 for probability - this calculates model weighted PE (Rosauer, in prep)
  # any value for abundance - this calculation is equivalent to BED (Cadotte & Davies 2010)
  
  #default value for presence
  if (is.na(presence)) {presence="presence"}
  
  # change to a phylobase phylo4 object
  if (class(tree) == "phylo") {tree <- phylo4(tree)}
  
  sites_x_branches <- data.frame(cbind(rep(0,nrow(sites_x_tips))))
  
  for (i in 1:nTips(tree)) {
    sites_x_branches[,i] <- sites_x_tips[,which(labels(tree)[i]==names(sites_x_tips))]
    names( sites_x_branches)[i] <- labels(tree)[i]
  }
  
  rm(sites_x_tips); gc()
  branch.labels <- as.character(labels(tree))
  branch.count <- length(labels(tree))
  
  # add names and occupancy columns for internal branches
  for (i in (nTips(tree)+1):branch.count) {
    branch.labels[i] <- paste("b",i,sep="")
    desc <- as.integer(descendants(tree,i, type="tips"))
    if (presence=="abundance") {
      branch_col <- as.numeric(apply(sites_x_branches[,desc],MARGIN=1,FUN=sum))
    } else if (presence=="presence") {
      branch_col <- as.numeric(apply(sites_x_branches[,desc],MARGIN=1,FUN=max))
    } else if (presence=="probability") {
      branch_col <- as.numeric(apply(sites_x_branches[,desc],MARGIN=1,FUN=parent.prob))
    }
    sites_x_branches[,i] <- branch_col
    names(sites_x_branches[i]) <- branch.labels[i]
    #cat(i,branch.labels[i],length(desc),"\n")
    gc(verbose=F)
  }
  
  #scale columns (branches) to sum to 1
  sites_x_branches <- apply(sites_x_branches,MARGIN=2,FUN=scale.to,1)
  
  #now scale branches to sum to their length
  branch.lengths <- as.numeric(edgeLength(tree,1:branch.count))
  branch.lengths[is.na(branch.lengths)] <- 0
  for (i in 1:length(branch.lengths)) {
    sites_x_branches[,i] <- sites_x_branches[,i] * branch.lengths[i]
  }
  
  PE.vec <- apply(sites_x_branches,MARGIN=1,FUN=sum,na.rm=T)
  PE <- data.frame(cbind(1:nrow(sites_x_branches),PE.vec))
  names(PE) <- c("site","PE")
  return(PE)
}


parent.prob <- function(probabilities) {
  # probabilities is a vector of values between 0 and 1
  # add code to check values are of correct type!
  parent.prob <- 1 - prod(1-probabilities)
  return(parent.prob)
}

scale.to <- function(vec,vec.sum) {
  #mat is a vector
  #this function rescales each the vector values to sum to 'vec.sum'
  vec.tot <- sum(vec,na.rm=TRUE)
  if (vec.tot > 0) {
    vec.out <- vec.sum*vec/vec.tot
  } else {
    vec.out <- rep(0,times=length(vec))  #should columns for tips / branches with no occurrences be removed?
  }
  return(vec.out)
}

