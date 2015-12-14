## this script calculates richness and endemism from modelled suitability surfaces
## it requires all the model grids to have the same extent

library(phylobase)
library(foreach)
library(plyr)
library(doParallel)

parent.prob <- function(probabilities) {
  # probabilities is a vector of values between 0 and 1
  # add code to check values are of correct type!
  parent.prob <- 1 - prod(1-probabilities)
  return(parent.prob)
}

parent.prob.groups <- function(probabilities,model.groups) {
  # probabilities is a vector of values between 0 and 1
  # values in the same model group are added,
  # across different groups, calculate the parent probability

  if (length(probabilities) != length(model.groups)) {
    cat("\nError in function parent.prob.groups.  The probabilities and groups vectors are of different lengths\n")
    return(NULL)
  }

  groups <- unique(model.groups)
  if (length(groups) == 1) {
    parent.prob <- min(sum(probabilities),1)
  } else {
    group.prob <- vector("numeric")
    i <- 1
    for (group in groups) {
      group.prob[i] <- min(sum(probabilities[which(model.groups==group)]),1)
      probabilities[which(model.groups==group)] <- 0
      i <- i+1
    }
    probabilities <- append(probabilities,group.prob)
    parent.prob <- 1 - prod(1-probabilities)
  }

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

calc_PE <- function(tree, sites_x_tips,presence=c("presence","abundance","probability")) {
  # add code to check that the values are correct for the presence type:
    # 0 or 1 for presence - this calculates PE (Rosauer et al 2009)
    # from 0 to 1 for probability - this calculates model weighted PE (Rosauer, in prep)
    # any value for abundance - this calculation is equivalent to BED (Cadotte & Davies 2010)

  # change to a phylobase phylo4 object
  if (class(tree) == "phylo") {tree <- phylo4(tree)}

  sites_x_branches <- data.frame(cbind(rep(0,nrow(sites_x_tips))))

  for (i in 1:nTips(tree)) {
    sites_x_branches[,i] <- sites_x_tips[,which(labels(tree)[i]==names(sites_x_tips))]
    names( sites_x_branches)[i] <- labels(tree)[i]
    cat(i,dim(sites_x_branches),"\n")
  }
  rm(sites_x_tips); gc()
  branch.labels <- labels(tree)
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
    cat(i,branch.labels[i],length(desc),"\n")
    gc(verbose=F)
  }

  #scale columns (branches) to sum to 1
  branch.lengths <- as.numeric(edgeLength(tree,1:branch.count))
  sites_x_branches <- apply(sites_x_branches,MARGIN=2,FUN=scale.to,1)

  sites_x_branches <- sites_x_branches[,1:branch.count] * branch.lengths
  PE.vec <- apply(sites_x_branches,MARGIN=1,FUN=sum,na.rm=T)

  PE <- data.frame(cbind(1:nrow(sites_x_branches),PE.vec))
  names(PE) <- c("site","PE")
  return(PE)
}

calc_PE_from_models <- function(tree, sites_x_tips, core_count=1) {
# tree must be a phylo4d tree containing only the tips for which there is a model of the same name
# the tree has one data column called model.groups

  # change to a phylobase phylo4 object
  if (class(tree) == "phylo") {tree <- phylo4(tree)}

  sites_x_branches <- data.frame(cbind(rep(0,nrow(sites_x_tips))))

  for (i in 1:nTips(tree)) {
    sites_x_branches[,i] <- sites_x_tips[,which(labels(tree)[i]==names(sites_x_tips))]
    names(sites_x_branches)[i] <- labels(tree)[i]
    cat(i,dim(sites_x_branches),"\n")
  }

  # calculate 'species' measures
  cat("\nCalculating species measures\n")

  SR.vec <- apply(sites_x_tips, MARGIN = 1, FUN = sum)
  cat("SR done\n")

  sites_x_tips    <- apply(sites_x_tips, MARGIN=2, FUN=scale.to, 1)
  WE.vec          <- apply(sites_x_tips, MARGIN=1, FUN=sum, na.rm=TRUE)
  #WE.vec          <- clusterApplyLB(cl=cl, x=sites_x_tips, margins=1, fun=sum, na.rm=TRUE)
  cat("WE done\n\n")

  rm(sites_x_tips); gc()
  branch.labels <- labels(tree)
  branch.count <- length(labels(tree))

  # add names and occupancy columns for internal branches
  for (i in (nTips(tree)+1):branch.count) {     # this can be a foreach loop
    branch.labels[i] <- paste("b",i,sep="")
    desc <- as.integer(descendants(tree,i, type="tips"))

    branch_col <- as.numeric(apply(sites_x_branches[,desc],MARGIN=1,FUN=parent.prob.groups,tipData(tree)[desc,1]))

    sites_x_branches[,i] <- branch_col
    names(sites_x_branches[i]) <- branch.labels[i]
    cat(i, "\t", branch.labels[i], "\t", length(desc), "\n")
    gc(verbose=F)
  }

  # calculate PD
  cat("\nCalculating PD\n")
  branch.lengths <- as.numeric(edgeLength(tree,1:branch.count))
  sites_x_branches_PD <- matrix(0,nrow=nrow(sites_x_branches),ncol=ncol(sites_x_branches))
  for (i in 1:ncol(sites_x_branches)) {
    sites_x_branches_PD[,i] <- sites_x_branches[,i] * branch.lengths[i]
  }

  PD.vec <- apply(sites_x_branches_PD,MARGIN=1,FUN=sum,na.rm=T)
  rm(sites_x_branches_PD)
  gc()

  cat("\nCalculating PE\n")

  #scale columns (branches) to sum to 1
  sites_x_branches <- apply(sites_x_branches,MARGIN=2,FUN=scale.to,1)

  #sites_x_branches <- sites_x_branches[,1:branch.count] * branch.lengths
  for (i in 1:branch.count) {
    sites_x_branches[,i] <- sites_x_branches[,i] * branch.lengths[i]
  }

  PE.vec <- apply(sites_x_branches,MARGIN=1,FUN=sum,na.rm=T)

  PE <- data.frame(cbind(1:nrow(sites_x_branches),PE.vec, PD.vec, WE.vec, SR.vec))
  names(PE) <- c("site","PE","PD","WE","SR")
  return(PE)
}
