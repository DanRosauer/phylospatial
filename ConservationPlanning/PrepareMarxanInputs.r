
################################################################################
## This code prepares inputs files to run conservation planning tool Marxan   ##
## To optimize selection of areas for conservation of phylogenetic diversity  ##
##                                                                            ##
## Dan Rosauer dan.rosauer@anu.edu.au                                         ##
##                                                                            ##
## Inputs are:                                                                ##
##  - a phylogeny                                                             ##
##  - occurrence in each planning unit, of species (or whatever the tree tips ##
##    are)                                                                    ##
##                                                                            ##
##  This code is at an early stage of development!                            ##
##  Feedback is welcome, but please be patient.                               ##
##  August 2016                                                               ##
##                                                                            ##
################################################################################

library(ape)
library(gdata)
library(phylobase)
library(dplyr)

get_branch_ranges = function(trees,
                             planningUnits,
                             occ,
                             output_dir,
                             amount_max = 0,
                             spf_multiplier = 1,
                             target_args=list(percent=0.1))
{

  ##################
  #   PARAMETERS   #
  ##################

  #   trees is an ape object of class phylo  or phylo_multi
  #
  #   planningUnits lists each planning unit as follows:
  #     PU_ID is an integer to uniquely identify the area
  #     cost is a relative cost of conserving the planning unit.
  #       Optionally cost can be 1 in all cases, if cost is not considered in the analysis
  #     status indicates whether the planning unit is available for protection, or already protected
  #       (see Marxan manual)
  #
  #   occ lists occurences in each planning unit as follows:
  #     TipName is a species or other entity. Only names which match a tip label in Tree_in, are used
  #     PU_ID matches a planning unit in planningUnits
  #     Amount is a number for the amount of the item found in the planning unit.  It could refer to abundance or
  #       area occupied by a species. Optionally the amount can be 1 in all cases, if amounts within planning units
  #       are not used.
  #
  #   amount_max is needed to manage the way that amount values are combined when a branch is represented by more than one
  #     tip in the same planning unit (PU).  If amount_max = 0 (default) if has no effect.  Otherwise the amount for a branch in
  #     a PU is capped at this value.
  #     For example, when calculating the amount in a PU for a branch which has several decendent tips in the PU:
  #     * if amount represents abundance, then the branch amount would be the sum of amounts of the tips
  #     * if amount represents the proportion of the PU occupied, then values could be capped at 1
  #     * a third option, not implemented, could cap amount at the area of a PU.
  #

  library(ape);

  StartTime    <- date()
  StartSysTime <- Sys.time()

  cat("\n\n************************************************\n")
  cat("Starting PhyloSpatial at: ", StartTime)
  cat  ("\n************************************************\n")

  ## check object validity before proceeding ##
  ## full validity checking not yet implemented, but the following is a start

  # validate the planningUnits object and add default values for missing columns
  if (! exists("planningUnits")) { stop("The argument planningUnits was not provided.") }
  planningUnits <- validate_planningUnits(planningUnits)

  # validate the occ object and add default values for missing columns
  if (! exists("occ")) { stop("The argument occ was not provided.") }
  occ <- validate_occ(occ)

  ### end of argument validation ###

  ##-- input data --
  # check if trees is a single tree or a set of trees
  if (class(trees) == "multiPhylo") {
    phy <- trees[[1]] #use the first tree for name matching etc
    number_of_trees <- length(trees)
  } else {
    phy <- trees
    number_of_trees <- 1
  }

  tips_spatial    <- as.character(unique(occ$TipName))
  tips_tree       <- phy$tip.label
  tips_tree_count <- length(length(tips_tree))

  ##-------- Tree calculations ------------------------------------------
  cat("\nDoing tree calculations\n")

  cat("\nStarting counts\n")
  # Counts
  cat("\n  Count of species or tips in the tree: ",tips_tree_count)
  cat("\n  Count of species in the spatial data: ",length(tips_spatial)) # count of species in the data
  cat("\n  Count of species in spatial data that are in tree: ",length(intersect(tips_tree, tips_spatial))) # count of species in data that are in tree
  tips_matching <- tips_spatial[is.element(tips_spatial, tips_tree)] # SpecIDs of species in data that are in tree
  cat("\n  Count of species in spatial data that are not in tree:", length(setdiff(tips_spatial, tips_tree))) # count of species in data that are not in tree
  unmatched_tips <- setdiff(tips_tree, tips_spatial)  # vector of species in tree that are not in data
  cat("\n  Count of species in tree that are not in spatial data",length(unmatched_tips),"\n")

  # trim the first tree to species that are in the spatial data, if required
  branch_count_orig <- length(phy$edge[,1])
  if (length(unmatched_tips) >= 1) {
    phy <- drop.tip(phy,unmatched_tips, trim.internal=TRUE)
    cat("\n  Tips and internal branches not represented in the spatial data have been removed.\n")
    cat("  Branches original: ", branch_count_orig,"\n")
    cat("  Branches remaining: ", length(phy$edge[,1]),"\n")
  } else {
    cat("\n  All", branch_count_orig ,"tree branches successfully matched to spatial data.\n")
  }

  # subset occ to species that are in the tree
  cat ("\n  Occurrence records loaded: ", nrow(occ),"\n")
  occ <- occ[occ$TipName %in% tips_matching,]
  cat ("\n  Occurrence records linked to tree were retained: ", nrow(occ),"\n")

  SpecList <- unique(occ$TipName)    # this way follows list and sequence as actually observed
  TotalSpecCount <- length(SpecList); TotalSpecCount

  # tree preparation outside of the loop
  reorder(phy, order = "cladewise")

  if (number_of_trees > 1) {
    trees_to_use <- 1:number_of_trees
  } else {
    trees_to_use <- 1
  }

  # i is the tree number, j is the iteration number
  run_count <- length(trees_to_use)

  spf.df      <- vector("double", run_count) # a vector to record the mean species SPF from each run

  #START LOOP THROUGH ALL TREES HERE
  for (j in 1:length(trees_to_use)) {

    i <- trees_to_use[j]

    cat("\n********************************************************")
    cat("\n  Starting iteration", j, "using tree ", i)
    cat(" at: ",date())
    cat("********************************************************\n")

    #use the same version of Occ and SpecAll for each iteration
    occ_working <- occ
    #SpecAll <- SpecAll_orig

    #trim unused tips and branches to species that are in the spatial data, if required
    if (j > 1) {
      phy <- trees[[i]]
      branch_count_orig <- length(phy$edge[,1])
      if (length(unmatched_tips) >= 1) {
        phy <- drop.tip(phy,unmatched_tips, trim.internal=TRUE)
        cat("  Branches original: ", branch_count_orig,"\n")
        cat("  Branches after trim to match spatial data: ", length(phy$edge[,1]),"\n")
      } else {
        cat("\n  All", branch_count_orig ,"tree branches successfully matched to spatial data.\n")
      }
    }

    phy4 <- phylo4(phy)
    count_nodes <- nrow(phy4@edge)

    # give names to internal branches (added to species names)
    node_labels <- paste("node_", 1:nNodes(phy4), sep="")
    nodeLabels(phy4) <- node_labels

    #################################################################
    cat("\nCalculating branch ranges for iteration",j,"\n")
    nodes <- getNode(phy4,1:count_nodes)
    BranchDone  <- rep(FALSE,count_nodes)
    BranchRange <- as.numeric(rep(0,count_nodes))
    BranchData  <- data.frame(cbind(names(nodes), BranchDone, BranchRange),stringsAsFactors=FALSE)
    BranchData$BranchDone <- as.logical((BranchDone))
    names(BranchData)[1] <- "BranchNames"
    rm(BranchDone,BranchRange)  # cleaning up

    tree_table <- as.data.frame(print(phy4))
    is.tip <- tree_table$node.type == "tip"

    for (m in 1:count_nodes) {

      if (is.tip[m]) {
        branch_tipname  <- tree_table$label[m]
        #BranchSpecID <- SpecMaster$SpecID[SpecMaster$taxon_name_tree == branch_tipname]  #SpecIDs of species descendant from this node
        branch_occ    <- occ_working[occ_working$TipName %in% branch_tipname,]  # Subset occ_working list to those species
        branch_PU   <- branch_occ[, 2:3]

      } else {

        ######### TESTING
        children <-  tree_table$node[tree_table$ancestor==m]
        children_done <- BranchData$BranchDone[children]

        if (all(children_done)) {
          children_occ  <- branch_occ_all[branch_occ_all$BranchID %in% children,]
          children_occ  <- summarise(group_by(children_occ[,2:3], PU_ID), Amount=sum(Amount))
          branch_PU   <- children_occ
          cat(".")  # remove once its working
        } else {

          ######### TESTING

          branch_tipname  <- names(descendants(phy4, m, type="tips"))
          branch_occ      <- occ_working[occ_working$TipName %in% branch_tipname,]  # Subset occ_working list to descendent tips
          branch_PU       <- summarise(group_by(branch_occ[,2:3], PU_ID), Amount=sum(Amount))
          #branch_PU   <- branch_PU[,-2]
        }

        # ensure that where multiple descendants of a branch occur in the same cell, the total amount does not sum to > amount_max
        if (amount_max > 0) {
          branch_PU$Amount[which(branch_PU$Amount > amount_max)] <- amount_max
        }
      }

      node <- as.integer(nodes[m])
      branch_PU     <- cbind(rep(node,nrow(branch_PU)),branch_PU)
      names(branch_PU) <- c("BranchID","PU_ID","Amount")

      if (m==1) {
        branch_occ_all <- branch_PU
      } else {
        branch_occ_all <- bind_rows(branch_occ_all, branch_PU)
        #branch_occ_all <- rbind(branch_occ_all, branch_PU)
      }

      BranchData$BranchRange[m] <- sum(branch_PU$Amount)
      BranchData$BranchDone[m] <- TRUE

      # progress output - remove once working
      if (m %% 100 == 0) {
        cat("\n",m,"of",count_nodes,"done for tree\t",i)
        minutes_elapsed <- as.double(round(difftime(Sys.time(), StartSysTime, units = "mins"), 2))
        cat("\n", tree_table$label[m],"\t",nrow(branch_PU),"\t",nrow(branch_occ_all), "\t", minutes_elapsed, "minutes\n")
      }
    }

    #name the output folder
    if (run_count > 1) {
      this_output_dir <- paste(output_dir, "_", i ,sep="")  # use tree number as the directory number
    } else {
      this_output_dir <- output_dir
    }

    gc()

    result <- marxan_inputs(write.dir = this_output_dir,
                           planningUnits = planningUnits,
                           tree = phy4,
                           branch_ranges = BranchData$BranchRange,
                           branch_occ = branch_occ_all,
                           target = target_args,
                           spf_multiplier)

#     if (length(spf.df) > 1) {
#       spf.df[j] <- result
#     } else {
#       spf.df    <- result
#     }

    try (rm(BranchData, branch_occ_all))
    gc() # garbage collect from memory

  }
}

############################################

marxan_inputs <- function(
                  tree,
                  branch_ranges,
                  branch_occ,
                  planningUnits,
                  target,  # see documentation for function get_target()
                  spf_multiplier,
                  write.dir)
  {

  cat("\n\nGenerating marxan outputs for",write.dir,"\n\n")

  cons.feature.count <- nrow(tree@edge)
  branch_ranges   <- as.numeric(branch_ranges)
  branch_names    <- as.character(labels(tree))
  branch_IDs      <- as.numeric(getNode(tree))

  branch_targets <- get_target_new(branch_ranges,
                                   target$type,
                                   target$percent,
                                   target$minTarget,
                                   target$maxTarget)

  branch_info     <- as.data.frame(print(tree))
  branch_lengths  <- branch_info$edge.length
  branch_lengths[which(is.na(branch_lengths))] <- 0
  specIDs <- 1:cons.feature.count
  spec <- data.frame(cbind(specIDs, branch_targets, round(branch_lengths * spf_multiplier,4), branch_names), stringsAsFactors=F)
  names(spec) <- c("species_id","target","spf","name")
  spec <- spec[-which(spec$spf==0),]

  # get the mean spf to define the spf for the species version
  species_count <- length(tipLabels(tree))
  mean_spf      <- mean(as.double(spec$spf[1:species_count]))

  pu_from_occ <- unique(branch_occ$PU_ID)
  pu_from_occ <- data.frame(PU_ID = sort(pu_from_occ))

  planningUnits      <- merge(pu_from_occ, planningUnits, by="PU_ID", all=F)
  pu.count     <- nrow(planningUnits)

  # the next 2 lines would create new PU_IDs starting from 1
  planningUnits <- data.frame(cbind(rep(1:pu.count), planningUnits))
  names(planningUnits) <- c("pu_id", "pu_name", "amount", "status")

  pu <- data.frame(planningUnits[,c(1,3,4)])
  names(pu) <- c("id","cost", "status")
  node_by_cell_SpID <- merge(spec, branch_occ, by.x="species_id", by.y="BranchID")

  puvspr2 <- node_by_cell_SpID[,c("species_id","PU_ID","Amount")]
  names(puvspr2) <- c("species","pu_name","amount")  # this data frame has the species number, but the PU name
  rm(node_by_cell_SpID)

  puvspr2 <- merge(planningUnits, puvspr2, by="pu_name")
  puvspr2 <-  puvspr2[,c("species","pu_id","amount.y")]

  names(puvspr2) <- c("species","pu","amount")
  puvspr2$species <- as.numeric(puvspr2$species)

  # now order files to use Marxan quick preparation method
  puvspr2 <- puvspr2[order(puvspr2$pu),]
  sporder <- puvspr2[order(puvspr2$species),]

  # ensure correct column names for maxent format
  names(spec)[1] <- "id"

  try(dir.create(write.dir),silent = TRUE)

  cat("\nAbout to write files to ",write.dir,"\n")

  setwd(write.dir)
  write.csv(spec,"spec.dat", row.names=F)
  write.csv(pu,"pu.dat", row.names=F)
  write.csv(puvspr2,"puvspr2.dat", row.names=F)
  write.csv(sporder,"sporder.dat", row.names=F)
  write.csv(planningUnits,"pu_id_lookup.csv", row.names=F)

  return(mean_tip_spf=mean_spf)
}


marxan_inputs_no_tree <- function(
  branch_lengths,
  branch_ranges,
  branch_names,
  branch_occ,
  planningUnits,
  only_occupied_PU = FALSE,
  target,  # see documentation for function get_target()
  spf_multiplier,
  write.dir)
{

  cat("\n\nGenerating marxan outputs for",write.dir,"\n\n")

  cons.feature.count <- length(branch_lengths)
  branch_ranges   <- as.numeric(branch_ranges)
  branch_names    <- as.character(branch_names)

  branch_targets <- get_target_new(branch_ranges,
                                   target$type,
                                   target$percent,
                                   target$minTarget,
                                   target$maxTarget,
                                   target$maxRangePercent)

  specIDs <- 1:cons.feature.count
  spec <- data.frame(cbind(specIDs, branch_targets, round(branch_lengths * spf_multiplier,4), branch_names), stringsAsFactors=F)
  names(spec) <- c("species_id","target","spf","name")
  spec_to_remove <- which(spec$spf==0)
  if (length(spec_to_remove) > 0) {
    spec <- spec[-which(spec$spf==0),]
  }
  rm(spec_to_remove)

  # get the mean spf to define the spf for the species version
  #species_count <- length(tipLabels(tree))
  #mean_spf      <- mean(as.double(spec$spf[1:species_count]))
  pu.count     <- nrow(planningUnits)
  pu_from_occ <- unique(branch_occ$PU_ID)
  pu_from_occ <- data.frame(PU_ID = sort(pu_from_occ))

  pu.occupied.count  <- length(intersect(pu_from_occ$PU_ID, planningUnits$PU_ID))
  pu.in.occ.not.planningUnits <- length(setdiff(pu_from_occ$PU_ID, planningUnits$PU_ID))

  cat("\nThere are", pu.count, "planning units, of which", pu.occupied.count, "are occupied.\n")
  if (pu.in.occ.not.planningUnits > 0) {
    cat("\nNote: there are", pu.in.occ.not.planningUnits, " occupied planning units which are not included in the supplied table of planning units.\n")
  }

  if (only_occupied_PU) {
  # this option includes only occupied planning units in pu.dat.  It is not a
  # good idea if bounds.dat includes all PUs and will not be changed to match.
    planningUnits      <- merge(pu_from_occ, planningUnits, by="PU_ID", all=F)
    pu.count     <- nrow(planningUnits)
  }

  pu <- planningUnits
  names(pu) <- c("id","cost", "status")

  node_by_cell_SpID <- merge(spec, branch_occ, by.x="name", by.y="TipName")

  puvspr2 <- node_by_cell_SpID[,c("species_id","PU_ID","Amount")]
  names(puvspr2) <- c("species","pu","amount")  # this data frame has the species number, but the PU name
  rm(node_by_cell_SpID)

    puvspr2$species <- as.numeric(puvspr2$species)

  # now order files to use Marxan quick preparation method
  puvspr2 <- puvspr2[order(puvspr2$pu),]
  sporder <- puvspr2[order(puvspr2$species),]

  # ensure correct column names for maxent format
  names(spec)[1] <- "id"

  try(dir.create(write.dir),silent = TRUE)

  cat("\nAbout to write files to ",write.dir,"\n")

  setwd(write.dir)
  write.csv(spec,"spec.dat", row.names=F)
  write.csv(pu,"pu.dat", row.names=F)
  write.csv(puvspr2,"puorder.dat", row.names=F)
  write.csv(sporder,"sporder.dat", row.names=F)
  write.csv(planningUnits,"pu_id_lookup.csv", row.names=F)

  return(1)
}


############################################

get_target <- function(range, option=1) {

# returns a target area, according to a function (this could be made more general with more parameters)
  # option 1: target = 50% if < 20 cells, 10 cells if > 20 cells
  # option 2: target = 25% if < 100 cells, then a flat 25.  Minimum target = 1
  # option 3: target = 100% if < 10 cells, 10 cells if < 100, 10% if < 2000, then flat at 200.

  if (!(option %in% c(1,2,3))) {
    option <- 1
    cat("\nDefaulting to target option 1\n")
  }

  n <- length(range)
  target <- vector(mode="numeric",length=n)

  for (i in 1:n) {

    if (option == 1) {
      target[i] <- min(range[i] * 0.5,10)

    } else if (option == 2) {
      target[i] <- min(range[i] * 0.25,25)
      target[i] <- max(target[i], 1)

    } else if (option ==3) {
      if (range[i] >= 100) {
        target[i] <- min(range[i],10)
      } else if (range[i] > 100) {
        target[i] <- min((range[i] / 10),200)
      }
    }

    # ensure that the target is never greater than the range
    if (target[i] > range[i]) {target[i] <- range[i]}

  }
  return(target)
}

############################################
######### RENAME AS get_target and delete original version, once working
get_target_new <- function(ranges,
                           type=c('percent'),
                           percent, minTarget=0,
                           maxTarget=0,
                           maxRangePercent=100) {

  # returns a target area, according to a function
  # type='percent':
  #   target = the specified percentage of range
  #   arguments minTarget and maxTarget are ignored
  #
  # type='bounded_percent':
  #   target = the specified percentage of range, but a minimum and maximum target may be specified
  #   target values below minTarget are set to the minTarget and those above maxTarget are set to maxTarget
  #   if minTarget or maxTarget are 0, then that limit is ignored, and the percentage is applied without a
  #   minimum or maximum limit

  if (!(tolower(type) %in% c('percent', 'bounded_percent'))) {
    type <- 'percent'
    cat("\nA vaild target type was not provided. Defaulting to a simple percentage target option.\n")
  }

  # check for out of range values of percent
  if (!(percent>0 & percent<=100)) {
    stop("Percent must be between 0 and 100.")
  }

  # check for missing values of minTarget and maxTarget
  if (!(minTarget > 0 | maxTarget > 0)) {
    type <- 'percent'
  }

  if (type == 'percent') {
    targ <- function(x) {return(x * (percent / 100))}
    target <- vapply(ranges, FUN=targ, FUN.VALUE = 1)
    rm(targ)
  }
    else if (type == 'bounded_percent')
  {

    n <- length(ranges)
    target <- vector(mode="numeric",length=n)

    for (i in 1:n) {

      if (minTarget > 0) {
        target[i] <- max(minTarget, ranges[i] * (percent / 100))
      } else {
        target[i] <- ranges[i] * (percent / 100)
      }

      # ensure that target doesn't exceed maxTarget
      if (target[i] > maxTarget & maxTarget > 0) { target[i] <- maxTarget }

      # ensure that target doesn't exceed maxRangePercent
      if (target[i] > (ranges[i] * (maxRangePercent/100))) {
          target[i] <- ranges[i] * (maxRangePercent/100)
      }
    }

  }

  return(target)
}

############################################

fixnames <- function (text_in) {
  #   if (grep("___",text_in)==1) {
  #     text_in <- paste("b",text_in,sep="")
  #   }
  text_out <- gsub('___',"",text_in)
  text_out <- gsub(" ","_",text_out)
  text_out <- paste("x_",text_out,sep="")
  return(text_out)
}

validate_planningUnits = function(planningUnits) {
  # this function checks the validity of the planningUnits object
  # in future validity of all objects should be checked

  # planningUnits$PU_ID is required, so if not valid, execution stops
  # if planningUnits$cost is not valid, a column of 1 is added with a warning
  # if planningUnits$status is not valid, a column of 0 is added with a warning
  if (class(planningUnits) == "data.frame") {
    if (class(planningUnits$PU_ID) == "integer") {
      if (length(planningUnits$PU_ID) == length(unique(planningUnits$PU_ID))) {

        # get to here if planningUnits is a valid data frame, and PU_ID is an integer
        # column with unique values

        # check if the cost column exists and is numeric. If not, give cost=1 to all PU
        # and warn
        if(! is.numeric(planningUnits$cost)) {
          costOne <- rep (1, rep=nrow(planningUnits))
          planningUnits$cost <- costOne
          warningText <- "planningUnits$cost was not provided or was not numeric, so cost=1 was assigned to each PU."
          cat("\nNote:", warningText, "\n\n")
          warning(warningText)
        }

        # check if the status column exists and is numeric. If not, give status=0 to all PU
        # and warn
        if (! is.numeric(planningUnits$status)) {
          statusZero <- rep (0, rep=nrow(planningUnits))
          planningUnits$status <- statusZero
          warningText <- "planningUnits$status was not provided or was not numeric, so status=0 (available) was assigned to each PU."
          cat("\nNote:", warningText, "\n\n")
          warning(warningText)
        }

      } else {stop("The planningUnits$PU_ID column contains duplicate values.") }
    } else { stop("The planningUnits$PU_ID column does not exist or is not of class integer.") }
  } else { stop("The argument planningUnits is not a valid data frame.") }

  return(planningUnits)
}

validate_occ = function(occ) {
  # this function checks the validity of the occ object

  # occ$TipName is required, so if not valid, execution stops
  # occ$PU_ID is required, so if not valid, execution stops
  # if occ$Amount is not valid, a column of 1 is added with a warning
  if (class(occ) == "data.frame") {
    if (class(Occ$TipName) %in% c("character", "factor", "integer")) {
      if (class(occ$PU_ID) == "integer") {

        # get to here if occ is a valid data frame, TipName is a value column
        # (character, factor or integer) and PU_ID is an integer column.

        # check if the amount column exists and is numeric. If not, give
        # amount=1 to each occurrence and warn.
        if(! is.numeric(occ$Amount)) {
          amountOne <- rep (1, rep=nrow(occ))
          occ$Amount <- amountOne
          warningText <- "occ$Amount was not provided or not numeric, so Amount=1 was assigned to each occurence."
          cat("\nNote:", warningText, "\n\n")
          warning(warningText)
        }

      } else {stop("The occ$PU_ID column does not exist or is not of class integer.") }
    } else { stop("The occ$TipName column does not exist or is not a valid class.") }
  } else { stop("The argument occ is not a valid data frame.") }

  return(occ)
}