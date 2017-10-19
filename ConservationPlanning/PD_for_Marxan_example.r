
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
##                                                                            ##
################################################################################

library(ape)
library(raster)

rm(list=ls())

##### PARAMETERS #####

# location for R code
source.dir <- "C:/Users/u3579238/Work/Software/dan-github/phylospatial-dev/Marxan/"

# location for the input data files
working.dir <- "C:/Users/u3579238/Work/Software/dan-github/phylospatial-dev/Marxan/Sample_Data_Australia/"

# location to write the output files formatted for Marxan
output_dir  <- "C:/Users/u3579238/Work/Software/dan-github/phylospatial-dev/Marxan/Sample_Data_Australia/Marxan_files"

# Inputs
Occ_filename            <- "Occ360x114Mammals_2010_Australia.csv"
PU_shapefile            <- "Australia_360x114.shp"
PU_ID_column            <- "Grid360ID"  # the column name in PU_shapefile which matches the PU_ID column in Occ_filename

phylogeny_filename      <- "Fritz.Resolved.Normal_first5_Aust.nwk"
land_cap                <- TRUE
number_of_trees         <- 1
spf_multiplier          <- 3

target_args <- list(type="bounded_percent",
                    percent=25,
                    minTarget=1,
                    maxTarget=25) # target 25% bounded by minimum 1, maximum 25

##### END OF PARAMETERS #####

StartTime <- date()
cat("\n\n************************************************\n")
cat("Starting Phylo for Marxan at: ", StartTime)
cat("\n************************************************\n")


setwd(source.dir)
source("PrepareMarxanInputs.r")

setwd(working.dir)

# Read and prepare occurrence file
cat("\n\nReading spatial files\n")
Occ   <- read.csv(Occ_filename)
names(Occ) <- c("TipName", "PU_ID", "Amount")

PU.shp          <- shapefile(PU_shapefile)
planning_units  <- PU.shp@data

# create the Marxan file PU.dat
PU_ID_col_index <- which(names(planning_units) == PU_ID_column)
planning_units  <-  data.frame(PU_ID = planning_units[, PU_ID_col_index], cost = planning_units$LandProp, status = 0)

cat("\nReading phylogeny files\n")

# try two different read tree functions
try(
  phy.orig <- read.nexus(file = phylogeny_filename), silent = TRUE
)
if(!exists("phy.orig")) {
  phy.orig <- read.tree(file = phylogeny_filename)
}
cat("\nTree loading finished\n")

if (class(phy.orig) == "multiPhylo") {
  count_all_trees <- length(phy.orig)
} else {
  count_all_trees <- 1
}

# sample from trees
if (count_all_trees > 1) {
  number_of_trees <- min(number_of_trees, count_all_trees)
  trees_to_use <- sample(1:count_all_trees, size=number_of_trees)
  phy <- phy.orig[trees_to_use]
} else {
  phy <- phy.orig
}

#### TEMPORARY FEEDBACK ####
cat("\nRows in Occ:",nrow(Occ),"\n")
cat("\nUnique rows in Occ:",nrow(unique(Occ)),"\n")
cat("\nUnique rows in Occ cols 1-2:",nrow(unique(Occ[,1:2])),"\n")
#### TEMPORARY FEEDBACK ####

target_args <- list(type="bounded_percent",
                    percent=25,
                    minTarget=1,
                    maxTarget=25) # target 25% bounded by minimum 1, maximum 25

get_branch_ranges(trees = phy,
                  planningUnits = planning_units,
                  occ = Occ,
                  output_dir,
                  amount_max = 1,
                  spf_multiplier,
                  target_args)

cat ("\n\n Started  :",StartTime,"\n")
cat ("\n Finished :",date(),"\n")
