#### for intraspecific lineages to generate lineage distribution models.
#### Copyright Dan Rosauer 2016         ####
#### Australian National University     ####
#### September 2012 - November 2016     ####
#### dan.rosauer@anu.edu.au             ####

########################################################################
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #

# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #

# You should have received a copy of the GNU General Public License    #
# along with this program.  If not, see <http://www.gnu.org/licenses/> #
########################################################################

## This script uses a set of species distribution models and a set of points for intraspecific lineages
## to generate lineage distribution models.

## It was rewritten to run in R, from an earlier version with used ArcGIS functions via Python.

## STEPS WHICH THE CODE DOES
## 1. import the points for the whole species
##
## 2. to bound the whole analysis, use euclidian distance to create a grid to define a boundary at a specified distance
##
## 3. load a species distribution model for the whole species, generated before running this script
##
## 4. loop through all of the lineages in the species
##    5a. generate a euc distance layer from sequenced locations for each lineage, bounded by the total species euclidean
##        distance layer from (2)
##     or
##    5b. generate a cost distance layer from sequenced locations for each lineage, using the maxent model to define the cost.
##        Cost = 1 - suitability
##
##    6a. generate a weight layer for each lineage as 1 / distance  from (5)
##     or
##    6b. generate a weight layer for each lineage as 1 / distance^n  from (5) where n = 2, 3 or 4
##
##    7.  set all weights below a threshold to 0, to reduce the effect of distant lineages
##
## 8. sum all of the lineage weight layers
##
## 9. divide each lineage weight layer by the sum of weights (8) so that the weights for each pixel
##    sum to 1 An option in this step, is to exclude lineages from a pixel where they have a low probablility
##    of occurring. Initial models found lineages predicted over a wide area beyond their primary range,
##    but with very low values.  To use this option set the parameters:
##          handle_minor = "threshold"
##          omit_minor_threshold = 0.1
##    as an example, 0.1 means that a lineage with less that 10% of the total of all potential lineages for
##    that species, in the pixel, would be omitted, with the model scaled across the lineages more likely
##    to occur in that pixel.
##    To not use this option, simply set
##          handle_minor = ""
##
## 10. multiply each lineage weight layer by the model likelihood so that the weights for each pixel sum to the original
##    SDM model value for that pixel.

library(raster)
library(sp)
library(stringr)
library(gdistance)

rm(list=ls())

############## START OF PARAMETERS ##############

genera              <- "Cryptoblepharus" # allows script to run for one or more genera

base_dir            <- "C:/Users/u3579238/work/AMT/Models/"
target_dir          <- paste(base_dir, "lineage_models/asc_Rtemp/", sep="")  # where the lineage model grids and working data go

buffer_dist         <- 2.5    # the buffer distance in map units (presumably decimal degrees)
additional_buffer   <- 0      # how much (as a proportion) the output grids should extend beyond the buffered points
grid_resolution     <- 0.01
spRef               <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

output_model_extent        <- extent(112.9, 153.64, -25, -9)    # the maximum extent for all lineage models
lineage_field_name  <- "lineage_from_mtDNA"   # the column for lineage name in the site data
distance_method     <- "model-cost"           # determines whether distance is calculated as euclidean or model-weighted cost distance
## so far, can be "euclidian" or "model-cost"
weight_function     <- "inverse_cube"    ## determines whether lineage weight is calculated as 1/distance or 1/(distance^2), or simply closest distance
## so far, can be "inverse" or "inverse_square" or "cost_allocation"
min_dist_value      <- grid_resolution/2  ## remove as a parameter, once working
min_SDM_value       <- 0.005  # values below this for the SDM are set to zero, to simplify calculation in areas of essentially unsuitable habitat
##   but keep current value for consistency in this study
min_weight_threshold <- 0.02        ## weights below this for any layer are set to 0.  If the value here is 0, then no threshold is applied
scale_to            <- "model"      ## determines whether lineage weights within a model group sum to the SDM value or to 1
## can be "model" or "one"

handle_minor        <- "threshold"  # if handle_minor = 'threshold', then lineages with < than the specified proportion of the lineage sum
omit_minor_threshold <- 0.1         # for that cell, are set to 0

skip_distance_layers <- FALSE       # skip creating the distance layers - they are already done.
                                    # THIS OPTION IS ONLY TO SAVE TIME DURING DEBUGGING

named_species       <- c("metallicus")          # a changeable list to allow for species in the dataset to be skipped
use_list            <- "do"           # specify what to do with the named species:
  #do - the named species (use_list="do")
  #skip - the named species (use_list="skip")
  #do all the species in the data and ignore the named species list (use_list="" or anything else);

lin_exclude_list    <- c()   # this list allows for skipping at the lineage level

############### END OF PARAMETERS ###############

cat ("\n\n*************************************** \n Lineage Distribution Estimation Tool ")
cat ("\n    Dan Rosauer \n    September 2012 - November 2016")
cat ("\n***************************************\n")

SDM_model_base     <- paste(base_dir, "species_models/maxent/", sep="")

time_start <- Sys.time()

# loop through each genus in the vector genera
for (genus in genera) {

  cat("\nGenus:", genus, "\n")

  # Load the sequence site data
  cat ("Loading the lineage locations...\t")
  lineage_site_filename <- paste(base_dir, "species_sites/", genus, "_sites.csv", sep="")

  all_lineages_sites <- read.csv(lineage_site_filename)
  orig_rowcount <- nrow(all_lineages_sites)

  # filter the records
  not_sequenced <- grep("not_sequenced", all_lineages_sites[, lineage_field_name])  # remove records with 'not_sequenced' in the lineage name
  if (length(not_sequenced >0)) {all_lineages_sites <- all_lineages_sites[-not_sequenced, ]}
  all_lineages_sites <- all_lineages_sites[which(all_lineages_sites$Use %in% c(-1,1)), ]    # remove records where 'use' <> -1 or 1
  all_lineages_sites <- all_lineages_sites[which(is.finite(all_lineages_sites[, "latitude"]) & is.finite(all_lineages_sites$longitude)), ]    # remove without a numeric latitude and longitude

  all_lineages_sites$model_group <- str_trim(all_lineages_sites$model_group)
  all_lineages_sites[, lineage_field_name] <- str_trim(all_lineages_sites[, lineage_field_name])

  groupLineageList  <- unique(all_lineages_sites[, 1:2])
  groups            <- unique(all_lineages_sites$model_group)
  groupLineages     <- unique(all_lineages_sites[, lineage_field_name])

  cat (orig_rowcount, "rows read,", nrow(all_lineages_sites), "valid records loaded\n")

  # turn lineage points into a SpatialPoints object
  allLineagePoints <- SpatialPointsDataFrame(all_lineages_sites[, c("longitude", "latitude")], data=all_lineages_sites, proj4string = CRS(spRef))
  allLineagePoints <- crop(allLineagePoints, output_model_extent)
  # python writes to file at this point, but prob not needed

  # restrict groups to particular species based on the named_species parameter
  if (use_list == "do") {
    groups <- intersect(groups, named_species)
  } else if (use_list == "skip") {
    groups <- setdiff(groups, named_species)
  }

  # print a list of model groups
  cat ("\nModel groups to do in", genus, "\n")
  for (group in groups) {
    cat("\t", group, "\n")
  }

  # start looping through the model groups for this genus
  for (group in groups) {
    group <- as.character(group)
    if (group == "0" | group == "") {next}

    cat ("\nStarting group", genus, group, "on", date(), "\n")

    # load the SDM model raster
    SDM_model <- paste(SDM_model_base, genus, "/", str_replace(group," ","_"), "_median.asc", sep="")
    SDM.ras <- raster(SDM_model)
    projection(SDM.ras) <- CRS(spRef)
    group_extent <- extent(SDM.ras)

    # get a list of the lineages in this group
    lineages <- groupLineageList[which(groupLineageList$model_group==group), lineage_field_name]
    if (length(lin_exclude_list) > 0) {
      lineages <- setdiff(lineages,lin_exclude_list)
    }

    cat ("\nLineages in", group, ":")
    for (lineage in lineages) {
      cat ("\n   ", lineage)
    }

    # # set the environment
    # env.snapRaster  = maxent_model
    # env.mask        = maxent_model
    # env.extent      = maxent_model

    if (length(lineages) > 1) { # proceed with lineage models if there are multiple
                                # lineages - otherwise just copy the SDM for the model group

      thisGroupPoints <- allLineagePoints[allLineagePoints$model_group == group, ]
      thisGroupPoints <- crop(thisGroupPoints, group_extent)

      # calculate a new extent
      group_points_extent <- extent(thisGroupPoints)
      buffer_ratio  <- 1 + additional_buffer
      extent_buffer <- buffer_dist * buffer_ratio

      # new extent is the same as points layer + a buffer, but where the extended buffer
      # goes beyond the extent of the maxent model, limit to the output model extent.
      xmin <- max(group_points_extent@xmin - extent_buffer, output_model_extent@xmin)
      ymin <- max(group_points_extent@ymin - extent_buffer, output_model_extent@ymin)
      xmax <- min(group_points_extent@xmax + extent_buffer, output_model_extent@xmax)
      ymax <- min(group_points_extent@ymax + extent_buffer, output_model_extent@ymax)
      thisGroupExtent <- extent(xmin, xmax, ymin, ymax)

      ### generate a weight grid for each lineage  START OF STEP 4

      #SDM.ras[SDM.ras < min_SDM_value] <- min_SDM_value
      if (distance_method == "model-cost") {
        cat("\n\nCreating a transition matrix based on the SDM for", group, "\n")

        model_cost.ras  <- -1 * log(SDM.ras)     # this is the original version of model cost
        model_trans.ras <- 1 / model_cost.ras  # but using the inverse here, to fit the accCost
          # function which is based on a transition matrix. For accCost() the cost of moving between
          # cells = 1 / transition value.
        trans     <- transition(model_trans.ras, transitionFunction=mean, directions = 8)
        trans     <- geoCorrection(trans, type="c", multpl=F)
      }

      cat ("\nLooping through the lineages in group", group, "to generate weight grids\n")
      count <- 0

      for (lineage in lineages) {

        count <- count +1

        thisLineagePoints <- thisGroupPoints[thisGroupPoints@data[, lineage_field_name] == lineage, ]

        # create a cost distance layer for the current lineage
        if (lineage == "0") {
          cat ("\nCreating distance layer for sequenced locations of unknown lineage")
        } else {
          cat ("\nCreating distance layer for lineage", lineage)
        }

        if (distance_method == "model-cost") {                                   ## STEP 5b
          # calculates the least cost distance to the nearest lineage point
          # the result is written directly to lineage_dist_gridname
          lin_dist.ras = accCost(trans, thisLineagePoints)

        } else {
          thisLineagePoints.ras <- rasterize(thisLineagePoints, SDM.ras, 1)
          lin_dist.ras <- distance(thisLineagePoints.ras) / 1000     ## STEP 5a
        }

        # change zero values to a very small non-zero value, to avoid nodata in division
        lin_dist.ras[lin_dist.ras < min_dist_value] <- min_dist_value

        if (weight_function == "inverse_square") {                 ## STEP 6b
          lin_weight.ras <- 1 / (lin_dist.ras ^ 2)
        } else if (weight_function == "inverse_cube") {
          lin_weight.ras <- 1 / (lin_dist.ras ^ 3)
        } else if (weight_function == "inverse_quad") {
          lin_weight.ras <- 1/(lin_dist.ras ^ 4)
        } else {
          lin_weight.ras <- 1/lin_dist.ras
        }

        # add the results to stack of distance layers for this group
        if (count == 1) {
          weight.stack <- stack(lin_weight.ras)
        } else {
          weight.stack <- stack(weight.stack, lin_weight.ras)
        }
        names(weight.stack)[count] <- lineage
      }

      weight_sum.ras <- sum(weight.stack)
      model.stack <- weight.stack / weight_sum.ras
      model.stack <- model.stack * SDM.ras
      names(model.stack) <- names(weight.stack)

    } else {
      SDM.ras[SDM.ras < min_SDM_value] <- min_SDM_value
      model.stack <- stack(SDM.ras)
      names(model.stack)[1] <- lineages[1]
    }

    # write results to file for this group
    model_names <- names(model.stack)

    cat("\n")

    for (i in 1:nlayers(model.stack)) {
      model.ras <- model.stack[[i]]
      model_filename <- paste("lin_model_", genus, "_", group, "_", model_names[i], ".asc", sep="")
      model_path <- paste(target_dir, model_filename, sep="")
      cat("writing model ascii for:", model_filename, "\n")
      writeRaster(model.ras, model_path, overwrite=T)
    }
    cat ("\nLineage models done for", genus, "\n")
    time_diff <- difftime(Sys.time(), time_start, units='mins')
    cat ("\nTime elapsed:", round(time_diff,2), "minutes\n")

  }

}
