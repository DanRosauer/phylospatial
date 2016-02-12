#### Dan Rosauer                        ####
#### Australian National University     ####
#### September 2012 - February 2016     ####
#### dan.rosauer@anu.edu.au             ####

## This script uses a set of species distribution models and a set of points for intraspecific lineages
## to generate lineage distribution models.

## It uses spatial functions in the program ArcGIS, so an ArcGIS Installation and license are required.
## A truly open source version (in R not Python) is envisaged, but not close to happening yet.

## STEPS WHICH THE CODE DOES
## 1. import the points for the whole species
##
## 2. to bound the whole analysis, use euclidian distance to create a grid to define a boundary at a set distance
##
## 3. load a species distribution model for the whole species, generated before running this script
##
## 4. loop through all of the lineages in the species
##    5a. generate a euc distance layer from sequenced locations for each lineage, bounded by the total species euclidean distance layer from (2)
##     or 
##    5b. generate a cost distance layer from sequenced locations for each lineage, using the maxent model to define the cost.  Cost = 1 - suitability
##    
##    6a. generate a weight layer for each lineage as 1 / distance  from (5)
##     or
##    6b. generate a weight layer for each lineage as 1 / distance^2  from (5)
##
##    7.  set all weights below a threshold to 0, to reduce the effect of distant lineages
##
## 8. sum all of the lineage weight layers
##
## 9. divide each lineage weight layer by the sum of weights (8) so that the weights for each pixel sum to 1
##      An option in this step, is to exclude lineages from a pixel where they have a low probablility of occurring.  Initial models found
##      lineages predicted over a wide area beyond their primary range, but with very low values.  To use this option set the parameters:
##          handle_minor = "threshold"
##          omit_minor_threshold = 0.1  - as an example, 0.1 means that a lineage with less that 10% of the total of all potential lineages
##                                      for that species, in the pixel, would be omitted, with the model scaled across the lineages more
##                                      likely to occur in that pixel
##      To not use this option, simply set
##          handle_minor = ""
##
## 10. multiply each lineage weight layer by the model likelihood so that the weights for each pixel sum to the model likelihood.

## All of the above could be nested within a loop that iterates through species.

############### PARAMETERS ###############
import arcpy
import sys, os, math, numpy, csv, os, string, datetime
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

sys.path.append("FILE PATH TO SpatialFunctions.py")
from SpatialFunctions import *

### PARAMETERS ###

genus_list = ["Demo"] # allows script to run for one or more genera

base_dir = "FILE PATH TO ROOT OF MODEL DIRECTORY STRUCTURE\\Models\\"
target_location = base_dir + "lineage_models\\"  # where the lineage model grids and working data
output_gdb_name = "lineage_models.gdb"
scratch_workspace = base_dir  #scratch workspace is used by ArcGIS for temporary files during analysis
scratch_gdb_name = "scratch.gdb"
export_asc = True
asc_target_location = base_dir + "lineage_models\\asc\\"  # assumes that output lineage models are wanted as ascii grids

buffer_dist = 2.5      # the buffer distance in map units (presumably decimal degrees)
additional_buffer = 0  ## how much (as a proportion) the output grids should extend beyond the buffered points
grid_resolution = 0.01
spRef = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"

Model_extent = arcpy.Extent(112.9,-25,153.64,-9)    # the maximum extent for all lineage models
Lineage_field_name = "lineage_from_mtDNA"           # the column for lineage name in the site data    
Distance_method = "model-cost"      ## determines whether distance is calculated as euclidean or model-weighted cost distance
                                    ## so far, can be "euclidian" or "model-cost"
Weight_function = "inverse_cube"    ## determines whether lineage weight is calculated as 1/distance or 1/(distance^2), or simply closest distance
                                    ## so far, can be "inverse" or "inverse_square" or "cost_allocation"
Min_dist_value = grid_resolution/2  ## remove as a parameter, once working
                                    ##   but keep current value for consistency in this study
Min_weight_threshold = 0.02         ## weights below this for any layer are set to 0.  If the value here is 0, then no threshold is applied
Scale_to = "model"                  ## determines whether lineage weights within a model group sum to the model suitability or to 1
                                    ## can be "model" or "one"

handle_minor = "threshold"          ## if handle_minor = 'threshold', then lineages with < than the specified proportion of the lineage sum 
omit_minor_threshold = 0.1          ## for that cell, are set to 0

skip_distance_layers = False         ## skip creating the distance layers - they are already done.  THIS OPTION IS ONLY TO SAVE TIME DURING DEBUGGING

# a changeable list to allow for species in the dataset to be skipped
named_species   = []
use_list        = ""  #specify whether to:
                        #do - the named species (use_list="do")
                        #skip - the named species (use_list="skip")
                        #do all the species in the data and ignore the named species list (use_list="" or anything else);

Lin_exclude_list = []   # this list allows for skipping at the lineage level

############### END OF PARAMETERS ###############

try:
    print "Lineage Distribution Estimation Tool"
    print "Dan Rosauer - October 2015\n"
    
    for genus in genus_list:
        print "\nGenus: " + genus + "\n"
    
        lineage_site_filename = base_dir + "species_sites\\" + genus + "_sites.csv"
        maxent_model_base = base_dir + "species_models\\maxent\\maxent_models.gdb"
    
        # create the scratch geodatabase if needed (to store temporary layers during analysis)
        scratch_gdb_path = scratch_workspace + scratch_gdb_name
        if not arcpy.Exists(scratch_gdb_path):
            arcpy.CreateFileGDB_management(scratch_workspace, scratch_gdb_name)
        
        # create the output geodatabase if needed
        if not os.path.exists(target_location):
            os.makedirs(target_location)
        target_location_ESRI = string.replace(target_location,"\\","/")
        if not arcpy.Exists(target_location_ESRI + output_gdb_name):        
            arcpy.CreateFileGDB_management(target_location_ESRI,output_gdb_name)
            print "Geodatabase created to store lineage models"
            print target_location_ESRI,output_gdb_name
            
        # remove the file schema.ini which causes problems for reading .csv files correctly
        if arcpy.Exists(base_dir + "species_sites\\" +'schema.ini'):
            os.remove(base_dir + "species_sites\\" +'schema.ini')
        
        # Load the sequence site data
        print "Loading the lineage locations...\t",
        
        # make each column being used, into a list
        # and also create lists of unique ModelGroups and Lineages
        with open(lineage_site_filename, 'rb') as csvfile:
            sequence_csv = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_MINIMAL)
            rownum = 0
            recnum = 0
            ModelGroup=[]
            Lineage = []
            Lat=[]
            Long=[]
            GroupList=[]
            GroupLineageList =[]
            usecol = -1
            
            for row in sequence_csv:
    
                if rownum == 0:
                    header = row
                    rownum += 1
                    # find column numbers
                    usecol = getColumnNumber(row, ['Use', 'use'])
                    lat_col  = getColumnNumber(row, ['latitude', 'Lat', 'lat', 'Latitude'])
                    long_col = getColumnNumber(row, ['longitude', 'Long', 'long', 'Longitude'])
                    model_group_col = getColumnNumber(row, ["model_group"])
                    lineage_col = getColumnNumber(row, [Lineage_field_name])
    
                else:
                    rownum += 1
                    if string.find(row[lineage_col], 'not_sequenced') >=0: continue # skip records where lineage name includes 'not_sequenced'
                    
                    if (usecol == -1 or row[usecol] == '1'):
                        try:        # if the lat or long can't be converted to a number, then skip that row
                            Lat.append(float(row[lat_col]))
                            try:
                                Long.append(float(row[long_col]))
                                # code gets to here for valid lat and long, so other steps can go here too
                                ModelGroup.append(row[model_group_col])
                                Lineage.append(row[lineage_col])
                                if row[model_group_col] not in GroupList:
                                    GroupList.append(row[model_group_col])
                                if row[0:4] not in GroupLineageList:
                                    GroupLineageList.append(row[0:4])
                                #print row[model_group_col] + "\t" + row[lineage_col] + "\t" + str(rownum)  # FOR DEBUGGING ONLY
                                recnum +=1
                            except:
                                pass
                        except:
                            pass
            # so now we have a list for each column, excluding where 'Use' is not 1 and rows with null coordinates
        print str(rownum-1) + " rows read, " + str(recnum) + " records loaded\n"
        
        # set the geoprocessing environment
        env.workspace  = target_location_ESRI + output_gdb_name
        env.scratchWorkspace = scratch_workspace + scratch_gdb_name
        env.extent = Model_extent
        
        # Make XY Event Layer
        points_layer = "sequenced_sites"
        layer_result = arcpy.MakeXYEventLayer_management(lineage_site_filename, header[long_col], header[lat_col], points_layer, spRef)
        points_layer = layer_result[0]
        print "Point layer generated for " + genus + "\n"
        
        # export the layer as a feature class
        outLocation = env.workspace
        outFeatureClass = genus + "_lin_points"
        
        try:
            arcpy.env.overwriteOutput=True
            arcpy.CopyFeatures_management(points_layer, outFeatureClass, "", "0", "0", "0")
        except:
            print "\n" + outFeatureClass + " could not be created."
            print arcpy.GetMessages()
        points_fc = outLocation + "/" + outFeatureClass
        
        # restrict the GroupList to particular species based on the names_species parameter
        if use_list == "do":
            GroupList = list(set(named_species).intersection(set(GroupList)))
        elif use_list == "skip":
            GroupList = list(set(GroupList).difference(set(named_species)))
            
        # print a list of model groups
        print "Model groups to do in", genus
        for group in GroupList:
            print "   ", group
        
        # Loop through the list of Model Groups
        for group in GroupList:
            if group != 0 and group != "":
            
                print "\nStarting group " + genus + " " + group + " at " + datetime.datetime.now().strftime("%I:%M %p") + "\n"
                
                maxent_model = maxent_model_base + "\\" + genus + "_" + string.replace(group," ","_")
                    
                ## get a list of the lineages in this group
                lineage_list=[]
                print "Lineages in " + group + ":"
                for row in GroupLineageList:
                    if row[lineage_col] not in lineage_list:
                        if row[model_group_col] == group and "," not in row[lineage_col] and "not_sequenced" not in row[lineage_col] and row[lineage_col] not in Lin_exclude_list:  # exclude particular lineage names and those with dodgy punctuation (fix in data later)
                            lineage = row[lineage_col]
                            lineage = lineage.strip()
                            lineage_list.append(lineage)
                            print "   ", row[lineage_col]
            
                # start a list of layers to delete at the end
                layers_to_delete = []
                
                # set the environment
                env.snapRaster  = maxent_model
                env.mask        = maxent_model
                env.extent      = maxent_model
                
                maxent_raster = arcpy.sa.Raster(maxent_model)
                    
                if len(lineage_list) > 1:  # proceed with lineage models if there are multiple lineages - otherwise just copy the SDM for the model group
                    
                    # define spatial data for this group
                    groupDefQuery = "[ModelGroup] = '" + group + "'"
                    points_layer.definitionQuery = groupDefQuery
        
                    maxent_extent = maxent_raster.extent
                    
                    # get the extent of the points for the model group
                    yrange = getFieldMinMax(points_fc,header[lat_col])
                    ymin   = yrange[0]
                    ymax   = yrange[1]
                    xrange = getFieldMinMax(points_fc,header[long_col])
                    xmin   = xrange[0]
                    xmax   = xrange[1]
            
                    buffer_ratio = (1 + additional_buffer)
                    extent_buffer = buffer_dist * buffer_ratio
                
                    # new extent is the same as points layer + a buffer
                    # but where the extended buffer goes beyond the extent of the maxent model, limit to the model extent (determined by the buffered environemnt grids)
                    xmin=math.floor(max([xmin - extent_buffer,maxent_extent.XMin]))
                    ymin=math.floor(max([ymin - extent_buffer,maxent_extent.YMin]))
                    xmax=math.ceil(min([xmax + extent_buffer,maxent_extent.XMax]))
                    ymax=math.ceil(min([ymax + extent_buffer,maxent_extent.YMax]))
                    env.extent = arcpy.Extent(xmin,ymin,xmax,ymax)
                
                    ### generate a weight grid for each lineage  START OF STEP 4
                    
                    ## get a selectable layer for the sequenced sites
                    lin_lyr = arcpy.MakeFeatureLayer_management(points_fc,"lineage_layer")[0]
                
                    print "\nLooping through the lineages in " + group + " to generate weight grids\n"
                    count = 0
                    
                    if Distance_method == "model-cost":
                        model_cost = -1 * arcpy.sa.Ln(maxent_raster)
        
                    for lineage in lineage_list:
                        # if lineage != 'planD':
                        #     continue
                        
                        count += 1
                    
                        where_clause = '"' + header[model_group_col] + '" = ' + "'" + group + "' and " + '"' + Lineage_field_name + '"' + " = " + "'" + lineage + "'"
                        arcpy.SelectLayerByAttribute_management(lin_lyr, "NEW_SELECTION", where_clause)
                        arcpy.CopyFeatures_management(lin_lyr, "lineage_points", "", "0", "0", "0")
                        layers_to_delete.append("lineage_points")
                    
                        # create a distance layer for the current lineage
                        if str(lineage) == "0":
                            print "Creating distance layer for sequenced locations of unknown lineage"
                        else:
                            print "Creating distance layer for lineage " + lineage
                        lineage_weight_gridname = "lin_wt_" + string.replace(group," ","_") +"_" + lineage
                        lineage_weight_gridname = lineage_weight_gridname.rstrip()
                        
                        if Distance_method == "model-cost":                                   ## STEP 5b
                            ## calculates the least cost distance to the nearest lineage point
                            ## the result is written directly to lineage_dist_gridname
                            lin_dist = arcpy.sa.PathDistance("lineage_points",model_cost)
                            # change zero values to a very small non-zero value, to avoid nodata in division
                            lin_dist = arcpy.sa.Con(lin_dist==0,0.0001,lin_dist)
                            
                            if Weight_function == "inverse_square":                 ## STEP 6b
                                lin_weight = 1/(lin_dist ** 2)
                            elif Weight_function == "inverse_cube":
                                lin_weight = 1/(lin_dist ** 3)
                            elif Weight_function == "inverse_quad":
                                lin_weight = 1/(lin_dist ** 4)
                            else:
                                lin_weight = 1/lin_dist                             ## STEP 6a       this comes 2nd, as it is the default, for any other values of Weight_function            
                    
                        else:
                            lin_dist = arcpy.sa.EucDistance(lin_lyr,"",grid_resolution)     ## STEP 5a
                            # change zero values to a very small non-zero value, to avoid nodata in division
                            lin_dist = arcpy.sa.Con(lin_dist > Min_dist_value,Min_dist_value,lin_dist)
                            
                            if Weight_function == "inverse_square":                 ## STEP 6b
                                lin_weight = 1/(lin_dist ** 2)
                            elif Weight_function == "inverse_cube":
                                lin_weight = 1/(lin_dist ** 3)
                            elif Weight_function == "inverse_quad":
                                lin_weight = 1/(lin_dist ** 4)
                            else:
                                lin_weight = 1/lin_dist                             ## STEP 6a       this comes 2nd, as it is the default, for any other values of Weight_function            
                       
                        # reset points definition to whole Model Group
                        arcpy.SelectLayerByAttribute_management(lin_lyr, "CLEAR_SELECTION")
                    
                        # apply a threshold to each weight grid                 ## STEP 7
                        if Min_weight_threshold > 0:
                            where_clause = '"VALUE" >= ' + str(Min_weight_threshold)
                            lin_weight = arcpy.sa.Con(lin_weight, lin_weight, 0, where_clause)
                            
                        # set NoData values to 0
                        lin_weight=arcpy.sa.Con(arcpy.sa.IsNull(lin_weight),0,lin_weight)
                    
                        lin_weight.save(lineage_weight_gridname)  ## this layer should be kept until the final weights are calculated
                        
                        #lineage_weight_gridname = arcpy.sa.Con(lineage_weight_gridname, 0, lineage_weight_gridname, lineage_weight_gridname = NULL)
    
                        layers_to_delete.append(lineage_weight_gridname)
                    
                        # calculate the sum of weights for each pixel to scale values later  ## STEP 8
                        if count == 1:
                            weight_sum = lin_weight
                        else:
                            weight_sum = weight_sum + lin_weight
                    
                    print "\nDistance layers done.\n"
                    
                    
                    # remove very small values for each pixel and rescale so the rest add up to the species model
                    env.mask = maxent_model
                   
                    if handle_minor == "threshold":
                        count = 0
                        for lineage in lineage_list:                                ## STEPS 9 and 10
                            count += 1
                            if str(lineage) != "0":  # lineage 0 is used to refer to sequenced locations without a named lineage.
                                print "Removing lineage values which account for less than " + (str(omit_minor_threshold*100)) + "% of pixel score. " + lineage
                                lineage_weight_gridname = "lin_wt_" + string.replace(group," ","_") +"_" + lineage 
                                lin_weight = arcpy.sa.Raster(lineage_weight_gridname)
                                
                                lin_proportion = lin_weight / weight_sum
                                where_clause = "VALUE > " + str(omit_minor_threshold)
                                lin_weight = arcpy.sa.Con(lin_proportion,lin_weight,0,where_clause)
                                lin_weight.save(lineage_weight_gridname)
                            if count == 1:
                                new_weight_sum = lin_weight
                            else:
                                new_weight_sum = new_weight_sum + lin_weight
                        weight_sum = new_weight_sum
                        layers_to_delete.append(lin_proportion)
                        layers_to_delete.append(new_weight_sum)
                        print "\nFinished removing low scores\n"
    
                    # calculate the scaled weight for each lineage, so they sum to a) 1 or b) the model suitability
                    for lineage in lineage_list:                                ## STEPS 9 and 10
                        if str(lineage) != "0":  #lineage 0 is used to refer to sequenced locations without a named lineage.
                            print "creating final scaled weight grid for lineage " + lineage
                            lineage_weight_gridname = "lin_wt_" + string.replace(group," ","_") +"_" + lineage 
                            lin_weight = arcpy.sa.Raster(lineage_weight_gridname)
                        
                            if Scale_to == "model":
                                lin_weight_scaled = (lin_weight / weight_sum) * maxent_raster
                            else:
                                lin_weight_scaled = (lin_weight / weight_sum)
                        
                            lineage_scaled_weight_name = "lin_model_" + genus + "_" + string.replace(group," ","_") +"_" + lineage
                            lin_weight_scaled.save(lineage_scaled_weight_name)  ## THIS is a final layer to keep!
                    
                            if export_asc:
                                asc_filename = asc_target_location + lineage_scaled_weight_name + ".asc"
                                arcpy.RasterToASCII_conversion(lin_weight_scaled, asc_filename)
                                arcpy.DefineProjection_management(asc_filename, spRef)
                    
                        else:   # lineage 0 is not an actual lineage, but is used to refer to sequenced locations without a named lineage.
                                # a final model is not created for lineage 0, but it should affect the values for the other lineages.
                            print "Sequenced locations of unnamed lineage were not modelled, but do affect values for other lineages."
                
                # where the model group has just a single lineage, simply copy the maxent model as the lineage model
                else:
                    lineage = lineage_list[0]
                    lineage_scaled_weight_name = "lin_model_" + genus + "_" + string.replace(group," ","_") +"_" + lineage
                    maxent_raster.save(lineage_scaled_weight_name)  # saving the maxent model for the model_group as the lineage model
                    print "created lineage grid for lineage " + lineage + " as a copy of the model for group: " + group
                    
                    if export_asc:
                        asc_filename = asc_target_location + lineage_scaled_weight_name + ".asc"
                        arcpy.RasterToASCII_conversion(maxent_raster, asc_filename)
                        arcpy.DefineProjection_management(asc_filename, spRef)                
                    
                print "\nAnalysis for " + group + " completed - now deleting temporary data."
    
                # and finally, delete temporary layers
                layers_to_delete = list(set(layers_to_delete)) # remove duplicates
                for layer in layers_to_delete:
                    try:
                        arcpy.Delete_management(layer)
                        print "   ", layer + " deleted"
                    except:
                        print "   ", layer + " NOT deleted"
            
            print "\n   **************************\n   * FINISHED Model Group :", group, "at", datetime.datetime.now().strftime("%I:%M %p on %B %d, %Y"), "*\n   **************************\n"
                    
        env.workspace = scratch_workspace + scratch_gdb_name
        datasets = arcpy.ListDatasets()
        for dataset in datasets:
            try:
                arcpy.Delete_management(dataset)
                print dataset + " deleted"
            except:
                print dataset + " NOT deleted"
                
        print "\nLineage models done for " + genus + "\n"
            
    print "\nLineage models done!\n" + datetime.datetime.now().strftime("%I:%M %p on %B %d, %Y") + "\n"

except:
    print "Unexpected error:", sys.exc_info()[0]
    raise
finally:
    os.system('pause')
