#### Dan Rosauer                        ####
#### Australian National University     ####
#### September 2012 - October 2015      ####
#### dan.rosauer@anu.edu.au             ####

# This script has some utility functions used in the spatial analyses

def getColumnNumber (header_row, column_name_options):
    # returns the index for the item which matches any of the alternative column names
    # normal use is to find which column has a given name (eg could be lat, Lat, latitude etc)
    # if there are multiple matches, just the first is returned
    
    col_index = -1
    if len(column_name_options) == 1:
        col_index = header_row.index(column_name_options[0])
    else:
        for name_option in column_name_options:
            col_index = header_row.index(name_option)
            if col_index >= 0:
                break
    return col_index


def getFieldValues(shapefile, target_field):
    ##returns the number of values in an attribute table
    import arcpy, sys, os

    rows=arcpy.SearchCursor(shapefile,"","",target_field)
    value_list =[]
    for row in rows:
        this_value = row.getValue(target_field)
        if this_value not in value_list:
            value_list.append(this_value)        
    return value_list
        
def getFieldMinMax(layer, target_field):
    ##returns the minimum and maximum values in a numeric field of an attribute table
    import arcpy.analysis

    #fields = [target_field]
    rows=arcpy.da.SearchCursor(layer,(target_field))
    #value_list =[row[0] for row in rows]
    
    value_list=[]
    for row in rows:
        #this_value = row.getValue(target_field)
        this_value = row[0]
        if this_value not in value_list:
            if this_value != None:
                value_list.append(this_value)
    minval=min(value_list)
    maxval=max(value_list)
    return [minval,maxval]
