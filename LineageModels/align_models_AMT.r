rm(list=ls())
library(raster)

#define directories
base.dir <- 'Your working location/Models/'   # modify to the base directory for your lineage modelling

input.dir    = paste(base.dir, 'lineage_models/asc/', sep ='')          # location of existing lineage distribution models
output.dir   = paste(base.dir, 'lineage_models/asc_aligned/', sep='')   # location to save aligned lineage distribution models
template_ext = paste(base.dir, 'species_models/maxent/Heteronotia/binoei_median.asc', sep='')   # an .asc grid with the extent to which all models will be cropped and aligned.

new_only <- TRUE  # if true, skip grids which are already in the output directory
file.pattern    <- '*.asc$'  #regex

setwd(input.dir)

input_files = list.files(path=input.dir, pattern=file.pattern, full.names=FALSE, recursive=FALSE, ignore.case=TRUE, include.dirs=FALSE)
output_files= list.files(path=output.dir, pattern=file.pattern, recursive=FALSE, ignore.case=TRUE, include.dirs=FALSE)

template.ras = raster(template_ext)
new_extent =      extent(template.ras)

raster_names <- ""

for (tfile in input_files) {
  filepath=paste(input.dir,tfile,sep='')
  outname = paste(output.dir,tfile,sep="")

  if ((!tfile %in% output_files) | (! new_only)) {
    grid.ras = raster(tfile)
    grid_ext.ras = extend(grid.ras,new_extent,value=0) # extend to the union of current grid and new extent
    grid_ext.ras = crop(grid_ext.ras,new_extent) # crop back to new extent
    writeRaster(grid_ext.ras,outname,overwrite=TRUE, NAflag=-9999)
    cat("\nExtended asc written for",tfile)

    # make a vector of the layer names
    if (raster_names == "") {
      raster_names <- outname
    } else {
      raster_names <- c(raster_names,outname)
    }

  } else {
    cat("\nSkipped",tfile)
  }
}

# now make a raster stack
lin.stack <- stack(raster_names)
writeRaster(lin.stack, "lin_models.stack")

setwd(output.dir)

maxval  <- stackApply(lin.stack,rep(1,nlayers(lin.stack)),fun=max,filename="max_val.asc")
sum     <- stackApply(lin.stack,rep(1,nlayers(lin.stack)),fun=sum,filename="sum.asc")
maxprop <- maxval / sum
writeRaster(maxprop,"max_prop.asc")

maxlin  <- which.max(lin.stack)
writeRaster(maxlin,"max_lin.asc")
stack_names <- data.frame(cbind(1:nlayers(lin.stack),names(lin.stack)))
names(stack_names) <- c("layer_num","layer_name")
write.csv(stack_names,"layer_names")
