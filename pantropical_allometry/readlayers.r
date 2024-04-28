##################################################
##### Manipulate raster and netCDF datasets ######
##### J Chave 15 Dec 2013                   ######
##### M Rejou-Mechain & J Chave 25 Mar 2014 ######
##### Written in R version 3.0.3            ######
##### 5feb2016 minor fix
##################################################

#=========================================================================
## In Chave et al. E and CWD are environmental variables
## CWD is long-term climatic water deficit (see Chave et al for a precise definition)
## TS is temperature seasonality extracted from the Worldclim dataset
## PS is the precipitation seasonality extracted from the Worldclim dataset
## E is an environmental stress variable defined in Eq (6b) of Chave et al as:
## E=(0,178*TS-0.938*CWD-6.61*PS)/1000
##
## The layers are provided at 2.5 arc-minute resolution, or ca. 5 km
## This 24 cells per degree, or 0.041666666666667 degree per cell
# For 'raster' beginners:
# raster starts at top left corner; Latitude is from +90 to -60, longitude from -180 to +180 
# thus, French Guiana is at (+5,-55), or 85 rows from top and (180-55)=125 cols from left.
# grain size is 0.041666666666667 degree, hence there are 24 cells per degree, or about 1 cell every 5 km
#=========================================================================

## filename is either CWD or E
## format 'nc' denotes a file in netCDF format  (default; much faster in linux-based environments)
## format 'bil' denotes a file in raster format
## default interpolation method used here is 'bilinear', where the exact value at the coordinates
## is interpolated along both the x and y axes. The other option is 'simple' (the cell value is retrieved)

retrieve_raster=function(filename,coord,plot=F,format="nc"){
  require("raster")
  if(format=='nc') require("ncdf4")
  zipurl <- "http://chave.ups-tlse.fr/pantropical_allometry"
  if(format=='nc') zipurl = paste(zipurl,"/",filename,".nc.zip",sep="")
  if(format=='bil') zipurl = paste(zipurl,"/",filename,".bil.zip",sep="")
  # Read the raster file in netCDF format from 
  # http://chave.ups-tlse.fr/pantropical_allometry.htm
  DEMzip <- download.file(zipurl, destfile = "zipdir")
  unzip("zipdir", exdir = "unzipdir")
  nam=paste("unzipdir/",filename,".",format,sep="") ## thxs to Rebbeca Senior @ Sheffield
  RAST <- raster(nam)
  # Check that the dataset has properly been imported and that 
  # your coordinates are correct
  # This step takes time
  if(plot==T){
    plot(RAST)
    points(coord,pch="x")
  }
  # Extract the raster value
  # coord=cbind(longitude,latitude)
  RASTval=extract(RAST,coord,method="bilinear")
  return(RASTval)
}

