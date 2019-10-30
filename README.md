# STAdmix
Spatio-temporal analysis of admixture components

These are scripts and commands used for the analysis in Racimo et al. 'A geostatistical approach to modelling human Holocene
migrations in Europe using ancient DNA' which is available in bioRxiv. The full list of analyses is in the file 'SpatioTemporal_Notes.txt'.

Below is a brief tutorial for performing spatio-temporal kriging of ancestry data:

First, we load the necessary functions and libraries in R:

```
R
source("STKrigFunc.R")
```

We then load the land feature data, which was obtained from https://www.naturalearthdata.com:

```
landfeatures <- "NaturalEarth/ne_50m_land/ne_50m_land.shp"
projectiontype <- "+init=epsg:4326"
land <- readOGR(landfeatures)
```

Then, we load the labels for all individuals (a fam file in plink format), the ancestry data (a Q matrix from Ohana or other admixture inference program), the labels for the ancestries (defined by the user) and meta-data with time and location information for each individual. The readkgfiles function combines all this data into a new variable called 'combined'. The option type fed into this function can be 'ohana' (if the ancestry Q matrix is in the format outputted by Ohana, Cheng et al. 2017) or 'admixture' (if the ancestry Q matrix is in the format outputted by the program Admixture, Alexander et al. 2009).

```
famfile <- "Data/ancient_modernEurope.fam"
ancestryfile <- "qpas_results/K4.Q.matrix"
ancall <- c("ANCE1","ANCE2","ANCE3","ANCE4")
timelocfile <- "Data/indstokeep.txt"
combined <- readkgfiles(famfile,ancestryfile,timelocfile,type="ohana",allowmissloc=TRUE,oldesttime=13000,minlat=35,maxlat=72,minlon=-20,maxlon=80)
```

# Create a spatial object
data <- combined
coordinates(data)=~LON+LAT
projection(data)=CRS(projectiontype)
data.UTM <- spTransform(data,CRS(projectiontype))

# Create spatial grid
#numgridpoints <- 5000
#numgridpoints <- 1000
numgridpoints <- 200
SpatialList <- CreateSpatialGrid(combined,landfeatures,numgridpoints,projectiontype)
dataSP <- SpatialList[[1]]; data.UTM <- SpatialList[[2]]; sp.grid.UTM <- SpatialList[[3]]; ReducedMap <- SpatialList[[4]]
allanc <- names(data.UTM)[-c(1,2)]
labanc <- c("ANC:\nNAFR","ANC:\nNEOL","ANC:\nHG","ANC:\nYAM")

# Create spatio-temporal grid
#oldesttime <- -10800; youngesttime <- 0; twindow <- 200; ntslots <- 55
#oldesttime <- -10800; youngesttime <- 0; twindow <- 600; ntslots <- 19
oldesttime <- -10800; youngesttime <- 0; twindow <- 200; ntslots <- 55
SPList <- CreateSpatioTemporalGrid(data.UTM, sp.grid.UTM, ntslots,oldesttime,youngesttime)
dataTM <- SPList[[1]]; grid.ST <- SPList[[2]]; rawtimegrid <- SPList[[3]]; tm.grid <- SPList[[4]]

# Matrix of neighbors
pointdists <- spDists(sp.grid.UTM,longlat=TRUE)
pointsnear <- t(apply(pointdists,1,function(allpoints){neighbortrue <- rank(allpoints) <= 10; return(as.numeric(neighbortrue))}))
pointsweights <- mat2listw(pointsnear)
plot(pointsweights,coordinates(sp.grid.UTM))

# Perform Kriging of ancestry data
AncestryKrigged <- list(); Ancestry <- list()
for(targetanc in allanc){
print(targetanc)
VarioList <- ComputeVariogram(data.UTM, dataSP, dataTM, targetanc, 50)
timeDF <- VarioList[[1]]; var <- VarioList[[2]]; finalVgm <- VarioList[[3]]; finalVgmMSE <- VarioList[[4]]
pred <- PerformSPKriging(timeDF,finalVgm,grid.ST)
attributes(pred)$data <- BoundKriging(attributes(pred)$data, 0, 1)
AncestryKrigged[[targetanc]] <- list(spobject=timeDF,variogram=var,varmodel=finalVgm,spgrid=sp.grid.UTM,tmgrid=tm.grid,rawtimegrid=rawtimegrid,pred=pred)
Ancestry[[targetanc]] <- pred
}


# VISUALIZATION - RAW DATA
land <- readOGR(landfeatures)
for( i in seq(1,length(allanc))){
toplot <- stplot(AncestryKrigged[[allanc[i]]]$spobject,colorkey=TRUE,main=labanc[i],number=10,mode="tp",sp.layout=list("sp.polygons",land)) 
toplot[[30]]$name <- as.numeric(as.POSIXct(toplot[[30]]$name,origin="1970-01-01"))
toplot
}

# VISUALIZATION - ST KRIGING
# Plot animation
for(targetanc in allanc){
for(i in seq(1,length(rawtimegrid))){
for(rep in seq(1,10)){
name <- rename((i-1)*10+rep)
png(name)
plot(AncestryKrigged[[targetanc]]$pred[,i],cex=1.5,pch=20,col=bpy.colors(21)[-1][cut(AncestryKrigged[[targetanc]]$pred[,i]@data[[1]],breaks=c(-0.2,seq(0.05,0.95,0.05),1.2),labels=FALSE)],main=as.character(AncestryKrigged[[targetanc]][["rawtimegrid"]])[i])
dev.off()
}}
my_command <- paste('convert -delay 3 *.png ',' ',targetanc,'_animation.gif',sep=""); system(my_command)
my_command <- 'rm *png'; system(my_command)
}
# Plot variogram
plot(AncestryKrigged[[targetanc]]$variogram,AncestryKrigged[[targetanc]]$varmodel,map=F,all=T)
attr(AncestryKrigged[[targetanc]]$varmodel,"MSE")
# Plot spatiotemporal plots
print(stplot(Ancestry[["ANCE1"]],cex=0.5, names.attr = as.character(AncestryKrigged[[targetanc]][["rawtimegrid"]]),main="NAF ancestry",colorkey=TRUE,ylim=c(35,70),xlim=c(-12,35)))
print(stplot(Ancestry[["ANCE2"]],cex=0.5, names.attr = as.character(AncestryKrigged[[targetanc]][["rawtimegrid"]]),main="NEOL ancestry",colorkey=TRUE,ylim=c(35,70),xlim=c(-12,35)))
print(stplot(Ancestry[["ANCE3"]],cex=0.5, names.attr = as.character(AncestryKrigged[[targetanc]][["rawtimegrid"]]),main="HG ancestry",colorkey=TRUE,ylim=c(35,70),xlim=c(-12,35)))
print(stplot(Ancestry[["ANCE4"]],cex=0.5, names.attr = as.character(AncestryKrigged[[targetanc]][["rawtimegrid"]]),main="YAM ancestry",colorkey=TRUE,ylim=c(35,70),xlim=c(-12,35)))
