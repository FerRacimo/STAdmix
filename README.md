# Spatio-temporal analysis of admixture components

These are scripts and commands used for the analysis in Racimo et al. 'A geostatistical approach to modelling human Holocene
migrations in Europe using ancient DNA' which is available in bioRxiv. The full list of analyses is in the file 'SpatioTemporal_Notes.txt'.

Below is a brief tutorial for performing spatio-temporal kriging of ancestry data. We will assume we have already performed admixture component inference on our data, using the program Ohana (Cheng et al. 2017). The resulting Q matrix is available in this repository, in the folder 'Data'.

First, we load the necessary functions and libraries in R:

```
R
source("STKrigFunc.R")
```

We then load the land feature data, which was obtained from https://www.naturalearthdata.com:

```
landfeatures <- "NaturalEarth/ne_50m_land/ne_50m_land.shp"
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

We now create a spatial object for our data, using the WGS84 (EPSG: 4326) coordinate reference system.

```
data <- combined
coordinates(data)=~LON+LAT
projectiontype <- "+init=epsg:4326"
projection(data)=CRS(projectiontype)
data.UTM <- spTransform(data,CRS(projectiontype))
```

We now create a spatial grid with 200 grid points.

```
numgridpoints <- 200
SpatialList <- CreateSpatialGrid(combined,landfeatures,numgridpoints,projectiontype)
dataSP <- SpatialList[[1]]; data.UTM <- SpatialList[[2]]; sp.grid.UTM <- SpatialList[[3]]; ReducedMap <- SpatialList[[4]]
allanc <- names(data.UTM)[-c(1,2)]
```

We then create a spatio-temporal grid spanning the last 10,800 years, in 200-year intervals.

```
oldesttime <- -10800; youngesttime <- 0; twindow <- 200; ntslots <- 55
SPList <- CreateSpatioTemporalGrid(data.UTM, sp.grid.UTM, ntslots,oldesttime,youngesttime)
dataTM <- SPList[[1]]; grid.ST <- SPList[[2]]; rawtimegrid <- SPList[[3]]; tm.grid <- SPList[[4]]
```


We are now ready to perform the spatio-temporal kriging of each of our ancestries. This is a two-step process. First, we must fit the spatio-temporal variogram (using the ComputeVariogram function, specifying 50 time lags of 60 years each). Then, we perform the actual kriging projection using the fitted variogram (using the PerformSPKriging function). The BoundKriging function ensures all our projected ancestry values are between 0 and 1.

```
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
```

For a given ancestry, we can visualize the empirical and fitted spatio-temporal variograms as follows:

```
targetanc <- "ANCE4"
plot(AncestryKrigged[[targetanc]]$variogram,AncestryKrigged[[targetanc]]$varmodel,map=F,all=T)
```

<img src="https://github.com/FerRacimo/STAdmix/blob/master/Variogram_YAM.png" height="300">

We can also visualize the kriging results as follows:

```
print(stplot(Ancestry[["ANCE1"]],cex=0.5, names.attr = as.character(AncestryKrigged[[targetanc]][["rawtimegrid"]]),main="ANCESTRY 1",colorkey=TRUE,ylim=c(35,70),xlim=c(-12,35)))
print(stplot(Ancestry[["ANCE2"]],cex=0.5, names.attr = as.character(AncestryKrigged[[targetanc]][["rawtimegrid"]]),main="ANCESTRY 2",colorkey=TRUE,ylim=c(35,70),xlim=c(-12,35)))
print(stplot(Ancestry[["ANCE3"]],cex=0.5, names.attr = as.character(AncestryKrigged[[targetanc]][["rawtimegrid"]]),main="ANCESTRY 3",colorkey=TRUE,ylim=c(35,70),xlim=c(-12,35)))
print(stplot(Ancestry[["ANCE4"]],cex=0.5, names.attr = as.character(AncestryKrigged[[targetanc]][["rawtimegrid"]]),main="ANCESTRY 4",colorkey=TRUE,ylim=c(35,70),xlim=c(-12,35)))
```

<img src="https://github.com/FerRacimo/STAdmix/blob/master/NEOL.png" height="300">
<img src="https://github.com/FerRacimo/STAdmix/blob/master/YAM.png" height="300">
<img src="https://github.com/FerRacimo/STAdmix/blob/master/HG.png" height="300">
<img src="https://github.com/FerRacimo/STAdmix/blob/master/NAF.png" height="300">

