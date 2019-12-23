library(tidyverse)
library(gstat)
library(sp)
library(spacetime)
library(spatialreg)
library(spdep)
library(lattice)
library(raster)
library(rgdal)
library(rgeos) 
library(ncdf4)
library(reshape)
library(xts)
library(igraph)
library(spTimer)
library(gridExtra)
library(parallel)
library(stringr)
library(geosphere)
library(lmodel2)
library(viridis)

#library(GWmodel)
#library(starma)
#library(SpatioTemporal)
#library(Matrix)
#library(plotrix)
#library(maps)

rename <- function(x){
  if (x < 10) {
    return(paste('000',x,'plot.png',sep=''))
  }
  if (x < 100 && x >= 10) {
    return(paste('00',x,'plot.png', sep=''))
  }
  if (x >= 100) {
    return(paste('0', x,'plot.png', sep=''))
  }
}


readkgfiles <- function(famfile,ancestryfile,timelocfile,type="ohana",allowmissloc=TRUE,oldesttime=13000,minlat=35,maxlat=72,minlon=-20,maxlon=80){

# Read fam file
fam <- read.table(famfile,header=FALSE)[,c(1,2)]
colnames(fam) <- c("IND","INDB")
# Read ancestry file - skip first line if using Ohana
if(type=="ohana"){ancestry <- read.table(ancestryfile,header=FALSE,skip=1)
} else if(type == "admixture"){ancestry <- read.table(ancestryfile,header=FALSE)}
colnames(ancestry) <- paste("ANCE",seq(1,dim(ancestry)[2]),sep="")
# Read time/lat/lon file
timeloc <- read.table(timelocfile,header=FALSE)
colnames(timeloc) <- c("IND","INDB","LAT","LON","TIME")
timeloc <- timeloc[which(timeloc[,1] %in% fam[,1]),]
famance <- data.frame(fam,ancestry)
combined <- inner_join(timeloc,famance,by="IND")
# Remove unwanted columns
combined$INDB.x <- NULL; combined$INDB.y <- NULL
# Filter for dates younger than 13,000 years ago
combined <- combined[which(combined$TIME < oldesttime),]
# Remove sites with missing locations
if(allowmissloc==FALSE){ combined <- combined[which(!is.na(combined$LAT) & !is.na(combined$LON)),] }
# Filter for sites inside square of interest
combined <- combined[which(combined$LAT > minlat),]
combined <- combined[which(combined$LAT < maxlat),]
combined <- combined[which(combined$LON > minlon),]
combined <- combined[which(combined$LON < maxlon),]
# Add nudge to space to avoid collisions
combined$LON <- combined$LON + rnorm(dim(combined)[1],0,0.1)
combined$LAT <- combined$LAT + rnorm(dim(combined)[1],0,0.1)
return(combined)

}






CreateSpatialGrid <- function(datatab,landfeatures,numsppoints,projectiontype,nmaxsp=Inf){

print("Creating spatial object...")

#Create a spatial object
data <- datatab
coordinates(data)=~LON+LAT
projection(data)=CRS(projectiontype)

#Transform into projection
data.UTM <- spTransform(data,CRS(projectiontype))
projCRS <- CRS(projectiontype)
dataSP <- SpatialPoints(data.UTM@coords,projCRS)

if(!is.na(landfeatures)){
	# Load spatial features
	print("Loading spatial features...")
	land <- readOGR(landfeatures)
	landProj <- spTransform(land, projCRS)	
	# Create spatial grid
	ReducedMap <- crop(landProj,data.UTM)

	if(!is.na(numsppoints)){
		sp.grid.UTM <- spsample(ReducedMap,n=numsppoints,type="regular")
	}
	else{
		dataSP <- dataSP[ReducedMap,]
		data.UTM <- data.UTM[ReducedMap,]
		sp.grid.UTM <- dataSP[ReducedMap,]
	}
} else {
       sp.grid.UTM <- dataSP
       ReducedMap <- NULL
}

sp.grid.UTM <- remove.duplicates(sp.grid.UTM)

return( list( dataSP, data.UTM, sp.grid.UTM , ReducedMap) )

}




CreateSpatioTemporalGrid <- function(data.UTM, sp.grid.UTM, numtimepoints, oldesttime = -Inf, youngesttime = 0){

# Years as seconds
dataTM <- as.POSIXlt(-data.UTM$TIME, origin="1970-01-01")
oldestdata <- max( oldesttime, min(-data.UTM$TIME) )
oldestdataTM <- as.POSIXlt(oldestdata, origin="1970-01-01")
youngestdata <- min( youngesttime, max(-data.UTM$TIME) )
youngestdataTM <- as.POSIXlt(youngestdata, origin="1970-01-01")

# Create temporal grid
if(!is.na(numtimepoints)){
	rawtimegrid <- round(seq(oldestdata,youngestdata,length.out=numtimepoints))
	tm.grid <- seq(as.POSIXct(oldestdataTM),as.POSIXct(youngestdataTM),length.out=numtimepoints)
} else {
       rawtimegrid <- sort(unique(-data.UTM$TIME))
       tm.grid <- as.POSIXct(rawtimegrid, origin="1970-01-01")
}

# Create spatiotemporal grid
grid.ST <- STF(sp.grid.UTM,tm.grid)

return( list( dataTM, grid.ST, rawtimegrid, tm.grid))

}


ComputeVariogram <- function(data.UTM, dataSP, dataTM, vartokrig, maxtlags, anistart=10, aniend=500, anistep=10,tunits="mins",width=NA){

# Select data
dataDF <- data.frame(data.UTM[[vartokrig]])
colnames(dataDF) <- "ANCE"

timeDF <- STIDF(dataSP,dataTM,data=dataDF)

print("Computing variogram...")

# Compute variogram - as default use "minutes" (years*60) for lag size
# From 60 years (="60 seconds = 1 minute") to 3000 years (="3000 seconds = 50 minutes")
if(is.na(width)){var <- variogramST(ANCE~1,data=timeDF,tunit=tunits,assumeRegular=F,na.omit=T,tlags=1:maxtlags) } else{
var <- variogramST(ANCE~1,data=timeDF,tunit=tunits,assumeRegular=F,na.omit=T,tlags=1:maxtlags,width=width)}

print("Estimating model for variogram...")

# Estimate model for variogram
pars.l <- c(sill.s = 0, range.s = 0, nugget.s = 0,sill.t = 0, range.t = 0, nugget.t = 0,sill.st = 0, range.st = 0, nugget.st = 0, anis = 0)
#pars.u <- c(sill.s = 1e7, range.s = 1e7, nugget.s = 1e7,sill.t = 1e7, range.t = 1e7, nugget.t = 1e7,sill.st = 1e7, range.st = 1e7, nugget.st = 1e7,anis = 1e7)
finalVgmMSE <- Inf
finalVgm <- NULL
for( anisotropy in seq(anistart,aniend,anistep)){
    try( {
    metric <- vgmST("metric", joint = vgm(psill=1,"Exp", range=5e3, nugget=1e1), stAni=anisotropy)
    metric_Vgm <- fit.StVariogram(var, metric, method="L-BFGS-B",lower=pars.l,tunit="mins")
    mse <- attr(metric_Vgm,"MSE")
    #print(paste("Anisotropy: ",anisotropy,"; MSE: ",mse,sep=""))
    if(mse < finalVgmMSE){
        finalVgmMSE <- mse
        finalVgm <- metric_Vgm
    }
    }, silent = TRUE)
}


return( list(timeDF, var,finalVgm,finalVgmMSE) )

}


PerformSPKriging <- function(timeDF,finalVgm,grid.ST,nmaxsp=Inf){
print("Performing spatio-temporal kriging")
pred <- krigeST(ANCE~1, data=timeDF, modelList=finalVgm, newdata=grid.ST,nmax=nmaxsp)
return(pred)

}


LoadncinSTData <- function(ncin, vartoget){

# First convert ncdf to same format at "combined"
lonvec <- ncvar_get(ncin,"lon")
latvec <- ncvar_get(ncin,"lat")
timevec <- ncvar_get(ncin,"time")
latlontime <- expand.grid(timevec,latvec,lonvec)[,c(2,3,1)]

envvartab <- c()
i <- 1
for(lon in lonvec){
envvartab <- rbind(envvartab, cbind(i, melt(t(ncvar_get( ncin, vartoget)[i,,]))[,c(2,1,3)]))
i <- i+1
}
envvartab[,c(1,2,3)] <- latlontime
envvartab <- cbind(seq(1,dim(envvartab)[1]), envvartab)
colnames(envvartab) <- c("IND","LAT","LON","TIME","VALUE")

return(envvartab)
}





RealifyValues <- function(presentdayfile, ReducedMap, data.UTM, vartoget, rastcols=100, rastrows=100){

# Present-day temperature
print("Loading present-day dataset...")
imported_raster=raster(presentdayfile)
presentday <- as(imported_raster, 'SpatialPointsDataFrame')
presentday <- spTransform(presentday,CRS(projectiontype))
presentday <- presentday[ReducedMap,]


# Raasterize presentday dataset
print("Rasterizing present-day dataset...")
rast <- raster()
extent(rast) <- extent(presentday)
ncol(rast) <- rastcols; nrow(rast) <- rastrows
presentdayrast = rasterize(presentday,rast,presentday[,1][[1]],fun=mean)
presentdayrast <- as(presentdayrast, 'SpatialPointsDataFrame')

# Find nearest points
print("Finding nearest points...")
realvalues <- unlist(sapply(seq(1,dim(data.UTM)[1]), function(pointidx){
alldist <- spDistsN1(presentday, data.UTM[pointidx,])
mindistidx <- which(alldist == min(alldist))[1]
basevalue <- presentday[mindistidx,][[1]]
realvalue <- basevalue + data.UTM[pointidx,3][[1]]
return(realvalue)
}))

# Replace values
print("Replacing values...")
data.UTM$VALUE <- realvalues
names(data.UTM)[which(names(data.UTM)=="VALUE")] <- vartoget
return(data.UTM)
}



ConvertPredToSpatioTemporal <- function(pred){

kriggedvar <- data.frame(pred)

# Observations
kriggedvar.obs <- data.frame(kriggedvar$var1.pred)
kriggedvar.obs <- cbind( kriggedvar$time, paste(kriggedvar$LAT, kriggedvar$LON,sep="_"), kriggedvar.obs)
kriggedvar.obs <- unique(kriggedvar.obs[,c(1,2,3)])
colnames(kriggedvar.obs) <- c("date","ID","obs")
kriggedvar.obs$date <- as.numeric(as.POSIXlt(kriggedvar.obs$date))
kriggedvar.obsdf <- c()
for(i in unique(kriggedvar.obs$ID)){ idx <- which(kriggedvar.obs$ID == i); newvec <- kriggedvar.obs$obs[idx]; kriggedvar.obsdf <- cbind(kriggedvar.obsdf,newvec) }
colnames(kriggedvar.obsdf) <- unique(kriggedvar.obs$ID); rownames(kriggedvar.obsdf) <- unique(kriggedvar.obs$date)

# Spatial covariates
kriggedvar.covars <- data.frame(cbind( paste(kriggedvar$LAT, kriggedvar$LON,sep="_"), kriggedvar$LAT, kriggedvar$LON ))
kriggedvar.covars <- unique(kriggedvar.covars[,c(1,2,3)])
colnames(kriggedvar.covars) <- c("ID","LAT","LON")

return(list(kriggedvar.obsdf,kriggedvar.covars))
}


ConvertEnvToSpatioTemporal <- function(EnvtimeDF){
dfenv <- data.frame(EnvtimeDF)
dfenv <- cbind( dfenv$time, paste(dfenv$LAT, dfenv$LON,sep="_"), dfenv$VALUE )
dfenv <- unique(dfenv[,c(1,2,3)])
dfenv <- data.frame(dfenv,stringsAsFactors = FALSE)
dfenv[,3] <- as.numeric(dfenv[,3])
colnames(dfenv) <- c("date","ID","obs")
dfenvdf <- c()
for(i in unique(dfenv$ID)){ idx <- which(dfenv$ID == i); newvec <- dfenv$obs[idx]; dfenvdf <- cbind(dfenvdf,newvec) }
colnames(dfenvdf) <- unique(dfenv$ID); rownames(dfenvdf) <- unique(dfenv$date)
return(dfenvdf)
}


# Load vegetation file
LoadVegTifFiles <- function(vegidxlist,vegfilelist,ReducedMap,sp.grid.UTM,projectiontype="+init=epsg:4326",rastrows=100,rastcols=100){
veglength <- length(vegidxlist)
finalveg <- NULL
for( idx in seq(1,veglength)){
i <- vegidxlist[idx]
vegfile <- vegfilelist[idx]
print(paste("Time slice: ",i,sep=""))
imported_raster=raster(vegfile)
vegcover <- as(imported_raster, 'SpatialPointsDataFrame')
vegcover <- spTransform(vegcover,CRS(projectiontype))
names(vegcover) <- "VALUE"
vegcover <- vegcover[ReducedMap,]
print(paste("Total points: ",length(vegcover),sep=""))
spplot(vegcover,cex=0.5)
# Rasterize dataset
rast <- raster()
extent(rast) <- extent(vegcover)
ncol(rast) <- rastcols; nrow(rast) <- rastrows
vegcoverrast = rasterize(vegcover,rast,vegcover[,1][[1]],fun=mean)
vegcoverrast <- as(vegcoverrast, 'SpatialPointsDataFrame')
vegcoverrast <- spTransform(vegcoverrast,CRS(projectiontype))
names(vegcoverrast) <- "VALUE"
print(paste("Points after rasterizing: ",length(vegcoverrast),sep=""))
# Median estimation (not kriging)
vegcovermed <- krige(VALUE~1, vegcoverrast, sp.grid.UTM, nmax=5,set = list(method = "med"))
vegcovermed <- vegcovermed[,1]; names(vegcovermed) <- "VALUE"
#vegcovermed$TIME <- rep(as.numeric(i),length(vegcovermed))
vegcovermed$IND <- length(vegcovermed)*(idx-1) + seq(1,length(vegcovermed))
print(paste("Points after median estimation: ",length(vegcovermed),sep=""))
if( suppressWarnings(is.null(finalveg))){
    #finalveg <- vegcovermed
    finalveg <- data.frame(VALUE=vegcovermed$VALUE,IND=vegcovermed$IND)
} else {
    #finalveg <- finalveg + vegcovermed
    finalveg <- rbind(finalveg,data.frame(VALUE=vegcovermed$VALUE,IND=vegcovermed$IND))
    }
print(" ")
}
#vegDF <- data.frame(finalveg$VALUE)
vegDF <- finalveg
vegTM <- xts(1:length(vegidxlist), as.POSIXlt(vegidxlist,origin="1970-01-01"))
row.names(sp.grid.UTM) = paste("point", 1:length(sp.grid.UTM), sep="")
#vegtimeDF <- as(STFDF(sp.grid.UTM,vegTM,data=vegDF),"STIDF")
vegtimeDF <- STFDF(sp.grid.UTM,vegTM,data=vegDF)
return(vegtimeDF)
}



ComputeDifferential <- function(testvar,grid,mode=1){
window <- length(grid)
totaltimes <- length(testvar)/window
totalchanges <- totaltimes - 1
changevar <- c()
if(mode == 1){
for(i in seq(1,totalchanges)){
    toadd <- testvar[seq( window*(i)+1,window*(i+1) )] - testvar[seq(window*(i-1)+1,window*(i))]
    changevar <- c(changevar, toadd) 
}
return(changevar)
} else if(mode == 2){
for(i in seq(1,totalchanges)){
    toadd <- testvar[seq( window*(i-1)+1,window*(i) )] - testvar[seq( window*(totalchanges)+1,window*(totalchanges+1) )]
    changevar <- c(changevar, toadd)
}
return(changevar)
} else if(mode == 0){
return(testvar)
} else if(mode == 3){
return(testvar[seq( 1,window*(totalchanges) )])
}
}


InterpolateTimes <- function(envvartab, rawtimegrid, mode=1){
newtab <- c()
testtimes <- -rawtimegrid
primtimes <- unique(envvartab$TIME)
for( rettime in testtimes){
rettab <- envvartab[which(envvartab$TIME == rettime),]
if( dim(rettab)[1] > 0 ){
    addtab <- rettab
    colnames(addtab) <- c("IND","LAT","LON","TIME","VALUE")
    newtab <- rbind(newtab,addtab)
} else{
    distances <- primtimes - rettime
    updistances <- distances
    downdistances <- distances
    updistances[which(updistances < 0)] <- NA
    downdistances[which(downdistances > 0)] <- NA
    if(sum(is.na(updistances)) == length(updistances)){upclosest <- NA
    } else{ upclosest <- primtimes[which(abs(updistances) == min(abs(updistances),na.rm=TRUE))]}
    if(sum(is.na(downdistances)) == length(downdistances)){downclosest <- NA
    } else{ downclosest <- primtimes[which(abs(downdistances) == min(abs(downdistances),na.rm=TRUE))]}
    twoclosest <- c(upclosest,downclosest)
    #twoclosest <- primtimes[abs(primtimes-rettime) %in% sort(abs(primtimes-rettime), partial=1:2)[1:2]]
    closest <- primtimes[abs(primtimes-rettime) %in% sort(abs(primtimes-rettime), partial=1)[1]]
    secondclosest <- twoclosest[which(twoclosest != closest)]
    if( !(NA %in% twoclosest) ){
        firststretch <- rettime - twoclosest[1]
    	secondstretch <- rettime - twoclosest[2]
	twostretch <- c(firststretch,secondstretch)
        totalt <- abs(twoclosest[2] - twoclosest[1])
    	firstprop <- abs(rettime - twoclosest[1])/totalt
    	secondprop <- abs(rettime - twoclosest[2])/totalt
    	firsttab <- envvartab[which(envvartab$TIME == twoclosest[1]),]
    	firstvalues <- firsttab[,5]
    	secondtab <- envvartab[which(envvartab$TIME == twoclosest[2]),]
    	secondvalues <- secondtab[,5]
    	firstvalues[which(is.na(firstvalues))] <- secondvalues[which(is.na(firstvalues))]
    	secondvalues[which(is.na(secondvalues))] <- firstvalues[which(is.na(secondvalues))]
    	if( rettime == 0 & mode == 1){ firstvalues <- 0; firstprop <- 1; secondvalues <- 0; secondprop <- 1 }
    	addtab <- cbind( firsttab[,seq(1,3)] , rettime, (firstprop*firstvalues + secondprop*secondvalues) )
	} else {
	closesttab <- envvartab[which(envvartab$TIME == closest),]
	closestvalues <- closesttab[,5]
        secondclosesttab <- envvartab[which(envvartab$TIME == secondclosest),]
        secondclosestvalues <- secondclosesttab[,5]
	closestvalues[which(is.na(closestvalues))] <- secondclosestvalues[which(is.na(closestvalues))]
	if( rettime == 0 & mode == 1){ closestvalues <- 0 }
	addtab <- cbind( closesttab[,seq(1,3)] , rettime, closestvalues )
	}
    colnames(addtab) <- c("IND","LAT","LON","TIME","VALUE")
    newtab <- rbind(newtab, addtab) 
}
}
return(newtab)
}




SpatialInterpolation <- function(rawtimegrid,RealEnvdata.UTM,sp.grid.UTM,medianmax=5){
timelength <- length(rawtimegrid)
finalenv <- NULL
for( idx in seq(1,timelength)){
i <- -rawtimegrid[idx]
envperiod <- RealEnvdata.UTM[which(RealEnvdata.UTM$TIME == i),]
names(envperiod)[3] <- "VALUE"
envperiod <- envperiod[which(!is.na(envperiod$VALUE)),]
envmed <- krige(VALUE~1, envperiod, sp.grid.UTM, nmax=medianmax,set = list(method = "med"))
envmed <- envmed[,1]
names(envmed)[1] <- "VALUE"
envmed$IND <- length(envmed)*(idx-1) + seq(1,length(envmed))
if( suppressWarnings(is.null(finalenv))){
    finalenv <- data.frame(VALUE=envmed$VALUE,IND=envmed$IND)
} else {
    finalenv <- rbind(finalenv,data.frame(VALUE=envmed$VALUE,IND=envmed$IND))
    }
}
envDF <- finalenv
envTM <- xts(1:length(rawtimegrid), as.POSIXlt(rawtimegrid,origin="1970-01-01"))
row.names(sp.grid.UTM) = paste("point", 1:length(sp.grid.UTM), sep="")
envtimeDF <- STFDF(sp.grid.UTM,envTM,data=envDF)
return(envtimeDF)
}



# Bound values of krigged variables
BoundKriging <- function(KriggedData,min,max){
leftbound <- which(KriggedData < min)
rightbound <- which(KriggedData > max)
KriggedData[leftbound,1] <- min
KriggedData[rightbound,1] <- max
return(KriggedData)
}


# Compute correlations between krigged variables
ComputeCors <- function(list1,list2,allvar1,allvar2,modecor,grid){
allvar1 <- sort(allvar1)
allvar2 <- sort(allvar2)
allcors <- c()
for(var1 in allvar1){
    corvec <- c()
    for(var2 in allvar2){
        vec1 <- attributes(list1[[var1]])$data[,1]
        vec2 <- attributes(list2[[var2]])$data[,1]
        vec1 <- ComputeDifferential(vec1,grid,modecor[1])
        vec2 <- ComputeDifferential(vec2,grid,modecor[2])
        newcor <- cor( vec1, vec2 )
	corvec <- c(corvec,newcor)
    }
    allcors <- rbind(allcors,corvec)
}
return(allcors)
}


# Plot correlations between krigged variables or coefficients of model

PlotCor <- function(allcors,allvar1,allvar2,title,labels1=NULL,labels2=NULL,legplot=TRUE,legtil="",legpos="topright",order=TRUE,vertexsize=25,asp=0.7,vertexcex=0.6){
if(order==TRUE){
	if(!is.null(labels1)){ order1 <- order(allvar1)} else{ order1 <- order(labels1)}
	if(!is.null(labels2)){ order2 <- order(allvar2)} else{ order2 <- order(labels2)}
	order2 <- order(allvar2)
	allvar1 <- allvar1[order1]
	allvar2 <- allvar2[order2]
	if(!is.null(labels1)){labels1 <- labels1[order1]} else{ labels1 <- allvar1}
	if(!is.null(labels2)){labels2 <- labels2[order2]} else{ labels2 <- allvar2}
	allcors <- allcors[order1,]
	allcors <- allcors[,order2]
}
allcors <- as.vector(t(allcors))
corcols <- rep(NA,length(allcors))
corcols[which(allcors < 0)] <- "blue"
corcols[which(allcors > 0)] <- "red"
widthscale <- 20
alphascale <- 8
corwidths <- abs(allcors)*widthscale
#alpha <- abs(allcors)*alphascale
alpha <- rep(1,length(allcors))
links <- matrix(rep(rep(1,length(allvar1)),length(allvar2)),ncol=length(allvar2))
net <- graph_from_incidence_matrix(links)
E(net)$color <- sapply(seq(1,length(allcors)), function(i){ adjustcolor( corcols[i], alpha[i] ) } )
E(net)$width <- corwidths
#net %>%add_layout_(as_bipartite()) %>%plot()
plot(net, vertex.label=c(labels1,labels2), vertex.size=vertexsize,layout=layout_as_bipartite,main=title,asp=asp,vertex.label.cex = vertexcex,vertex.color = adjustcolor("yellow", alpha.f = 1))
legvalues <- c(0.5,0.25,0.1,-0.1,-0.25,-0.5)
legcols <- legvalues
legcols[which(legcols<0)] <- "blue"
legcols[which(legcols!="blue")] <- "red"
if(legplot){legend(legpos,as.character(legvalues), col=legcols,lwd=abs(legvalues)*widthscale,title=legtil)}
}



# Process means from spTimer model
ProcessMeans <- function(test,zero=TRUE){
means <- test$parameter$Mean
if(zero==TRUE){means[which(sign(test$parameter$Low2.5p) != sign(test$parameter$Up97.5p))] <- 0}
means <- means[c(-1,-length(means),-length(means)+1,-length(means)+2)]
return(means)
}




# Create dataset for linear model
DataForModel <- function(AncestryKrigged, Fyfe, Climate, allanc, allfyfevar, allclim, rawtimegrid, type, standardize=TRUE){
if( type == "all"){
    # Actual variables
    lon <- rep(data.frame(AncestryKrigged[[1]]$spgrid)[,1],length(AncestryKrigged[[1]]$tmgrid))
    lat <- rep(data.frame(AncestryKrigged[[1]]$spgrid)[,2],length(AncestryKrigged[[1]]$tmgrid))
    time <- unlist(lapply(rawtimegrid, function(time){rep(time,length(AncestryKrigged[[1]]$spgrid))}))
    dataall <- cbind(lon,lat,time)
    for( ancvar in allanc){ dataall <- cbind(dataall, attributes(AncestryKrigged[[ancvar]]$pred)$data[,1]) }
    for( vegvar in allfyfevar){dataall <- cbind(dataall, attributes(Fyfe[[vegvar]])$data[,1])}
    for( climvar in allclim){dataall <- cbind(dataall,attributes(Climate[[climvar]])$data[,1])}
    dataall.sc <- scale(dataall)
    colnames(dataall) <- c("lon","lat","time",allanc,allfyfevar,allclim)
    dataall <- data.frame(dataall)
    colnames(dataall.sc) <- c("lon","lat","time",allanc,allfyfevar,allclim)
    dataall.sc <- data.frame(dataall.sc)
    if(standardize == TRUE){ return(dataall.sc) } else{ return(dataall) }
} else if( type == "diff"){
    # Variable differences
    lon <- rep(data.frame(AncestryKrigged[[1]]$spgrid)[,1],length(AncestryKrigged[[1]]$tmgrid)-1)
    lat <- rep(data.frame(AncestryKrigged[[1]]$spgrid)[,2],length(AncestryKrigged[[1]]$tmgrid)-1)
    time <- unlist(lapply(rawtimegrid[-length(rawtimegrid)], function(time){rep(time,length(AncestryKrigged[[1]]$spgrid))}))
    datadiff <- cbind(lon,lat,time)
    for( ancvar in allanc){ datadiff <- cbind(datadiff, ComputeDifferential(attributes(AncestryKrigged[[ancvar]]$pred)$data[,1],sp.grid.UTM,1)) }
    for( vegvar in allfyfevar){ datadiff <- cbind(datadiff, ComputeDifferential(attributes(Fyfe[[vegvar]])$data[,1],sp.grid.UTM,1)) }
    for( climvar in allclim){ datadiff <- cbind(datadiff, ComputeDifferential(attributes(Climate[[climvar]])$data[,1],sp.grid.UTM,1)) }
    datadiff.sc <- scale(datadiff)
    colnames(datadiff) <- c("lon","lat","time",allanc,allfyfevar,allclim)
    datadiff <- data.frame(datadiff)
    colnames(datadiff.sc) <- c("lon","lat","time",allanc,allfyfevar,allclim)
    datadiff.sc <- data.frame(datadiff.sc)
    if(standardize == TRUE){ return(datadiff.sc) } else{ return(datadiff) }
} else if( type == "anom"){
    # Variable anomalies
    lon <- rep(data.frame(AncestryKrigged[[1]]$spgrid)[,1],length(AncestryKrigged[[1]]$tmgrid)-1)
    lat <- rep(data.frame(AncestryKrigged[[1]]$spgrid)[,2],length(AncestryKrigged[[1]]$tmgrid)-1)
    time <- unlist(lapply(rawtimegrid[-length(rawtimegrid)], function(time){rep(time,length(AncestryKrigged[[1]]$spgrid))}))
    dataanom <- cbind(lon,lat,time)
    for( ancvar in allanc){ dataanom <- cbind(dataanom, ComputeDifferential(attributes(AncestryKrigged[[ancvar]]$pred)$data[,1],sp.grid.UTM,2)) }
    for( vegvar in allfyfevar){ dataanom <- cbind(dataanom, ComputeDifferential(attributes(Fyfe[[vegvar]])$data[,1],sp.grid.UTM,2)) }
    for( climvar in allclim){ dataanom <- cbind(dataanom, ComputeDifferential(attributes(Climate[[climvar]])$data[,1],sp.grid.UTM,2)) }
    dataanom.sc <- scale(dataanom)
    colnames(dataanom) <- c("lon","lat","time",allanc,allfyfevar,allclim)
    dataanom <- data.frame(dataanom)
    colnames(dataanom.sc) <- c("lon","lat","time",allanc,allfyfevar,allclim)
    dataanom.sc <- data.frame(dataanom.sc)
    if(standardize == TRUE){ return(dataanom.sc) } else{ return(dataanom) }
}
}



# Load climate data - PaleoClim
LoadPaleoClim <- function(paleoclimfol,worldclimfol,sp.grid.UTM,rawtimegrid,projectiontype){
totalclimnum <- 19
allpaleoclim <- paste("bio_",seq(1,totalclimnum),sep="")
labpaleoclim <- paste("bio_",seq(1,totalclimnum),sep="")
paleoclimtimes <- as.character(c(10013,6263,2250,0))
# Present-day climate files
PaleoClimFiles <- list()
for( i in seq(1,totalclimnum)){
PaleoClimFiles[[paste("bio_",i,sep="")]] <- list()
PaleoClimFiles[[paste("bio_",i,sep="")]][[paleoclimtimes[1]]] <- paste(paleoclimfol,"/EH_v1_10m/bio_",i,".tif",sep="")
PaleoClimFiles[[paste("bio_",i,sep="")]][[paleoclimtimes[2]]] <- paste(paleoclimfol,"/MH_v1_10m/bio_",i,".tif",sep="")
PaleoClimFiles[[paste("bio_",i,sep="")]][[paleoclimtimes[3]]] <- paste(paleoclimfol,"/LH_v1_10m/bio_",i,".tif",sep="")
PaleoClimFiles[[paste("bio_",i,sep="")]][[paleoclimtimes[4]]] <- paste(worldclimfol,"/wc2.0_bio_10m_",str_pad(i, 2, pad = "0"),".tif",sep="")
}
PaleoClim <- list()
for( bioclim in seq(1,totalclimnum)){
print(bioclim)
for(j in seq(1,length(paleoclimtimes))){
imported_raster=raster(PaleoClimFiles[[paste("bio_",bioclim,sep="")]][[paleoclimtimes[j]]])
vegcover <- as(imported_raster, 'SpatialPointsDataFrame')
vegcover <- spTransform(vegcover,CRS(projectiontype))
names(vegcover) <- "VALUE"
vegcover <- vegcover[ReducedMap,]
if(bioclim %in% c(1,2,5,6,7,8,9,10,11) ){ if(j < 4) {vegcover$VALUE <- vegcover$VALUE/10}}
vegcovermed <- krige(VALUE~1, vegcover, sp.grid.UTM, nmax=5,set = list(method = "med"))
vegcovermed <- vegcovermed[,1]
vegcovermed <- cbind(vegcovermed,paleoclimtimes[j])
names(vegcovermed) <- c("VALUE","TIME")
#print(vegcovermed)
if(j==1){finalmed <- vegcovermed} else{finalmed <- rbind(finalmed,vegcovermed)}
}
finalmed <- data.frame(finalmed)
finalmed <- as.data.frame(cbind(seq(1,dim(finalmed)[1]),finalmed[,4],finalmed[,3],as.numeric(as.character(finalmed[,2])),finalmed[,1]))
colnames(finalmed) <- c("IND","LAT","LON","TIME","VALUE")
finalmed <- InterpolateTimes(finalmed, rawtimegrid, mode=2)
finalmed$IND <- seq(1,dim(finalmed)[1])
envDF <- finalmed[,c(5,1)]
envTM <- xts(1:length(rawtimegrid), as.POSIXlt(rawtimegrid,origin="1970-01-01"))
finalpaleoclim <- STFDF(sp.grid.UTM,envTM,data=envDF)
PaleoClim[[paste("bio_",bioclim,sep="")]] <- finalpaleoclim
}
return(PaleoClim)
}

# Find nearest point in a ReacTran grid
FindNearestPointInGrid <- function(lat,lon,time,grid2D,times){
    gridlon <- grid2D$x.mid
    gridlat <- grid2D$y.mid
    newlat <- sapply(lat, function(latnum){gridlat[which.min(abs(gridlat  - latnum))]})
    newlon <- sapply(lon, function(lonnum){gridlon[which.min(abs(gridlon  - lonnum))]})
    newtime <- sapply(time, function(timenum){times[which.min(abs(times  - timenum))]})
    newtab <- cbind(newlat,newlon,newtime)
    colnames(newtab) <- c("LAT","LON","TIME")
    return(newtab)
}


