# load MCMC results 
# create index of clusters for the data needed for making predictions for each sub tif files
# since we used a tapered covariate model, not all data points are needed for kriging for predicting every location;
# in order to save computation memory, we only load in the clusters of data which are relevant


load("mcmc_results_SOC.RData")

rm(list=ls()[!ls()=="locations.est.list"])
#load defined R functions
sph.cor.func <- function(d, phi){
	cor.value <- ((1-3/2*d/phi+1/2*(d/phi)^3)*(d<phi))
	return(cor.value)
}

loccords.jc <- function(predloc, locations.est.cluster){
	loc.dif <- cbind(locations.est.cluster[,1]-predloc[1], locations.est.cluster[,2]-predloc[2])
	dist.dif <- sqrt(rowSums(loc.dif^2))
	return(dist.dif)
}

loccords.jc.outer <- function(predloc, locations.est.cluster){
	loc.dif <- (outer( predloc[,1], locations.est.cluster[,1], "-"))^2+(outer( predloc[,2], locations.est.cluster[,2], "-"))^2
	dist.dif <- sqrt(loc.dif)
	return(dist.dif)
}


taper.range <- 40
library(rgdal)

setwd("/data3/JCdata/tif_data/")

path <- getwd()

file.create("datalist.keptcol.txt")

for(nm in list.files(path, pattern = "\\.tif$")[]){
		kept.chol <- NA
		tiffile <- nm
		tif.info <- GDALinfo(tiffile, silent=TRUE)
		totalcols <- tif.info[2]
		totalrows <- tif.info[1]
		origin.x <- tif.info[4]
		res.x <- tif.info[6]
		res.y <- tif.info[7]
		if(attr(tif.info, "ysign")==-1){
		#origin of data array should be the left upper corner
			origin.y <- tif.info[5] + totalrows*res.y
		}else{
		origin.y <- tif.info[5]	
		}
		evi.mask <- readGDAL(tiffile, band=5)@data$band1
		if(sum(!is.na(evi.mask))>0){
		predict.gridx <- seq(origin.x, origin.x + (totalcols-1)*res.x, by=res.x)
		predict.gridy <- seq(origin.y - (totalrows-1)*res.y, origin.y, by=res.y)
		predict.grid <- cbind(x=rep(predict.gridx, totalrows), y=sort(rep(predict.gridy, totalcols), decreasing=TRUE))
		predict.grid <- predict.grid/1000		
		predict.grid <- predict.grid[!is.na(evi.mask),]
		if(sum(!is.na(evi.mask))==1){
			predict.grid <- t(as.matrix(predict.grid))
		}
				
		
	dist.toest <- mapply(loccords.jc.outer, list(predict.grid), locations.est.list)
	dist.toest.taper.index <- (mapply(sum, mapply("<",dist.toest, taper.range))>0)*(1:length(locations.est.list))	
	kept.chol <- c(kept.chol, (dist.toest.taper.index[-1]-1)[dist.toest.taper.index[-1]>0])
		}
	kept.chol <- unique(kept.chol[-1])
	write.table(paste(nm, paste(as.character(kept.chol), collapse=" "), sep=","), file="datalist.keptcol_SOC.txt", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
	cat(paste(nm, paste(as.character(kept.chol), collapse=" "), sep=","), "\n")
}



