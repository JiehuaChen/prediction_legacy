library(proj4)
library(nlme)

afhori <- read.table("Horizon.csv", header=T, sep=",", stringsAsFactors=FALSE)
afprof <- read.table("Profile.csv", header=T, sep=",")
covprof <- read.table("legacy_profile_withcov.csv", header=T, sep=",")

afprof.withcov <- cbind(afprof, covprof)
aflegacy <- merge(afhori, afprof.withcov, by="SPID")
aflegacy <-aflegacy[aflegacy$EVI_mask_1==1, ]

index.na <- (1:dim(aflegacy)[1])*is.na(rowSums(aflegacy[, c("SOC", "Bot", "bio1", "bio12", "CTI_1K", "ELEV_1K", "EVIM_1K", "M13RB1ALT", "M13RB2ALT", "M13RB3ALT", "M13RB7ALT", "NPP_Mean_1", "RELIEF_1K", "SLOPE_1K", "lstday", "lstnight")]))
aflegacy <- aflegacy[-index.na, ]

#band names:bio1,bio12,CTI_1K,ELEV_1K,EVI_mask_1K,EVIM_1K,M13RB1ALT,M13RB2ALT,M13RB3ALT,M13RB7ALT,MASK,NPP_Mean,RELIEF,SLOPE_1K
covar.selected.map <- c(1, 2, 3, 4, 6, 7, 12, 13, 15, 16)
covar.selected.legacyprof <- c(1, 2, 3, 4, 6, 7, 11, 12, 14, 15)
#mapnames <- read.table("mapnames.txt", stringsAsFactors=FALSE, header=FALSE)
mean.covar <- colMeans(aflegacy[, 19+covar.selected.legacyprof])
sd.covar <- sqrt(diag(var(aflegacy[, 19+covar.selected.legacyprof])))


###prediction set up
library(rgdal)
tif.info <- GDALinfo("predictgrid.tif", silent=TRUE)
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

predict.gridx <- seq(origin.x, origin.x + (totalcols-1)*res.x, by=res.x)
predict.gridy <- seq(origin.y, origin.y + (totalrows-1)*res.y, by=res.y)
predict.grid <- cbind(x=rep(predict.gridx, totalrows), y=sort(rep(predict.gridy, totalcols), decreasing=TRUE))
predict.grid <- predict.grid/1000		


pred.post.20 <- matrix(NA, dim(locations.pred)[1],  length(phi.est))
pred.post.110 <- matrix(NA,  dim(locations.pred)[1],  length(phi.est))

log.depth.pred1 <- log(20)
log.depth.pred2 <- log(110)

#prepare MCMC results
P <- dim(X)[2]
K <- dim(locations.est)[1]
mcmc.results.onechain <- t(mcmc.results[, 1, ])
mean.coef.est <- mcmc.results.onechain[(1:P), ]
alpha.est <- mcmc.results.onechain[((P+1):(P+K)), ]
sigma2.est <- mcmc.results.onechain[(P+K+1), ]
sigma2.alpha.est <- mcmc.results.onechain[(P+K+2), ]
phi.est <- mcmc.results.onechain[(P+K+3),]

sph.cor.func <- function(d, phi){
	cor.value <- ((1-3/2*d/phi+1/2*(d/phi)^3)*(d<phi))
	return(cor.value)
}

dist.toest <- mapply(loccoords, list(t(as.matrix(predict.grid))), locations.est.list)
dist.toest.taper <- (mapply(sum, mapply("<",dist.toest, 100)))*(1:length(locations.est.list))
sph.toest <-  mapply(sph.cor.func, dist.toest[dist.toest.taper], list(100))

dist.toest.rep <- mapply(matrix, mapply(rep, dist.toest[dist.toest.taper], list(length(phi.est))), nrow=length(phi.est), byrow=TRUE)
dist.toest.rep <- mapply("t", dist.toest.rep)
cor.toest <- lapply(lapply(lapply(dist.toest.rep, "-"), "%*%", diag(1/phi.est)), "exp")
cor.toest <- mapply("*", cor.toest, sph.toest)


sources("invcor_list.R")
inv.cor.est <- vector("list", length=length(phi.est))

for(i in 1:length(phi.est)){
	inv.cor.est[[i]] <- incor(d.site, phi.est[i], sph.cor20.site, sizes.noise, sizes.site, size.points)
	print(i)
}


residuals.Y <- matrix(rep(Y, length(phi.est)), N, length(phi.est))- X%*%mean.coef.est

predict.value <- matrix(NA ,dim(locations.pred)[1],  length(phi.est))

for(i in 1:(dim(predict.grid)[1])){
	ptm <- proc.time()
	print(i)
	if(!is.nan(evi.mask[i])){
			rowoffset <- floor(i/totalcols)
			coloffset <- i-i*rowoffset
			predcov.values <- readGDAL("predictgrid.tif", band=covar.selected.map, offset=c(rowoffset, coloffset), region.dim=c(1,1))@data
			predcov.values.scaled <- t((t(predcov.values)-mean.covar)*(1/sd.covar))
			predcov.values.whole <- cbind(1,log.depth.pred1, predcov.values.scaled, predcov.values.scaled*log.depth.pred1)
			pred.fixed <- (predcov.values.whole%*%mean.coef.est)
			dist.toest <- loccoords(t(as.matrix(predict.grid[i, ])),(locations.est))
			if(sum(dist.toest<100)==0){
				pred.random <- 0
			}else{
				taper.sph <-  ((1-3/2*dist.toest/100+1/2*(dist.toest/100)^3)*(dist.toest<100))
				dist.toest.scaled <- lapply(lapply(lapply(phi.est, solve), "*", dist.toest),"-")
				exp.cor <- lapply(dist.toest.scaled, exp)
				exp.cor.tap <- mapply("*", exp.cor, taper.sph)
				pred.random <- do.call("rbind", mapply("%*%", inv.cor.est, exp.cor.tap))
				}	
				
			pred.post.20[i, ] <- exp(pred.fixed + pred.random + rnorm(0, sqrt(sigma2.est)))
			
			predcov.values.whole <- cbind(1,log.depth.pred2, predcov.values.scaled, predcov.values.scaled*log.depth.pred2)
			pred.fixed <- (predcov.values.whole%*%mean.coef.est)	
			pred.post.110[(i+1), ] <- exp(pred.fixed + pred.random+rnorm(0, sqrt(sigma2.est)))
			print(paste("processing time", i))
			print(proc.time()-ptm)
		
	}
}
		
	
	