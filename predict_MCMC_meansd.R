###prediction set up
#!/usr/bin/Rscript

library(rgdal)

		setwd("../..")
		tiffile <- "predictgrid.tif"
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
	
		setwd("Legacy_Data/predict_SOC")
		meansd.read <- read.table("meansd.txt",header=TRUE)
		mean.covar <- meansd.read[1, ]
		sd.covar <- meansd.read[2, ]
		
		load("mcmc_results_SOC.RData")
		aflegacy.lambertcord <- project(cbind(aflegacy$Lon, aflegacy$Lat), "+proj=laea +datum=WGS84 +lat_0=5 +lon_0=20")
		aflegacy.lambertcord <- aflegacy.lambertcord/1000

		#ordered unique locations
		locations.est <- cbind(aggregate(aflegacy.lambertcord[,1], by=list(PID = PIDn), mean)[,2], aggregate(aflegacy.lambertcord[,2], by=list(PID = PIDn), mean)[,2])
		locations.est.list  <-  vector("list", (length(d.site)+1))
		sizes.noise <- table.dbscan[1]
		sizes.site <- lapply(d.site, dim)
		locations.est.list[[1]] <- locations.est[1:sizes.noise, ]
		k <- sizes.noise+1
		for(i in 1:(length(d.site))){
			locations.est.list[[(i+1)]] <-locations.est[(k:(k+sizes.site[[i]][1]-1)), ]
			k <- k+(sizes.site[[i]])[1]
		}
		
		source("invcor_list_pred.R")
		depth.pred <- read.table("depthpred.txt")
				
		#prepare MCMC results
		P <- dim(X)[2]
		K <- dim(locations.est)[1]
		mcmc.results.onechain <- t(mcmc.results[, 1, ])
		mean.coef.est <- mcmc.results.onechain[(1:P), ]
		alpha.est <- mcmc.results.onechain[((P+1):(P+K)), ]
		alpha.est.list <- vector("list", (length(d.site)+1))
		alpha.est.list[[1]] <- alpha.est[1:sizes.noise, ]
		k <- sizes.noise+1
		for(i in 1:(length(d.site))){
			alpha.est.list[[i+1]] <- alpha.est[(k:(k+sizes.site[[i]][1]-1)), ]
			k <- k+(sizes.site[[i]])[1]
		}

		sigma2.est <- mcmc.results.onechain[(P+K+1), ]
		sigma2.alpha.est <- mcmc.results.onechain[(P+K+2), ]
		phi.est <- mcmc.results.onechain[(P+K+3),]
		
		covar.selected.map <-  c(1, 2, 3, 4, 6, 7, 12, 13, 15, 16)
		setwd("../..")
		
	#predict for the first point
		i <- 1
			ptm <- proc.time()
			pixeln.x <- ((predict.grid[i,1]*1000-origin.x)/res.x)
			pixeln.y <- ((origin.y-predict.grid[i,2]*1000)/res.y)
			predcov.values <- readGDAL(tiffile, band=covar.selected.map, offset=c(pixeln.y, pixeln.x), region.dim=c(1,1))@data
			predcov.values.scaled <- t((t(predcov.values)-mean.covar)*(1/sd.covar))
			dist.toest <- mapply(loccords.jc, list(predict.grid[i,]), locations.est.list)
			dist.toest.taper.index <- (mapply(sum, mapply("<",dist.toest, 100))>0)*(1:length(locations.est.list))	
			if(sum(dist.toest.taper.index)==0){
				pred.random <- 0
			}else{
			kg.term.values <- kg.term(dist.toest, d.site, phi.est, sph.cor20.site, dist.toest.taper.index)
			alpha.selected <- do.call("rbind", alpha.est.list[dist.toest.taper.index])
			kg.mean <- diag(t(alpha.selected)%*%kg.term.values[[1]])
			kg.var <- sigma2.alpha.est - sigma2.alpha.est*kg.term.values[[2]]
			pred.random <- rnorm(length(phi.est), kg.mean, kg.var)
			}
			for(dp in 1:dim(depth.pred)[1]){
				log.depth.pred1 <- log(depth.pred[dp,1])
				predcov.values.whole <- cbind(1, t(predcov.values.scaled), t(predcov.values.scaled*log.depth.pred1), log.depth.pred1)
				pred.fixed <- (predcov.values.whole%*%mean.coef.est)
				pred.post.draw <- exp(pred.fixed + pred.random + rnorm(length(phi.est), 0, sqrt(sigma2.est)))
				pred.mean <- mean(pred.post.draw)
				pred.sd <- apply(log(pred.post.draw), 1, sd)
				write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.mean), paste("depth", depth.pred[dp,1], "SOC.mean.txt", sep="."), quote=FALSE, row.names=FALSE, col.names=FALSE)
				write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.sd), paste("depth", depth.pred[dp,1], "SOC.sd.txt", sep="."), quote=FALSE, row.names=FALSE, col.names=FALSE)
			}						
		print(paste("processing time", i))
		print(proc.time()-ptm)
		


		for(i in 2:(dim(predict.grid)[1])){
			ptm <- proc.time()
			pixeln.x <- ((predict.grid[i,1]*1000-origin.x)/res.x)
			pixeln.y <- ((origin.y-predict.grid[i,2]*1000)/res.y)
			predcov.values <- readGDAL(tiffile, band=covar.selected.map, offset=c(pixeln.y, pixeln.x), region.dim=c(1,1))@data
			predcov.values.scaled <- t((t(predcov.values)-mean.covar)*(1/sd.covar))
			dist.toest <- mapply(loccords.jc, list(predict.grid[i,]), locations.est.list)
			dist.toest.taper.index <- (mapply(sum, mapply("<",dist.toest, 100))>0)*(1:length(locations.est.list))	
			if(sum(dist.toest.taper.index)==0){
				pred.random <- 0
			}else{
			kg.term.values <- kg.term(dist.toest, d.site, phi.est, sph.cor20.site, dist.toest.taper.index)
			alpha.selected <- do.call("rbind", alpha.est.list[dist.toest.taper.index])
			kg.mean <- diag(t(alpha.selected)%*%kg.term.values[[1]])
			kg.var <- sigma2.alpha.est - sigma2.alpha.est*kg.term.values[[2]]
			pred.random <- rnorm(length(phi.est), kg.mean, kg.var)
			}
			for(dp in 1:dim(depth.pred)[1]){
				log.depth.pred1 <- log(depth.pred[dp,1])
				predcov.values.whole <- cbind(1, t(predcov.values.scaled), t(predcov.values.scaled*log.depth.pred1), log.depth.pred1)
				pred.fixed <- (predcov.values.whole%*%mean.coef.est)
				pred.post.draw <- exp(pred.fixed + pred.random + rnorm(length(phi.est), 0, sqrt(sigma2.est)))
				pred.mean <- mean(pred.post.draw)
				pred.sd <- apply(log(pred.post.draw), 1, sd)
				write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.mean), paste("depth", depth.pred[dp,1], "SOC.mean.txt", sep="."), append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
				write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.sd), paste("depth", depth.pred[dp,1], "SOC.sd.txt", sep="."), append=TRUE,quote=FALSE, row.names=FALSE, col.names=FALSE)
			}							
		print(paste("processing time", i))
		print(proc.time()-ptm)
		}
	}		
	