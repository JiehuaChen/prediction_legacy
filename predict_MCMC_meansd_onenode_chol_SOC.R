#!/usr/bin/Rscript
###prediction set up

f <- file("stdin",open = "r")

key <- 0
options(warn=2)
taper.range <- 40
heartbeat_threshold <- 300

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
		
kg.term <- function(dist.toest, phi.est, chol.est.list, dist.toest.taper.index, taper.range){
	kg.selected <- rep(0, length(phi.est))
	kg.var <- rep(0, length(phi.est))
	if(dist.toest.taper.index[1]>0){
		sph.toest <-  sph.cor.func(dist.toest[[dist.toest.taper.index[1]]], taper.range)
		kg.selected <- rbind(kg.selected, mapply("*", lapply(lapply(1/phi.est, "*", (-dist.toest[[1]])), exp), list(sph.toest)))
		kg.var <- colSums(kg.selected^2)
	}
	
	if(sum(dist.toest.taper.index[-1]>0)){
		dist.toest.selected <- dist.toest[dist.toest.taper.index[-1]]
		dist.toest.taper.index[-1] <- ifelse(dist.toest.taper.index[-1]>0, dist.toest.taper.index[-1]-1, 0)

		chol.selected <- chol.est.list[dist.toest.taper.index[-1]] 
		s <- sum(dist.toest.taper.index>0&dist.toest.taper.index!=1)

		for(i in 1:s){
			cor.toest <-  lapply(lapply(1/phi.est,"*", (-dist.toest.selected[[i]])), exp)
			sph.toest <-  sph.cor.func(dist.toest.selected[[i]], taper.range)
			cor.toest.taper <- lapply(cor.toest,"*", sph.toest)
			cor.toest.kg <- mapply(backsolve, chol.selected[[i]], mapply(l=lapply(chol.selected[[i]], t), forwardsolve, x=cor.toest.taper, SIMPLIFY=FALSE))
			kg.selected <- rbind(kg.selected, cor.toest.kg)	
			kg.var <- kg.var + diag(t(mapply("*",  cor.toest, list(sph.toest)))%*%cor.toest.kg)
		}	
	}
	kg.selected <- kg.selected[-1, ]
	return(list(kg.selected, kg.var))
}

load("mcmc_results_SOC.RData")
covar.selected.est <- c( "bio1", "bio12", "CTI_1K", "ELEV_1K", "EVIM_1K", "M13RB1ALT", "NPP_Mean_1", "RELIEF_1K", "lstday", "lstnight")
mean.covar <- colMeans(aflegacy[, (names(aflegacy)%in%covar.selected.est)])
sd.covar <- sqrt(diag(var(aflegacy[, (names(aflegacy)%in%covar.selected.est)])))
	
#depth data		
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
	
#layer selected
	covar.selected.map <-  c(1, 2, 3, 4, 6, 7, 12, 13, 15, 16)

#load rgdal
	library(rgdal)
	
while(length(line <- readLines(f,n=1, warn=FALSE)) > 0) {
	line.split <- strsplit(line, split=",")
	system(paste("hadoop dfs -get s3n://afsis.legacy.prediction/data/dataForEachNode/", line.split[[1]][1], " .", sep=""))
	system(paste("hadoop dfs -get s3n://afsis.legacy.prediction/data/dataForEachNode/", line.split[[1]][2], " .", sep=""))
	if(length(grep("\\.tif$", line.split[[1]][1]))>0){	
		startt <- proc.time()
		tiffile <- line.split[[1]][1]
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
		block_count <- 0
	if(sum(!is.na(evi.mask))>0){
		predict.gridx <- seq(origin.x, origin.x + (totalcols-1)*res.x, by=res.x)
		predict.gridy <- seq(origin.y - (totalrows-1)*res.y, origin.y, by=res.y)
		predict.grid <- cbind(x=rep(predict.gridx, totalrows), y=sort(rep(predict.gridy, totalcols), decreasing=TRUE))
		predict.grid <- predict.grid/1000		
		predict.grid <- predict.grid[!is.na(evi.mask),]
			
		block_count <- dim(predict.grid)[1]
				
		load(line.split[[1]][2])
							
		#predict for the first point	
		i <- 1	
		ptm <- proc.time()
		pixeln.x <- ((predict.grid[i,1]*1000-origin.x)/res.x)
		pixeln.y <- ((origin.y-predict.grid[i,2]*1000)/res.y)
		predcov.values <- readGDAL(tiffile, band=covar.selected.map, offset=c(pixeln.y, pixeln.x), region.dim=c(1,1), silent=TRUE)@data
		predcov.values.scaled <- t((t(predcov.values)-mean.covar)*(1/sd.covar))
		# only predict the locations with reasonable covariates
		if(sum(abs(predcov.values.scaled)>10)==0){
			dist.toest <- mapply(loccords.jc, list(predict.grid[i,]), locations.est.list)
			dist.toest.taper.index <- (mapply(sum, mapply("<",dist.toest, taper.range))>0)*(1:length(locations.est.list))	
			if(sum(dist.toest.taper.index)==0){
				pred.random <- rnorm(length(phi.est), 0, sqrt(sigma2.alpha.est))
			}else{
			kg.term.values <- kg.term(dist.toest, phi.est, chol.est.list.used, dist.toest.taper.index, taper.range)
			alpha.selected <- do.call("rbind", alpha.est.list[dist.toest.taper.index])
			kg.mean <- diag(t(alpha.selected)%*%kg.term.values[[1]])
			kg.var <- sigma2.alpha.est - sigma2.alpha.est*kg.term.values[[2]]
			pred.random <- rnorm(length(phi.est), kg.mean, sqrt(kg.var))
			}
			for(dp in 1:dim(depth.pred)[1]){
				log.depth.pred1 <- log(depth.pred[dp,1])
				predcov.values.whole <- cbind(1, (predcov.values.scaled), (predcov.values.scaled*log.depth.pred1), log.depth.pred1)
				pred.fixed <- (predcov.values.whole%*%mean.coef.est)
				pred.post.draw <- exp(pred.fixed + pred.random + rnorm(length(phi.est), 0, sqrt(sigma2.est)))
				pred.mean <- mean(pred.post.draw)
				pred.sd <- apply(log(pred.post.draw), 1, sd)
				write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.mean), paste(line.split[[1]][1],"depth", depth.pred[dp,1], "SOC.mean.txt", sep="."), quote=FALSE, row.names=FALSE, col.names=FALSE)
				write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.sd), paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.sd.txt", sep="."), quote=FALSE, row.names=FALSE, col.names=FALSE)	
				}
		}else{
			for(dp in 1:dim(depth.pred)[1]){
				pred.mean <- NA
				pred.sd <- NA
				write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.mean), paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.mean.txt", sep="."), quote=FALSE, row.names=FALSE, col.names=FALSE)
				write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.sd), paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.sd.txt", sep="."), quote=FALSE, row.names=FALSE, col.names=FALSE)		
			}
		}							
		write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), as.matrix(proc.time()-ptm)[3]), paste(line.split[[1]][1], "SOC.processtime.txt", sep="."), quote=FALSE, row.names=FALSE, col.names=FALSE)		
		if((proc.time()-startt)[3]>heartbeat_threshold){
			write((proc.time()-startt)[3], stderr())
		}
		
		for(i in 2:(dim(predict.grid)[1])){
			ptm <- proc.time()
			pixeln.x <- ((predict.grid[i,1]*1000-origin.x)/res.x)
			pixeln.y <- ((origin.y-predict.grid[i,2]*1000)/res.y)
			predcov.values <- readGDAL(tiffile, band=covar.selected.map, offset=c(pixeln.y, pixeln.x), region.dim=c(1,1), silent=TRUE)@data
			predcov.values.scaled <- t((t(predcov.values)-mean.covar)*(1/sd.covar))
			if(sum(abs(predcov.values.scaled)>10)==0){
				dist.toest <- mapply(loccords.jc, list(predict.grid[i,]), locations.est.list)
				dist.toest.taper.index <- (mapply(sum, mapply("<",dist.toest, taper.range))>0)*(1:length(locations.est.list))	
				if(sum(dist.toest.taper.index)==0){
					pred.random <- rnorm(length(phi.est), 0, sqrt(sigma2.alpha.est))
				}else{
				kg.term.values <- kg.term(dist.toest, phi.est, chol.est.list.used, dist.toest.taper.index, taper.range)
				alpha.selected <- do.call("rbind", alpha.est.list[dist.toest.taper.index])
				kg.mean <- diag(t(alpha.selected)%*%kg.term.values[[1]])
				kg.var <- sigma2.alpha.est - sigma2.alpha.est*kg.term.values[[2]]
				pred.random <- rnorm(length(phi.est), kg.mean, sqrt(kg.var))
				}
				for(dp in 1:dim(depth.pred)[1]){
					log.depth.pred1 <- log(depth.pred[dp,1])
					predcov.values.whole <- cbind(1, (predcov.values.scaled), (predcov.values.scaled*log.depth.pred1), log.depth.pred1)
					pred.fixed <- (predcov.values.whole%*%mean.coef.est)
					pred.post.draw <- exp(pred.fixed + pred.random + rnorm(length(phi.est), 0, sqrt(sigma2.est)))
					pred.mean <- mean(pred.post.draw)
					pred.sd <- apply(log(pred.post.draw), 1, sd)
					write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.mean), paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.mean.txt", sep="."), append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
					write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.sd), paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.sd.txt", sep="."), append=TRUE,quote=FALSE, row.names=FALSE, col.names=FALSE)
				}
			}else{
				for(dp in 1:dim(depth.pred)[1]){
					pred.mean <- NA
					pred.sd <- NA
					write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.mean), paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.mean.txt", sep="."), append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
					write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.sd), paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.sd.txt", sep="."), append=TRUE,quote=FALSE, row.names=FALSE, col.names=FALSE)
					}
				}								
			write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), as.matrix(proc.time()-ptm)[3]), paste(line.split[[1]][1], "SOC.processtime.txt", sep="."), append=TRUE,quote=FALSE, row.names=FALSE, col.names=FALSE)	
			if((proc.time()-startt)[3]>heartbeat_threshold){
				write((proc.time()-startt)[3], stderr())
			}
		}
		for(dp in 1:dim(depth.pred)[1]){
			system(paste("hadoop dfs -put", paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.mean.txt", sep="."), "s3n://afsis.legacy.prediction/results/"))
			system(paste("hadoop dfs -put", paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.sd.txt", sep="."), "s3n://afsis.legacy.prediction/results/"))
		}	
		system(paste("hadoop dfs -put", paste(line, "SOC.processtime.txt", sep="."), "s3n://afsis.legacy.prediction/results/"))
		cat("LongValueSum:1\t" , block_count)
		rm(list="chol.est.list.used")
		}		
	}
	system(paste("rm", line.split[[1]][1]))
	system(paste("rm", line.split[[1]][2]))		
}

close(f)
