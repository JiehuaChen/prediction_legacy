#!/usr/bin/Rscript
###prediction set up



key <- 0
options(warn=2)
taper.range <- 40
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
		s <- sum(dist.toest.taper.index>0&dist.toest.taper.index!=1)			
		dist.toest.taper.index[-1] <- ifelse(dist.toest.taper.index[-1]>0, dist.toest.taper.index[-1]-1, 0)
		chol.selected <- chol.est.list[dist.toest.taper.index[-1]] 
	
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

#load in calibrated data
load("mcmc_results_Clay.RData")


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



		
for(dp in 1:dim(depth.pred)[1]){
	file.create( paste(line.split[[1]][1],"depth", depth.pred[dp,1], "SOC.mean.txt", sep="."))
	file.create(paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.sd.txt", sep="."))
	file.create(paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.quantiles.txt", sep="."))
}
file.create(paste(line.split[[1]][1], "SOC.processtime.txt", sep="."))

for(i in 1:(dim(predict.grid)[1])){
	ptm <- proc.time()
	predcov.values <- covariatesdata[i, covar.selected.map]
	predcov.values.scaled <- t((t(predcov.values)-mean.covar)*(1/sd.covar))
	if(sum(abs(predcov.values.scaled)>100)==0){
		dist.toest <- mapply(loccords.jc, list(predict.grid[i,]), locations.est.list)
		dist.toest.taper.index <- (mapply(sum, mapply("<",dist.toest, taper.range))>0)*(1:length(locations.est.list))	
		if(sum(dist.toest.taper.index)==0){
			pred.random <- rnorm(length(phi.est), 0, sqrt(sigma2.alpha.est))
		}else{
		kg.term.values <- kg.term(dist.toest, phi.est, chol.est.list.used.combined, dist.toest.taper.index, taper.range)
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
			pred.quantiles <- exp(quantile(log(pred.post.draw), probs=c(0.025, 0.975)))
			write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.mean), paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.mean.txt", sep="."), append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
			write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.sd), paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.sd.txt", sep="."), append=TRUE,quote=FALSE, row.names=FALSE, col.names=FALSE)
			write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.quantiles), paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.quantiles.txt", sep="."), append=TRUE,quote=FALSE, row.names=FALSE, col.names=FALSE)
		}
	}else{
		for(dp in 1:dim(depth.pred)[1]){
			pred.mean <- NA
			pred.sd <- NA
			write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.mean), paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.mean.txt", sep="."), append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
			write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.sd), paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.sd.txt", sep="."), append=TRUE,quote=FALSE, row.names=FALSE, col.names=FALSE)
			write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), pred.quantiles), paste(line.split[[1]][1], "depth", depth.pred[dp,1], "SOC.quantiles.txt", sep="."), append=TRUE,quote=FALSE, row.names=FALSE, col.names=FALSE)
			}
		}								
	write.table(cbind(t(as.matrix(predict.grid[i, ]*1000)), as.matrix(proc.time()-ptm)[3]), paste(line.split[[1]][1], "SOC.processtime.txt", sep="."), append=TRUE,quote=FALSE, row.names=FALSE, col.names=FALSE)	
	if((proc.time()-startt)[3]>heartbeat_threshold){
		write((proc.time()-startt)[3], stderr())
	}
}

		

