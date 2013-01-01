sph.cor.func <- function(d, phi){
	cor.value <- ((1-3/2*d/phi+1/2*(d/phi)^3)*(d<phi))
	return(cor.value)
}

loccords.jc <- function(predloc, locations.est.cluster){
	loc.dif <- cbind(locations.est.cluster[,1]-predloc[1], locations.est.cluster[,2]-predloc[2])
	dist.dif <- sqrt(rowSums(loc.dif^2))
	return(dist.dif)
}
		
kg.term <- function(dist.toest, phi.est, chol.est.list, dist.toest.taper.index){
	kg.selected <- rep(0, length(phi.est))
	kg.var <- rep(0, length(phi.est))
	if(dist.toest.taper.index[1]>0){
		sph.toest <-  sph.cor.func(dist.toest[[dist.toest.taper.index[1]]], 100)
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
			sph.toest <-  sph.cor.func(dist.toest.selected[[i]], 100)
			cor.toest.taper <- lapply(cor.toest,"*", sph.toest)
			cor.toest.kg <- mapply(backsolve, chol.selected[[i]], mapply(l=lapply(chol.selected[[i]], t), forwardsolve, x=cor.toest.taper, SIMPLIFY=FALSE))
			kg.selected <- rbind(kg.selected, cor.toest.kg)	
			kg.var <- kg.var + diag(t(mapply("*",  cor.toest, list(sph.toest)))%*%cor.toest.kg)
		}	
	}
	kg.selected <- kg.selected[-1, ]
	return(list(kg.selected, kg.var))
}


