# generate chol.est.list: the cholesky decomposition list of the MCMC covariance matrices

load("mcmc_results_SOC.RData")
print(proc.time()-ptm)

mcmc.results.onechain <- t(mcmc.results[,1,])
P <- dim(X)[2]
K <- length(unique(PIDn))
phi.est <- mcmc.results.onechain[(P+K+3), ]



chol.est.list <- vector("list", length(d.site))

for(i in 1:length(d.site)){
	d.site.scaled  <- lapply(1/phi.est,"*", (-d.site[i][[1]]))
	cor <- lapply(d.site.scaled, exp)
	cor.taper <- lapply(cor, "*", (sph.cor20.site[i][[1]]))
	chol.est.list[[i]] <- lapply(cor.taper, chol)
	print(i)	
}

save(chol.est.list, file="chol_est_list.RData")




