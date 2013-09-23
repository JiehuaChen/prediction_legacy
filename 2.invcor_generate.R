load("mcmc_results_exbase.RData")

print(proc.time()-ptm)
mcmc.results.onechain <- t(mcmc.results[,1,])
P <- dim(X)[2]
K <- length(unique(PIDn))
phi.est <- mcmc.results.onechain[(P+K+3), ]

# for(i in 1:length(d.site)){
	# d.site.scaled  <- lapply(1/phi.est,"*", (-d.site[i][[1]]))
	# cor <- lapply(d.site.scaled, exp)
	# cor.taper <- lapply(cor, "*", (sph.cor20.site[i][[1]]))
	# invcor <- lapply(cor.taper, solve)
	# save("invcor", file=paste("invcor", i, "RData", sep="."))
	# print(i)	
# }

chol.est.list <- vector("list", length(d.site))

for(i in 1:length(d.site)){
	d.site.scaled  <- lapply(1/phi.est,"*", (-d.site[i][[1]]))
	cor <- lapply(d.site.scaled, exp)
	cor.taper <- lapply(cor, "*", (sph.cor20.site[i][[1]]))
	chol.est.list[[i]] <- lapply(cor.taper, chol)
	print(i)	
}

save(chol.est.list, file="chol_est_list_exbase.RData")

# for(j in 1:length(phi.est)){
	# d.site.scaled  <- lapply(lapply(d.site, "-"), "/", phi.est[j]) 
	# cor <- lapply(d.site.scaled, exp)
	# cor.taper <- mapply("*", cor, sph.cor20.site)
	# chol.est <- lapply(cor.taper,chol)
	# save("chol.est", file=paste("cholcov", "draw=",j,"RData", sep="."))	
	# print(j)	
# }


