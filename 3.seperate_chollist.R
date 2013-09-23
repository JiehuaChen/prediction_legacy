load("chol_est_list_exbase.RData")

 
for(i in 1:length(chol.est.list)){
	chol.est.list.used <-chol.est.list[[i]]
	save(chol.est.list.used, file=paste("chol_est_list_exbase", i, "RData", sep="."))
	print(i)
}
