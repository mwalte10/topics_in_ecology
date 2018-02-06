setwd("~/Desktop/outputs/test_vac_parm")
out.list = list()
FOI.list = list()
FOI.h.list = list()
FOI.l.list = list()

for(i in 1:4){
  load(paste('output_', i, '.RData', sep = ''))
  out.list[[i]] = out_last[2:673]
  load(paste('FOI_', i, '.RData', sep = ''))
  FOI.list[[i]] = FOI
  load(paste('FOI.h_', i, '.RData', sep = ''))
  FOI.h.list[[i]] = FOI_h
  load(paste('FOI.l_', i, '.RData', sep = ''))
  FOI.l.list[[i]] = FOI_l
}

FOI_mat <- matrix(NA, nrow = 4, ncol = 3)
for(i in 1:4){
  FOI <- as.vector(unlist(FOI.list[i])) 
  FOI_mat[i,1] <- min((!is.nan(FOI)))
  FOI_mat[i,2] <- mean((!is.nan(FOI)))
  FOI_mat[i,3] <- max((!is.nan(FOI)))
}

FOI_gen <- as.vector(unlist(FOI.list[1])) 
for(i in 2:4){
  FOI <- as.vector(unlist(FOI.list[i])) 
  plot(FOI_gen, log = "y", type = "l", col = rainbow(4)[1])
  lines(FOI, col = rainbow(4)[i])
}

FOI.h_mat <- matrix(NA, nrow = 4, ncol = 3)
for(i in 1:4){
  FOI.h <- as.vector(unlist(FOI.h.list[i])) 
  FOI.h_mat[i,1] <- min((!is.nan(FOI.h)))
  FOI.h_mat[i,2] <- mean((!is.nan(FOI.h)))
  FOI.h_mat[i,3] <- max((!is.nan(FOI.h)))
}
FOI.h_gen <- as.vector(unlist(FOI.h.list[1])) 
for(i in 2:4){
  FOI <- as.vector(unlist(FOI.h.list[i])) 
  plot(FOI.h_gen, log = "y", type = "l", col = rainbow(4)[1])
  lines(FOI, col = rainbow(4)[i])
}


FOI.l_mat <- matrix(NA, nrow = 4, ncol = 3)
for(i in 1:4){
  FOI <- as.vector(unlist(FOI.l.list[i])) 
  FOI.l_mat[i,1] <- min((!is.nan(FOI)))
  FOI.l_mat[i,2] <- mean((!is.nan(FOI)))
  FOI.l_mat[i,3] <- max((!is.nan(FOI)))
}
FOI.l_gen <- as.vector(unlist(FOI.l.list[1])) 
for(i in 2:4){
  FOI <- as.vector(unlist(FOI.l.list[i])) 
  plot(FOI.l_gen, log = "y", type = "l", col = rainbow(4)[1])
  lines(FOI, col = rainbow(4)[i])
}

out_gen <- as.vector(unlist(out.list[1]))
for(i in 2:21){
  out <- as.vector(unlist(out.list[i]))
  plot(out_gen, type = "l", col = rainbow(21)[1])
  lines(out, col = rainbow(21)[i])
}

FOI <- as.vector(unlist(FOI.list[1])) 
out <- as.vector(unlist(out.list[1]))

setwd("~/Desktop/outputs/test_vac_parm_short")
out.list_s = list()
FOI.list_s = list()
FOI.h.list_s = list()
FOI.l.list_s = list()

for(i in 1:4){
  load(paste('output_', i, '.RData', sep = ''))
  out.list_s[[i]] = out_last[2:673]
  load(paste('FOI_', i, '.RData', sep = ''))
  FOI.list_s[[i]] = FOI
  load(paste('FOI.h_', i, '.RData', sep = ''))
  FOI.h.list_s[[i]] = FOI_h
  load(paste('FOI.l_', i, '.RData', sep = ''))
  FOI.l.list_s[[i]] = FOI_l
}

FOI_mat_s <- matrix(NA, nrow = 4, ncol = 3)
for(i in 1:4){
  FOI <- as.vector(unlist(FOI.list_s[i])) 
  FOI_mat_s[i,1] <- min((!is.nan(FOI)))
  FOI_mat_s[i,2] <- mean((!is.nan(FOI)))
  FOI_mat_s[i,3] <- max((!is.nan(FOI)))
}

FOI.h_mat_s <- matrix(NA, nrow = 4, ncol = 3)
for(i in 1:4){
  FOI.h <- as.vector(unlist(FOI.h.list_s[i])) 
  FOI.h_mat_s[i,1] <- min((!is.nan(FOI.h)))
  FOI.h_mat_s[i,2] <- mean((!is.nan(FOI.h)))
  FOI.h_mat_s[i,3] <- max((!is.nan(FOI.h)))
}

FOI.l_mat_s <- matrix(NA, nrow = 4, ncol = 3)
for(i in 1:4){
  FOI <- as.vector(unlist(FOI.l.list_s[i])) 
  FOI.l_mat_s[i,1] <- min((!is.nan(FOI)))
  FOI.l_mat_s[i,2] <- mean((!is.nan(FOI)))
  FOI.l_mat_s[i,3] <- max((!is.nan(FOI)))
}

FOI <- as.vector(unlist(FOI.list_s[1])) 
plot(FOI, type = "l")

out <- as.vector(unlist(out.list_s[1])) 
plot(out, type = "l")


