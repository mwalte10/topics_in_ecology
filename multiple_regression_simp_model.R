library(fields)
sp9.h <- c()
sp9 <- c()
sp9.l <- c()
# setwd('~/Desktop/topics_in_ecology/new_calib/travel')
setwd('~/Desktop/topics_in_ecology/new_calib/no_travel')
#setwd('~/Desktop/topics_in_ecology/new_calib/skew_l')
#setwd('~/Desktop/topics_in_ecology/training_beta/skew_l')

for(length in 1:100){
  
    tryCatch({
      load(paste('sp9_', length, '.RData', sep = ''))
      sp9.h[length] <- sp9.vec[2]
      sp9.l[length] <- sp9.vec[3]
      sp9[length] <- sp9.vec[1]
    }, error=function(e){})

}



beta_h.vec <- seq(0,0.2, length.out = 10)
beta_l.vec <- seq(0.2,0.6, length.out = 10)
beta.mat <- as.vector(sapply(beta_h.vec, function(i)rep(i,10)))
beta.mat <- cbind(beta.mat, rep(beta_l.vec, 10))
parms <- cbind(beta.mat, rep(0.05,100), sp9, sp9.h, sp9.l)
parm_mat <- as.data.frame(parms)
travel_levels <- unique(parm_mat[,3])


outputs.list <- list()
for(i in 1:length(travel_levels)){
  parms <- parm_mat.list[[i]]
  sp9 <- parms[,4]
  sp9.h <- parms[,5]
  sp9.l <- parms[,6]
  beta_h <- parms[,1]
  beta_l <- parms[,2]
  
  travel_1 <- data.frame(cbind(beta_h, beta_l, sp9))
  names(travel_1) <- c('beta_h', 'beta_l', 'sp9')
  
  spline_sp9 = Tps(x = as.matrix(travel_1[c("beta_h","beta_l")]),
                   Y = matrix(travel_1["sp9"][,1]),
                   scale.type = "unscaled",lambda = 0)
  x <- seq(0, 0.3, length.out = 50)
  list.x <- list()
  for(j in 1:length(x)){
    list.x[[j]] <- rep(x[j], 50)
  }
  list.x <- unlist(list.x)
  y <- seq(0.4, 0.9, length.out = 50)
  new.data <- data.frame(beta_h = list.x, beta_l = rep(y, 50))
  sp9.predict <- predict(spline_sp9,  new.data)
  # surface(object = spline_sp9,extrap = TRUE,type = 'C', zlim = c(0.45, 0.55),
  #         main = "SP9 at various beta_l and beta_h")
  
  restrict <- which(sp9.predict > 0.4 & sp9.predict < 0.6)
  sp9.predict <- sp9.predict[which(sp9.predict > 0.4 & sp9.predict < 0.6)]
  new.data <- new.data[restrict,]
  
  ##predict sp9.h
  travel_1 <- data.frame(cbind(beta_h, beta_l, sp9.h))
  names(travel_1) <- c('beta_h', 'beta_l', 'sp9.h')
  spline_sp9.h = Tps(x = as.matrix(travel_1[c("beta_h","beta_l")]),
                     Y = matrix(travel_1["sp9.h"][,1]),
                     scale.type = "unscaled",lambda = 0)
  sp9.h.predict <- predict(spline_sp9.h,  new.data)
  # surface(object = spline_sp9.h,extrap = TRUE,type = 'C', zlim = range(sp9.h.predict),
  #         main = "SP9 at various beta_l and beta_h")
  
  ##predict sp9.l
  travel_1 <- data.frame(cbind(beta_h, beta_l, sp9.l))
  names(travel_1) <- c('beta_h', 'beta_l', 'sp9.l')
  spline_sp9.l = Tps(x = as.matrix(travel_1[c("beta_h","beta_l")]),
                     Y = matrix(travel_1["sp9.l"][,1]),
                     scale.type = "unscaled",lambda = 0)
  sp9.l.predict <- predict(spline_sp9.l,  new.data)
  # surface(object = spline_sp9.l,extrap = TRUE,type = 'C', zlim = range(sp9.l.predict),
  #         main = "SP9 at various beta_l and beta_h")
  
  sums.mat <- cbind((sp9.predict - 0.5), (sp9.h.predict - 0.2), (sp9.l.predict - 0.8))
  sums.mat <- data.frame(abs(sums.mat))
  sums <- sums.mat[,1] + sums.mat[,2] + sums.mat[,3]
  output <- c(new.data[which.min(sums),], sp9.predict[which.min(sums)], 
              sp9.l.predict[which.min(sums),], sp9.h.predict[which.min(sums),])
  output <- unlist(output)
  names(output) <- c('beta_h', 'beta_l', 'sp9', 'sp9.l', 'sp9.h')
  output <- as.vector(output)
  outputs.list[[i]] <- output
}

###plot the different SP9s
{
sp9.mat <- do.call(rbind, outputs.list)
colnames(sp9.mat) <- c('beta_h', 'beta_l', 'sp9', 'sp9.l', 'sp9.h')

par(mfrow = c(1,2))

par(mar = c(5, 5.5, 4, 2))
plot(x = rev(1- travel_levels), y = rev(sp9.mat[,1]), pch = 16, ylim  = c(0,1),
     ylab = "Transmission \nCoefficient", xlab = "Travel Parameter",yaxt = 'none', yaxs = 'i',
     main = "Adjusted\ntransmission\ncoefficient", col = 'blue')
axis(2, seq(0,1, by = 0.05), at = seq(0,1, by = 0.05), las = 2, cex.axis = 0.8)
points(x = rev(1- travel_levels), y = rev(sp9.mat[,2]), pch = 16, col = 'red')
legend("bottomright", legend = c(expression(beta['Low \ntransmission']), expression(beta['High \ntransmission'])), 
       col = c('blue', 'red'), pch = rep(16,3), cex = 0.8)

plot(x = rev(1- travel_levels), y = rev(sp9.mat[,3]), pch = 16, ylim  = c(0,1),
     ylab = "SP9", xlab = "Travel Parameter", yaxt = 'none', 
     main = "SP9 resulting from \nadjusted transmission coefficients",  yaxs = "i")
axis(2, seq(0,1, by = 0.05), at = seq(0,1, by = 0.05), las = 2, cex.axis = 0.8)
points(x = rev(1- travel_levels), y = rev(sp9.mat[,4]), pch = 16, col = 'red')
points(x = rev(1- travel_levels), y = rev(sp9.mat[,5]), pch = 16, col = 'blue')
abline(h = mean(sp9.mat[,4]), col = 'red')
abline(h = mean(sp9.mat[,3]), col = 'black')
abline(h = mean(sp9.mat[,5]), col = 'blue')
legend("bottomright", legend = c('SP9', expression('SP9'['Low \ntransmission']), expression('SP9'['High \ntransmission'])), 
       col = c('black', 'blue', 'red'), pch = rep(16,3), cex = 0.8)
text(x = 0.0015, y = 0.27, 'SP9 = 0.2481', col = 'blue', cex = 0.7)
text(x = 0.0015, y = 0.57, 'SP9 = 0.5488', cex = 0.7)
text(x = 0.0015, y = 0.87, 'SP9 = 0.8494', col = 'red', cex = 0.7)

mtext(expression(bold('A')), side=1, line=2.3, at= - 0.0775, cex = 2)
mtext(expression(bold('B')), side=1, line=2.3, at= - 0.007, cex = 2)

}






beta_h <- c()
beta_l <- c()
for(i in 1:1){
  beta_h[i] <- output[1]
  beta_l[i] <- output[2]
}
parms.mat <- cbind(beta_h, beta_l, travel_levels)

parms.list <- list()
for(i in 1:20){
  parms.list[[i]] <- parms.mat
}
parms.mat <- do.call(rbind, parms.list)
vac <- seq(0, 1, length.out = 20)
vac.add <- cbind(vac, rev(vac))
colnames(vac.add) <- c('vac_h', 'vac_l')




new.parms.mat <- cbind(parms.mat, vac.add)
setwd('~/Desktop/topics_in_ecology/')
save(new.parms.mat, file = "parms.mat.simp_model_SKEW_L.RData")





