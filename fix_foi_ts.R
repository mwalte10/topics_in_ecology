############################
#CRC info
############################
args = commandArgs(TRUE)
input = as.numeric(args[1])

FOI_list.h.vac <- list()
FOI_list.l.vac <- list()
FOI_list.h <- list()
FOI_list.l <- list()



  load(paste('FOI_output_list_',input, '.RData', sep = ''))
  w <- c()
  x <- c()
  y <- c()
  z <- c()
  for(i in 1:length(FOI_output[[1]])){
    index <- c(1:i)
    w[i] <- sum(diffinv(FOI_output[[1]][index])) / length(index)
    x[i] <- sum(diffinv(FOI_output[[2]][index])) / length(index)
    y[i] <- sum(diffinv(FOI_output[[3]][index])) / length(index)
    z[i] <- sum(diffinv(FOI_output[[4]][index])) / length(index)
  }
  FOI_list.h.vac <- w
  FOI_list.l.vac <- x
  FOI_list.h <- y
  FOI_list.l <- z



save(FOI_list.h.vac, file = paste('FOI_list.h.vac_', input, '.RData', sep = ''))
save(FOI_list.l.vac, file = paste('FOI_list.l.vac_', input, '.RData', sep = ''))
save(FOI_list.h, file = paste('FOI_list.h_', input, '.RData', sep = ''))
save(FOI_list.l, file = paste('FOI_list.l_', input, '.RData', sep = ''))
