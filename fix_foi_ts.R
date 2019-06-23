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
    w <- sum(diffinv(FOI_output[[1]][index])) / length(index)
    x <- sum(diffinv(FOI_output[[2]][index])) / length(index)
    y <- sum(diffinv(FOI_output[[3]][index])) / length(index)
    z <- sum(diffinv(FOI_output[[4]][index])) / length(index)
  }
  FOI_list.h.vac <- w
  FOI_list.l.vac <- x
  FOI_list.h <- y
  FOI_list.l <- z



save(FOI_list.h.vac, file = 'FOI_list.h.vac.RData')
save(FOI_list.l.vac, file = 'FOI_list.l.vac.RData')
save(FOI_list.h, file = 'FOI_list.h.RData')
save(FOI_list.l, file = 'FOI_list.l.RData')
