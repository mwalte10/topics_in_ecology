setwd("~/Desktop/outputs/vaccination")

vac_h <- rep(seq(0.1, 0.9, length.out = 20), 20)
vac_l <- rep(NA, 400)
for(i in 1:20){
  j <- 20* (i - 1)
  vac_l[(1:20) + j] <- rep(vac_h[i], 20)
  remove(i)
  remove(j)
}
vac <- cbind(vac_h, vac_l)

# cases = list()
# cases_l = c()
# cases_h = c()
infections = c()
infections.l = c()
infections.h = c()

for(i in 1:312){
  #load(file = paste('cases.ti_', i, '.RData', sep = ''))
 #  load(file = paste('cases.l.ti_', i, '.RData', sep = ''))
 #  load(file = paste('cases.h.ti_', i, '.RData', sep = ''))
 #  cases_l[i] = cases.l[length(cases.l)]
 #  cases_h[i] = cases.h[length(cases.h)]
 # cases[[i]] = cases
 load(file = paste('track.infected.ti_', i, '.RData', sep = ''))
 load(file = paste('track.l.ti_', i, '.RData', sep = ''))
 load(file = paste('track.h.ti_', i, '.RData', sep = ''))
 infections[i] <- track_infected[length(track_infected)]
 infections.l[i] <- track_l[length(track_l)]
 infections.h[i] <- track_h[length(track_h)]
}
remove(cases.h) ; remove(cases.l) ; remove(track_infected) ; remove(track_l) ; remove(track_h)

load('cases.nv_1.RData'); cases.nv <- cases[length(cases)]; remove(cases)
load('cases.nv.h_1.RData'); cases.h.nv <- cases.h[length(cases.h)]; remove(cases.h)
load('cases.nv.l_1.RData'); cases.l.nv <- cases.l[length(cases.l)]; remove(cases.l)
load('track.infected.nv_1.RData'); infections.nv <- track_infected[length(track_infected)]; remove(track_infected)
load('track.nv.l_1.RData'); infections.nv.l <- track_l[length(track_l)]; remove(track_l)
load('track.nv.h_1.RData'); infections.nv.h <- track_h[length(track_h)]; remove(track_h)

infections_averted <- c()
infections_averted.l <- c()
infections_averted.h <- c()
for(i in 1:312){
  infections.t <- infections[i]
  infections_averted[i] <- (((infections.nv - infections.t) / infections.nv) * 100)
  infections.th <- infections.h[i]
  infections_averted.h[i] <- (((infections.nv.h - infections.th) / infections.nv.h) * 100)
  infections.tl <- infections.l[i]
  infections_averted.l[i] <- (((infections.nv.l - infections.tl) / infections.nv.l) * 100)
}
infections_averted <- c(infections_averted, rep(NA, 88))
infections.mat <- matrix(NA, nrow = 20, ncol = 20, byrow = TRUE)
for(i in 1:20){
  j <- i - 1
  infections.mat[i,] <- infections_averted[(1:20) + j * (1:20)]
}

row.names(infections.mat) <- seq(0.1, 0.9, length.out = 20)
colnames(infections.mat) <- seq(0.1, 0.9, length.out = 20)

infections.mat.l <- matrix(infections_averted.l, nrow = 20, ncol = 20, byrow = FALSE)
#rows are vaccination in pop.h
row.names(infections.mat.l) <- seq(0.1, 0.9, length.out = 20)
#columns are vaccination in pop.l
colnames(infections.mat.l) <- seq(0.1, 0.9, length.out = 20)

infections.mat.h <- matrix(infections_averted.h, nrow = 20, ncol = 20, byrow = TRUE)
#rows are vaccination in pop.h
row.names(infections.mat.h) <- seq(0.1, 0.9, length.out = 20)
#columns are vaccination in pop.l
colnames(infections.mat.h) <- seq(0.1, 0.9, length.out = 20)

cases_averted <- rep(NA, 400)
cases_averted.l <- rep(NA, 312)
cases_averted.h <- rep(NA, 312)
for(i in 1:312){
  # cases.v <- cases[[i]]
  # cases_averted[i] <- (((cases.nv - cases.v) / cases.nv) * 100)[length(cases.v)]
  cases.l.v <- cases_l[i]
  cases_averted.l[i] <- (((cases.l.nv - cases.l.v) / cases.l.nv) * 100)
  cases.h.v <- cases_h[i]
  cases_averted.h[i] <- (((cases.h.nv - cases.h.v) / cases.h.nv) * 100)
}
remove(cases.l.v) ; remove(cases.h.v)

cases.mat <- matrix(cases_averted, nrow = 20, ncol = 20, byrow = TRUE)
row.names(cases.mat) <- seq(0.1, 0.9, length.out = 20)
colnames(cases.mat) <- seq(0.1, 0.9, length.out = 20)

cases.mat.l <- matrix(NA, nrow = 20, ncol = 20)
cases_averted.l <- c(cases_averted.l, rep(NA, 88))
for(i in 1:20){
  j <- i - 1
  cases.mat.l[i,] <- cases_averted.l[(1:20) + j * (1:20)]
}
#rows are vaccination in pop.h
row.names(cases.mat.l) <- seq(0.1, 0.9, length.out = 20)
#columns are vaccination in pop.l
colnames(cases.mat.l) <- seq(0.1, 0.9, length.out = 20)

cases_averted.h <- c(cases_averted.h, rep(NA, 88))
cases.mat.h <- matrix(cases_averted.h, nrow = 20, ncol = 20, byrow = TRUE)
#rows are vaccination in pop.h
row.names(cases.mat.h) <- seq(0.1, 0.9, length.out = 20)
#columns are vaccination in pop.l
colnames(cases.mat.h) <- seq(0.1, 0.9, length.out = 20)

vac <- seq(0.1, 0.9, length.out = 20)
{par(mfrow = c(2,3))
  hist3D(x = vac, y = vac, z = infections.mat.h,
         xlab = "Beta_h Vaccination", ylab = "Beta_l Vaccination", z = "Infections.h Averted",
         main = "Infections.h Averted")
  hist3D(x = vac, y = vac, z = infections.mat, 
         xlab = "Beta_h Vaccination", ylab = "Beta_l Vaccination", z = "Infections Averted",
         main = "Infections Averted")
  hist3D(x = vac, y = vac, z = infecionts.mat.l,
         xlab = "Beta_h Vaccination", ylab = "Beta_l Vaccination", z = "Infections.l Averted",
         main = "Infections.l Averted")
  hist3D(x = vac, y = vac, z = cases.mat.h,
         xlab = "Beta_h Vaccination", ylab = "Beta_l Vaccination", z = "Cases.h Averted",
         main = "Cases.h Averted")
  hist3D(x = vac, y = vac, z = cases.mat,
         xlab = "Beta_h Vaccination", ylab = "Beta_l Vaccination", z = "Cases Averted",
         main = "Cases Averted")
  hist3D(x = vac, y = vac, z = cases.mat.l,
         xlab = "Beta_h Vaccination", ylab = "Beta_l Vaccination", z = "Cases.l Averted",
         main = "Cases.l Averted")
}



