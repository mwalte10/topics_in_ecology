############################
#CRC info
############################
args = commandArgs(TRUE)
input = as.numeric(args[1])
beta.vec = seq(from = 0, to = 1, by = .1)
beta = beta.vec[input]

library(deSolve)

############################
#Initial condidtions and parameters
############################
percentage_vec <- c(rep(1.8,5), rep(1.8,5), rep(1.84,5), rep(1.86,5),
                    17.2, 15, 12.7, 8.8, 5.5, 2.8, 1.3, 0.2)
percentage_vec <- percentage_vec / 100
initial_conditions <- as.data.frame(matrix(NA, nrow = 28*4, ncol = 3))
initial_conditions[,2] <- rep(1:28,4)
for(i in 1:4){
  x <- i - 1
  initial_conditions[x*(28) + (1:28),1] <- rep(i,28)
}
initial_conditions[,3] <- rep(percentage_vec,4)

susceptible_init_h <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 112))
susceptible_init_h[,3] <- c(initial_conditions[1:28,3] * 6 *10^6, initial_conditions[29:56,3] * 0, 
                            initial_conditions[57:84,3] * 0, initial_conditions[85:112,3] * 0)
susceptible_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 112))
susceptible_init_l[,3] <- c(initial_conditions[1:28,3] * 6 *10^6, initial_conditions[29:56,3] * 0, 
                            initial_conditions[57:84,3] * 0, initial_conditions[85:112,3] * 0)

infected_init_h <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 112))
infected_init_h[,3] <- c(initial_conditions[1:28,3] * 1, initial_conditions[29:56,3] * 0, 
                         initial_conditions[57:84,3] * 0, initial_conditions[85:112,3] * 0)
infected_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 112))
infected_init_l[,3] <- c(initial_conditions[1:28,3] * 1, initial_conditions[29:56,3] * 0, 
                         initial_conditions[57:84,3] * 0, initial_conditions[85:112,3] * 0)

recovered_init_h <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 112))
recovered_init_h[,3] <- c(initial_conditions[1:28,3] * 0, initial_conditions[29:56,3] * 0, 
                          initial_conditions[57:84,3] * 0, initial_conditions[85:112,3] * 0)
recovered_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 112))
recovered_init_l[,3] <- c(initial_conditions[1:28,3] * 0, initial_conditions[29:56,3] * 0, 
                          initial_conditions[57:84,3] * 0, initial_conditions[85:112,3] * 0)


susceptible_h <- function(exposure){
  x <- exposure - 1
  susceptible <- c(susceptible_init_h[(x * 28) + 1:28,3])
}
susceptible_total_h <- sum(susceptible_h(1) + susceptible_h(2) + susceptible_h(3) + susceptible_h(4))
susceptible_l <- function(exposure){
  x <- exposure - 1
  susceptible <- c(susceptible_init_l[(x * 28) + 1:28,3])
}
susceptible_total_l <- sum(susceptible_l(1) + susceptible_l(2) + susceptible_l(3) + susceptible_l(4))

infected_h <- function(exposure){
  x <- exposure - 1
  infected <- c(infected_init_h[(x * 28) + 1:28,3])
} 
infected_total_h <- sum(infected_h(1) + infected_h(2) + infected_h(3) + infected_h(4))
infected_l <- function(exposure){
  x <- exposure - 1
  infected <- c(infected_init_l[(x * 28) + 1:28,3])
} 
infected_total_l <- sum(infected_l(1) + infected_l(2) + infected_l(3) + infected_l(4))

recovered_h <- function(exposure){
  x <- exposure - 1
  recovered <- c(recovered_init_h[(x * 28) + 1:28,3])
} 
recovered_total_h <- sum(recovered_h(1) + recovered_h(2) + recovered_h(3) + recovered_h(4))
recovered_l <- function(exposure){
  x <- exposure - 1
  recovered <- c(recovered_init_l[(x * 28) + 1:28,3])
} 
recovered_total_l <- sum(recovered_l(1) + recovered_l(2) + recovered_l(3) + recovered_l(4))

population_h <- sum(susceptible_total_h + infected_total_h + recovered_total_h)
population_l <- sum(susceptible_total_l + infected_total_l + recovered_total_l)


native <- c(rep(1, 7), rep(0.86, 9), rep(0.842, 4), 0.814, 0.7676, 0.7784, rep(0.809, 5))

parms <- c(beta_h = beta,
           beta_l = beta,
           gamma = 1/4,
           sigma = 1/(365 * 1.2),
           mu = 19 / (1000 * 365),
           delta = c((0.013 / 365), rep(0.001 / 365, 4), rep(0, 10), rep(0.001, 5), 
                     0.002 / 365, 0.0025 / 365, 
                     0.004 / 365, 0.0085 / 365,
                     0.0175 / 365, 0.0425 / 365,
                     0.114 / 365, 0.151 / 365),
           age_window = c(rep(1, 21), rep(10, 7)),
           native = c(rep(1, 7), rep(0.86, 9), rep(0.842, 4), 0.814, 0.7676, 0.7784, rep(0.809, 5)),
           travel <- 1 - native)

############################
#MODEL
############################
model <- function(t, y, parms){
  
  beta_h <- parms[1]
  beta_l <- parms[2]
  gamma <- parms[3]
  sigma <- parms[4]
  mu <- parms[5]
  delta <- parms[6:33]
  age_window <- parms[34:61]
  native <- parms[62:89]
  travel <- parms[90:117]
  
  #Gen pop
  S1_h <- y[1:28]
  I1_h <- y[29:56]
  R1_h <- y[57:84]
  S2_h <- y[85:112]
  I2_h <- y[113:140]
  R2_h <- y[141:168]
  S3_h <- y[169:196]
  I3_h <- y[197:224]
  R3_h <- y[225:252]
  S4_h <- y[253:280]
  I4_h <- y[281:308]
  R4_h <- y[309:336]
  
  
  S1_l <- y[(1:28) + 336]
  I1_l <- y[(29:56) + 336]
  R1_l <- y[(57:84) + 336]
  S2_l <- y[(85:112) + 336]
  I2_l <- y[(113:140) + 336]
  R2_l <- y[(141:168) + 336]
  S3_l <- y[(169:196) + 336]
  I3_l <- y[(197:224) + 336]
  R3_l <- y[(225:252) + 336]
  S4_l <- y[(253:280) + 336]
  I4_l <- y[(281:308) + 336]
  R4_l <- y[(309:336) + 336]
  
  
  ###Tracking infected individuals
  I_total <- y[673]
  
  
  infected_total_h <- sum(sum(I1_h), sum(I2_h), sum(I3_h), sum(I4_h))
  population_h <- y[1:336]
  pop_h <- sum(population_h)
  infected_total_l <- sum(sum(I1_l), sum(I2_l), sum(I3_l), sum(I4_h))
  population_l <- y[337:672]
  pop_l <- sum(population_l)
  
  effective_population_h <- population_h + rep(travel,12) * population_l
  effective_population_h <- sum(effective_population_h)
  effective_population_l <- population_l + rep(travel,12) * population_h
  effective_population_l <- sum(effective_population_l)
  
  #first infection
  dS1_h <-  
    mu * c(pop_h, rep(0,27)) +
    c(0, 1/365 / head(age_window, -1) * head(S1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S1_h, -1), 0) -
    native * S1_h * (beta_h * infected_total_h / effective_population_h) -
    travel * S1_h * (beta_l * infected_total_l / effective_population_l) -
    S1_h * delta 
  dS1_l <-  
    mu * c(pop_l, rep(0,27)) +
    c(0, 1/365 / head(age_window, -1) * head(S1_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S1_l, -1), 0) -
    travel * S1_l * (beta_h * infected_total_h / effective_population_h) -
    native * S1_l * (beta_l * infected_total_l / effective_population_l) -
    S1_l * delta 
  
  dI1_h <- 
    native * S1_h * (beta_h * infected_total_h / effective_population_h) +
    travel * S1_h * (beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I1_h, -1), 0) -
    I1_h * gamma -
    I1_h * delta
  dI1_l <- 
    travel * S1_l * (beta_h * infected_total_h / effective_population_h) +
    native * S1_l * (beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I1_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I1_l, -1), 0) -
    I1_l * gamma -
    I1_l * delta
  
  I_total <- native * S1_h * (beta_h * infected_total_h / effective_population_h) +
    travel * S1_h * (beta_l * infected_total_l / effective_population_l) +
    travel * S1_l * (beta_h * infected_total_h / effective_population_h) +
    native * S1_l * (beta_l * infected_total_l / effective_population_l)
  
  dR1_h <-
    I1_h * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(R1_h, -1), 0) -
    R1_h * sigma -
    R1_h * delta 
  dR1_l <-
    I1_l * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R1_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R1_l, -1), 0) -
    R1_l * sigma -
    R1_l * delta 
  
  #second infection
  dS2_h <- 
    R1_h * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_h, -1), 0) -
    native * S2_h * (0.75 * beta_h * infected_total_h / effective_population_h) -
    travel * S2_h * (0.75 * beta_l * infected_total_l / effective_population_l) -
    S2_h * delta  
  dS2_l <- 
    R1_l * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_l, -1), 0) -
    travel * S2_l * (0.75 * beta_h * infected_total_h / effective_population_h) -
    native * S2_l * (0.75 * beta_l * infected_total_l / effective_population_l) -
    S2_l * delta 
  
  dI2_h <- 
    native * S2_h * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S2_h * (0.75 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I2_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_h, -1), 0) -
    I2_h * gamma -
    I2_h * delta
  dI2_l <- 
    travel * S2_l * (0.75 * beta_h * infected_total_h / effective_population_h) +
    native * S2_l * (0.75 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_l, -1), 0) -
    I2_l * gamma -
    I2_l * delta
  
  I_total <- 
    native * S2_h * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S2_h * (0.75 * beta_l * infected_total_l / effective_population_l) +
    travel * S2_l * (0.75 * beta_h * infected_total_h / effective_population_h) +
    native * S2_l * (0.75 * beta_l * infected_total_l / effective_population_l) 
  
  dR2_h <-  
    I2_h * gamma +  
    c(0, 1/365 / head(age_window, -1) * head(R2_h, -1)) -
    c(1/365 / head(age_window, -1) * head(R2_h, -1), 0) -
    R2_h * sigma -
    R2_h * delta 
  dR2_l <-  
    I2_l * gamma +  
    c(0, 1/365 / head(age_window, -1) * head(R2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R2_l, -1), 0) -
    R2_l * sigma -
    R2_l * delta 
  
  #third infection
  dS3_h <- 
    R2_h * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_h, -1), 0) -
    native * S3_h * (0.5 * beta_h * infected_total_h / effective_population_h) -
    travel * S3_h * (0.5 * beta_l * infected_total_l / effective_population_l) -
    S3_h * delta
  dS3_l <- 
    R2_l * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_l, -1), 0) -
    travel * S3_l * (0.5 * beta_h * infected_total_h / effective_population_h) -
    native * S3_l * (0.5 * beta_l * infected_total_l / effective_population_l) -
    S3_l * delta 
  
  dI3_h <- 
    native * S3_h * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_h * (0.5 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_h, -1), 0) -
    I3_h * gamma -
    I3_h * delta
  dI3_l <- 
    travel * S3_l * (0.5 * beta_h * infected_total_h / effective_population_h) +
    native * S3_l * (0.5 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_l, -1), 0) -
    I3_l * gamma -
    I3_l * delta
  
  I_total <- 
    native * S3_h * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_h * (0.5 * beta_l * infected_total_l / effective_population_l) +
    travel * S3_l * (0.5 * beta_h * infected_total_h / effective_population_h) +
    native * S3_l * (0.5 * beta_l * infected_total_l / effective_population_l) 
  
  dR3_h <- 
    I3_h * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_h, -1), 0) -
    R3_h * sigma -
    R3_h * delta 
  dR3_l <- 
    I3_l * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_l, -1), 0) -
    R3_l * sigma -
    R3_l * delta 
  
  #fourth infection
  dS4_h <- 
    R3_h * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S4_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_h, -1), 0) -
    native * S4_h * (0.25 * beta_h * infected_total_h / effective_population_h) -
    travel * S4_h * (0.25 * beta_l * infected_total_l / effective_population_l) -
    S4_h * delta 
  dS4_l <- 
    R3_l * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_l, -1), 0) -
    travel * S4_l * (0.25 * beta_h * infected_total_h / effective_population_h) -
    native * S4_l * (0.25 * beta_l * infected_total_l / effective_population_l) -
    S4_l * delta 
  
  dI4_h <- 
    native * S4_h * (0.25 * beta_h * infected_total_h / effective_population_h) +
    travel * S4_h * (0.25 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I4_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_h, -1), 0) -
    I4_h * gamma -
    I4_h * delta
  dI4_l <- 
    travel * S4_l * (0.25 * beta_h * infected_total_h / effective_population_h) +
    native * S4_l * (0.25 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_l, -1), 0) -
    I4_l * gamma -
    I4_l * delta
  
  I_total <- 
    native * S4_h * (0.25 * beta_h * infected_total_h / effective_population_h) +
    travel * S4_h * (0.25 * beta_l * infected_total_l / effective_population_l) +
    travel * S4_l * (0.25 * beta_h * infected_total_h / effective_population_h) +
    native * S4_l * (0.25 * beta_l * infected_total_l / effective_population_l) 
  
  dR4_h <- 
    I4_h * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R4_h, -1)) -
    c(1/365 / head(age_window, -1) * head(R4_h, -1), 0) -
    R4_h * delta 
  dR4_l <- 
    I4_l * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R4_l, -1), 0) -
    R4_l * delta 

  I_total <- sum(I_total)
  
  list(c(dS1_h, dI1_h, dR1_h,
         dS2_h, dI2_h, dR2_h,
         dS3_h, dI3_h, dR3_h,
         dS4_h, dI4_h, dR4_h,
         dS1_l, dI1_l, dR1_l,
         dS2_l, dI2_l, dR2_l,
         dS3_l, dI3_l, dR3_l,
         dS4_l, dI4_l, dR4_l,
         I_total))
  
}

############################
#RUN MODEL
############################
y_init <- c(susceptible_h(1), infected_h(1), recovered_h(1),
            susceptible_h(2), infected_h(2), recovered_h(2),
            susceptible_h(3), infected_h(3), recovered_h(3),
            susceptible_h(4), infected_h(4), recovered_h(4),
            susceptible_l(1), infected_l(1), recovered_l(1),
            susceptible_l(2), infected_l(2), recovered_l(2),
            susceptible_l(3), infected_l(3), recovered_l(3),
            susceptible_l(4), infected_l(4), recovered_l(4),
            0)
years = 50
times <- seq(from = 0, to = 365 * years, by = .1)
out <- ode(times = times, y = y_init, func = model, parms = parms)


############################
#OUTPUT
############################
out_last <- out[nrow(out),]
save(out_last, file = paste('output.nv_', input, '.RData', sep = ''))


track_infected <- out[,ncol(out)]
save(track_infected, file = paste('track.infected.nv_', input, '.RData', sep = ''))


###############################
#plot proportions SIR over time 
###############################

#png(filename = paste('prop.h.SIR.nv_', input, '.png', sep = ''))
#ts <- matrix(NA, nrow = 12, ncol = years * 365 * 10)
#for(j in 1:12){
#  for(i in 1:(years* 365 * 10)){
#    x <- j - 1
#    ts[j,i] <- sum(out[i,((x * 28) + 2):((x * 28) + 29)]) / sum(out[i,2:337])
#  }
#}
#rownames(ts) <- c("dS1_h", "dI1_h", "dR1_h",
#                  "dS2_h", "dI2_h", "dR2_h",
#                  "dS3_h", "dI3_h", "dR3_h",
#                  "dS4_h", "dI4_h", "dR4_h")

#par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
########
##set colors
########
#colors <- c(rgb(0, 191/255, 255/255, alpha = 0.25), rgb(0, 191/255, 255/255, alpha = 0.5), rgb(0, 191/255, 255/255, alpha = 1),
#            rgb(0, 0, 255/255, alpha = 0.25), rgb(0, 0, 255/255, alpha = 0.5), rgb(0, 0, 255/255, alpha = 1),
#            rgb(0, 139/255, 0, alpha = 0.25), rgb(0, 139/255, 0, alpha = 0.5), rgb(0, 139/255, 0, alpha = 1),
#            rgb(104/255, 34/255, 139/255, alpha = 0.25), rgb(104/255, 34/255, 139/255, alpha = 0.5), rgb(104/255, 34/255, 139/255, alpha = 1))
#densities <- c(rep(100, 12), rep(10,10))
#colors_c <- c(rgb(0, 191/255, 255/255, alpha = 0.25), rgb(0, 191/255, 255/255, alpha = 0.5), rgb(0, 191/255, 255/255, alpha = 1),
#              rgb(0, 0, 255/255, alpha = 0.25), rgb(0, 0, 255/255, alpha = 0.5), rgb(0, 0, 255/255, alpha = 1),
#              rgb(0, 139/255, 0, alpha = 0.25), rgb(0, 139/255, 0, alpha = 0.5), rgb(0, 139/255, 0, alpha = 1),
#              rgb(104/255, 34/255, 139/255, alpha = 0.25), rgb(104/255, 34/255, 139/255, alpha = 0.5), rgb(104/255, 34/255, 139/255, alpha = 1))
#colors_f <- c(colors_c, colors_c[3:length(colors_c)])    
#######
#{par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#  barplot(ts.new, col = colors_f, 
#          border=colors, space=0.04, main ="Proportion of Each Category, High SES", 
#          xlab = "Timestep", angle = angles)
#  legend("topright",  inset=c(-0.1,0),
#         c("S1", "I1", "R1", 
#           "S2", "I2", "R2",
#           "S3", "I3", "R3",
#           "S4", "I4", "R4"), 
#         fill=colors_f[1:12], horiz=FALSE, cex=0.8)
#  segments(x0 = (years_vac * 10 * 365), x1 = (years_vac * 10 * 365), y0 = 0, y1 = 1)
#}
#dev.off()

#png(filename = paste('prop.l.SIR.nv_', input, '.png', sep = ''))
#ts_l <- matrix(NA, nrow = 12, ncol = years * 365 * 10)
#for(j in 1:22){
#  for(i in 1:(years* 365 * 10)){
#    x <- j - 1
#    ts_l[j,i] <- sum(out[i,((x * 28) + 2+ 336):((x * 28) + 29 + 336)]) / sum(out[i,337:673])
#  }
#}
#rownames(ts_l) <- c("dS1_l", "dI1_l", "dR1_l",
#                    "dS2_l", "dI2_l", "dR2_l",
#                    "dS3_l", "dI3_l", "dR3_l",
#                    "dS4_l", "dI4_l", "dR4_l")
#{par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#  barplot(ts_l, col = topo.colors(12), 
#          border=colors, space=0.04, main ="Proportion of Each Category, Low SES",
#          xlab = "Timestep")
#  legend("topright",  inset=c(-0.1,0),
#         c("S1", "I1", "R1", 
#           "S2", "I2", "R2",
#           "S3", "I3", "R3",
#           "S4", "I4", "R4"), 
#         fill=colors[1:12], horiz=FALSE, cex=0.8)
#  segments(x0 = (years_vac * 10 * 365), x1 = (years_vac * 10 * 365), y0 = 0, y1 = 1)}
#dev.off()

#####################
#SP9 over time
#####################
nines <- 11 + c(1, 29, 57,
                85, 113, 141,
                169, 197, 225, 
                253, 281, 309)
#sp9 <- function(ts, low){
#  x <- out[ts, nines[1] + low] / sum(out[ts, (nines + low)])
#  return(1 - x)
#} 

sp9.h <- rep(NA, years * 10 * 365)
for(i in 1:(years * 10 *365)){
  sp9.h[i] <- 1 - (out[i, nines[1]] / sum(out[i, nines]))
}
save(sp9.h, file = paste('sp9.h.nv_', input, '.RData', sep = ''))



sp9.l <- rep(NA, years * 10 * 365)
for(i in 1:(years * 10 *365)){
  sp9.l[i] <- 1 - (out[i, (nines[1] + 336)] / sum(out[i, (nines + 336)]))
}
save(sp9.l, file = paste('sp9.l.nv_', input, '.RData', sep = ''))



