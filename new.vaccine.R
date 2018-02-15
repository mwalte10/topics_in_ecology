############################
#CRC info
############################
args = commandArgs(TRUE)
input = as.numeric(args[1])
vac.vec = seq(from = 0, to = 1, by = .2)
vac = vac.vec[input]

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

parms <- c(beta_h = 0.3447,
           beta_l = 0.4529,
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
           travel <- 1 - native,
           vac_h = vac,
           vac_l = vac)

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
  vac_h <- c(rep(0,8), ifelse((t>(365*30)), parms[118] / 365, 0), rep(0,19))
  vac_l <- c(rep(0,8), ifelse((t>(365*30)), parms[119] / 365, 0), rep(0,19))
  
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

  #Vaccinated
  R1_h.v <- y[337:364]
  S2_h.v <- y[365:392]
  I2_h.v <- y[393:420]
  R2_h.v <- y[421:448]
  S3_h.v <- y[449:476]
  I3_h.v <- y[477:504]
  R3_h.v <- y[505:532]
  S4_h.v <- y[533:560]
  I4_h.v <- y[561:588]
  R4_h.v <- y[589:616]
  
  S1_l <- y[617:644]
  I1_l <- y[645:672]
  R1_l <- y[673:700]
  S2_l <- y[701:728]
  I2_l <- y[729:756]
  R2_l <- y[757:784]
  S3_l <- y[785:812]
  I3_l <- y[813:840]
  R3_l <- y[841:868]
  S4_l <- y[869:896]
  I4_l <- y[897:924]
  R4_l <- y[925:952]
  
  R1_l.v <- y[953:980]
  S2_l.v <- y[981:1008]
  I2_l.v <- y[1009:1036]
  R2_l.v <- y[1037:1064]
  S3_l.v <- y[1065:1092]
  I3_l.v <- y[1093:1120]
  R3_l.v <- y[1121:1148]
  S4_l.v <- y[1149:1176]
  I4_l.v <- y[1177:1204]
  R4_l.v <- y[1205:1232]

  
  infected_total_h <- sum(sum(I1_h), sum(I2_h), sum(I3_h), sum(I4_h), sum(I2_h.v), sum(I3_h.v), sum(I4_h.v))
  population_h <- y[1:616]
  pop_h <- sum(population_h)
  infected_total_l <- sum(sum(I1_l), sum(I2_l), sum(I3_l), sum(I4_l), sum(I2_l.v), sum(I3_l.v), sum(I4_l.v))
  population_l <- y[617:1232]
  pop_l <- sum(population_l)
  
  effective_population_h <- population_h + rep(travel,22) * population_l
  effective_population_h <- sum(effective_population_h)
  effective_population_l <- population_l + rep(travel,22) * population_h
  effective_population_l <- sum(effective_population_l)
  
  #first infection
  dS1_h <-  
    mu * c(pop_h, rep(0,27)) +
    c(0, 1/365 / head(age_window, -1) * head(S1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S1_h, -1), 0) -
    native * S1_h * (beta_h * infected_total_h / effective_population_h) -
    travel * S1_h * (beta_l * infected_total_l / effective_population_l) -
    S1_h * delta -
    S1_h * vac_h
  dS1_l <-  
    mu * c(pop_l, rep(0,27)) +
    c(0, 1/365 / head(age_window, -1) * head(S1_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S1_l, -1), 0) -
    travel * S1_l * (beta_h * infected_total_h / effective_population_h) -
    native * S1_l * (beta_l * infected_total_l / effective_population_l) -
    S1_l * delta - 
    S1_l * vac_l
  
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
  
  dR1_h <-
    I1_h * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(R1_h, -1), 0) -
    R1_h * sigma -
    R1_h * delta - 
    R1_h * vac_h
  dR1_h.v <- 
    S1_h * vac_h +
    c(0, 1/365 / head(age_window, -1) * head(R1_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R1_h.v, -1), 0) -
    R1_h.v * delta -
    R1_h.v * sigma
  dR1_l <-
    I1_l * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R1_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R1_l, -1), 0) -
    R1_l * sigma -
    R1_l * delta +
    S1_l * vac_l - 
    R1_l * vac_l
  dR1_l.v <- 
    S1_l * vac_l +
    c(0, 1/365 / head(age_window, -1) * head(R1_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R1_l.v, -1), 0) -
    R1_l.v * delta -
    R1_l.v * sigma
  
  #second infection
  dS2_h <- 
    R1_h * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_h, -1), 0) -
    native * S2_h * (0.75 * beta_h * infected_total_h / effective_population_h) -
    travel * S2_h * (0.75 * beta_l * infected_total_l / effective_population_l) -
    S2_h * delta - 
    S2_h * vac_h
  dS2_h.v <-
    R1_h.v * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_h.v, -1), 0) -
    native * S2_h.v * (0.75 * beta_h * infected_total_h / effective_population_h) -
    travel * S2_h.v * (0.75 * beta_l * infected_total_l / effective_population_l) -
    S2_h.v * delta
  dS2_l <- 
    R1_l * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_l, -1), 0) -
    travel * S2_l * (0.75 * beta_h * infected_total_h / effective_population_h) -
    native * S2_l * (0.75 * beta_l * infected_total_l / effective_population_l) -
    S2_l * delta - 
    S2_l * vac_l
  dS2_l.v <-
    R1_l.v * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_l.v, -1), 0) -
    native * S2_l.v * (0.75 * beta_h * infected_total_h / effective_population_h) -
    travel * S2_l.v * (0.75 * beta_l * infected_total_l / effective_population_l) -
    S2_l.v * delta
    
  dI2_h <- 
    native * S2_h * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S2_h * (0.75 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I2_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_h, -1), 0) -
    I2_h * gamma -
    I2_h * delta
  dI2_h.v <-
    native * S2_h.v * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S2_h.v * (0.75 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I2_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_h.v, -1), 0) -
    I2_h.v * gamma -
    I2_h.v * delta
  dI2_l <- 
    travel * S2_l * (0.75 * beta_h * infected_total_h / effective_population_h) +
    native * S2_l * (0.75 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_l, -1), 0) -
    I2_l * gamma -
    I2_l * delta
  dI2_l.v <-
    native * S2_l.v * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S2_l.v * (0.75 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I2_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_l.v, -1), 0) -
    I2_l.v * gamma -
    I2_l.v * delta
  
  dR2_h <-  
    I2_h * gamma +  
    c(0, 1/365 / head(age_window, -1) * head(R2_h, -1)) -
    c(1/365 / head(age_window, -1) * head(R2_h, -1), 0) -
    R2_h * sigma -
    R2_h * delta -
    R2_h * vac_h
  dR2_h.v <-
    I2_h.v * gamma +  
    c(0, 1/365 / head(age_window, -1) * head(R2_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R2_h.v, -1), 0) -
    R2_h.v * sigma -
    R2_h.v * delta +
    R1_h.v * vac_h +
    S2_h.v * vac_h
  dR2_l <-  
    I2_l * gamma +  
    c(0, 1/365 / head(age_window, -1) * head(R2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R2_l, -1), 0) -
    R2_l * sigma -
    R2_l * delta + 
    R1_l * vac_l +
    S2_l * vac_l -
    R2_l * vac_l
  dR2_l.v <-
    I2_l.v * gamma +  
    c(0, 1/365 / head(age_window, -1) * head(R2_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R2_l.v, -1), 0) -
    R2_l.v * sigma -
    R2_l.v * delta +
    R1_l.v * vac_l +
    S2_l.v * vac_l
  
  #third infection
  dS3_h <- 
    R2_h * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_h, -1), 0) -
    native * S3_h * (0.5 * beta_h * infected_total_h / effective_population_h) -
    travel * S3_h * (0.5 * beta_l * infected_total_l / effective_population_l) -
    S3_h * delta -
    S3_h * vac_h
  dS3_h.v <-
    R2_h.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_h.v, -1), 0) -
    native * S3_h.v * (0.5 * beta_h * infected_total_h / effective_population_h) -
    travel * S3_h.v * (0.5 * beta_l * infected_total_l / effective_population_l) -
    S3_h.v * delta 
  dS3_l <- 
    R2_l * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_l, -1), 0) -
    travel * S3_l * (0.5 * beta_h * infected_total_h / effective_population_h) -
    native * S3_l * (0.5 * beta_l * infected_total_l / effective_population_l) -
    S3_l * delta -
    S3_l * vac_l
  dS3_l.v <-
    R2_l.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_l.v, -1), 0) -
    native * S3_l.v * (0.5 * beta_h * infected_total_h / effective_population_h) -
    travel * S3_l.v * (0.5 * beta_l * infected_total_l / effective_population_l) -
    S3_l.v * delta
  
  dI3_h <- 
    native * S3_h * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_h * (0.5 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_h, -1), 0) -
    I3_h * gamma -
    I3_h * delta
  dI3_h.v <-
    native * S3_h.v * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_h.v * (0.5 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_h.v, -1), 0) -
    I3_h.v * gamma -
    I3_h.v * delta
  dI3_l <- 
    travel * S3_l * (0.5 * beta_h * infected_total_h / effective_population_h) +
    native * S3_l * (0.5 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_l, -1), 0) -
    I3_l * gamma -
    I3_l * delta
  dI3_l.v <-
    native * S3_l.v * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_l.v * (0.5 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I3_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_l.v, -1), 0) -
    I3_l.v * gamma -
    I3_l.v * delta
  
  dR3_h <- 
    I3_h * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_h, -1), 0) -
    R3_h * sigma -
    R3_h * delta -
    R3_h * vac_h
  dR3_h.v <-
    I3_h.v * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_h.v, -1), 0) -
    R3_h.v * sigma + 
    R2_h.v * vac_h +
    S3_h.v * vac_h
  dR3_l <- 
    I3_l * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_l, -1), 0) -
    R3_l * sigma -
    R3_l * delta - 
    R3_l * vac_l
  dR3_l.v <-
    I3_l.v * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_l.v, -1), 0) -
    R3_l.v * sigma + 
    R2_l.v * vac_l +
    S3_l.v * vac_l
  
  #fourth infection
  dS4_h <- 
    R3_h * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S4_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_h, -1), 0) -
    native * S4_h * (0.25 * beta_h * infected_total_h / effective_population_h) -
    travel * S4_h * (0.25 * beta_l * infected_total_l / effective_population_l) -
    S4_h * delta - 
    S4_h * vac_h
  dS4_h.v <-
    R3_h.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S4_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_h.v, -1), 0) -
    native * S4_h.v * (0.5 * beta_h * infected_total_h / effective_population_h) -
    travel * S4_h.v * (0.5 * beta_l * infected_total_l / effective_population_l) -
    S4_h.v * delta 
  dS4_l <- 
    R3_l * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_l, -1), 0) -
    travel * S4_l * (0.25 * beta_h * infected_total_h / effective_population_h) -
    native * S4_l * (0.25 * beta_l * infected_total_l / effective_population_l) -
    S4_l * delta - 
    S4_l * vac_l
  dS4_l.v <-
    R3_l.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S4_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_l.v, -1), 0) -
    native * S4_l.v * (0.5 * beta_h * infected_total_h / effective_population_h) -
    travel * S4_l.v * (0.5 * beta_l * infected_total_l / effective_population_l) -
    S4_l.v * delta 
  
  dI4_h <- 
    native * S4_h * (0.25 * beta_h * infected_total_h / effective_population_h) +
    travel * S4_h * (0.25 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I4_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_h, -1), 0) -
    I4_h * gamma -
    I4_h * delta
  dI4_h.v <-
    native * S4_h.v * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S4_h.v * (0.5 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I4_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_h.v, -1), 0) -
    I4_h.v * gamma -
    I4_h.v * delta
  dI4_l <- 
    travel * S4_l * (0.25 * beta_h * infected_total_h / effective_population_h) +
    native * S4_l * (0.25 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_l, -1), 0) -
    I4_l * gamma -
    I4_l * delta
  dI4_l.v <-
    native * S4_l.v * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S4_l.v * (0.5 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I4_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_l.v, -1), 0) -
    I4_l.v * gamma -
    I4_l.v * delta
  
  dR4_h <- 
    I4_h * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R4_h, -1)) -
    c(1/365 / head(age_window, -1) * head(R4_h, -1), 0) -
    R4_h * delta +
    R3_h * vac_h +
    S4_h * vac_h
  dR4_h.v <-
    I4_h.v * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_h.v, -1), 0) -
    R4_h.v * delta + 
    R3_h.v * vac_h +
    S4_h.v * vac_h
  dR4_l <- 
    I4_l * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R4_l, -1), 0) -
    R4_l * delta + 
    R3_l * vac_l +
    S4_l * vac_l
  dR4_l.v <-
    I4_l.v * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_l.v, -1), 0) -
    R4_l.v * delta + 
    R3_l.v * vac_l +
    S4_l.v * vac_l
  
  
  list(c(dS1_h, dI1_h, dR1_h,
         dS2_h, dI2_h, dR2_h,
         dS3_h, dI3_h, dR3_h,
         dS4_h, dI4_h, dR4_h,
         dR1_h.v,
         dS2_h.v, dI2_h.v, dR2_h.v,
         dS3_h.v, dI3_h.v, dR3_h.v,
         dS4_h.v, dI4_h.v, dR4_h.v,
         dS1_l, dI1_l, dR1_l,
         dS2_l, dI2_l, dR2_l,
         dS3_l, dI3_l, dR3_l,
         dS4_l, dI4_l, dR4_l,
         dR1_l.v,
         dS2_l.v, dI2_l.v, dR2_l.v,
         dS3_l.v, dI3_l.v, dR3_l.v,
         dS4_l.v, dI4_l.v, dR4_l.v))
  
}

############################
#RUN MODEL
############################
y_init <- c(susceptible_h(1), infected_h(1), recovered_h(1),
            susceptible_h(2), infected_h(2), recovered_h(2),
            susceptible_h(3), infected_h(3), recovered_h(3),
            susceptible_h(4), infected_h(4), recovered_h(4),
            rep(0, 280),
            susceptible_l(1), infected_l(1), recovered_l(1),
            susceptible_l(2), infected_l(2), recovered_l(2),
            susceptible_l(3), infected_l(3), recovered_l(3),
            susceptible_l(4), infected_l(4), recovered_l(4),
            rep(0, 280))
years = 50
times <- seq(from = 0, to = 365 * years, by = .1)
out <- ode(times = times, y = y_init, func = model, parms = parms)

############################
#PROCESSING
############################
#png(filename = paste('infections_', input, '.png', sep = ''))
#par(mfrow = c(2,2))
#infected_names <- c("Infected 1", "Infected 2", "Infected 3", "Infected 4")
#for(i in 0:3){
#  x <- 2 + 3 * i
#  y <- x + 12
#  legend_y <- max(rowSums(out[,(1+(x-1)*28+1):(1+(x)*28)]))
#  {plot(out[,1],rowSums(out[,(1+(x-1)*28+1):(1+(x)*28)]),type='l', 
#        col = "red", lwd = 2, ylab = '',
#        main = infected_names[i + 1])
#    lines(out[,1],rowSums(out[,(1+(y-1)*28+1):(1+(y)*28)]),type='l', lwd = 2)
#    legend(8000, (2/3) * legend_y, legend=c("High SES", "Low SES"),
#           col=c("red", "black"), lwd = c(2,2), cex=0.8)
#  }}
#dev.off()


############################
#OUTPUT
############################
out_last <- out[nrow(out),]
save(out_last, file = paste('output_', input, '.RData', sep = ''))


row_names <- c("I1", "I2", "I3", "I4", "I2.v", "I3.v", "I4.v")
infected.h <- matrix(NA, nrow = 7, ncol = years * 365 * 10)
for(i in 1:4){
  for(j in 1:(years * 365 * 10)){
    x <- i - 1 
    y <- sum(out[j,((30:57)+ 84 * x)]) - sum(out[(j-1),((30:57)+ 84 * x)])
    infected.h[i,j] <- ifelse(y>0, y, 0)
  }
}
for(i in 1:3){
  for(j in 1:(years * 365 * 10)){
    x <- i - 1
    y <- sum(out[j,((393:420)+ 84 * x)]) - sum(out[(j-1),((393:420)+ 84 * x)]) 
    infected_h[i + 4,j] <- ifelse(y>0, y, 0)
  }
}
row.names(infected.h) <- row_names

infected.l <- matrix(NA, nrow = 7, ncol = years * 365 * 10)
for(i in 1:4){
  for(j in 1:(years * 365 * 10)){
    x <- i - 1 
    y <- sum(out[j,((646:673)+ 84 * x)]) - sum(out[(j -1),((646:673)+ 84 * x)])
    infected.l[i,j] <- ifelse(y>0, y, 0)
  }
}
for(i in 1:3){
  for(j in 1:(years * 365 * 10)){
    x <- i - 1
    y <- sum(out[j,((1010:1037)+ 84 * x)]) - sum(out[(j-1),((1010:1037)+ 84 * x)]) 
    infected.l[i + 4,j] <- ifelse(y>0, y, 0)
  }
}
row.names(infected.l) <- row_names

save(infected.l, file = paste('infected.l_', input, '.RData', sep = ''))
save(infected.h, file = paste('infected.h_', input, '.RData', sep = ''))

#################
#prop infected over time
#infected <- matrix(NA, nrow = 4, ncol = years * 365 * 10)
#for(i in 1:4){
#  for(j in 1:(years*365*10)){
#    x <- i - 1 
#    infected[i,j] <- sum(out[j,((30:57)+ 84 * x)]) / sum(out[j,2:337])
#  }
#}
#png(filename = paste('infected.prop.h_', input, '.png', sep = ''))
#{plot(infected[1,], type = 'l', main = "Proportion of infections over time, high SES")
  lines(infected[2,], type = 'l', col = "red")
  lines(infected[3,], type = 'l', col = "blue")
  lines(infected[4,], type = 'l', col = "darkgreen")
  legend(((7/8) * years * 3650), 0.04, legend=c("I1", "I2", "I3", "I4"),
                   col=c("black", "red", "blue", "darkgreen"), lwd = c(2,2), cex=0.8)
#}
#dev.off()
#save(infected, file = paste('infected.h_', input, '.RData', sep = ''))

#infected_l <- matrix(NA, nrow = 4, ncol = years * 365 * 10)
#for(i in 1:4){
#  for(j in 1:(years*365*10)){
#    x <- i - 1 
#    infected_l[i,j] <- sum(out[j,((366:393)+ 84 * x)]) / sum(out[j,2:337])
#  }
#}
#png(filename = paste('infected.prop.l_', input, '.png', sep = ''))
#{plot(infected_l[1,], type = 'l', main = "Proportion of infections over time, low SES")
#  lines(infected_l[2,], type = 'l', col = "red")
#  lines(infected_l[3,], type = 'l', col = "blue")
#  lines(infected_l[4,], type = 'l', col = "darkgreen")
#  legend(((7/8) * years * 3650), 0.04, legend=c("I1", "I2", "I3", "I4"),
#         col=c("black", "red", "blue", "darkgreen"), lwd = c(2,2), cex=0.8)
#}
#dev.off()
#save(infected_l, file = paste('infected.l_', input, '.RData', sep = ''))


#plot infected
#png(filename = paste('infected_', input, '.png', sep = ''))
#par(mfrow = c(2,4))
#infected_names <- c("Infected 1, High", "Infected 2, High", 
#                    "Infected 3, High", "Infected 4, High",
#                    "Infected 1, Low", "Infected 2, Low", 
#                    "Infected 3, Low", "Infected 4, Low")
#for(i in 0:7){
#  x <- 2 + 3 * i
#  plot(out[,1],rowSums(out[,(1+(x-1)*28+1):(1+(x)*28)]),type='l', 
#       col = "black",
#       main = infected_names[i + 1])
#}
#dev.off()

#Plot susceptible
#png(filename = paste('susceptible_', input, '.png', sep = ''))
#par(mfrow = c(2,4))
#susceptible_names <- c("Susceptible 1, High", "Susceptible 2, High", 
#                       "Susceptible 3, High", "Susceptible 4, High",
#                       "Susceptible 1, Low", "Susceptible 2, Low", 
#                       "Susceptible 3, Low", "Susceptible 4, Low")
#for(i in 0:7){
#  x <- 1 + 3 * i
#  plot(out[,1],rowSums(out[,(1+(x-1)*28+1):(1+(x)*28)]),type='l', 
#         col = "black", main = susceptible_names[i + 1])
#}
#dev.off()

#plot recovered
#png(filename = paste('recovered_', input, '.png', sep = ''))
#par(mfrow = c(2,4))
#recovered_names <- c("Recovered 1, High", "Recovered 2, High", 
#                     "Recovered 3, High", "Recovered 4, High",
#                     "Recovered 1, Low", "Recovered 2, Low", 
#                     "Recovered 3, Low", "Recovered 4, Low")
#for(i in 0:7){
#  x <- 3 + 3 * i
#  plot(out[,1],rowSums(out[,(1+(x-1)*28+1):(1+(x)*28)]),type='l', 
#         col = "black", main = recovered_names [i + 1])
#}
#dev.off()

###############################
#plot proportions SIR over time 
###############################

png(filename = paste('prop.h.SIR_', input, '.png', sep = ''))
ts <- matrix(NA, nrow = 22, ncol = years * 365 * 10)
for(j in 1:22){
  for(i in 1:(years* 365 * 10)){
    x <- j - 1
    ts[j,i] <- sum(out[i,((x * 28) + 2):((x * 28) + 29)]) / sum(out[i,2:617])
  }
}

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
########
##set colors
########
colors <- c(rgb(0, 191/255, 255/255, alpha = 0.25), rgb(0, 191/255, 255/255, alpha = 0.5), rgb(0, 191/255, 255/255, alpha = 1),
            rgb(0, 0, 255/255, alpha = 0.25), rgb(0, 0, 255/255, alpha = 0.5), rgb(0, 0, 255/255, alpha = 1),
            rgb(0, 139/255, 0, alpha = 0.25), rgb(0, 139/255, 0, alpha = 0.5), rgb(0, 139/255, 0, alpha = 1),
            rgb(104/255, 34/255, 139/255, alpha = 0.25), rgb(104/255, 34/255, 139/255, alpha = 0.5), rgb(104/255, 34/255, 139/255, alpha = 1),
           gray.colors(10))
             #rgb(0, 191/255, 255/255, alpha = 1),
            #rgb(0, 0, 255/255, alpha = 0.25), rgb(0, 0, 255/255, alpha = 0.5), rgb(0, 0, 255/255, alpha = 1),
            #rgb(0, 139/255, 0, alpha = 0.25), rgb(0, 139/255, 0, alpha = 0.5), rgb(0, 139/255, 0, alpha = 1),
            #rgb(104/255, 34/255, 139/255, alpha = 0.25), rgb(104/255, 34/255, 139/255, alpha = 0.5), rgb(104/255, 34/255, 139/255, alpha = 1))
#######
barplot(ts, col = colors, 
        border=colors, space=0.04, main ="Proportion of Each Category, High SES", 
        xlab = "Timestep")
legend("topright",  inset=c(-0.1,0),
       c("S1", "I1", "R1", 
         "S2", "I2", "R2",
         "S3", "I3", "R3",
         "S4", "I4", "R4"), 
       fill=colors[1:12], horiz=FALSE, cex=0.8)
legend("topright",  inset=c(-0.25,0.083),
       c("R1.v", 
         "S2.v", "I2.v", "R2.v",
         "S3.v", "I3.v", "R3.v",
         "S4.v", "I4.v", "R4.v"), 
       fill=colors[13:22], horiz=FALSE, cex=0.8)
rgdev.off()

png(filename = paste('prop.l.SIR_', input, '.png', sep = ''))
ts_l <- matrix(NA, nrow = 22, ncol = years * 365 * 10)
for(j in 1:22){
  for(i in 1:(years* 365 * 10)){
    x <- j - 1
    ts_l[j,i] <- sum(out[i,((x * 28) + 2+ 616):((x * 28) + 29 + 616)]) / sum(out[i,618:1233])
  }
}
{par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(ts_l, col = topo.colors(12), 
        border=colors, space=0.04, main ="Proportion of Each Category, Low SES",
        xlab = "Timestep")
legend("topright",  inset=c(-0.1,0),
       c("S1", "I1", "R1", 
         "S2", "I2", "R2",
         "S3", "I3", "R3",
         "S4", "I4", "R4"), 
       fill=colors[1:12], horiz=FALSE, cex=0.8)
legend("topright",  inset=c(-0.25,0.083),
       c("R1.v", 
         "S2.v", "I2.v", "R2.v",
         "S3.v", "I3.v", "R3.v",
         "S4.v", "I4.v", "R4.v"), 
       fill=colors[13:22], horiz=FALSE, cex=0.8)}
dev.off()
