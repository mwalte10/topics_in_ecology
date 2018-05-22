############################
#CRC info
############################
args = commandArgs(TRUE)
input = as.numeric(args[1])
beta_h_more <- 0.3684211 
beta_l_more <- 0.95
beta_h.more <- beta_h_more[input]
beta_l.more <- beta_l_more[input]

beta_h_less <- 0.4710526 
beta_l_less <- 0.8815789
beta_h.less <- beta_h_less[input]
beta_l.less <- beta_l_less[input]


library(deSolve)

############################
#Initial conidtions and parameters
############################
percentage_vec <- c(rep(1.8,5), rep(1.8,5), rep(1.84,5), rep(1.86,5),
                    17.2, 15, 12.7, 8.8, 5.5, 2.8, 1.3, 0.20)
percentage_vec <- percentage_vec / 100
initial_conditions <- as.data.frame(matrix(NA, nrow = 28*4, ncol = 3))
initial_conditions[,2] <- rep(1:28,4)
for(i in 1:4){
  x <- i - 1
  initial_conditions[x*(28) + (1:28),1] <- rep(i,28)
}
initial_conditions[,3] <- rep(percentage_vec,4)

susceptible_init_h <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 112))
susceptible_init_h[,3] <- c(initial_conditions[1:28,3] * (6 * 10^6), initial_conditions[29:56,3] * 0, 
                            initial_conditions[57:84,3] * 0, initial_conditions[85:112,3] * 0)
susceptible_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 112))
susceptible_init_l[,3] <- c(initial_conditions[1:28,3] * (6 * 10^6), initial_conditions[29:56,3] * 0, 
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
travel = 1 - native
travel_more <- 2 * travel
travel_less <- 0.5 * travel

parms_more <- c(beta_h = beta_h.more,
                beta_l = beta_l.more,
                gamma = 1/4,
                sigma = 1/(365 * 1.2),
                mu = 19 / (1000 * 365),
                delta = c((0.013 / 365), rep(0.001 / 365, 4), rep(0, 10), rep(0.001, 5),
                          0.002 / 365, 0.0025 / 365,
                          0.004 / 365, 0.0085 / 365,
                          0.0175 / 365, 0.0425 / 365,
                          0.114 / 365, 0.151 / 365),
                age_window = c(rep(1, 21), rep(10, 7)),
                native_more = c(rep(1, 7), rep(0.86, 9), rep(0.842, 4), 0.814, 0.7676, 0.7784, rep(0.809, 5)),
                travel <- travel_more,
                vac_h = 0,
                vac_l = 0,
                sens = 0.8,
                spec = 0.95)
parms_less <- c(beta_h = beta_h.less,
                beta_l = beta_l.less,
                gamma = 1/4,
                sigma = 1/(365 * 1.2),
                mu = 19 / (1000 * 365),
                delta = c((0.013 / 365), rep(0.001 / 365, 4), rep(0, 10), rep(0.001, 5),
                          0.002 / 365, 0.0025 / 365,
                          0.004 / 365, 0.0085 / 365,
                          0.0175 / 365, 0.0425 / 365,
                          0.114 / 365, 0.151 / 365),
                age_window = c(rep(1, 21), rep(10, 7)),
                native_less = c(rep(1, 7), rep(0.86, 9), rep(0.842, 4), 0.814, 0.7676, 0.7784, rep(0.809, 5)),
                travel <- travel_less,
                vac_h = 0,
                vac_l = 0,
                sens = 0.8,
                spec = 0.95)

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
  t_vac.h <- ifelse((t>(365*(years_vac))), parms[118] / 365, 0)
  t_vac.l <- ifelse((t>(365*(years_vac))), parms[119] / 365, 0)
  vac_h <- c(rep(0,8), t_vac.h, rep(0,19))
  vac_l <- c(rep(0,8), t_vac.l, rep(0,19))
  sens <- parms[120]
  spec <- parms[121]
  
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
  R1_h.v <- y[(1:28) + 336]
  S2_h.v <- y[(29:56) + 336]
  I2_h.v <- y[(57:84) + 336]
  R2_h.v <- y[(85:112) + 336]
  S3_h.v <- y[(113:140) + 336]
  I3_h.v <- y[(141:168) + 336]
  R3_h.v <- y[(169:196) + 336]
  S4_h.v <- y[(197:224) + 336]
  I4_h.v <- y[(225:252) + 336]
  R4_h.v <- y[(253:280) + 336]
  
  S1_l <- y[(1:28) + 616]
  I1_l <- y[(29:56) + 616]
  R1_l <- y[(57:84) + 616]
  S2_l <- y[(85:112) + 616]
  I2_l <- y[(113:140) + 616]
  R2_l <- y[(141:168) + 616]
  S3_l <- y[(169:196) + 616]
  I3_l <- y[(197:224) + 616]
  R3_l <- y[(225:252) + 616]
  S4_l <- y[(253:280) + 616]
  I4_l <- y[(281:308) + 616]
  R4_l <- y[(309:336) + 616]
  
  R1_l.v <- y[(1:28) + 952]
  S2_l.v <- y[(29:56) + 952]
  I2_l.v <- y[(57:84) + 952]
  R2_l.v <- y[(85:112) + 952]
  S3_l.v <- y[(113:140) + 952]
  I3_l.v <- y[(141:168) + 952]
  R3_l.v <- y[(169:196) + 952]
  S4_l.v <- y[(197:224) + 952]
  I4_l.v <- y[(225:252) + 952]
  R4_l.v <- y[(253:280) + 952]
  
  ###Tracking infected individuals
  I_total <- y[1233]
  I_l <- y[1234]
  ###Tracking cases
  cases <- y[1235]
  cases.l <- y[1236]
  
  
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
    S1_h * vac_h * (1 - spec)
  dS1_l <-  
    mu * c(pop_l, rep(0,27)) +
    c(0, 1/365 / head(age_window, -1) * head(S1_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S1_l, -1), 0) -
    travel * S1_l * (beta_h * infected_total_h / effective_population_h) -
    native * S1_l * (beta_l * infected_total_l / effective_population_l) -
    S1_l * delta - 
    S1_l * vac_l * (1 - spec)
  
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
    R1_h * vac_h * sens
  dR1_h.v <- 
    S1_h * vac_h * (1 - spec) +
    c(0, 1/365 / head(age_window, -1) * head(R1_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R1_h.v, -1), 0) -
    R1_h.v * delta -
    R1_h.v * sigma
  dR1_l <-
    I1_l * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R1_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R1_l, -1), 0) -
    R1_l * sigma -
    R1_l * delta - 
    R1_l * vac_l * sens
  dR1_l.v <- 
    S1_l * vac_l * (1 - spec) +
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
    S2_h * vac_h * sens
  dS2_h.v <-
    R1_h.v * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_h.v, -1), 0) -
    native * S2_h.v * (beta_h * infected_total_h / effective_population_h) -
    travel * S2_h.v * (beta_l * infected_total_l / effective_population_l) -
    S2_h.v * delta
  dS2_l <- 
    R1_l * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_l, -1), 0) -
    travel * S2_l * (0.75 * beta_h * infected_total_h / effective_population_h) -
    native * S2_l * (0.75 * beta_l * infected_total_l / effective_population_l) -
    S2_l * delta - 
    S2_l * vac_l * sens
  dS2_l.v <-
    R1_l.v * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_l.v, -1), 0) -
    native * S2_l.v * (beta_h * infected_total_h / effective_population_h) -
    travel * S2_l.v * (beta_l * infected_total_l / effective_population_l) -
    S2_l.v * delta
  
  dI2_h <- 
    native * S2_h * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S2_h * (0.75 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I2_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_h, -1), 0) -
    I2_h * gamma -
    I2_h * delta
  dI2_h.v <-
    native * S2_h.v * (beta_h * infected_total_h / effective_population_h) +
    travel * S2_h.v * (beta_l * infected_total_l / effective_population_l) +
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
    native * S2_l.v * (beta_h * infected_total_h / effective_population_h) +
    travel * S2_l.v * (beta_l * infected_total_l / effective_population_l) +
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
    R2_h * vac_h * sens
  dR2_h.v <-
    I2_h.v * gamma +  
    c(0, 1/365 / head(age_window, -1) * head(R2_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R2_h.v, -1), 0) -
    R2_h.v * sigma -
    R2_h.v * delta +
    R1_h * vac_h * sens +
    S2_h * vac_h * sens
  dR2_l <-  
    I2_l * gamma +  
    c(0, 1/365 / head(age_window, -1) * head(R2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R2_l, -1), 0) -
    R2_l * sigma -
    R2_l * delta -
    R2_l * vac_l * sens
  dR2_l.v <-
    I2_l.v * gamma +  
    c(0, 1/365 / head(age_window, -1) * head(R2_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R2_l.v, -1), 0) -
    R2_l.v * sigma -
    R2_l.v * delta +
    R1_l * vac_l * sens +
    S2_l * vac_l * sens
  
  #third infection
  dS3_h <- 
    R2_h * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_h, -1), 0) -
    native * S3_h * (0.5 * beta_h * infected_total_h / effective_population_h) -
    travel * S3_h * (0.5 * beta_l * infected_total_l / effective_population_l) -
    S3_h * delta -
    S3_h * vac_h * sens
  dS3_h.v <-
    R2_h.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_h.v, -1), 0) -
    native * S3_h.v * (0.75 * beta_h * infected_total_h / effective_population_h) -
    travel * S3_h.v * (0.75 * beta_l * infected_total_l / effective_population_l) -
    S3_h.v * delta 
  dS3_l <- 
    R2_l * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_l, -1), 0) -
    travel * S3_l * (0.5 * beta_h * infected_total_h / effective_population_h) -
    native * S3_l * (0.5 * beta_l * infected_total_l / effective_population_l) -
    S3_l * delta -
    S3_l * vac_l * sens
  dS3_l.v <-
    R2_l.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_l.v, -1), 0) -
    native * S3_l.v * (0.75 * beta_h * infected_total_h / effective_population_h) -
    travel * S3_l.v * (0.75 * beta_l * infected_total_l / effective_population_l) -
    S3_l.v * delta
  
  dI3_h <- 
    native * S3_h * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_h * (0.5 * beta_l * infected_total_l / effective_population_l) +
    c(0, 1/365 / head(age_window, -1) * head(I3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_h, -1), 0) -
    I3_h * gamma -
    I3_h * delta
  dI3_h.v <-
    native * S3_h.v * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_h.v * (0.75 * beta_l * infected_total_l / effective_population_l) +
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
    native * S3_l.v * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_l.v * (0.75 * beta_l * infected_total_l / effective_population_l) +
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
    R3_h * vac_h * sens
  dR3_h.v <-
    I3_h.v * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_h.v, -1), 0) -
    R3_h.v * sigma + 
    R2_h * vac_h * sens +
    S3_h * vac_h * sens - 
    R3_h.v * delta
  dR3_l <- 
    I3_l * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_l, -1), 0) -
    R3_l * sigma -
    R3_l * delta - 
    R3_l * vac_l * sens
  dR3_l.v <-
    I3_l.v * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_l.v, -1), 0) -
    R3_l.v * sigma + 
    R2_l * vac_l * sens +
    S3_l * vac_l * sens -
    R3_l.v * delta
  
  #fourth infection
  dS4_h <- 
    R3_h * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S4_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_h, -1), 0) -
    native * S4_h * (0.25 * beta_h * infected_total_h / effective_population_h) -
    travel * S4_h * (0.25 * beta_l * infected_total_l / effective_population_l) -
    S4_h * delta - 
    S4_h * vac_h * sens
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
    S4_l * vac_l * sens
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
    R4_h * delta 
  dR4_h.v <-
    I4_h.v * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R4_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R4_h.v, -1), 0) -
    R4_h.v * delta + 
    R3_h * vac_h * sens +
    S4_h * vac_h * sens
  dR4_l <- 
    I4_l * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R4_l, -1), 0) -
    R4_l * delta 
  dR4_l.v <-
    I4_l.v * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R4_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R4_l.v, -1), 0) -
    R4_l.v * delta + 
    R3_l * vac_l * sens+
    S4_l * vac_l * sens
  
  I_total <-
    native * S1_h * (beta_h * infected_total_h / effective_population_h) +
    travel * S1_h * (beta_l * infected_total_l / effective_population_l) +
    travel * S1_l * (beta_h * infected_total_h / effective_population_h) +
    native * S1_l * (beta_l * infected_total_l / effective_population_l) +
    native * S2_h * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S2_h * (0.75 * beta_l * infected_total_l / effective_population_l) +
    native * S2_h.v * (beta_h * infected_total_h / effective_population_h) +
    travel * S2_h.v * (beta_l * infected_total_l / effective_population_l) +
    travel * S2_l * (0.75 * beta_h * infected_total_h / effective_population_h) +
    native * S2_l * (0.75 * beta_l * infected_total_l / effective_population_l) +
    native * S2_l.v * (beta_h * infected_total_h / effective_population_h) +
    travel * S2_l.v * (beta_l * infected_total_l / effective_population_l) +
    native * S3_h * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_h * (0.5 * beta_l * infected_total_l / effective_population_l) +
    native * S3_h.v * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_h.v * (0.75 * beta_l * infected_total_l / effective_population_l) +
    travel * S3_l * (0.5 * beta_h * infected_total_h / effective_population_h) +
    native * S3_l * (0.5 * beta_l * infected_total_l / effective_population_l) +
    native * S3_l.v * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_l.v * (0.75 * beta_l * infected_total_l / effective_population_l) +
    native * S4_h * (0.25 * beta_h * infected_total_h / effective_population_h) +
    travel * S4_h * (0.25 * beta_l * infected_total_l / effective_population_l) +
    native * S4_h.v * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S4_h.v * (0.5 * beta_l * infected_total_l / effective_population_l) +
    travel * S4_l * (0.25 * beta_h * infected_total_h / effective_population_h) +
    native * S4_l * (0.25 * beta_l * infected_total_l / effective_population_l) +
    native * S4_l.v * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S4_l.v * (0.5 * beta_l * infected_total_l / effective_population_l)
  I_total <- sum(I_total)
  
  I_l <-
    travel * S1_l * (beta_h * infected_total_h / effective_population_h) +
    native * S1_l * (beta_l * infected_total_l / effective_population_l) +
    travel * S2_l * (0.75 * beta_h * infected_total_h / effective_population_h) +
    native * S2_l * (0.75 * beta_l * infected_total_l / effective_population_l) +
    native * S2_l.v * (beta_h * infected_total_h / effective_population_h) +
    travel * S2_l.v * (beta_l * infected_total_l / effective_population_l) +
    travel * S3_l * (0.5 * beta_h * infected_total_h / effective_population_h) +
    native * S3_l * (0.5 * beta_l * infected_total_l / effective_population_l) +
    native * S3_l.v * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_l.v * (0.75 * beta_l * infected_total_l / effective_population_l) +
    travel * S4_l * (0.25 * beta_h * infected_total_h / effective_population_h) +
    native * S4_l * (0.25 * beta_l * infected_total_l / effective_population_l) +
    native * S4_l.v * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S4_l.v * (0.5 * beta_l * infected_total_l / effective_population_l)
  I_l <- sum(I_l)
  
  primary <-
    native * S1_h * (beta_h * infected_total_h / effective_population_h) +
    travel * S1_h * (beta_l * infected_total_l / effective_population_l) +
    travel * S1_l * (beta_h * infected_total_h / effective_population_h) +
    native * S1_l * (beta_l * infected_total_l / effective_population_l) +
    native * S2_h.v * (beta_h * infected_total_h / effective_population_h) +
    travel * S2_h.v * (beta_l * infected_total_l / effective_population_l) +
    native * S2_l.v * (beta_h * infected_total_h / effective_population_h) +
    travel * S2_l.v * (beta_l * infected_total_l / effective_population_l)
  primary <- sum(primary) * 0.18
  
  secondary <-
    native * S2_h * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S2_h * (0.75 * beta_l * infected_total_l / effective_population_l) +
    travel * S2_l * (0.75 * beta_h * infected_total_h / effective_population_h) +
    native * S2_l * (0.75 * beta_l * infected_total_l / effective_population_l) +
    native * S3_h.v * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_h.v * (0.75 * beta_l * infected_total_l / effective_population_l) +
    native * S3_l.v * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_l.v * (0.75 * beta_l * infected_total_l / effective_population_l)
  secondary <- sum(secondary) * 0.41
  
  tertiary <-
    native * S3_h * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_h * (0.5 * beta_l * infected_total_l / effective_population_l) +
    travel * S3_l * (0.5 * beta_h * infected_total_h / effective_population_h) +
    native * S3_l * (0.5 * beta_l * infected_total_l / effective_population_l) +
    native * S4_h.v * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S4_h.v * (0.5 * beta_l * infected_total_l / effective_population_l) +
    native * S4_l.v * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S4_l.v * (0.5 * beta_l * infected_total_l / effective_population_l)
  tertiary <- sum(tertiary) * 0.063425
  
  quaternary <-
    native * S4_h * (0.25 * beta_h * infected_total_h / effective_population_h) +
    travel * S4_h * (0.25 * beta_l * infected_total_l / effective_population_l) +
    travel * S4_l * (0.25 * beta_h * infected_total_h / effective_population_h) +
    native * S4_l * (0.25 * beta_l * infected_total_l / effective_population_l)
  quaternary <- sum(quaternary) * 0.063425
  
  cases <- primary + secondary + tertiary + quaternary
  
  primary.l <-
    travel * S1_l * (beta_h * infected_total_h / effective_population_h) +
    native * S1_l * (beta_l * infected_total_l / effective_population_l) +
    native * S2_l.v * (beta_h * infected_total_h / effective_population_h) +
    travel * S2_l.v * (beta_l * infected_total_l / effective_population_l)
  primary.l <- sum(primary.l) * 0.18
  
  secondary.l <-
    travel * S2_l * (0.75 * beta_h * infected_total_h / effective_population_h) +
    native * S2_l * (0.75 * beta_l * infected_total_l / effective_population_l) +
    native * S3_l.v * (0.75 * beta_h * infected_total_h / effective_population_h) +
    travel * S3_l.v * (0.75 * beta_l * infected_total_l / effective_population_l)
  secondary.l <- sum(secondary.l) * 0.41
  
  tertiary.l <-
    travel * S3_l * (0.5 * beta_h * infected_total_h / effective_population_h) +
    native * S3_l * (0.5 * beta_l * infected_total_l / effective_population_l) +
    native * S4_l.v * (0.5 * beta_h * infected_total_h / effective_population_h) +
    travel * S4_l.v * (0.5 * beta_l * infected_total_l / effective_population_l)
  tertiary.l <- sum(tertiary.l) * 0.063425
  
  quaternary.l <-
    travel * S4_l * (0.25 * beta_h * infected_total_h / effective_population_h) +
    native * S4_l * (0.25 * beta_l * infected_total_l / effective_population_l)
  quaternary.l <- sum(quaternary.l) * 0.063425
  
  cases.l <- primary.l + secondary.l  + tertiary.l + quaternary.l
  
  
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
         dS4_l.v, dI4_l.v, dR4_l.v,
         I_total, I_l, cases, cases.l))
  
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
            rep(0, 280), 0, 0, 0, 0)
years = 80
years_vac = 30
times <- seq(from = 0, to = 365 * years, by = .1)


###############
#RECORD NULL RESPONSES
###############

out.more.null <- ode(times = times, y = y_init, func = model, parms = parms_more)
out.less.null <- ode(times = times, y = y_init, func = model, parms = parms_less)

{track_infected.more <- out.more.null[,(ncol(out.more.null) - 3)]
  track_infected.more <- track_infected.more[length(track_infected.more)]
  track_l.more <- out.more.null[,(ncol(out.more.null) - 2)]
  track_l.more <- track_l.more[length(track_l.more)]
  track_h.more <- track_infected.more - track_l.more
  cases.more <- out.more.null[,(ncol(out.more.null) - 1)]
  cases.more <- cases.more[length(cases.more)]
  cases.l.more <- out.more.null[,(ncol(out.more.null))]
  cases.l.more <- cases.l.more[length(cases.l.more)]
  cases.h.more <- cases.more - cases.l.more}

output_more <- c(track_h.more, track_infected.more, track_l.more, cases.h.more, cases.more, cases.l.more)
save(output_more, file = "output_more.RData")

{track_infected.less <- out.less.null[,(ncol(out.less.null) - 3)]
  track_infected.less <- track_infected.less[length(track_infected.less)]
  track_l.less <- out.less.null[,(ncol(out.less.null) - 2)]
  track_l.less <- track_l.less[length(track_l.less)]
  track_h.less <- track_infected.less - track_l.less
  cases.less <- out.less.null[,(ncol(out.less.null) - 1)]
  cases.less <- cases.less[length(cases.less)]
  cases.l.less <- out.less.null[,(ncol(out.less.null))]
  cases.l.less <- cases.l.less[length(cases.l.less)]
  cases.h.less <- cases.less - cases.l.less}

output_less <- c(track_h.less, track_infected.less, track_l.less, cases.h.less, cases.less, cases.l.less)
save(output_less, file = "output_less.RData")




