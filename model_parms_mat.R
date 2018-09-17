#track_i.R
############################
#CRC info
############################
args = commandArgs(TRUE)
input = as.numeric(args[1])


##########################
#IF DOING VAC COVERAGE
##########################
load('parameters.RData')
beta_h <- new.parms.mat[input, 1]
beta_l <- new.parms.mat[input, 2]
native <- rep(new.parms.mat[input, 3], 28)
vac_h <- new.parms.mat[input, 4]
vac_l <- new.parms.mat[input, 5]


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
susceptible_init_h[,3] <- c(initial_conditions[1:28,3] * (0.5 * 12 * 10^6), initial_conditions[29:56,3] * 0, 
                            initial_conditions[57:84,3] * 0, initial_conditions[85:112,3] * 0)
susceptible_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 112))
susceptible_init_l[,3] <- c(initial_conditions[1:28,3] * (0.5 * 12 * 10^6), initial_conditions[29:56,3] * 0, 
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


parms <- c(beta_h = beta_h,
           beta_l = beta_l,
           gamma = 1/4,
           sigma = 1/(365 * 1.2),
           mu = 19 / (1000 * 365),
           delta = c((0.013 / 365), rep(0.001 / 365, 4), rep(0, 10), rep(0.001, 5),
                     0.002 / 365, 0.0025 / 365,
                     0.004 / 365, 0.0085 / 365,
                     0.0175 / 365, 0.0425 / 365,
                     0.114 / 365, 0.151 / 365),
           age_window = c(rep(1, 21), rep(10, 7)),
           native = native,
           travel = 1 - native,
           vac_h = vac_h,
           vac_l = vac_l,
           sens = 0.85,
           spec = 0.95)
parms_null <- c(beta_h = beta_h,
                beta_l = beta_l,
                gamma = 1/4,
                sigma = 1/(365 * 1.2),
                mu = 19 / (1000 * 365),
                delta = c((0.013 / 365), rep(0.001 / 365, 4), rep(0, 10), rep(0.001, 5),
                          0.002 / 365, 0.0025 / 365,
                          0.004 / 365, 0.0085 / 365,
                          0.0175 / 365, 0.0425 / 365,
                          0.114 / 365, 0.151 / 365),
                age_window = c(rep(1, 21), rep(10, 7)),
                native = native,
                travel = 1 - native,
                vac_h = 0,
                vac_l = 0,
                sens = 0.85,
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
  #paste("I_", seq(1,28), sep = '') <- y[1234:1261]
  
  infecteds_h <- cbind((I1_h), (I2_h), (I3_h), (I4_h), (I2_h.v), (I3_h.v), (I4_h.v))
  infected_total_h <- rowSums(infecteds_h)
  population_h <- cbind(S1_h, I1_h, R1_h,
                        S2_h, I2_h, R2_h,
                        S3_h, I3_h, R3_h,
                        S4_h, I4_h, R4_h,
                        R1_h.v,
                        S2_h.v, I2_h.v, R2_h.v,
                        S3_h.v, I3_h.v, R3_h.v,
                        S4_h.v, I4_h.v, R4_h.v)
  pop_h <- rowSums(population_h)
  birth_pop_h <- sum(pop_h)
  infecteds_l <- cbind((I1_l), (I2_l), (I3_l), (I4_l), (I2_l.v), (I3_l.v), (I4_l.v))
  infected_total_l <- rowSums(infecteds_l)
  population_l <- cbind(S1_l, I1_l, R1_l,
                 S2_l, I2_l, R2_l,
                 S3_l, I3_l, R3_l,
                 S4_l, I4_l, R4_l,
                 R1_l.v,
                 S2_l.v, I2_l.v, R2_l.v,
                 S3_l.v, I3_l.v, R3_l.v,
                 S4_l.v, I4_l.v, R4_l.v)
  pop_l <- rowSums(population_l)
  birth_pop_l <- sum(pop_l)
  # effective_population_h <- population_h + rep(travel,22) * population_l
  # effective_population_h <- sum(effective_population_h)
  # effective_population_l <- population_l + rep(travel,22) * population_h
  # effective_population_l <- sum(effective_population_l)
  
  #first infection
  dS1_h <-  
    mu * c(birth_pop_h, rep(0,27)) +
    c(0, 1/365 / head(age_window, -1) * head(S1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S1_h, -1), 0) -
    native * S1_h * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    travel * S1_h * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    S1_h * delta -
    S1_h * vac_h * (1 - spec)
  dS1_l <-  
    mu * c(birth_pop_l, rep(0,27)) +
    c(0, 1/365 / head(age_window, -1) * head(S1_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S1_l, -1), 0)  -
    native * S1_l * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    travel * S1_l * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    S1_l * delta - 
    S1_l * vac_l * (1 - spec)
  
  dI1_h <- 
    native * S1_h * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S1_h * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I1_h, -1), 0) -
    I1_h * gamma -
    I1_h * delta
  dI1_l <- 
    native * S1_l * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S1_l * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
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
    native * S2_h * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    travel * S2_h * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    S2_h * delta - 
    S2_h * vac_h * sens
  dS2_h.v <-
    R1_h.v * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_h.v, -1), 0) -
    native * S2_h.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    travel * S2_h.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    S2_h.v * delta
  dS2_l <- 
    R1_l * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_l, -1), 0) -
    native * S2_l * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    travel * S2_l * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    S2_l * delta - 
    S2_l * vac_l * sens
  dS2_l.v <-
    R1_l.v * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_l.v, -1), 0) -
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    S2_l.v * delta
  
  dI2_h <- 
    native * S2_h * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I2_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_h, -1), 0) -
    I2_h * gamma -
    I2_h * delta
  dI2_h.v <-
    native * S2_h.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I2_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_h.v, -1), 0) -
    I2_h.v * gamma -
    I2_h.v * delta
  dI2_l <- 
    native * S2_l * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_l, -1), 0) -
    I2_l * gamma -
    I2_l * delta
  dI2_l.v <-
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
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
    native * S3_h * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    travel * S3_h * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    S3_h * delta -
    S3_h * vac_h * sens
  dS3_h.v <-
    R2_h.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_h.v, -1), 0) -
    native * S3_h.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    travel * S3_h.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    S3_h.v * delta 
  dS3_l <- 
    R2_l * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_l, -1), 0) -
    native * S3_l * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    travel * S3_l * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    S3_l * delta -
    S3_l * vac_l * sens
  dS3_l.v <-
    R2_l.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_l.v, -1), 0) -
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    S3_l.v * delta
  
  dI3_h <- 
    native * S3_h * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_h, -1), 0) -
    I3_h * gamma -
    I3_h * delta
  dI3_h.v <-
    native * S3_h.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_h.v, -1), 0) -
    I3_h.v * gamma -
    I3_h.v * delta
  dI3_l <- 
    native * S3_l * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_l, -1), 0) -
    I3_l * gamma -
    I3_l * delta
  dI3_l.v <-
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
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
    native * S4_h * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    travel * S4_h * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    S4_h * delta - 
    S4_h * vac_h * sens
  dS4_h.v <-
    R3_h.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S4_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_h.v, -1), 0) -
    native * S4_h.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    travel * S4_h.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    S4_h.v * delta 
  dS4_l <- 
    R3_l * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_l, -1), 0) -
    native * S4_l * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    travel * S4_l * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    S4_l * delta - 
    S4_l * vac_l * sens
  dS4_l.v <-
    R3_l.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S4_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_l.v, -1), 0) -
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    S4_l.v * delta 
  
  dI4_h <- 
    native * S4_h * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I4_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_h, -1), 0) -
    I4_h * gamma -
    I4_h * delta
  dI4_h.v <-
    native * S4_h.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I4_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_h.v, -1), 0) -
    I4_h.v * gamma -
    I4_h.v * delta
  dI4_l <- 
    native * S4_l * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_l, -1), 0) -
    I4_l * gamma -
    I4_l * delta
  dI4_l.v <-
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
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
    R3_l * vac_l * sens +
    S4_l * vac_l * sens
  
  I_tot <-
    native * S2_h * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_h.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_l * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
    native * S1_h * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S1_h * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S1_l * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S1_l * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S3_h * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_h.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_l * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_h * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_h.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_l * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  I_total <- sum(I_tot)
  
  I_primary_tot <- 
    native * S1_h * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S1_h * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S1_l * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S1_l * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  I_primary_tot <- sum(I_primary_tot)
  
  I_secondary_tot <-
    native * S2_h * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_h.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_l * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  I_secondary_tot <- sum(I_secondary_tot)
  
  I_secondary_vac <- 
    native * S2_h.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  I_secondary_vac <- sum(I_secondary_vac)
    
  
  I_post_sec_tot <- 
    native * S3_h * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_h.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_l * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_h * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_h.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_l * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  I_post_sec_tot <- sum(I_post_sec_tot)
  
  I_post_sec_vac <- 
    native * S3_h.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_h.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  I_post_sec_vac <- sum(I_post_sec_vac)
    
    
 
   I_l <-
    native * S1_l * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S1_l * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S2_l * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S3_l * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_l * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  I_l <- sum(I_l)
  
  I_l_vac <- 
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  I_l_vac <- sum(I_l_vac)
    
  
  I_l_primary_tot <- 
    native * S1_l * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S1_l * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  I_l_primary_tot <- sum(I_l_primary_tot)
  
  I_l_sec_tot <-
    native * S2_l * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  I_l_sec_tot <- sum(I_l_sec_tot)
  
  I_l_sec_vac <- 
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  I_l_sec_vac <- sum(I_l_sec_vac)
  
  I_l_post_sec_tot <- 
    native * S3_l * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_l * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  I_l_post_sec_tot <- sum(I_l_post_sec_tot)
  
  I_l_post_sec_vac <- 
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  I_l_post_sec_vac <- sum(I_l_post_sec_vac)
    
  
  primary <-
    native * S1_h * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S1_h * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S1_l * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S1_l * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  primary <- sum(primary) * 0.18
  
  ##primary cases incidence
{  primary_tot.cases <- 
    native * S1_h * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S1_h * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S1_l * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S1_l * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  primary_tot.cases <- sum(primary_tot.cases) * 0.18
  
  primary_tot.ncases <- 
    native * S1_h * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S1_h * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S1_l * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S1_l * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  primary_tot.ncases <- sum(primary_tot.ncases) * 0.82
  
  primary_cases.l <-
    travel * S1_h * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S1_l * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) 
  primary_cases.l <- sum(primary_cases.l) * 0.18
  
  primary_ncases.l <- 
    travel * S1_h * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S1_l * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) 
  primary_ncases.l <- sum(primary_ncases.l) * 0.82
  
  primary_cases.h <- 
    native * S1_h * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S1_l * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  primary_cases.h <- sum(primary_cases.h) * 0.18
  
  primary_ncases.h <- 
    native * S1_h * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S1_l * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  primary_ncases.h <- sum(primary_ncases.h) * 0.82}
  
  
  secondary <-
    native * S2_h * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_l * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S2_h.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  secondary <- sum(secondary) * 0.41
  
  #vaccinated case incidence
{  secondary_vac_tot.cases <- 
    native * S2_h.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  secondary_vac_tot.cases <- sum(secondary_vac_tot.cases) * 0.41
  
  secondary_vac_tot.ncases <- 
    native * S2_h.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  secondary_vac_tot.ncases <- sum(secondary_vac_tot.ncases) * 0.59
  
  secondary_vac_cases.l <- 
    travel * S2_h.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) 
  secondary_vac_cases.l <- sum(secondary_vac_cases.l) * 0.41
  
  secondary_vac_ncases.l <- 
    travel * S2_h.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) 
  secondary_vac_ncases.l <- sum(secondary_vac_ncases.l) * 0.59
  
  secondary_vac_cases.h <- 
    native * S2_h.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  secondary_vac_cases.h <- sum(secondary_vac_cases.h) * 0.41
  
  secondary_vac_ncases.h <- 
    native * S2_h.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  secondary_vac_ncases.h <- sum(secondary_vac_ncases.h) * 0.59}
  
  #nonvaccinated case incidence
{  secondary_tot.cases <-
    native * S2_h * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_l * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  secondary_tot.cases <- sum(secondary_tot.cases) * 0.41
  
  secondary_tot.ncases <-
    native * S2_h * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_l * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  secondary_tot.ncases <- sum(secondary_tot.ncases) * 0.59
  
  secondary_cases.l <-
    travel * S2_h * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_l * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) 
  secondary_cases.l <- sum(secondary_cases.l) * 0.41
  
  secondary_ncases.l <-
    travel * S2_h * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_l * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) 
  secondary_ncases.l <- sum(secondary_ncases.l) * 0.59
  
  secondary_cases.h <-
    native * S2_h * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_l * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  secondary_cases.h <- sum(secondary_cases.h) * 0.41
  
  secondary_ncases.h <-
    native * S2_h * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_l * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  secondary_ncases.h <- sum(secondary_ncases.h) * 0.59}
  
  
  tertiary <-
    native * S3_h * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_l * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S3_h.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  tertiary <- sum(tertiary) * 0.063425
  
  
  ##post sec vaccinated case incidence
 { postsec_vac_tot.cases <- 
    native * S3_h.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_h.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  postsec_vac_tot.cases <- sum(postsec_vac_tot.cases) * 0.063425
  
  postsec_vac_tot.ncases <- 
    native * S3_h.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_h.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  postsec_vac_tot.ncases <- sum(postsec_vac_tot.ncases) * 0.936575
  
  postsec_vac_cases.l <- 
    travel * S3_h.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_h.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) 
  postsec_vac_cases.l <- sum(postsec_vac_cases.l) * 0.063425
  
  postsec_vac_ncases.l <- 
    travel * S3_h.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_h.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) 
  postsec_vac_ncases.l <- sum(postsec_vac_ncases.l) * 0.936575
  
  postsec_vac_cases.h <- 
    native * S3_h.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_h.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  postsec_vac_cases.h <- sum(postsec_vac_cases.h) * 0.063425
  
  postsec_vac_ncases.h <- 
    native * S3_h.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_h.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  postsec_vac_ncases.h <- sum(postsec_vac_ncases.h) * 0.936575}
  
  ##post sec nonvaccinated case incidence
{  postsec_tot_cases <- 
    native * S3_h * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_l * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_h * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_l * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  postsec_tot_cases <- sum(postsec_tot_cases) * 0.063425
  
  postsec_tot_ncases <- 
    native * S3_h * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_l * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_h * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_l * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  postsec_tot_ncases <- sum(postsec_tot_ncases) * 0.936575
  
  postsec_cases.l <- 
    travel * S3_h * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_l * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_h * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_l * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) 
  postsec_cases.l <- sum(postsec_cases.l) * 0.063425
  
  postsec_ncases.l <- 
    travel * S3_h * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_l * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_h * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_l * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) 
  postsec_ncases.l <- sum(postsec_ncases.l) * 0.936575
  
  postsec_cases.h <- 
    native * S3_h * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_l * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_h * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_l * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  postsec_cases.h <- sum(postsec_cases.h) * 0.063425
  
  postsec_ncases.h <- 
    native * S3_h * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_l * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_h * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_l * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  postsec_ncases.h <- sum(postsec_ncases.h) * 0.936575
  
}
{  
  quaternary <-
    native * S4_h * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_l * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_h.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  quaternary <- sum(quaternary) * 0.063425
  
  cases <- primary + secondary + tertiary + quaternary
  post_sec <- tertiary + quaternary

  primary.l <-
    native * S1_l * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S1_l * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)  
  primary.l <- sum(primary.l) * 0.18
  
  secondary.l <-
    native * S2_l * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  secondary.l <- sum(secondary.l) * 0.41
  
  tertiary.l <-
    native * S3_l * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
    
  tertiary.l <- sum(tertiary.l) * 0.063425
  
  quaternary.l <-
    native * S4_l * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) + 
    native * S4_l.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  quaternary.l <- sum(quaternary.l) * 0.063425
  
  cases.l <- primary.l + secondary.l  + tertiary.l + quaternary.l
  post_sec.l <- tertiary.l + quaternary.l}

    
  
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
         I_secondary_vac, I_post_sec_vac,
         I_l_sec_vac, I_l_post_sec_vac,
         
         #1238 & 1239
         primary_tot.cases, primary_tot.ncases, 
         #1240, 1241
         primary_cases.l, primary_ncases.l,
         #1242, 1243
         primary_cases.h, primary_ncases.h,
         
         #1244, 1245
         secondary_vac_tot.cases, secondary_vac_tot.ncases,
         #1246, 1247
         secondary_vac_cases.l, secondary_vac_ncases.l,
         #1248, 1249
         secondary_vac_cases.h, secondary_vac_ncases.h,
         #1250, 1251
         secondary_tot.cases, secondary_tot.ncases,
         #1252, 1253
         secondary_cases.l, secondary_ncases.l,
         #1254, 1255
         secondary_cases.h, secondary_ncases.h,
         
         #1256, 1257
         postsec_vac_tot.cases, postsec_vac_tot.ncases,
         #1258, 1259
         postsec_vac_cases.l, postsec_vac_ncases.l,
         #1260, 1261
         postsec_vac_cases.h, postsec_vac_ncases.h,
         #1262, 1263
         postsec_tot_cases, postsec_tot_ncases,
         #1264, 1265
         postsec_cases.l, postsec_ncases.l,
         #1266, 
         postsec_cases.h, postsec_ncases.h,
         
         I_total, I_l, cases, cases.l,
         I_primary_tot, I_secondary_tot, I_post_sec_tot,
         I_l_primary_tot, I_l_sec_tot, I_l_post_sec_tot,
         primary, secondary, post_sec,
         primary.l, secondary.l, post_sec.l
         ))

  
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
            rep(0, 280), rep(0, 30),
            rep(0,20))
years = 50
years_vac = 30
times <- seq(from = 0, to = 365 * years, by = .1)
out <- ode(times = times, y = y_init, func = model, parms = parms)
out_null <- ode(times = times, y = y_init, func = model, parms = parms_null)

out_last <- out[nrow(out),(2:(ncol(out) - 16))]
out_last.null <- out_null[nrow(out_null),(2:(ncol(out_null) - 16))]
save(out_last, file = paste('out_last_', input, '.RData', sep = ''))
save(out_last.null, file = paste('out_last.null_', input, '.RData', sep = ''))


##incidence calcs 
{I_prim_tot <- diff(out[,1242])
I_sec_tot <- diff(out[,1243])
I_psec_tot <- diff(out[,1244])
I_prim_l <- diff(out[,1245])
I_sec_l <- diff(out[,1246])
I_psec_l <- diff(out[,1247])
I_prim_h <- I_prim_tot - I_prim_l
I_sec_h <- I_sec_tot - I_sec_l
I_psec_h <- I_psec_tot - I_psec_l
I_sec_vac <- diff(out[,1234])
I_psec_vac <- diff(out[,1235])
I_sec_vac_l <- diff(out[,1236])
I_psec_vac_l <- diff(out[,1237])
I_sec_vac_h <- I_sec_vac - I_sec_vac_l
I_psec_vac_h <- I_psec_vac - I_psec_vac_l}

sv.h <- sum(diff(out[,1248]), diff(out[,1260]))
sn.h <- sum(diff(out[(3650 * years_vac + 1):nrow(out),1242]), diff(out[(3650 * years_vac + 1):nrow(out),1254]), diff(out[(3650 * years_vac + 1):nrow(out),1266]))
hv.h <- sum(diff(out[,1249]), diff(out[,1261]))
hn.h <- sum(diff(out[(3650 * years_vac + 1):nrow(out),1243]), diff(out[(3650 * years_vac + 1):nrow(out),1255]), diff(out[(3650 * years_vac + 1):nrow(out),1267]))

or.h <- (sv.h / sn.h) / (hv.h / hn.h)
rr.h <- (sv.h / (sv.h + hv.h)) / (sn.h / (sn.h + hn.h))

sv.l <- sum(diff(out[,1246]), diff(out[,1258]))
sn.l <- sum(diff(out[(3650 * years_vac + 1):nrow(out),1240]), diff(out[(3650 * years_vac + 1):nrow(out),1252]), diff(out[(3650 * years_vac + 1):nrow(out),1264]))
hv.l <- sum(diff(out[,1247]), diff(out[,1259]))
hn.l <- sum(diff(out[(3650 * years_vac + 1):nrow(out),1241]), diff(out[(3650 * years_vac + 1):nrow(out),1253]), diff(out[(3650 * years_vac + 1):nrow(out),1265]))

or.l <- (sv.l / sn.l) / (hv.l / hn.l)
rr.l <- (sv.l / (sv.l + hv.l)) / (sn.l / (sn.l + hn.l))

rr_or_vec <- c(or.h, rr.h, or.l, rr.l)
save(rr_or_vec, file = paste('rr_or_', input, '.RData', sep = ''))

incidence_mat <- out[,1238:1267]
save(incidence_mat, file = paste('incidence_counts_', input, '.RData', sep = ''))



###issue with the non-sec, non-vac individuals (no way that these are being counted currently)
##this is only a problem for the vaccination individuals bc they never experience a primary inf/case
##might also need to track the number of s1 - s2 nv individuals for this? 
##summing all of these would give total number of people experiencing that infection in that time period




#Averted calcs
{
{
  #out_last <- out[nrow(out),(2:(ncol(out) - 4))]
  track_infected <- out[,(ncol(out) - 15)]
  track_l <- out[,(ncol(out) - 14)]
  track_h <- track_infected - track_l
  track_infected <- track_infected[length(track_infected)]
  track_l <- track_l[length(track_l)]
  track_h <- track_infected - track_l
  cases <- out[,(ncol(out) - 13)]
  cases <- cases[length(cases)]
  cases.l <- out[,(ncol(out) -12)]
  cases.l <- cases.l[length(cases.l)]
  cases.h <- cases - cases.l
}

##Null outputs
{
  track_infected.null <- out_null[,(ncol(out_null) - 15)]
  track_infected.null <- track_infected.null[length(track_infected.null)]
  track_l.null <- out_null[,(ncol(out_null) - 14)]
  track_l.null <- track_l.null[length(track_l.null)]
  track_h.null <- track_infected.null - track_l.null

  cases.null <- out_null[,(ncol(out_null) - 13)]
  cases.null <- cases.null[length(cases.null)]
  cases.l.null <- out_null[,(ncol(out_null) -12)]
  cases.l.null <- cases.l.null[length(cases.l.null)]
  cases.h.null <- cases.null - cases.l.null
  }

#Averted calculations, NEED TO CHECK THIS
{
  infections_averted <- (((track_infected.null - track_infected) / track_infected.null) * 100)
  infections_averted.h <- (((track_h.null - track_h) / track_h.null) * 100)
  infections_averted.l <- (((track_l.null - track_l) / track_l.null) * 100)
  cases_averted <- (((cases.null - cases) / cases.null) * 100)
  cases_averted.h <- (((cases.h.null - cases.h) / cases.h.null) * 100)
  cases_averted.l <- (((cases.l.null - cases.l) / cases.l.null) * 100)
  output <- cbind(infections_averted.h, infections_averted, infections_averted.l,
                  cases_averted.h, cases_averted, cases_averted.l)
}
  save(output, file = paste('output_', input, '.RData', sep = ''))
}

#Secondary averted calculations, CHECK THIS
track_sec_infected <- out[,(ncol(out) - 10)]
track_sec_l <- out[,(ncol(out) - 7)]
track_sec_h <- track_sec_infected - track_sec_l
track_sec_infected <- track_sec_infected[length(track_sec_infected)]
track_sec_l <- track_sec_l[length(track_sec_l)]
track_sec_h <- track_sec_infected - track_sec_l
cases_sec <- out[,(ncol(out) - 4)]
cases_sec <- cases_sec[length(cases_sec)]
cases_sec.l <- out[,(ncol(out) - 1)]
cases_sec.l <- cases_sec.l[length(cases_sec.l)]
cases_sec.h <- cases_sec - cases_sec.l

track_sec_infected.null <- out_null[,(ncol(out_null) - 10)]
track_sec_infected.null <- track_sec_infected.null[length(track_sec_infected.null)]
track_sec_l.null <- out_null[,(ncol(out_null) - 7)]
track_sec_l.null <- track_sec_l.null[length(track_sec_l.null)]
track_sec_h.null <- track_sec_infected.null - track_sec_l.null

cases.sec_null <- out_null[,(ncol(out_null) - 4)]
cases.sec_null <- cases.sec_null[length(cases.sec_null)]
cases.sec_l.null <- out_null[,(ncol(out_null) - 1)]
cases.sec_l.null <- cases.sec_l.null[length(cases.sec_l.null)]
cases.sec_h.null <- cases.sec_null - cases.sec_l.null

infections.sec_averted <- (((track_sec_infected.null - track_sec_infected) / track_sec_infected.null) * 100)
infections.sec_averted.h <- (((track_sec_h.null - track_sec_h) / track_sec_h.null) * 100)
infections.sec_averted.l <- (((track_sec_l.null - track_sec_l) / track_sec_l.null) * 100)
cases.sec_averted <- (((cases.sec_null - cases_sec) / cases.sec_null) * 100)
cases.sec_averted.h <- (((cases.sec_h.null - cases_sec.h) / cases.sec_h.null) * 100)
cases.sec_averted.l <- (((cases.sec_l.null - cases_sec.l) / cases.sec_l.null) * 100)
output.sec <- cbind(infections.sec_averted.h, infections.sec_averted, infections.sec_averted.l,
                cases.sec_averted.h, cases.sec_averted, cases.sec_averted.l)
save(output.sec, file = paste('output.sec_', input, '.RData', sep = ''))

#prop cases calculations
primary.cases <- out[,(ncol(out) - 5)]
primary.cases <- primary.cases[length(primary.cases)]
secondary.cases <- out[,(ncol(out) - 4)]
secondary.cases <- secondary.cases[length(secondary.cases)]
postsecondary.cases <- out[,(ncol(out) - 3)]
postsecondary.cases <- postsecondary.cases[length(postsecondary.cases)]


primary.l.cases <- out[,(ncol(out) - 2)]
primary.l.cases <- primary.l.cases[length(primary.l.cases)]
secondary.l.cases <- out[,(ncol(out) - 1)]
secondary.l.cases <- secondary.l.cases[length(secondary.l.cases)]
postsecondary.l.cases <- out[,(ncol(out))]
postsecondary.l.cases <- postsecondary.l.cases[length(postsecondary.l.cases)]

primary.h.cases <- primary.cases - primary.l.cases
secondary.h.cases <- secondary.cases - secondary.l.cases
postsecondary.h.cases <- postsecondary.cases - postsecondary.l.cases

prop.cases.tot <- c(primary.cases, secondary.cases, postsecondary.cases,
                    primary.l.cases, secondary.l.cases, postsecondary.l.cases,
                    primary.h.cases, secondary.h.cases, postsecondary.h.cases) / c(rep(cases, 3), rep(cases.l, 3), rep(cases.h, 3))
names(prop.cases.tot) <- c("Primary Cases", "Secondary Cases", "Postsecondary Cases",
                           "Primary Cases, High Transmission", "Secondary Cases, High Transmission", "Postsecondary Cases, High Transmission",
                           "Primary Cases, Low Transmission", "Secondary Cases, Low Transmission", "Postsecondary Cases, Low Transmission")
save(prop.cases.tot, file = paste('prop.cases_', input, '.RData', sep = ''))







