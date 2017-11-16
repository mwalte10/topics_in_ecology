############################
#CRC info
############################
# args = commandArgs(TRUE)
# input = as.numeric(args[1])
df = expand.grid(beta_l = seq(0.1, 2, length.out = 20),
				 beta_h = seq(0.1, 2, length.out = 20))
beta_l = df[input,1]
beta_h = df[input,2]

############################
#setup
############################
# install.packages('deSolve', dependencies = TRUE, repos='http://cran.us.r-project.org')
library(deSolve)


############################
#Initial condidtions and parameters
############################
percentage_vec <- c(rep(1.92,5), rep(1.84,5), rep(1.8,5), rep(1.74,5),
                    17.5, 15, 12.1, 8.9, 5.7, 3, 1.4, 0.2)
percentage_vec <- percentage_vec / 100
initial_conditions <- as.data.frame(matrix(NA, nrow = 28*4, ncol = 3))
initial_conditions[,2] <- rep(1:28,4)
for(i in 1:4){
  x <- i - 1
  initial_conditions[x*(28) + (1:28),1] <- rep(i,28)
}
initial_conditions[,3] <- rep(percentage_vec,4)

susceptible_init <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 112))
susceptible_init[,3] <- c(initial_conditions[1:28,3] * 1 *10^6, initial_conditions[29:56,3] * 0, 
                          initial_conditions[57:84,3] * 0, initial_conditions[85:112,3] * 0)

infected_init <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 112))
infected_init[,3] <- c(initial_conditions[1:28,3] * 1, initial_conditions[29:56,3] * 0, 
                       initial_conditions[57:84,3] * 0, initial_conditions[85:112,3] * 0)

recovered_init <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 112))
recovered_init[,3] <- c(initial_conditions[1:28,3] * 0, initial_conditions[29:56,3] * 0, 
                        initial_conditions[57:84,3] * 0, initial_conditions[85:112,3] * 0)


susceptible <- function(exposure){
  x <- exposure - 1
  susceptible <- c(susceptible_init[(x * 28) + 1:28,3])
}
susceptible_total <- sum(susceptible(1) + susceptible(2) + susceptible(3) + susceptible(4))

infected <- function(exposure){
  x <- exposure - 1
  infected <- c(infected_init[(x * 28) + 1:28,3])
} 
infected_total <- sum(infected(1) + infected(2) + infected(3) + infected(4))

recovered <- function(exposure){
  x <- exposure - 1
  recovered <- c(recovered_init[(x * 28) + 1:28,3])
} 
recovered_total <- sum(recovered(1) + recovered(2) + recovered(3) + recovered(4))

population <- sum(susceptible_total + infected_total + recovered_total)

native <- c(rep(1, 7), rep(0.86, 9), rep(0.842, 4), 0.814, 0.7676, 0.7784, rep(0.809, 5))

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
  
  S1_l <- y[337:364]
  I1_l <- y[365:392]
  R1_l <- y[393:420]
  S2_l <- y[421:448]
  I2_l <- y[449:476]
  R2_l <- y[477:504]
  S3_l <- y[505:532]
  I3_l <- y[533:560]
  R3_l <- y[561:588]
  S4_l <- y[589:616]
  I4_l <- y[617:644]
  R4_l <- y[645:672]
  
  infected_total_h <- sum(sum(I1_h), sum(I2_h), sum(I3_h), sum(I4_h))
  population_h <- y[1:336]
  pop_h <- sum(population_h)
  infected_total_l <- sum(sum(I1_l), sum(I2_l), sum(I3_l), sum(I4_l))
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
  
  dR1_h <-
    I1_h * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(R1_h, -1), 0) -
    R1_h * sigma -
    R1_h * delta
  dR1_l <-
    I1_h * gamma +
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
  
  list(c(dS1_h, dI1_h, dR1_h,
         dS2_h, dI2_h, dR2_h,
         dS3_h, dI3_h, dR3_h,
         dS4_h, dI4_h, dR4_h,
         dS1_l, dI1_l, dR1_l,
         dS2_l, dI2_l, dR2_l,
         dS3_l, dI3_l, dR3_l,
         dS4_l, dI4_l, dR4_l))
  
}

############################
#RUN MODEL
############################
y_init <- c(susceptible(1), infected(1), recovered(1),
            susceptible(2), infected(2), recovered(2),
            susceptible(3), infected(3), recovered(3),
            susceptible(4), infected(4), recovered(4))
times <- seq(from = 0, to = 365 * 100, by = .1)
out <- ode(times = times, y = y_init, func = model, parms = parms)


############################
#OUTPUT
############################
save(out, file = paste('output.movement_', input, '.RData', sep = ''))