############################
#CRC info
############################
args = commandArgs(TRUE)
input = as.numeric(args[1])
beta.vec = seq(0.01, 0.5, length.out = 100)
beta = beta.vec[input]

############################
#setup
############################
#install.packages('deSolve', dependencies = TRUE, repos='http://cran.us.r-project.org')
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

parms <- c(beta = beta,
           gamma = 2.5e-1,
           sigma = 1/(365 * 1.2),
           mu = 19 / (1000 * 365),
           delta = c((0.013 / 365), rep(0.001 / 365, 4), rep(0, 10), rep(0.001, 5), 
                     0.002 / 365, 0.0025 / 365, 
                     0.004 / 365, 0.0085 / 365,
                     0.0175 / 365, 0.0425 / 365,
                     0.114 / 365, 0.151 / 365),
           age_window = c(rep(1, 21), rep(10, 7)))

############################
#MODEL
############################
model <- function(t, y, parms){
  
  beta <- parms[1]
  gamma <- parms[2]
  sigma <- parms[3]
  mu <- parms[4]
  delta <- parms[5:32]
  age_window <- parms[33:60]
  
  S1 <- y[1:28]
  I1 <- y[29:56]
  R1 <- y[57:84]
  S2 <- y[85:112]
  I2 <- y[113:140]
  R2 <- y[141:168]
  S3 <- y[169:196]
  I3 <- y[197:224]
  R3 <- y[225:252]
  S4 <- y[253:280]
  I4 <- y[281:308]
  R4 <- y[309:336]
  
  infected_total <- sum(sum(I1), sum(I2), sum(I3), sum(I4))
  population <- sum(y)
  
  #first infection
  dS1 <-  
    mu * c(population, rep(0,27)) + 
    c(0, 1/365 / head(age_window, -1) * head(S1, -1)) -
    c(1/365 / head(age_window, -1) * head(S1, -1), 0) -
    S1 * (beta * infected_total / population) -
    S1 * delta 
  
  dI1 <- 
    S1 * (beta * infected_total / population)+
    c(0, 1/365 / head(age_window, -1) * head(I1, -1)) -
    c(1/365 / head(age_window, -1) * head(I1, -1), 0) -
    I1 * gamma - 
    I1 * delta
  
  dR1 <-
    I1  * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R1, -1)) -
    c(1/365 / head(age_window, -1) * head(R1, -1), 0) -
    R1 * sigma -
    R1 * delta
  
  #second infection
  dS2 <- 
    R1 * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2, -1)) -
    c(1/365 / head(age_window, -1) * head(S2, -1), 0) -
    S2 * (3/4 * beta * infected_total / population) -
    S2 * delta
  
  dI2 <- 
    S2 * (3/4 * beta * infected_total / population) +
    c(0, 1/365 / head(age_window, -1) * head(I2, -1)) -
    c(1/365 / head(age_window, -1) * head(I2, -1), 0) -
    I2 * gamma -
    I2 * delta
  
  dR2 <-  
    I2 * gamma +  
    c(0, 1/365 / head(age_window, -1) * head(R2, -1)) -
    c(1/365 / head(age_window, -1) * head(R2, -1), 0) -
    R2 * sigma -
    R2 * delta
  
  #third infection
  dS3 <- 
    R2 * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3, -1)) -
    c(1/365 / head(age_window, -1) * head(S3, -1), 0) -
    S3 * (2/4 * beta * infected_total / population) -
    S3 * delta
  
  dI3 <- 
    S3 * (2/4 * beta * infected_total / population) +
    c(0, 1/365 / head(age_window, -1) * head(I3, -1)) -
    c(1/365 / head(age_window, -1) * head(I3, -1), 0) -
    I3 * gamma -
    I3 * delta
  
  dR3 <- 
    I3 * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3, -1)) -
    c(1/365 / head(age_window, -1) * head(R3, -1), 0) -
    R3 * sigma -
    R3 * delta
  
  #fourth infection
  dS4 <- 
    R3 * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S4, -1)) -
    c(1/365 / head(age_window, -1) * head(S4, -1), 0) -
    S4 * (1/4 * beta * infected_total / population) -
    S4 * delta 
  
  dI4 <- 
    S4 * (1/4 * beta * infected_total / population) +
    c(0, 1/365 / head(age_window, -1) * head(I4, -1)) -
    c(1/365 / head(age_window, -1) * head(I4, -1), 0) -
    I4 * gamma -
    I4 * delta
  
  dR4 <- 
    I4 * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R4, -1)) -
    c(1/365 / head(age_window, -1) * head(R4, -1), 0) -
    R4 * delta 
  
  list(c(dS1, dI1, dR1,
         dS2, dI2, dR2,
         dS3, dI3, dR3,
         dS4, dI4, dR4))
  
  
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
out <- out[nrow(out),]

############################
#OUTPUT
############################
save(out, file = paste('output_', input, '.RData', sep = ''))

