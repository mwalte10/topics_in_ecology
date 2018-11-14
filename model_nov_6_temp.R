#track_i.R
############################
#CRC info
############################
args = commandArgs(TRUE)
input = as.numeric(args[1])


##########################
#IF DOING VAC COVERAGE
##########################
redo <- c(3, 56, 57, 111, 199, 201)
# load('parms.mat.new.RData')
# beta_h <- new.parms.mat[input, 1]
# beta_l <- new.parms.mat[input, 2]
# native <- rep(new.parms.mat[input, 3], 100)
# vac_h <- new.parms.mat[input, 4]
# vac_l <- new.parms.mat[input, 5]
# load('pop.RData')
# load('birth.RData')
# load('death.RData')

load('parms.mat.new.RData')
beta_h <- new.parms.mat[redo[input], 1]
beta_l <- new.parms.mat[redo[input], 2]
native <- rep(new.parms.mat[redo[input], 3], 100)
vac_h <- new.parms.mat[redo[input], 4]
vac_l <- new.parms.mat[redo[input], 5]
load('pop.RData')
load('birth.RData')
load('death.RData')



library(deSolve)

############################
#Initial conidtions and parameters
############################
initial_conditions <- as.data.frame(matrix(NA, nrow = 100*4, ncol = 3))
initial_conditions[,2] <- rep(0:99,4)
for(i in 1:4){
  x <- i - 1
  initial_conditions[x*(100) + (1:100),1] <- rep(i,100)
}
initial_conditions[,3] <- rep(pop,4)

susceptible_init_h <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 400))
susceptible_init_h[,3] <- c(initial_conditions[1:100,3] * 0.5, initial_conditions[101:200,3] * 0, 
                            initial_conditions[201:300,3] * 0, initial_conditions[301:400,3] * 0)
susceptible_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 400))
susceptible_init_l[,3] <- c(initial_conditions[1:100,3] * 0.5, initial_conditions[101:200,3] * 0, 
                            initial_conditions[201:300,3] * 0, initial_conditions[301:400,3] * 0)

infected_init_h <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 400))
infected_init_h[,3] <- c(rep(0.01, 100), initial_conditions[101:200,3] * 0, 
                         initial_conditions[201:300,3] * 0, initial_conditions[301:400,3] * 0)
infected_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 400))
infected_init_l[,3] <- c(rep(0.01, 100), initial_conditions[101:200,3] * 0, 
                         initial_conditions[201:300,3] * 0, initial_conditions[301:400,3] * 0)

recovered_init_h <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 400))
recovered_init_h[,3] <- c(initial_conditions[1:100,3] * 0, initial_conditions[101:200,3] * 0, 
                          initial_conditions[201:300,3] * 0, initial_conditions[301:400,3] * 0)
recovered_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 400))
recovered_init_l[,3] <- c(initial_conditions[1:100,3] * 0, initial_conditions[101:200,3] * 0, 
                          initial_conditions[201:300,3] * 0, initial_conditions[301:400,3] * 0)


susceptible_h <- function(exposure){
  x <- exposure - 1
  susceptible <- c(susceptible_init_h[(x * 100) + 1:100,3])
}
susceptible_total_h <- sum(susceptible_h(1) + susceptible_h(2) + susceptible_h(3) + susceptible_h(4))
susceptible_l <- function(exposure){
  x <- exposure - 1
  susceptible <- c(susceptible_init_l[(x * 100) + 1:100,3])
}
susceptible_total_l <- sum(susceptible_l(1) + susceptible_l(2) + susceptible_l(3) + susceptible_l(4))

infected_h <- function(exposure){
  x <- exposure - 1
  infected <- c(infected_init_h[(x * 100) + 1:100,3])
} 
infected_total_h <- sum(infected_h(1) + infected_h(2) + infected_h(3) + infected_h(4))
infected_l <- function(exposure){
  x <- exposure - 1
  infected <- c(infected_init_l[(x * 100) + 1:100,3])
} 
infected_total_l <- sum(infected_l(1) + infected_l(2) + infected_l(3) + infected_l(4))

recovered_h <- function(exposure){
  x <- exposure - 1
  recovered <- c(recovered_init_h[(x * 100) + 1:100,3])
} 
recovered_total_h <- sum(recovered_h(1) + recovered_h(2) + recovered_h(3) + recovered_h(4))
recovered_l <- function(exposure){
  x <- exposure - 1
  recovered <- c(recovered_init_l[(x * 100) + 1:100,3])
} 
recovered_total_l <- sum(recovered_l(1) + recovered_l(2) + recovered_l(3) + recovered_l(4))

population_h <- sum(susceptible_total_h + infected_total_h + recovered_total_h)
population_l <- sum(susceptible_total_l + infected_total_l + recovered_total_l)



parms <- list(beta_h = beta_h,
           beta_l = beta_l,
           gamma = 1/4,
           sigma = 1/(365 * 1.2),
           mu = birth,
           delta = death,
           age_window = rep(1, 100),
           native = native,
           travel = 1 - native,
           vac_h = vac_h,
           vac_l = vac_l,
           sens = 0.85,
           spec = 0.95)
parms_null <- list(beta_h = beta_h,
              beta_l = beta_l,
              gamma = 1/4,
              sigma = 1/(365 * 1.2),
              mu = birth,
              delta = death,
              age_window = rep(1, 100),
              native = native,
              travel = 1 - native,
              vac_h = 0,
              vac_l = 0,
              sens = 0.85,
              spec = 0.95)
parms_notest <- list(beta_h = beta_h,
              beta_l = beta_l,
              gamma = 1/4,
              sigma = 1/(365 * 1.2),
              mu = birth,
              delta = death,
              age_window = rep(1, 100),
              native = native,
              travel = 1 - native,
              vac_h = vac_h,
              vac_l = vac_l,
              sens = 0,
              spec = 0)
parms_notest_null <- list(beta_h = beta_h,
                     beta_l = beta_l,
                     gamma = 1/4,
                     sigma = 1/(365 * 1.2),
                     mu = birth,
                     delta = death,
                     age_window = rep(1, 100),
                     native = native,
                     travel = 1 - native,
                     vac_h = 0,
                     vac_l = 0,
                     sens = 0,
                     spec = 0)

years = 50
years_vac = 30
times <- seq(from = 0, to = 365 * years, by = .1)
times <- times[1:(length(times) - 1)]

############################
#MODEL
############################
model <- function(t, y, parms){
  
  beta_h <- parms[[1]]
  beta_l <- parms[[2]]
  gamma <- parms[[3]]
  sigma <- parms[[4]]
  mu <- parms[[5]][floor(t / 365) + 1]
  delta <- parms[[6]][[floor(t / 365) + 1]]
  age_window <- parms[[7]]
  native <- parms[[8]]
  travel <- parms[[9]]
  t_vac.h <- ifelse((t>(365*(years_vac))), parms[[10]] / 365, 0)
  t_vac.l <- ifelse((t>(365*(years_vac))), parms[[11]] / 365, 0)
  vac_h <- c(rep(0,8), t_vac.h, rep(0,91))
  vac_l <- c(rep(0,8), t_vac.l, rep(0,91))
  sens <- parms[[12]]
  spec <- parms[[13]]
  
  #Gen pop
  S1_h <- y[which(names(y)== 'sh1')]
  I1_h <- y[which(names(y)== 'ih1')]
  R1_h <- y[which(names(y)== 'rh1')]
  S2_h <- y[which(names(y)== 'sh2')]
  I2_h <- y[which(names(y)== 'ih2')]
  R2_h <- y[which(names(y)== 'rh2')]
  S3_h <- y[which(names(y)== 'sh3')]
  I3_h <- y[which(names(y)== 'ih3')]
  R3_h <- y[which(names(y)== 'rh3')]
  S4_h <- y[which(names(y)== 'sh4')]
  I4_h <- y[which(names(y)== 'ih4')]
  R4_h <- y[which(names(y)== 'rh4')]
  
  #Vaccinated
  R1_h.v <- y[which(names(y)== 'rh1.v')]
  S2_h.v <- y[which(names(y)== 'sh2.v')]
  I2_h.v <- y[which(names(y)== 'ih2.v')]
  R2_h.v <- y[which(names(y)== 'rh2.v')]
  S3_h.v <- y[which(names(y)== 'sh3.v')]
  I3_h.v <- y[which(names(y)== 'ih3.v')]
  R3_h.v <- y[which(names(y)== 'rh3.v')]
  S4_h.v <- y[which(names(y)== 'sh4.v')]
  I4_h.v <- y[which(names(y)== 'ih4.v')]
  R4_h.v <- y[which(names(y)== 'rh4.v')]
  
  #Gen pop
  S1_l <- y[which(names(y)== 'sl1')]
  I1_l <- y[which(names(y)== 'il1')]
  R1_l <- y[which(names(y)== 'rl1')]
  S2_l <- y[which(names(y)== 'sl2')]
  I2_l <- y[which(names(y)== 'il2')]
  R2_l <- y[which(names(y)== 'rl2')]
  S3_l <- y[which(names(y)== 'sl3')]
  I3_l <- y[which(names(y)== 'il3')]
  R3_l <- y[which(names(y)== 'rl3')]
  S4_l <- y[which(names(y)== 'sl4')]
  I4_l <- y[which(names(y)== 'il4')]
  R4_l <- y[which(names(y)== 'rl4')]
  
  #Vaccinated
  R1_l.v <- y[which(names(y)== 'rl1.v')]
  S2_l.v <- y[which(names(y)== 'sl2.v')]
  I2_l.v <- y[which(names(y)== 'il2.v')]
  R2_l.v <- y[which(names(y)== 'rl2.v')]
  S3_l.v <- y[which(names(y)== 'sl3.v')]
  I3_l.v <- y[which(names(y)== 'il3.v')]
  R3_l.v <- y[which(names(y)== 'rl3.v')]
  S4_l.v <- y[which(names(y)== 'sl4.v')]
  I4_l.v <- y[which(names(y)== 'il4.v')]
  R4_l.v <- y[which(names(y)== 'rl4.v')]
  
  ###Tracking infected individuals

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
    mu * c(birth_pop_h, rep(0,99)) +
    c(0, 1/365 / head(age_window, -1) * head(S1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S1_h, -1), 0) -
    native * S1_h * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) -
    travel * S1_h * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) -
    S1_h * delta -
    S1_h * vac_h * (1 - spec)
  dS1_l <-  
    mu * c(birth_pop_l, rep(0,99)) +
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
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
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
  
  I_h <- 
    native * S2_h * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S2_h.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_h.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S1_h * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S1_h * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_h * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S3_h.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_h.v * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_h * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h * beta_l * 0.25 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    native * S4_h.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_h.v * beta_l * 0.5 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) 
  I_h <- sum(I_h)
  
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
  primary_inf <- primary/ 0.18

  
  
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
  secondary_inf <- secondary / 0.41
  
  
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
  post_sec_inf <- post_sec / 0.063425
  
  primary.h <- 
  native * S1_h * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
  travel * S1_l * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  primary.h <- sum(primary.h) * 0.18
  
  secondary.h <- 
    native * S2_h * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_l * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S2_h.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  secondary.h <- sum(secondary.h) * 0.41
  
  post_sec.h <- 
    native * S3_h * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_l * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S3_h.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S3_l.v * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
    native * S4_h * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_l * beta_h * 0.25 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S4_h.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    travel * S4_l.v * beta_h * 0.5 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)
  post_sec.h <- sum(post_sec.h) * 0.063425
  
  cases.h <- primary.h + secondary.h + post_sec.h

  primary.l <-
    native * S1_l * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S1_l * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l)  
  primary.l <- sum(primary.l) * 0.18
  primary.l.inf <- primary.l / 0.18
  
  secondary.l <-
    native * S2_l * beta_l * 0.75 * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) +
    native * S2_l.v * beta_l * (native * infected_total_l / pop_l + travel * infected_total_h / pop_h) +
    travel * S2_l.v * beta_h * (native * infected_total_h / pop_h + travel * infected_total_l / pop_l) 
  secondary.l <- sum(secondary.l) * 0.41
  secondary.l.inf <- secondary.l / 0.41
  
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
  post_sec.l <- tertiary.l + quaternary.l
  post_sec.l.inf <- post_sec.l / 0.063425
  }

    
  
  list(c(
         dS1_h, dI1_h, dR1_h,
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
         
         primary_tot.cases, primary_tot.ncases, 
         primary_cases.l, primary_ncases.l,
         primary_cases.h, primary_ncases.h,
         
         secondary_vac_tot.cases, secondary_vac_tot.ncases,
         secondary_vac_cases.l, secondary_vac_ncases.l,
         secondary_vac_cases.h, secondary_vac_ncases.h,
         secondary_tot.cases, secondary_tot.ncases,
         secondary_cases.l, secondary_ncases.l,
         secondary_cases.h, secondary_ncases.h,
         
         postsec_vac_tot.cases, postsec_vac_tot.ncases,
         postsec_vac_cases.l, postsec_vac_ncases.l,
         postsec_vac_cases.h, postsec_vac_ncases.h,
         postsec_tot_cases, postsec_tot_ncases,
         postsec_cases.l, postsec_ncases.l,
         postsec_cases.h, postsec_ncases.h,
         
         cases.h,
         
         I_total, I_l, cases, cases.l,
         I_primary_tot, I_secondary_tot, I_post_sec_tot,
         I_l_primary_tot, I_l_sec_tot, I_l_post_sec_tot,
         primary, secondary, post_sec,
         primary.l, secondary.l, post_sec.l, 
         primary_inf, secondary_inf, post_sec_inf, 
         primary.l.inf, secondary.l.inf, post_sec.l.inf,
         I_h
         ))

  
}

############################
#RUN MODEL
############################
y_init <- c(susceptible_h(1), infected_h(1), recovered_h(1),
            susceptible_h(2), infected_h(2), recovered_h(2),
            susceptible_h(3), infected_h(3), recovered_h(3),
            susceptible_h(4), infected_h(4), recovered_h(4),
            rep(0, 100), rep(0, 100), rep(0, 100), rep(0, 100), rep(0, 100),
            rep(0, 100), rep(0, 100), rep(0, 100), rep(0, 100), rep(0, 100),
            susceptible_l(1), infected_l(1), recovered_l(1),
            susceptible_l(2), infected_l(2), recovered_l(2),
            susceptible_l(3), infected_l(3), recovered_l(3),
            susceptible_l(4), infected_l(4), recovered_l(4),
            rep(0, 100), rep(0, 100), rep(0, 100), rep(0, 100), rep(0, 100),
            rep(0, 100), rep(0, 100), rep(0, 100), rep(0, 100), rep(0, 100),
            0,
            rep(0, 30),
            rep(0,27))
names(y_init) <- c(rep('sh1', 100), rep('ih1', 100), rep('rh1', 100),
                   rep('sh2', 100), rep('ih2', 100), rep('rh2', 100),
                   rep('sh3', 100), rep('ih3', 100), rep('rh3', 100),
                   rep('sh4', 100), rep('ih4', 100), rep('rh4', 100),
                   rep('rh1.v', 100),
                   rep('sh2.v', 100), rep('ih2.v', 100), rep('rh2.v', 100),
                   rep('sh3.v', 100), rep('ih3.v', 100), rep('rh3.v', 100),
                   rep('sh4.v', 100), rep('ih4.v', 100), rep('rh4.v', 100),
                   rep('sl1', 100), rep('il1', 100), rep('rl1', 100),
                   rep('sl2', 100), rep('il2', 100), rep('rl2', 100),
                   rep('sl3', 100), rep('il3', 100), rep('rl3', 100),
                   rep('sl4', 100), rep('il4', 100), rep('rl4', 100),
                   rep('rl1.v', 100),
                   rep('sl2.v', 100), rep('il2.v', 100), rep('rl2.v', 100),
                   rep('sl3.v', 100), rep('il3.v', 100), rep('rl3.v', 100),
                   rep('sl4.v', 100), rep('il4.v', 100), rep('rl4.v', 100),
                   'i_sec_vac', 'i_psec_vac',
                   'il_sec_vac', 'il_psec_vac',
                   'prim_tot.cases', 'prim_tot.ncases',
                   'prim_cases.l', 'prim_ncases.l',
                   'prim_cases.h', 'prim_ncases.h',
                   
                   'sec_tot.cases.v', 'sec_tot.ncases.v',
                   'sec_vac.cases.l', 'sec_vac.ncases.l',
                   'sec_vac.cases.h', 'sec_vac.ncases.h',
                   'sec_tot.cases', 'sec_tot.ncases',
                   'sec_cases.l', 'sec_ncases.l',
                   'sec_cases.h', 'sec_ncases.h',
                   
                   'psec_vac_tot.cases', 'psec_vac_tot.ncases',
                   'psec_vac.cases.l', 'psec_vac.ncases.l',
                   'psec_vac.cases.h', 'psec_vac.ncases.h',
                   'psec_tot.cases', 'psec_tot.ncases',
                   'psec_cases.l', 'psec_ncases.l',
                   'psec_cases.h', 'psec_ncases.h',
                   
                   'cases.h',
                   
                   'i_total', 'il', 'cases', 'cases.l',
                   'i_prim_tot', 'i_sec_tot', 'i_post_sec_tot',
                   'il_prim_tot', 'il_sec_tot', 'il_psec_tot',
                   'prim', 'sec', 'psec',
                   'prim.l', 'sec.l', 'psec.l',
                   'prim_inf', 'sec_inf', 'psec_inf',
                   'prim.l.inf', 'sec.l.inf', 'psec.l.inf',
                   'ih')
years = 50
years_vac = 30
out <- ode(times = times, y = y_init, func = model, parms = parms)
out_null <- ode(times = times, y = y_init, func = model, parms = parms_null)
out_notest <-  ode(times = times, y = y_init, func = model, parms = parms_notest)
out_notest_null <-  ode(times = times, y = y_init, func = model, parms = parms_notest_null)



# out_last <- out[nrow(out),(2:(ncol(out)))]
# out_last.null <- out_null[nrow(out_null),(2:(ncol(out_null)))]
# save(out_last, file = paste('out_last_', input, '.RData', sep = ''))
# save(out_last.null, file = paste('out_last.null_', input, '.RData', sep = ''))


##incidence calcs 

out <- out[,2:ncol(out)]
out_null <- out_null[,2:ncol(out_null)]
# 
orrr.h_calc <- function(out_mat){
  sv.h <- sum(diff(out_mat[,which(colnames(out_mat) == 'sec_vac.cases.h')]), diff(out_mat[,which(colnames(out_mat) == 'psec_vac.cases.h')]))
  sn.h <- sum(diff(out_mat[(3650 * years_vac + 1):nrow(out_mat),which(colnames(out_mat) == 'sec_cases.h')]), 
              diff(out_mat[(3650 * years_vac + 1):nrow(out_mat),which(colnames(out_mat) == 'prim_cases.h')]), 
              diff(out_mat[(3650 * years_vac + 1):nrow(out_mat),which(colnames(out_mat) == 'psec_cases.h')]))
  hv.h <- sum(diff(out_mat[,which(colnames(out_mat) == 'sec_vac.ncases.h')]), diff(out_mat[,which(colnames(out_mat) == 'psec_vac.ncases.h')]))
  hn.h <- sum(diff(out_mat[(3650 * years_vac + 1):nrow(out_mat),which(colnames(out_mat) == 'prim_ncases.h')]), 
              diff(out_mat[(3650 * years_vac + 1):nrow(out_mat),which(colnames(out_mat) == 'sec_ncases.h')]), 
              diff(out_mat[(3650 * years_vac + 1):nrow(out_mat),which(colnames(out_mat) == 'psec_ncases.h')]))
  rr.h <- (sv.h / (sv.h + hv.h)) / (sn.h / (sn.h + hn.h))
  or.h <- (sv.h / sn.h) / (hv.h / hn.h)
  return(c(or.h, rr.h))
}

orrr.l_calc <- function(out_mat){
  sv.l <- sum(diff(out_mat[,which(colnames(out_mat) == 'sec_vac.cases.l')]), diff(out_mat[,which(colnames(out_mat) == 'psec_vac.cases.l')]))
  sn.l <- sum(diff(out_mat[(3650 * years_vac + 1):nrow(out_mat),which(colnames(out_mat) == 'sec_cases.l')]), 
              diff(out_mat[(3650 * years_vac + 1):nrow(out_mat),which(colnames(out_mat) == 'prim_cases.l')]), 
              diff(out_mat[(3650 * years_vac + 1):nrow(out_mat),which(colnames(out_mat) == 'psec_cases.l')]))
  hv.l <- sum(diff(out_mat[,which(colnames(out_mat) == 'sec_vac.ncases.l')]), diff(out_mat[,which(colnames(out_mat) == 'psec_vac.ncases.l')]))
  hn.l <- sum(diff(out_mat[(3650 * years_vac + 1):nrow(out_mat),which(colnames(out_mat) == 'prim_ncases.l')]), 
              diff(out_mat[(3650 * years_vac + 1):nrow(out_mat),which(colnames(out_mat) == 'sec_ncases.l')]), 
              diff(out_mat[(3650 * years_vac + 1):nrow(out_mat),which(colnames(out_mat) == 'psec_ncases.l')]))
  rr.l <- (sv.l / (sv.l + hv.l)) / (sn.l / (sn.l + hn.l))
  or.l <- (sv.l / sn.l) / (hv.l / hn.l)
  return(c(or.l, rr.l))
}

rr_or_vec <- c(orrr.h_calc(out), orrr.l_calc(out))
save(rr_or_vec, file = paste('rr_or_', redo[input], '.RData', sep = ''))

rr_or_vac.notest <- c(orrr.h_calc(out_notest), orrr.l_calc(out_notest))
save(rr_or_vac.notest, file = paste('rr_or.notest_', redo[input], '.RData', sep = ''))


indexing <- c((3650 * years_vac + 1):(nrow(out) - 1))
# # 
cases_averted.func <- function(out_mat, out_mat_null){
  track_infected <- sum(diff(out_mat[indexing, which(colnames(out_mat) == 'i_total')]))
  track_l <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'il')]))
  track_h <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'ih')]))
  cases <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'cases')]))
  cases.l <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'cases.l')]))
  cases.h <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'cases.h')]))
  
  track_infected.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'i_total')]))
  track_l.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'il')]))
  track_h.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'ih')]))
  cases.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'cases')]))
  cases.l.null <- sum(diff(out_mat_null[indexing,which(colnames(out) == 'cases.l')]))
  cases.h.null <- sum(diff(out_mat_null[indexing,which(colnames(out) == 'cases.h')]))
  
  
  infections_averted <- ((track_infected.null - track_infected) / track_infected.null) * 100
  infections_averted.h <- ((track_h.null - track_h) / track_h.null) * 100
  infections_averted.l <- ((track_l.null - track_l) / track_l.null) * 100
  cases_averted <- ((cases.null - cases) / cases.null) * 100
  cases_averted.h <- ((cases.h.null - cases.h) / cases.h.null) * 100
  cases_averted.l <- ((cases.l.null - cases.l) / cases.l.null) * 100
  output <- c(infections_averted.h, infections_averted, infections_averted.l,
              cases_averted.h, cases_averted, cases_averted.l)
  
  return(output)
}

# #

output <- cases_averted.func(out, out_null)
save(output, file = paste('output_', redo[input], '.RData', sep = ''))

output.notest <- cases_averted.func(out_notest, out_notest_null)
save(output.notest, file = paste('output.notest_', redo[input], '.RData', sep = ''))

 
##time series, do for input 120

# infections_h.save.vac <- out[,c(29:56, 113:140, 197:224, 281:308, 393:420, 477:504, 561:588)]
# infections_l.save.vac <- 

  # infections_h.vac <- rowSums(out[,c(29:56, 113:140, 197:224, 281:308, 393:420, 477:504, 561:588)], na.rm=TRUE)
  # infections_l.vac <- rowSums(out[,c(645:672, 729:756, 813:840, 897:924, 1009:1036, 1093:1120, 1177:1204)], na.rm=TRUE)
  # 
  # 
  # 
  # infections_h <- rowSums(out_null[,c(29:56, 113:140, 197:224, 281:308, 393:420, 477:504, 561:588)], na.rm=TRUE)
  # infections_l <- rowSums(out_null[,c(645:672, 729:756, 813:840, 897:924, 1009:1036, 1093:1120, 1177:1204)], na.rm=TRUE)
  # 
  # 
  # ts_inf <- list(infections_h.vac, infections_l.vac, infections_h, infections_l)
  # # names(ts_inf) <- c('high ses, vac', 'low ses, vac',
  # #                    'high ses, nvac', 'low ses, nvac')
  # 
  # save(ts_inf, file = paste('ts_inf_', x[input], '.RData', sep = ''))



{

# 
# #prop cases calculations

# primary.cases <- out[nrow(out),1278]
# secondary.cases <- out[nrow(out),1279]
# postsecondary.cases <- out[nrow(out),1280]
# cases <- primary.cases + secondary.cases + postsecondary.cases
# 
# primary.l.cases <- out[nrow(out),1281]
# secondary.l.cases <- out[nrow(out),1282]
# postsecondary.l.cases <- out[nrow(out),1283]
# cases.l <- primary.l.cases + secondary.l.cases + postsecondary.l.cases
# 
# primary.h.cases <- primary.cases - primary.l.cases
# secondary.h.cases <- secondary.cases - secondary.l.cases
# postsecondary.h.cases <- postsecondary.cases - postsecondary.l.cases
# cases.h <- cases - cases.l
# 
# prop.cases.tot <- c(primary.cases / cases, secondary.cases / cases, postsecondary.cases / cases,
#                     primary.l.cases / cases.l, secondary.l.cases / cases.l, postsecondary.l.cases / cases.l,
#                     primary.h.cases / cases.h, secondary.h.cases / cases.h, postsecondary.h.cases / cases.h) 
# names(prop.cases.tot) <- c("Primary Cases", "Secondary Cases", "Postsecondary Cases",
#                            "Primary Cases, High Transmission", "Secondary Cases, High Transmission", "Postsecondary Cases, High Transmission",
#                            "Primary Cases, Low Transmission", "Secondary Cases, Low Transmission", "Postsecondary Cases, Low Transmission")
# save(prop.cases.tot, file = paste('prop.cases_', input, '.RData', sep = ''))
# 
# ####Infections
# primary.inf <- out[nrow(out),1278]
# secondary.inf <- out[nrow(out),1279]
# postsecondary.inf<- out[nrow(out),1280]
# inf <- primary.inf + secondary.inf + postsecondary.inf
# 
# primary.l.inf<- out[nrow(out),1281]
# secondary.l.inf <- out[nrow(out),1282]
# postsecondary.l.inf <- out[nrow(out),1283]
# inf.l <- primary.l.inf + secondary.l.inf + postsecondary.l.inf
# 
# primary.h.inf <- primary.inf - primary.l.inf
# secondary.h.inf <- secondary.inf - secondary.l.inf
# postsecondary.h.inf <- postsecondary.inf - postsecondary.l.inf
# inf.h <- inf - inf.l
# 
# prop.inf.tot <- c(primary.inf / inf, secondary.inf / inf, postsecondary.inf / inf,
#                     primary.l.inf / inf.l, secondary.l.inf / inf.l, postsecondary.l.inf / inf.l,
#                     primary.h.inf / inf.h, secondary.h.inf / inf.h, postsecondary.h.inf / inf.h) 
# names(prop.inf.tot) <- c("Primary Infections", "Secondary Infections", "Postsecondary Infections",
#                            "Primary Infections, High Transmission", "Secondary Infections, High Transmission", "Postsecondary Infections, High Transmission",
#                            "Primary Infections, Low Transmission", "Secondary Infections, Low Transmission", "Postsecondary Infections, Low Transmission")
# save(prop.inf.tot, file = paste('prop.inf_', input, '.RData', sep = ''))
# 
# 
# primary.cases <- out_null[nrow(out),1278]
# secondary.cases <- out_null[nrow(out),1279]
# postsecondary.cases <- out_null[nrow(out),1280]
# cases <- primary.cases + secondary.cases + postsecondary.cases
# 
# primary.l.cases <- out_null[nrow(out),1281]
# secondary.l.cases <- out_null[nrow(out),1282]
# postsecondary.l.cases <- out_null[nrow(out),1283]
# cases.l <- primary.l.cases + secondary.l.cases + postsecondary.l.cases
# 
# primary.h.cases <- primary.cases - primary.l.cases
# secondary.h.cases <- secondary.cases - secondary.l.cases
# postsecondary.h.cases <- postsecondary.cases - postsecondary.l.cases
# cases.h <- cases - cases.l
# 
# prop.cases.tot <- c(primary.cases / cases, secondary.cases / cases, postsecondary.cases / cases,
#                     primary.l.cases / cases.l, secondary.l.cases / cases.l, postsecondary.l.cases / cases.l,
#                     primary.h.cases / cases.h, secondary.h.cases / cases.h, postsecondary.h.cases / cases.h) 
# names(prop.cases.tot) <- c("Primary Cases", "Secondary Cases", "Postsecondary Cases",
#                            "Primary Cases, High Transmission", "Secondary Cases, High Transmission", "Postsecondary Cases, High Transmission",
#                            "Primary Cases, Low Transmission", "Secondary Cases, Low Transmission", "Postsecondary Cases, Low Transmission")
# save(prop.cases.tot, file = paste('prop.cases.null_', input, '.RData', sep = ''))
# 
# ####Infections
# primary.inf <- out_null[nrow(out),1278]
# secondary.inf <- out_null[nrow(out),1279]
# postsecondary.inf<- out_null[nrow(out),1280]
# inf <- primary.inf + secondary.inf + postsecondary.inf
# 
# primary.l.inf<- out_null[nrow(out),1281]
# secondary.l.inf <- out_null[nrow(out),1282]
# postsecondary.l.inf <- out_null[nrow(out),1283]
# inf.l <- primary.l.inf + secondary.l.inf + postsecondary.l.inf
# 
# primary.h.inf <- primary.inf - primary.l.inf
# secondary.h.inf <- secondary.inf - secondary.l.inf
# postsecondary.h.inf <- postsecondary.inf - postsecondary.l.inf
# inf.h <- inf - inf.l
# 
# prop.inf.tot <- c(primary.inf / inf, secondary.inf / inf, postsecondary.inf / inf,
#                   primary.l.inf / inf.l, secondary.l.inf / inf.l, postsecondary.l.inf / inf.l,
#                   primary.h.inf / inf.h, secondary.h.inf / inf.h, postsecondary.h.inf / inf.h) 
# names(prop.inf.tot) <- c("Primary Infections", "Secondary Infections", "Postsecondary Infections",
#                          "Primary Infections, High Transmission", "Secondary Infections, High Transmission", "Postsecondary Infections, High Transmission",
#                          "Primary Infections, Low Transmission", "Secondary Infections, Low Transmission", "Postsecondary Infections, Low Transmission")
# save(prop.inf.tot, file = paste('prop.inf.null_', input, '.RData', sep = ''))


}

































