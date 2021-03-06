############################
#CRC info
############################
args = commandArgs(TRUE)
input = as.numeric(args[1])


##########################
#IF DOING VAC COVERAGE
##########################

beta_h <- 0.15
beta_l <- 0.25
native <- 1 
travel <- 0
vac_h <- 0.5
vac_l <- 0.5
spec <- 0.85


# load('pop_1950.RData')
load('birth_1950.RData')
birth <- birth[51:80]
load('death_1950.RData')
death <- death[51:80]
# source("formulas_for_model.R")


hopkins <- c(0.53, 1, 0.115)
hopkins_inverse <- 1 - hopkins


library(deSolve)





############################
#Initial conidtions and parameters
###########################
load('last_row_1.RData')
y_init <- c(last_row[1:3520])


parms.h <- list(beta_h = beta_h,
                beta_l = beta_l,
                gamma = 1/4,
                sigma = 1/(365 * 1.2),
                mu = birth,
                delta = death,
                age_window = rep(1, 80),
                native = native,
                travel = 1 - native,
                vac_h = vac_h,
                vac_l = vac_l,
                vac_h.test = vac_h,
                vac_l.test = vac_l, 
                sens = 0.85,
                spec = spec,
                hopkins,
                hopkins_inverse)
parms_null.h <- list(beta_h = beta_h,
                     beta_l = beta_l,
                     gamma = 1/4,
                     sigma = 1/(365 * 1.2),
                     mu = birth,
                     delta = death,
                     age_window = rep(1, 80),
                     native = native,
                     travel = 1 - native,
                     vac_h = 0,
                     vac_l = 0,
                     vac_h.test = vac_h,
                     vac_l.test = vac_l, 
                     sens = 0.85,
                     spec = spec,
                     hopkins,
                     hopkins_inverse)

years = 30
years_vac = 0
times <- seq(from = 0, to = 365 * years, by = 0.1)
times <- times[1:(length(times) - 1)]

############################
#MODEL
############################
model <- function(t, y, parms, null){
  
  #############################
  ##SET UP PARAMETERS 
  #############################
  beta_h <- parms[[1]]
  beta_l <- parms[[2]]
  gamma <- parms[[3]]
  sigma <- parms[[4]]
  mu <- parms[[5]][floor(t / 365) + 1]
  delta <- parms[[6]][[floor(t / 365) + 1]]
  age_window <- parms[[7]]
  native <- parms[[8]]
  native_h <- parms[[8]]
  native_l <- parms[[8]]
  travel <- parms[[9]]
  travel_l <- parms[[9]]
  travel_h <- parms[[9]]
  t_vac.h <- ifelse((t>(365*(years_vac))), parms[[10]] / 365, 0)
  t_vac.l <- ifelse((t>(365*(years_vac))), parms[[11]] / 365, 0)
  vac_h <- c(rep(0,9), t_vac.h, rep(0,70))
  vac_l <- c(rep(0,9), t_vac.l, rep(0,70))
  
  t_vac.h.test <- ifelse((t>(365*(years_vac))), parms[[12]] / 365, 0)
  t_vac.l.test <- ifelse((t>(365*(years_vac))), parms[[13]] / 365, 0)
  vac_h.test <- c(rep(0,9), t_vac.h.test, rep(0,70))
  vac_l.test <- c(rep(0,9), t_vac.l.test, rep(0,70))

  sens <- parms[[14]]
  spec <- parms[[15]]
  inf <- parms[[16]]
  ninf <- parms[[17]]
  
#############################
##SET UP STATE VARIABLES
#############################
  #GEN POP
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
  
  #############################
  ##ESTABLISH SYMPTOMATIC INDIVIDUALS
  #############################
  infecteds_h.prim.c <- I1_h * inf[1]
  infecteds_h.sec.c <- cbind(I2_h, I2_h.v) * inf[2]
  infecteds_h.psec.c <- cbind(I3_h, I4_h, I3_h.v, I4_h.v) * inf[3]
  sym_inf_h <- cbind(infecteds_h.prim.c, infecteds_h.sec.c, infecteds_h.psec.c)
  sym_inf_h <- sum(sym_inf_h)
  
  infecteds_l.prim.c <- I1_l * inf[1]
  infecteds_l.sec.c <- cbind(I2_l, I2_l.v) * inf[2]
  infecteds_l.psec.c <- cbind(I3_l, I4_l, I3_l.v, I4_l.v) * inf[3]
  sym_inf_l <- cbind(infecteds_l.prim.c, infecteds_l.sec.c, infecteds_l.psec.c)
  sym_inf_l <- sum(sym_inf_l)
  
  
  #############################
  ##ESTABLISH ASYMPTOMATIC INDIVIDUALS
  #############################
  infecteds_h.prim <- I1_h * ninf[1]
  infecteds_h.sec <- cbind(I2_h, I2_h.v) * ninf[2]
  infecteds_h.psec <- cbind(I3_h, I4_h, I3_h.v, I4_h.v) * ninf[3]
  inf_h <- cbind(infecteds_h.prim, infecteds_h.sec, infecteds_h.psec)
  inf_h <- sum(inf_h)
  
  infecteds_l.prim <- I1_l * ninf[1]
  infecteds_l.sec <- cbind(I2_l, I2_l.v) * ninf[2]
  infecteds_l.psec <- cbind(I3_l, I4_l, I3_l.v, I4_l.v) * ninf[3]
  inf_l <- cbind(infecteds_l.prim, infecteds_l.sec, infecteds_l.psec)
  inf_l <- sum(inf_l)
  
  
  #############################
  ##ESTABLISH POPULATION
  #############################
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
  pop_h <- sum(pop_h)

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
  pop_l <- sum(pop_l)

  
  
  
  #############################
  ##FIRST INFECTION
  #############################
  dS1_h <- 
    #BIRTH 
    mu * c(birth_pop_h, rep(0,79)) +
    #AGE IN
    c(0, 1/365 / head(age_window, -1) * head(S1_h, -1)) -
    #AGE OUT
    c(1/365 / head(age_window, -1) * head(S1_h, -1), 0) -
    #INFECTIONS
    native * S1_h * 2 * beta_h * ((sym_inf_h) / pop_h) -
    travel * S1_h * 2 * beta_l * ((sym_inf_l) / pop_l) -
    native * S1_h * beta_h * ((inf_h) / pop_h) -
    travel * S1_h * beta_l * ((inf_l) / pop_l) -
    #DEATH
    S1_h * delta -
    #VACCINATION
    S1_h * vac_h * (1 - spec)
  dS1_l <-  
    mu * c(birth_pop_l, rep(0,79)) +
    c(0, 1/365 / head(age_window, -1) * head(S1_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S1_l, -1), 0)  -
    native * S1_l * 2 * beta_l * ((sym_inf_l) / pop_l) -
    travel * S1_l * 2 * beta_h * ((sym_inf_h)/ pop_h) -
    native * S1_l * beta_l * ((inf_l) / pop_l) -
    travel * S1_l * beta_h * ((inf_h)/ pop_h) -
    S1_l * delta - 
    S1_l * vac_l * (1 - spec)
  
  dI1_h <- 
    native * S1_h * 2 * beta_h * ((sym_inf_h) / pop_h) +
    travel * S1_h * 2 * beta_l * ((sym_inf_l) / pop_l) +
    native * S1_h * beta_h * ((inf_h) / pop_h) +
    travel * S1_h * beta_l * ((inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I1_h, -1), 0) -
    #RECOVERY
    I1_h * gamma -
    I1_h * delta
  dI1_l <- 
    native * S1_l * 2 * beta_l * ((sym_inf_l) / pop_l) +
    travel * S1_l * 2 * beta_h * ((sym_inf_h)/ pop_h) +
    native * S1_l * beta_l * ((inf_l) / pop_l) +
    travel * S1_l * beta_h * ((inf_h)/ pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I1_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I1_l, -1), 0) -
    I1_l * gamma -
    I1_l * delta
  
  dR1_h <-
    I1_h * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(R1_h, -1), 0) -
    #LOSS OF CROSS PROTECTIVE IMMUNITY
    R1_h * sigma -
    R1_h * delta - 
    #VACCINATION
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
  
  #############################
  ##SECOND INFECTION
  #############################
  dS2_h <- 
    R1_h * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_h, -1), 0) -
    native * S2_h * 2 * 0.75 * beta_h * ((sym_inf_h) / pop_h) -
    travel * S2_h * 2 * 0.75 * beta_l * ((sym_inf_l) / pop_l) -
    native * S2_h * 0.75 * beta_h * ((inf_h) / pop_h) -
    travel * S2_h * 0.75 * beta_l * ((inf_l) / pop_l) -
    S2_h * delta - 
    S2_h * vac_h * sens
  dS2_h.v <-
    R1_h.v * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_h.v, -1), 0) -
    native * S2_h.v * 2 * beta_h * ((sym_inf_h) / pop_h) -
    travel * S2_h.v * 2  * beta_l * ((sym_inf_l) / pop_l) -
    native * S2_h.v  * beta_h * ((inf_h) / pop_h) -
    travel * S2_h.v * beta_l * ((inf_l) / pop_l) -
    S2_h.v * delta
  dS2_l <- 
    R1_l * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_l, -1), 0) -
    native * S2_l * beta_l * 2 * 0.75 * ((sym_inf_l) / pop_l) -
    travel * S2_l * beta_h * 2 *  0.75 * ((sym_inf_h)  / pop_h) -
    native * S2_l * beta_l * 0.75 * ((inf_l) / pop_l) -
    travel * S2_l * beta_h * 0.75 * ((inf_h)  / pop_h) -
    S2_l * delta - 
    S2_l * vac_l * sens
  dS2_l.v <-
    R1_l.v * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_l.v, -1), 0) -
    native * S2_l.v * beta_l * 2 * ((sym_inf_l) / pop_l) -
    travel * S2_l.v * beta_h * 2  * ((sym_inf_h)  / pop_h) -
    native * S2_l.v * beta_l  * ((inf_l) / pop_l) -
    travel * S2_l.v * beta_h * ((inf_h)  / pop_h) -
    S2_l.v * delta
  
  dI2_h <- 
    native * S2_h * 2 * 0.75 * beta_h * ((sym_inf_h) / pop_h) +
    travel * S2_h * 2 * 0.75 * beta_l * ((sym_inf_l) / pop_l) +
    native * S2_h * 0.75 * beta_h * ((inf_h) / pop_h) +
    travel * S2_h * 0.75 * beta_l * ((inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I2_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_h, -1), 0) -
    I2_h * gamma -
    I2_h * delta
  dI2_h.v <-
    native * S2_h.v * 2 * beta_h * ((sym_inf_h) / pop_h) +
    travel * S2_h.v * 2  * beta_l * ((sym_inf_l) / pop_l) +
    native * S2_h.v * beta_h * ((inf_h) / pop_h) +
    travel * S2_h.v * beta_l * ((inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I2_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_h.v, -1), 0) -
    I2_h.v * gamma -
    I2_h.v * delta
  dI2_l <- 
    native * S2_l * beta_l * 2 * 0.75 * ((sym_inf_l) / pop_l) +
    travel * S2_l * beta_h * 2 *  0.75 * ((sym_inf_h)  / pop_h) +
    native * S2_l * beta_l * 0.75 * ((inf_l) / pop_l) +
    travel * S2_l * beta_h * 0.75 * ((inf_h)  / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_l, -1), 0) -
    I2_l * gamma -
    I2_l * delta
  dI2_l.v <-
    native * S2_l.v * beta_l * 2  * ((sym_inf_l) / pop_l) +
    travel * S2_l.v * beta_h * 2  * ((sym_inf_h)  / pop_h) +
    native * S2_l.v * beta_l  * ((inf_l) / pop_l) +
    travel * S2_l.v * beta_h  * ((inf_h)  / pop_h) +
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
  
  #############################
  ##THIRD INFECTION
  #############################
  dS3_h <- 
    R2_h * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_h, -1), 0) -
    native * S3_h * 2 * 0.5 * beta_h * ((sym_inf_h) / pop_h) -
    travel * S3_h * 2 * 0.5 * beta_l * ((sym_inf_l) / pop_l) -
    native * S3_h * 0.5 * beta_h * ((inf_h) / pop_h) -
    travel * S3_h * 0.5 * beta_l * ((inf_l) / pop_l) -
    S3_h * delta -
    S3_h * vac_h * sens
  dS3_h.v <-
    R2_h.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_h.v, -1), 0) -
    native * S3_h.v * 2 * 0.75 * beta_h * ((sym_inf_h) / pop_h) -
    travel * S3_h.v * 2 * 0.75 * beta_l * ((sym_inf_l) / pop_l) -
    native * S3_h.v * 0.75 * beta_h * ((inf_h) / pop_h) -
    travel * S3_h.v * 0.75 * beta_l * ((inf_l) / pop_l) -
    S3_h.v * delta 
  dS3_l <- 
    R2_l * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_l, -1), 0) -
    native * S3_l * beta_l * 2 * 0.5 * ((sym_inf_l) / pop_l ) -
    travel * S3_l * beta_h * 2 *  0.5 * ((sym_inf_h)  / pop_h ) -
    native * S3_l * beta_l * 0.5 * ((inf_l) / pop_l) -
    travel * S3_l * beta_h * 0.5 * ((inf_h)  / pop_h) -
    S3_l * delta -
    S3_l * vac_l * sens
  dS3_l.v <-
    R2_l.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_l.v, -1), 0) -
    native * S3_l.v * beta_l * 2 * 0.75 * ((sym_inf_l) / pop_l ) -
    travel * S3_l.v * beta_h * 2 *  0.75 * ((sym_inf_h)  / pop_h ) -
    native * S3_l.v * beta_l * 0.75 * ((inf_l) / pop_l) -
    travel * S3_l.v * beta_h * 0.75 * ((inf_h)  / pop_h) -
    S3_l.v * delta
  
  dI3_h <- 
    native * S3_h * 2 * 0.5 * beta_h * ((sym_inf_h) / pop_h) +
    travel * S3_h * 2 * 0.5 * beta_l * ((sym_inf_l) / pop_l) +
    native * S3_h * 0.5 * beta_h * ((inf_h) / pop_h) +
    travel * S3_h * 0.5 * beta_l * ((inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_h, -1), 0) -
    I3_h * gamma -
    I3_h * delta
  dI3_h.v <-
    native * S3_h.v * 2 * 0.75 * beta_h * ((sym_inf_h) / pop_h) +
    travel * S3_h.v * 2 * 0.75 * beta_l * ((sym_inf_l) / pop_l) +
    native * S3_h.v * 0.75 * beta_h * ((inf_h) / pop_h) +
    travel * S3_h.v * 0.75 * beta_l * ((inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_h.v, -1), 0) -
    I3_h.v * gamma -
    I3_h.v * delta
  dI3_l <- 
    native * S3_l * beta_l * 2 * 0.5 * ((sym_inf_l) / pop_l ) +
    travel * S3_l * beta_h * 2 *  0.5 * ((sym_inf_h)  / pop_h ) +
    native * S3_l * beta_l * 0.5 * ((inf_l) / pop_l) +
    travel * S3_l * beta_h * 0.5 * ((inf_h)  / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_l, -1), 0) -
    I3_l * gamma -
    I3_l * delta
  dI3_l.v <-
    native * S3_l.v * beta_l * 2 * 0.75 * ((sym_inf_l) / pop_l ) +
    travel * S3_l.v * beta_h * 2 *  0.75 * ((sym_inf_h)  / pop_h ) +
    native * S3_l.v * beta_l * 0.75 * ((inf_l) / pop_l) +
    travel * S3_l.v * beta_h * 0.75 * ((inf_h)  / pop_h) +
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
  
  #############################
  ##FOURTH INFECTION
  #############################
  dS4_h <- 
    R3_h * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S4_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_h, -1), 0) -
    native * S4_h * 2 * 0.25 * beta_h * ((sym_inf_h) / pop_h) -
    travel * S4_h * 2 * 0.25 * beta_l * ((sym_inf_l) / pop_l ) -
    native * S4_h * 0.25 * beta_h * ((inf_h) / pop_h ) -
    travel * S4_h * 0.25 * beta_l * ((inf_l) / pop_l) -
    S4_h * delta - 
    S4_h * vac_h * sens
  dS4_h.v <-
    R3_h.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S4_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_h.v, -1), 0) -
    native * S4_h.v * 2 * 0.5 * beta_h * ((sym_inf_h) / pop_h) -
    travel * S4_h.v * 2 * 0.5 * beta_l * ((sym_inf_l) / pop_l ) -
    native * S4_h.v * 0.5 * beta_h * ((inf_h) / pop_h ) -
    travel * S4_h.v * 0.5 * beta_l * ((inf_l) / pop_l) -
    S4_h.v * delta 
  dS4_l <- 
    R3_l * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_l, -1), 0) -
    native * S4_l * beta_l * 2 * 0.25 * ((sym_inf_l) / pop_l) -
    travel * S4_l * beta_h * 2 *  0.25 * ((sym_inf_h)  / pop_h) -
    native * S4_l * beta_l * 0.25 * ((inf_l) / pop_l) -
    travel * S4_l * beta_h * 0.25 * ((inf_h)  / pop_h) -
    S4_l * delta - 
    S4_l * vac_l * sens
  dS4_l.v <-
    R3_l.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S4_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_l.v, -1), 0) -
    native * S4_l.v * beta_l * 2 * 0.5 * ((sym_inf_l) / pop_l) -
    travel * S4_l.v * beta_h * 2 *  0.5 * ((sym_inf_h)  / pop_h) -
    native * S4_l.v * beta_l * 0.5 * ((inf_l) / pop_l) -
    travel * S4_l.v * beta_h * 0.5 * ((inf_h)  / pop_h) -
    S4_l.v * delta 
  
  dI4_h <- 
    native * S4_h * 2 * 0.25 * beta_h * ((sym_inf_h) / pop_h) +
    travel * S4_h * 2 * 0.25 * beta_l * ((sym_inf_l) / pop_l ) +
    native * S4_h * 0.25 * beta_h * ((inf_h) / pop_h ) +
    travel * S4_h * 0.25 * beta_l * ((inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I4_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_h, -1), 0) -
    I4_h * gamma -
    I4_h * delta
  dI4_h.v <-
    native * S4_h.v * 2 * 0.5 * beta_h * ((sym_inf_h) / pop_h) +
    travel * S4_h.v * 2 * 0.5 * beta_l * ((sym_inf_l) / pop_l ) +
    native * S4_h.v * 0.5 * beta_h * ((inf_h) / pop_h ) +
    travel * S4_h.v * 0.5 * beta_l * ((inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I4_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_h.v, -1), 0) -
    I4_h.v * gamma -
    I4_h.v * delta
  dI4_l <- 
    native * S4_l * beta_l * 2 * 0.25 * ((sym_inf_l) / pop_l) +
    travel * S4_l * beta_h * 2 *  0.25 * ((sym_inf_h)  / pop_h) +
    native * S4_l * beta_l * 0.25 * ((inf_l) / pop_l) +
    travel * S4_l * beta_h * 0.25 * ((inf_h)  / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_l, -1), 0) -
    I4_l * gamma -
    I4_l * delta
  dI4_l.v <-
    native * S4_l.v * beta_l * 2 * 0.5 * ((sym_inf_l) / pop_l) +
    travel * S4_l.v * beta_h * 2 *  0.5 * ((sym_inf_h)  / pop_h) +
    native * S4_l.v * beta_l * 0.5 * ((inf_l) / pop_l) +
    travel * S4_l.v * beta_h * 0.5 * ((inf_h)  / pop_h) +
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
    dS4_l.v, dI4_l.v, dR4_l.v))
  
}

############################
#RUN MODEL
############################


names(y_init) <- c(rep('sh1', 80), rep('ih1', 80), rep('rh1', 80),
                   rep('sh2', 80), rep('ih2', 80), rep('rh2', 80),
                   rep('sh3', 80), rep('ih3', 80), rep('rh3', 80),
                   rep('sh4', 80), rep('ih4', 80), rep('rh4', 80),
                   rep('rh1.v', 80),
                   rep('sh2.v', 80), rep('ih2.v', 80), rep('rh2.v', 80),
                   rep('sh3.v', 80), rep('ih3.v', 80), rep('rh3.v', 80),
                   rep('sh4.v', 80), rep('ih4.v', 80), rep('rh4.v', 80),
                   rep('sl1', 80), rep('il1', 80), rep('rl1', 80),
                   rep('sl2', 80), rep('il2', 80), rep('rl2', 80),
                   rep('sl3', 80), rep('il3', 80), rep('rl3', 80),
                   rep('sl4', 80), rep('il4', 80), rep('rl4', 80),
                   rep('rl1.v', 80),
                   rep('sl2.v', 80), rep('il2.v', 80), rep('rl2.v', 80),
                   rep('sl3.v', 80), rep('il3.v', 80), rep('rl3.v', 80),
                   rep('sl4.v', 80), rep('il4.v', 80), rep('rl4.v', 80))


#run intervention model
out.h <- ode(times = times, y = y_init, func = model, parms = parms.h)
#run null model 
out_null.h <- ode(times = times, y = y_init, func = model, parms = parms_null.h)

###remove the time column
out.h <- out.h[,2:ncol(out.h)]
out_null.h <- out_null.h[,2:ncol(out_null.h)]
