############################
#CRC info
############################
args = commandArgs(TRUE)
input = as.numeric(args[1])


##########################
#IF DOING VAC COVERAGE
##########################
load("parms.mat.simp_model_NO_TRAVEL.RData")



# missing.vec <- c(6,7,33,34,36,45,48,50,55,102)
# input <- missing.vec[input]



# spec <- new.parms.mat[input, 6]


validation_parms.test <- c(0.1632653, 0.2122449, 0.2938776 , 0.4163265, 0.7020408)
beta_h <- validation_parms.test[input]
beta_l <- validation_parms.test[input]
native <- 1 - new.parms.mat[input,3]
travel <- new.parms.mat[input,3]
vac_h <- new.parms.mat[input,4]
vac_l <- new.parms.mat[input,5]



load('pop_1950.RData')
load('birth_1950.RData')
load('death_1950.RData')


hopkins <- c(0.53, 1, 0.115)
hopkins_inverse <- 1 - hopkins


library(deSolve)




############################
#Initial conidtions and parameters
###########################

initial_conditions <- as.data.frame(matrix(NA, nrow = 80*4, ncol = 3))
initial_conditions[,2] <- rep(0:79,4)
for(i in 1:4){
  x <- i - 1
  initial_conditions[x*(80) + (1:80),1] <- rep(i,80)
}
initial_conditions[,3] <- rep(pop,4)

source("formulas_for_model.R")

susceptible_init_h <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 320))
susceptible_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 320))
susceptible_init_h[,3] <- fill_suscep(prop.h.1 = 0.5 , prop.h.2=0,  prop.h.3=0,  prop.h.4=0, prop.l.1 = 0.5, prop.l.2=0,  prop.l.3=0, prop.l.4=0)[[1]]
susceptible_init_l[,3] <-  fill_suscep(prop.h.1 = 0.5 , prop.h.2=0,  prop.h.3=0,  prop.h.4=0, prop.l.1 = 0.5, prop.l.2=0,  prop.l.3=0, prop.l.4=0)[[2]]

infected_init_h <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 320))
infected_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 320))
infected_init_h[,3] <- fill_inf(1,1)[[1]]
infected_init_l[,3] <- fill_inf(1,1)[[2]]


load('recovered_init_h.RData')
load('recovered_init_l.RData')


susceptible_total_h <- sum(susceptible_h(1) + susceptible_h(2) + susceptible_h(3) + susceptible_h(4))
susceptible_total_l <- sum(susceptible_l(1) + susceptible_l(2) + susceptible_l(3) + susceptible_l(4))


infected_total_h <- sum(infected_h(1) + infected_h(2) + infected_h(3) + infected_h(4))
infected_total_l <- sum(infected_l(1) + infected_l(2) + infected_l(3) + infected_l(4))

 
recovered_total_h <- sum(recovered_h(1) + recovered_h(2) + recovered_h(3) + recovered_h(4))
recovered_total_l <- sum(recovered_l(1) + recovered_l(2) + recovered_l(3) + recovered_l(4))

population_h <- sum(susceptible_total_h + infected_total_h + recovered_total_h)
population_l <- sum(susceptible_total_l + infected_total_l + recovered_total_l)



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
                sens = 0.85,
                spec = 0.95,
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
                     sens = 0.85,
                     spec = 0.95,
                     hopkins,
                     hopkins_inverse)



years = 50
years_vac = 30
times <- seq(from = 0, to = 365 * years, by = 0.1)
times <- times[1:(length(times) - 1)]

############################
#MODEL
############################
model <- function(t, y, parms){
  
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
  travel <- parms[[9]]
  t_vac.h <- ifelse((t>(365*(years_vac))), parms[[10]] / 365, 0)
  t_vac.l <- ifelse((t>(365*(years_vac))), parms[[11]] / 365, 0)
  vac_h <- c(rep(0,8), t_vac.h, rep(0,71))
  vac_l <- c(rep(0,8), t_vac.l, rep(0,71))
  sens <- parms[[12]]
  spec <- parms[[13]]
  inf <- parms[[14]]
  ninf <- parms[[15]]
  
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
    native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) -
    travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) -
    travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    #DEATH
    S1_h * delta -
    #VACCINATION
    S1_h * vac_h * (1 - spec)
  dS1_l <-  
    mu * c(birth_pop_l, rep(0,79)) +
    c(0, 1/365 / head(age_window, -1) * head(S1_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S1_l, -1), 0)  -
    native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) -
    native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l) -
    S1_l * delta - 
    S1_l * vac_l * (1 - spec)
  
  dI1_h <- 
    native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I1_h, -1), 0) -
    #RECOVERY
    I1_h * gamma -
    I1_h * delta
  dI1_l <- 
    native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
    native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l) +
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
    native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) -
    travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) -
    travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    S2_h * delta - 
    S2_h * vac_h * sens
  dS2_h.v <-
    R1_h.v * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_h.v, -1), 0) -
    native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) -
    travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) -
    travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    S2_h.v * delta
  dS2_l <- 
    R1_l * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_l, -1), 0) -
    native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) -
    travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) -
    native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) -
    travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) -
    S2_l * delta - 
    S2_l * vac_l * sens
  dS2_l.v <-
    R1_l.v * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_l.v, -1), 0) -
    native * S2_l.v * beta_l * 2 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) -
    travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) -
    native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) -
    travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) -
    S2_l.v * delta
  
  dI2_h <- 
    native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I2_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_h, -1), 0) -
    I2_h * gamma -
    I2_h * delta
  dI2_h.v <-
    native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I2_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_h.v, -1), 0) -
    I2_h.v * gamma -
    I2_h.v * delta
  dI2_l <- 
    native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_l, -1), 0) -
    I2_l * gamma -
    I2_l * delta
  dI2_l.v <-
    native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
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
    native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) -
    travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) -
    travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    S3_h * delta -
    S3_h * vac_h * sens
  dS3_h.v <-
    R2_h.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_h.v, -1), 0) -
    native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) -
    travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) -
    travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    S3_h.v * delta 
  dS3_l <- 
    R2_l * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_l, -1), 0) -
    native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) -
    travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) -
    native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) -
    travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) -
    S3_l * delta -
    S3_l * vac_l * sens
  dS3_l.v <-
    R2_l.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_l.v, -1), 0) -
    native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) -
    travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) -
    native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) -
    travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) -
    S3_l.v * delta
  
  dI3_h <- 
    native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_h, -1), 0) -
    I3_h * gamma -
    I3_h * delta
  dI3_h.v <-
    native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_h.v, -1), 0) -
    I3_h.v * gamma -
    I3_h.v * delta
  dI3_l <- 
    native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_l, -1), 0) -
    I3_l * gamma -
    I3_l * delta
  dI3_l.v <-
    native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
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
    native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) -
    travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) -
    travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    S4_h * delta - 
    S4_h * vac_h * sens
  dS4_h.v <-
    R3_h.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S4_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_h.v, -1), 0) -
    native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) -
    travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) -
    travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    S4_h.v * delta 
  dS4_l <- 
    R3_l * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_l, -1), 0) -
    native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) -
    travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) -
    native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) -
    travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) -
    S4_l * delta - 
    S4_l * vac_l * sens
  dS4_l.v <-
    R3_l.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S4_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_l.v, -1), 0) -
    native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) -
    travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) -
    native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) -
    travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) -
    S4_l.v * delta 
  
  dI4_h <- 
    native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I4_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_h, -1), 0) -
    I4_h * gamma -
    I4_h * delta
  dI4_h.v <-
    native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I4_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_h.v, -1), 0) -
    I4_h.v * gamma -
    I4_h.v * delta
  dI4_l <- 
    native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_l, -1), 0) -
    I4_l * gamma -
    I4_l * delta
  dI4_l.v <-
    native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
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
  
 
  
  #############################
  ##CUMULATIVE PRIMARY CASES
  #############################
  { 
    
    primary_cases.l <-
      native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l) 
    primary_cases.l.v <- sum(primary_cases.l[9]) * inf[1]
    primary_cases.l <- sum(primary_cases.l) * inf[1]
    
    
    primary_cases.h <- 
      native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    primary_cases.h.v <- sum(primary_cases.h[9]) * inf[1]
    primary_cases.h <- sum(primary_cases.h) * inf[1]
  }

  #############################
  ##CUMULATIVE SECONDARY CASES
  #############################
  {  
    secondary_vac_cases.l <- 
      native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l)  
    secondary_vac_cases.l <- sum(secondary_vac_cases.l) * inf[2]
    
    
    secondary_vac_cases.h <- 
      native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h)
    secondary_vac_cases.h <- sum(secondary_vac_cases.h) * inf[2] }

  #############################
  ##CUMULATIVE SECONDARY CASES, UNVACCINATED
  #############################
  {  
    
    secondary_cases.l <-
      travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) 
    secondary_cases.l.v <- sum(secondary_cases.l[9]) * inf[2]
     secondary_cases.l <- sum(secondary_cases.l) * inf[2]
    

    
    secondary_cases.h <-
      native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    secondary_cases.h.v <- sum(secondary_cases.h[9]) * inf[2]
    secondary_cases.h <- sum(secondary_cases.h) * inf[2]
    
  }
  
  
  #############################
  ##CUMULATIVE POST-SECONDARY CASES, VACCINATED
  #############################
{
    postsec_vac_cases.l <- 
      travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h)
    postsec_vac_cases.l <- sum(postsec_vac_cases.l) * inf[3]
  
    
    postsec_vac_cases.h <- 
      native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    postsec_vac_cases.h <- sum(postsec_vac_cases.h) * inf[3]
}
  
  #############################
  ##CUMULATIVE POST-SECONDARY CASES, UNVACCINATED
  #############################
{    postsec_cases.l <- 
      travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) 
    postsec_cases.l.v <- sum(postsec_cases.l[9]) * inf[3]
    postsec_cases.l <- sum(postsec_cases.l) * inf[3]
    
    
    postsec_cases.h <- 
      native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    postsec_cases.h.v <- sum(postsec_cases.h[9]) * inf[3]
    postsec_cases.h <- sum(postsec_cases.h) * inf[3]
 
  }
  
 
  #############################
  ##FOI
  ############################# 
{  FOI_h.travel <-
    sum(native * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
          travel * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
          native * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
          travel * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h))
  
  FOI_l.travel <-
    sum(native * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
          travel * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
          native * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
          travel * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l))
  }
  
  
  #############################
  ##Vaccination coverages
  ############################# 
  ####inappropriate vaccinations
  vac_1.h <-  sum(S1_h * vac_h * (1 - spec)) 
  vac_1.l <-  sum(S1_l * vac_l * (1 - spec)) 
  
  pop_1.h.9 <- sum(S1_h[9])
  pop_1.l.9 <- sum(S1_l[9])
  
  
  ##appropriate vaccinations, averts secondary infections
  vac_2.h <- sum(R1_h * vac_h * sens)  + sum(S2_h * vac_h * sens) 
  vac_2.l <- sum(R1_l * vac_l * sens )  + sum(S2_l * vac_l * sens) 
  
  pop_2.h.9 <- sum(S2_h[9] + R1_h[9])
  pop_2.l.9 <- sum(S2_l[9] + R1_l[9])
  
  ####averts tertiary infects
  vac_3.h <- sum(R2_h * vac_h * sens) + sum(S3_h * vac_h * sens) 
  vac_3.l <- sum(R2_l * vac_l * sens) + sum(S3_l * vac_l * sens) 
  
  pop_3.h.9 <- sum(S3_h[9] + R2_h[9])
  pop_3.l.9 <- sum(S3_l[9] + R2_l[9])
  
  ##### averts quat infects
  vac_4.h <- sum(R3_h * vac_h * sens) + sum(S4_h * vac_h * sens)
  vac_4.l <- sum(R3_l * vac_l * sens) + sum(S4_l * vac_l * sens) 
  
  pop_4.h.9 <- sum(S4_h[9] + R3_h[9])
  pop_4.l.9 <- sum(S4_l[9] + R3_l[9])

  
  
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
    

    
    primary_cases.l, 
    primary_cases.h, 
    
    secondary_vac_cases.l,
    secondary_vac_cases.h, 
    secondary_cases.l, 
    secondary_cases.h, 
    
    postsec_vac_cases.l, 
    postsec_vac_cases.h, 
    postsec_cases.l,
    postsec_cases.h, 
    

    FOI_h.travel, FOI_l.travel,
    
    birth_pop_h, birth_pop_l,
    
    primary_cases.l.v,  primary_cases.h.v,
    secondary_cases.l.v,  secondary_cases.h.v,
    postsec_cases.l.v, postsec_cases.h.v,
    
    vac_1.h,  vac_1.l,
    vac_2.h,  vac_2.l,
    vac_3.h,  vac_3.l,
    vac_4.h,  vac_4.l,
    
    pop_1.h.9, pop_1.l.9,
    pop_2.h.9, pop_2.l.9,
    pop_3.h.9, pop_3.l.9,
    pop_4.h.9, pop_4.l.9
    
  ))
  
}

############################
#RUN MODEL
############################
y_init <- c(susceptible_h(1), infected_h(1), recovered_h(1),
            susceptible_h(2), infected_h(2), recovered_h(2),
            susceptible_h(3), infected_h(3), recovered_h(3),
            susceptible_h(4), infected_h(4), recovered_h(4),
            rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80),
            rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80),
            susceptible_l(1), infected_l(1), recovered_l(1),
            susceptible_l(2), infected_l(2), recovered_l(2),
            susceptible_l(3), infected_l(3), recovered_l(3),
            susceptible_l(4), infected_l(4), recovered_l(4),
            rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80),
            rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80),
            rep(0, 12), rep(0,2), rep(0,6), 
            rep(0, 8), rep(0, 8))
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
                   rep('sl4.v', 80), rep('il4.v', 80), rep('rl4.v', 80),

                   'prim_cases.l', 
                   'prim_cases.h', 
                   
                   'sec_vac.cases.l', 
                   'sec_vac.cases.h', 
                   'sec_cases.l',
                   'sec_cases.h', 
                   
                   'psec_vac.cases.l', 
                   'psec_vac.cases.h', 
                   'psec_cases.l', 
                   'psec_cases.h',
                   

  
                   'FOI_h.travel', 'FOI_l.travel',
                   'pop_h', 'pop_l',
                   'vac_eleg_p.l', 'vac_eleg_p.h',
                   'vac_eleg_s.l', 'vac_eleg_s.h',
                   'vac_elege_ps.l', 'vac_eleg.ps.h',
                   
                   rep('vac_1.h', 1),  rep('vac_1.l', 1),
                   rep('vac_2.h', 1), rep('vac_2.l', 1),
                   rep('vac_3.h', 1),  rep('vac_3.l', 1),
                   rep('vac_4.h', 1),  rep('vac_4.l', 1),
                   
                   rep('pop_1.h', 1),  rep('pop_1.l', 1),
                   rep('pop_2.h', 1), rep('pop_2.l', 1),
                   rep('pop_3.h', 1),  rep('pop_3.l', 1),
                   rep('pop_4.h', 1),  rep('pop_4.l', 1)
                   
)


#run intervention model
out.h <- ode(times = times, y = y_init, func = model, parms = parms.h)
#run null model 
out_null.h <- ode(times = times, y = y_init, func = model, parms = parms_null.h)

###remove the time column
out.h <- out.h[,2:ncol(out.h)]
out_null.h <- out_null.h[,2:ncol(out_null.h)]



###################
#CASES AVERTED
###################

  
  ####time series coverage
  {
    timepoint_year <- years
    index <- (timepoint_year * 3650 )

    vac_h.1 <- out.h[index,which(colnames(out.h) == 'vac_1.h')] / out.h[index,which(colnames(out.h) == 'pop_1.h')]
    vac_h.2 <- out.h[index,which(colnames(out.h) == 'vac_2.h')] / out.h[index,which(colnames(out.h) == 'pop_2.h')]
    vac_h.3 <- sum(out.h[index,which(colnames(out.h) == 'vac_3.h')] + out.h[index,which(colnames(out.h) == 'vac_4.h')]) /
      sum(out.h[index,which(colnames(out.h) == 'pop_3.h')] + out.h[index,which(colnames(out.h) == 'pop_4.h')])

    cov_h <- sum(out.h[index,which(colnames(out.h) == 'vac_1.h')] +
                   out.h[index,which(colnames(out.h) == 'vac_2.h')] +
                   out.h[index,which(colnames(out.h) == 'vac_3.h')] +
                   out.h[index,which(colnames(out.h) == 'vac_4.h')]) / sum(out.h[index,which(colnames(out.h) == 'pop_1.h')] +
                                                                             out.h[index,which(colnames(out.h) == 'pop_2.h')] +
                                                                             out.h[index,which(colnames(out.h) == 'pop_3.h')] +
                                                                             out.h[index,which(colnames(out.h) == 'pop_4.h')])
    coverage_h <- c(vac_h.1, vac_h.2, vac_h.3)



    vac_l.1 <- out.h[index,which(colnames(out.h) == 'vac_1.l')] / out.h[index,which(colnames(out.h) == 'pop_1.l')]
    vac_l.2 <- out.h[index,which(colnames(out.h) == 'vac_2.l')] / out.h[index,which(colnames(out.h) == 'pop_2.l')]
    vac_l.3 <- sum(out.h[index,which(colnames(out.h) == 'vac_3.l')] + out.h[index,which(colnames(out.h) == 'vac_4.l')]) /
      sum(out.h[index,which(colnames(out.h) == 'pop_3.l')] + out.h[index,which(colnames(out.h) == 'pop_4.l')])

    cov_l <- sum(out.h[index,which(colnames(out.h) == 'vac_1.l')] +
                   out.h[index,which(colnames(out.h) == 'vac_2.l')] +
                   out.h[index,which(colnames(out.h) == 'vac_3.l')] +
                   out.h[index,which(colnames(out.h) == 'vac_4.l')]) / sum(out.h[index,which(colnames(out.h) == 'pop_1.l')] +
                                                                             out.h[index,which(colnames(out.h) == 'pop_2.l')] +
                                                                             out.h[index,which(colnames(out.h) == 'pop_3.l')] +
                                                                             out.h[index,which(colnames(out.h) == 'pop_4.l')])
    coverage_l <- c(vac_l.1, vac_l.2, vac_l.3)

    coverage <- c(cov_h, cov_h)
    names(coverage) <- c('h', 'l')
    save(coverage, file = paste('new.cov_', input, '.RData', sep = ''))

    cases.output.vec.h  <- cases_averted.func(out_mat = out.h, out_mat_null = out_null.h, timepoint_year = years, cases = 1)
#   # infections.output.vec.h  <- cases_averted.func(out_mat = out.h, out_mat_null = out_null.h, timepoint_year = years, cases = 0)
#   
#   
#   
save(cases.output.vec.h, file = paste('cases_averted_', input, '.RData', sep = ''))
#   # save(infections.output.vec.h, file = paste('infections_averted_', input, '.RData', sep = ''))
#   
 }



# last_row <- out.h[nrow(out.h),]
# save(last_row, file = paste('last_row_', input, '.RData', sep =''))

######check seroprevalence levels

# sp9.vec <- seroprevalence_fun(out_mat = out.h, age = 9)
# save(sp9.vec, file = paste('sp9_', input, '.RData', sep = ''))

######check FOI
# foi_h <- FOI_h.fun(years = years)
# foi_l <- FOI_l.fun(years = years)
# foi <- c(foi_h, foi_l)
# names(foi) <- c('h', 'l')
# save(foi, file = paste('foi_', input, '.RData', sep = ''))
# 
# 
