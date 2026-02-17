library(EpiEstim)
library(nimble)
library(bayesplot)
library(coda)
library(viridis)
library(dplyr)
library(gridExtra)
library(scales)
library(socialmixr)
library(rstan)
library(deSolve)
library(reshape2)
library(ggplot2)
library(shinySIR)
library(rootSolve)
library(R0) 

## Inserting Data ##

data <- read.csv("https://opendata.ecdc.europa.eu/covid19/nationalcasedeath_eueea_daily_ei/csv", na.strings = "", fileEncoding = "UTF-8-BOM")

## Taking the subset of Czech Republic ##

cz=data.frame(data[data$countriesAndTerritories=="Czechia",])
N=cz[1,10]  #population of Czech Republic
## manipulate Data to the desired subset ##
cz$dateRep=as.Date(cz$dateRep,format = "%d/%m/%Y")
cz=cz[,-c(2,3,4,7,8,9,10,11)]
cz=cz[-c(1:954),]
cz=cz[-c(49:16),]
cz$day=seq(nrow(cz),1)
cz=cz[order(cz$day),]

## 1st Question ##

# plot with gradient bars + connecting line 
ggplot(cz, aes(x = dateRep, y = cases)) +
  # gradient bars
  geom_col(aes(fill = cases), color = "white", width = 0.8, show.legend = FALSE) +
  scale_fill_viridis_c(option = "C", direction = -1, end = 0.9) +
  
  # connect the tops of the bars with a line
  geom_line(aes(group = 1), color = "#2A3F5F", size = 1) +
  
  
  
  # custom x/y scales
  scale_x_date(date_breaks = "1 week", date_labels = "%b %d, %Y",
               expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(labels = label_comma(),
                     expand = expansion(mult = c(0, 0.05))) +
  
  # labels
  labs(
    title    = "Epidemic Curve of Cases Over 15 days",
    subtitle = "Daily counts for the first 15 days of the pandemic in Czech Republic",
    x        = "Date of Onset",
    y        = "Number of Cases"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 12, color = "grey40", hjust = 0.5),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12, color = "grey20"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "#F7F7F7", color = NA)
  )

## 2nd Question ##

#we assume serial interval is 4 days based on early day covid-19 researches
#Link available on Word 
t_start <- seq(2, nrow(cz) - 6)

t_end <- t_start + 6

Rt <- estimate_R(incid = cz$cases,
                 method = "parametric_si",
                 config = make_config(list(
                   mean_si = 4, 
                   std_si = 2.8, 
                   t_start = t_start,
                   t_end = t_end)))

plot(Rt, legend = FALSE) 
plot(Rt, "R")

## 3rd Question ##

## 3a ## 

initial_state_values <- c(S = N-10, # the whole population is susceptible to infection, except for 10 people
                          E = 10,            # infected but not infectious
                          I = 0,          # the epidemic starts with 0 infected people
                          R = 0,            # no prior immunity
                          C = 0)          # cumulative number of infections

R0=round(Rt$R$`Mean(R)`[1],2)

gamma <- 1/7
beta <- R0 * gamma

parameters <- c(beta = beta,   
                sigma = 1/4,  
                gamma = gamma)

times <- seq(from = 0, to = 120) # 4 months/ 120 days 

seir_mod = function(time, state, parms) {
  
  with(as.list(c(state, parms)), {
    
    N <- S + E + I + R   # total population size
    lambda <- beta * I/N # define lambda
    
    # define differential equations
    dS <- - lambda * S
    dE <- lambda * S - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I 
    dC <- sigma * E             # keep track of cumulative incidence (all those who enter Ι)
    res = c(dS, dE, dI, dR, dC) # in the same order as the input state variables
    return(list(res))
  })
}

out <- ode(y = initial_state_values,
           times = times, 
           func = seir_mod, 
           parms = parameters)
out <- as.data.frame(out)



## 3b ##

ggplot(out, aes(x = time, y = E)) +
  geom_line(size = 1.2, color = "#D55E00") +
  labs(
    title = " Number of Exposed (E) over First 4 Months",
    x     = "Time (days since start)",
    y     = "E (currently exposed)"
  ) +
  theme_minimal(base_size = 14)

## 3c ##
  
  out$inf_per <- out$C/N*100
  inf_per<- out$C[nrow(out)] / N * 100
 
  cat ("The percentage of the population that got infected in the first 4 months of the pandemic was ",round(inf_per,2),"%")

  # Cumulative incidence
ggplot(out , aes(x = time, y = C)) +
  geom_line(size = 1.2, color = "#0072B2") +
  labs(
    title = "Cumulative Number of Infections (C)",
    x     = "Time (days since start)",
    y     = "Cumulative incidence"
  ) +
  theme_minimal(base_size = 14)

# % of population infected
ggplot(out, aes(x = time, y = inf_per)) +
  geom_line(size = 1.2, color = "#009E73") +
  labs(
    title = "% of Population Infected over First 4 Months",
    x     = "Time (days since start)",
    y     = "Percent infected (%)"
  ) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  theme_minimal(base_size = 14)

## 3d ## 

equil <- runsteady(y = c(S = N-10, E = 10,I=0, R = 0,C=0),  
                   times = c(0, 1E5), func = seir_mod, parms = parameters)

comma(round(equil$y, digits = 0))


## 4th Question ##

parametersSD <- c(beta = beta,   # infection rate (/day)
                  sigma = 1/4,   # rate from E to I = 1/duration of latent phase (/day)
                  gamma = gamma, # recovery rate = 1/duration of infectiousness (/day) 
                  f = 0.4,      # reduction in beta due to social distancing
                  tSD = 60)      # time at which social distancing is implemented

# sequence of timesteps to solve the model at 
times <- seq(from = 0, to = 120) # from 0 to 120 days with 10 time-increments per day

# SΕIR with social distancing model function
# takes as input arguments (in the following order): time, state and parameters
seir_SD_mod = function(time, state, parms) {
  
  with(as.list(c(state, parms)), {
    
    betaSD <- beta*(1-f)     # reduced infection rate
    N <- S + E + I + R       # total population size
    lambda <- beta * I/N     # define lambda with no measures
    lambdaSD <- betaSD * I/N # define lambda with social distancing 
    
    # define differential equations
    if (time <= tSD){  # if time <= tSD days: no measures
      dS <- - lambda * S
      dE <- lambda * S - sigma * E
      dI <- sigma * E - gamma * I
      dR <- gamma * I 
      dC <- sigma * E
    }
    else{              # if time > tSD days: social distancing
      dS = - lambdaSD * S
      dE = lambdaSD * S - sigma * E
      dI = sigma * E -  gamma * I
      dR = gamma * I
      dC = sigma * E
    }
    res = c(dS, dE, dI, dR, dC) 
    return(list(res))
  })
}

# solve the differential equations using the ode integration algorithm 
outSD <- ode(y = initial_state_values,
             times = times, 
             func = seir_SD_mod, 
             parms = parametersSD)
outSD <- as.data.frame(outSD)


ggplot() +
  # Baseline: solid line
  geom_line(
    data    = out,
    aes(x = times, y = E ,color = "Baseline"),
    size    = 1.2
  ) +
  # Social distancing: dashed line
  geom_line(
    data    = outSD,
    aes(x = times, y = E, color = "Social distancing"),
    size    = 1.2,
    linetype = "dashed"
  ) +
  
  # Manual color scale for consistency & accessibility
  scale_color_manual(
    name   = "Scenario",
    values = c(
      "Baseline"           = "#1f78b4",
      "Social distancing"  = "#33a02c"
    )
  ) +
  
  # Axis labels and title
  labs(
    title = "Daily Infected: Baseline vs Social Distancing",
    x     = "Time (days)",
    y     = "Daily number of new infections"
  ) +
  
  # Minimalist theme with legend at bottom
  theme_minimal(base_size = 14) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold"),
    legend.text        = element_text(size = 12),
    plot.title         = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text          = element_text(color = "grey20"),
    panel.grid.major   = element_line(color = "grey90"),
    panel.grid.minor   = element_blank()
  )

## 5th Question ##

## 5a ##


# where adults are adults between 18-64 years old 
# where seniors are adults 65+ years old

Cmat <- matrix(
  c(18,  5,  1,   # children’ contacts:   [with child, adult, senior]
    4,   6,  2,   # adults’ contacts:     [with child, adult, senior]
    1,   2,  2),  # seniors’ contacts:    [with child, adult, senior]
  nrow = 3, byrow = TRUE,
  dimnames = list(
    c("child","adult","senior"),
    c("child","adult","senior")
  )
)
matrix_plot(Cmat)


## 5b ##

initial_state_values <- c(S1 = N*0.2-5,    # the whole subpopulation of children is susceptible to infection, except for 5 child
                          E1 = 5,          # the epidemic starts with 1 infected child
                          I1 = 0,       
                          R1 = 0,          # no prior immunity
                          S2 = N*0.5-5,    # the whole subpopulation of adults is susceptible to infection, except for 5 adult 
                          E2 = 5,          # the epidemic starts with 5 infected  adults
                          I2 = 0,
                          R2 = 0,        # no prior immunity
                          S3= N*0.3,       # the whole subpopulation of seniors, is susceptible to infection
                          E3= 0,           #the epidemic starts with 0 infected seniors
                          I3= 0,
                          R3= 0)          # no prior immunity  

# parameters describing the transition rates
parameters <- c(p = 0.05,    # the probability of infection per contact is 5%
                c11 = 18,    # daily number of contacts that children make with each other
                c12 = 5,     # daily number of contacts that children make with adults 
                c13 = 1,     # daily number of contacts that children make with seniors
                c21 = 4,     # daily number of contacts that adults make with children
                c22 = 6,     # daily number of contacts that adults make with each other
                c23 = 2,     # daily number of contacts that adults make with seniors 
                c31 = 1,     # daiy number of contacts that seniors make with children
                c32 = 2,     # daily number of contacts that seniors make with adults
                c33 = 2,     # daily number of contacts that seniors make with each other
                sigma = 1/4, # rate from E to I = 1/duration of latent phase (/day)
                gamma = 1/7) # recovery rate = 1/duration of infectiousness (/day)

# sequence of timesteps to solve the model at 
times <- seq(from = 0, to = 120, by = 1) # from 0 to 120 days (4 months) in daily intervals

# age-structured SEIR model function
seir_mod_2grps = function(time, state, parms) {
  
  with(as.list(c(state, parms)), {
    
    N1 <- S1 + E1 + I1 + R1                       # total number of children 
    N2 <- S2 + E2 + I2 + R2                       # total number of 18-64 years old adults 
    N3 =  S3 + E3 + I3 + R3                       # total number of 65+ years old adults 
    
    lambda1 <- p * c11 * I1/N1 + p * c12 * I2/N2 + p * c13 * I3/N3 # force of infection for children
    lambda2 <- p * c21 * I1/N1 + p * c22 * I2/N2 + p * c23 * I3/N3  # force of infection for 18-64 years old adults
    lambda3 =  p * c31 * I1/N1 + p * c32 * I2/N2 + p * c33 * I3/N3 # force of infection for 65+ years old adults
    
    # define differential equations
    # Children
    dS1 <- - lambda1 * S1
    dE1 <- lambda1 * S1 - sigma * E1 
    dI1 <- sigma * E1 - gamma * I1 
    dR1 <- gamma * I1
    # Adults
    dS2 <- - lambda2 * S2
    dE2 <- lambda2 * S2 - sigma * E2
    dI2 <- sigma * E2 - gamma * I2
    dR2 <- gamma * I2
    # Seniors
    dS3 <- - lambda3 * S3
    dE3 <- lambda3 * S3 - sigma * E3
    dI3 <- sigma * E3 - gamma * I3
    dR3 <- gamma * I3
    res <- c(dS1, dE1, dI1, dR1, dS2, dE2, dI2, dR2,dS3,dE3,dI3,dR3) 
    return(list(res))
  })
}

# solve the differential equations 
out <- ode(y = initial_state_values,
           times = times, 
           func = seir_mod_2grps, 
           parms = parameters)
out <- as.data.frame(out)


ggplot(out, aes(x = time)) +
  geom_line(aes(y = I1, color = "Children"),       size = 1.2) +
  geom_line(aes(y = I2, color = "Adults"), size = 1.2, ) +
  geom_line(aes(y = I3, color = "Seniors"),  size = 1.2, ) +
  
  scale_color_manual(
    name   = "Age group",
    values = c(
      "Children" = "#1f78b4",
      "Adults"  = "#33a02c",
      "Seniors"   = "#e31a1c" )) +
  
  scale_x_continuous(
    breaks = seq(0, max(out$time), by = 14),
    expand = expansion(mult = c(0, 0))) +
  
  labs(
    title = "SEIR Model: Infectious over Time by Age Group",
    x     = "Day",
    y     = "Number Infectious") +
  theme_bw(base_size = 14) +
  theme(
    legend.position     = "bottom",
    legend.key.width    = unit(2, "line"),
    plot.title          = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x         = element_text(angle = 45, hjust = 1),
    panel.grid.major    = element_line(color = "grey90"),
    panel.grid.minor    = element_blank())


## 6th Question ## 



cz1=data.frame(data[data$countriesAndTerritories=="Czechia",])

cz1$dateRep=as.Date(cz1$dateRep,"%d/%m/%Y")

cz1=cz1[cz1$dateRep>="2020-09-01" & cz1$dateRep<= "2021-06-20",]

N=cz1[1,10] 

cz1=cz1[,-c(2,3,4,7,8,9,10,11)]

cz1$day=seq(nrow(cz1),1)

cz1=cz1[order(cz1$day),]
rownames(cz1)=cz1$day




cases    <- cz1$cases
lag_inf  <- 7
removals <- c(rep(0, lag_inf), head(cases, -lag_inf))

# 2) Bundle into NIMBLE lists
n_obs       <- length(cases)

data_list  <- list(
  new_cases    = cases,
  new_removals = removals
)
const_list <- list(
  n_obs = n_obs,
  n_pop = N,
  I0    = 1            # assume 1 initial infected
)

# 3) Define the chain‐binomial S→I model
sir_code <- nimbleCode({
  # prior on transmission rate
  beta ~ dlnorm(0, 5)
  
  # compute initial susceptibles
  S0 <- n_pop - I0
  
  # time t = 1
  p[1]         <- 1 - exp(-beta * I0 / n_pop)
  new_cases[1] ~ dbin(p[1], S0)
  S[1]         <- S0 - new_cases[1]
  I[1]         <- I0 + new_cases[1] - new_removals[1]
  
  # times t = 2 ... n_obs
  for(t in 2:n_obs) {
    p[t]         <- 1 - exp(-beta * I[t-1] / n_pop)
    new_cases[t] ~ dbin(p[t],    S[t-1])
    S[t]         <- S[t-1] - new_cases[t]
    I[t]         <- I[t-1] + new_cases[t] - new_removals[t]
  }
})

# 4) Build & compile the model
inits_list <- list(beta = runif(1, 0.01, 0.1))

nim_mod   <- nimbleModel(
  code      = sir_code,
  data      = data_list,
  constants = const_list,
  inits     = inits_list
)
c_nim_mod <- compileNimble(nim_mod)

# 5) Configure, build & compile MCMC
mcmc_conf  <- configureMCMC(nim_mod, monitors = "beta")
mcmc_run   <- buildMCMC(mcmc_conf)
c_mcmc_run <- compileNimble(mcmc_run, project = nim_mod)

# 6) Run MCMC
set.seed(1234)
samples    <- runMCMC(
  c_mcmc_run,
  niter             = 15000,
  nburnin           = 3000,
  nchains           = 3,
  samplesAsCodaMCMC = TRUE
)

# 7) Summaries & diagnostics
print(summary(samples), digits = 3)

mcmc_trace(samples, pars = "beta") + ggtitle("Trace Plot: β")
mcmc_dens( samples, pars = "beta") + ggtitle("Density Plot: β")
mcmc_acf(  samples, pars = "beta") + ggtitle("ACF Plot: β")



## 7th Question ##

ggplot(cz1, aes(x = dateRep, y = cases)) +
  geom_col(fill = "#0072B2", color = "white", width = 1) +
  labs(
    title = "Current Infected Individuals Over Time",
    x     = "Date",
    y     = "Number Infectious"
  ) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%d/%m/%y",
    expand      = expansion(add = c(0, 0))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

idx <- cz1$cases > 15000
cz1[idx, c("dateRep","cases","day")]


# 0) Pick your single change‐point (day index)
tau1_fixed <- 40

# 1) Bundle data & constants
data_list <- list(
  new_cases = cases ,
  new_removals =removals 
)
const_list <- list(
  n_obs = n_obs,
  n_pop = N,
  I0    = 1,
  tau1  = tau1_fixed
)

# 2) Piecewise‐beta chain‐binomial in nimbleCode
singleCP_code <- nimbleCode({
  # priors
  beta1 ~ dlnorm(0, 5)
  beta2 ~ dlnorm(0, 5)
  
  # initial states
  S0 <- n_pop - 1
  I0 <- 1
  
  # day t=1
  ind1[1]    <- step(tau1 - 1 + 0.5)     # 1 if 1 ≤ tau1
  beta_t[1]  <- beta1*ind1[1] + beta2*(1-ind1[1])
  p[1]       <- 1 - exp(-beta_t[1] * I0 / n_pop)
  new_cases[1]    ~ dbin(p[1], S0)
  S[1]       <- S0 - new_cases[1]
  I[1]       <- I0 + new_cases[1] - new_removals[1]
  
  # days t=2…n_obs
  for(t in 2:n_obs) {
    ind1[t]   <- step(tau1 - t + 0.5)
    beta_t[t] <- beta1*ind1[t] + beta2*(1-ind1[t])
    p[t]      <- 1 - exp(-beta_t[t] * I[t-1] / n_pop)
    
    new_cases[t]    ~ dbin(p[t], S[t-1])
    S[t]            <- S[t-1] - new_cases[t]
    I[t]            <- I[t-1] + new_cases[t] - new_removals[t]
  }
})

# 3) Compile & run MCMC
inits1 <- function() list(beta1 = runif(1,0.01,0.1),
                          beta2 = runif(1,0.01,0.1))

mod1   <- nimbleModel(singleCP_code,
                      data      = data_list,
                      constants = const_list,
                      inits     = inits1())
cmod1  <- compileNimble(mod1)

conf1  <- configureMCMC(mod1, monitors = c("beta1","beta2"))
rmcmc1 <- buildMCMC(conf1)
cmcmc1 <- compileNimble(rmcmc1, project = mod1)

set.seed(2025)
samples1 <- runMCMC(cmcmc1,
                    niter             = 10000,
                    nburnin           = 2000,
                    nchains           = 3,
                    samplesAsCodaMCMC = TRUE)

# --- 4) summarize & plot -------------------------------
print(summary(samples1), digits = 3)

# trace & density
mcmc_trace(samples1, pars = c("beta1","beta2")) +
  ggtitle("Trace Plots for b1 (pre-τ) & b2 (post-τ)")

mcmc_dens(samples1, pars = c("beta1","beta2")) +
  ggtitle("Posterior Densities for b1 (pre-τ) & b2 (post-τ)")

mcmc_acf(samples1, pars = c("beta1","beta2")) + ggtitle("ACF: β")






## Question 8 ##





idx <- cz1$cases<5000 & cz1$dateRep > "2020-11-04"
cz1[idx, c("dateRep","cases","day")]

# --- 0) specify your fixed change‐points ----------
tau1_fixed <- 40# first switch at day 40
tau2_fixed <- 94   # second switch at day 94
tau3_fixed <- 130  # third switch at day 130
tau4_fixed <- 174  # fourth switch at day 174
tau5_fixed <- 196 # fifth switch at day 196

# --- 1) bundle data & constants -------------------
data_list <- list(
  new_cases    = cz1$cases,
  new_removals = cz1$deaths
)
const_list <- list(
  n_obs = length(cz1$cases),
  n_pop = N,
  tau1  = tau1_fixed,
  tau2  = tau2_fixed,
  tau3  = tau3_fixed,
  tau4  = tau4_fixed,
  tau5  = tau5_fixed
)

# --- 2) define the multi‐phase chain‐binomial model ---
multiTauCode <- nimbleCode({
  # priors for each β in its regime
  beta1 ~ dlnorm(0, 5)  # t ≤ tau1
  beta2 ~ dlnorm(0, 5)  # tau1 < t ≤ tau2
  beta3 ~ dlnorm(0, 5)  # tau2 < t ≤ tau3
  beta4 ~ dlnorm(0, 5)  # tau3 < t ≤ tau4
  beta5 ~ dlnorm(0, 5)  # tau4 < t ≤ tau5
  beta6 ~ dlnorm(0, 5)  # t >  tau5
  
  # initial conditions
  S0 <- n_pop - 1
  I0 <- 1
  
  # day 1
  ind1[1] <- step(tau1 - 1 + 0.5)
  ind2[1] <- step(tau2 - 1 + 0.5) - ind1[1]
  ind3[1] <- step(tau3 - 1 + 0.5) - step(tau2 - 1 + 0.5)
  ind4[1] <- step(tau4 - 1 + 0.5) - step(tau3 - 1 + 0.5)
  ind5[1] <- step(tau5 - 1 + 0.5) - step(tau4 - 1 + 0.5)
  ind6[1] <- 1 - step(tau5 - 1 + 0.5)
  beta_t[1] <- beta1*ind1[1] +
    beta2*ind2[1] +
    beta3*ind3[1] +
    beta4*ind4[1] +
    beta5*ind5[1] +
    beta6*ind6[1]
  p[1] <- 1 - exp(-beta_t[1] * I0 / n_pop)
  new_cases[1] ~ dbin(p[1], S0)
  S[1] <- S0 - new_cases[1]
  I[1] <- I0 + new_cases[1] - new_removals[1]
  
  # days 2…n_obs
  for(t in 2:n_obs) {
    ind1[t] <- step(tau1 -  t + 0.5)
    ind2[t] <- step(tau2 -  t + 0.5) - ind1[t]
    ind3[t] <- step(tau3 -  t + 0.5) - step(tau2 -  t + 0.5)
    ind4[t] <- step(tau4 -  t + 0.5) - step(tau3 -  t + 0.5)
    ind5[t] <- step(tau5 -  t + 0.5) - step(tau4 -  t + 0.5)
    ind6[t] <- 1 - step(tau5 -  t + 0.5)
    
    beta_t[t] <- beta1*ind1[t] +
      beta2*ind2[t] +
      beta3*ind3[t] +
      beta4*ind4[t] +
      beta5*ind5[t] +
      beta6*ind6[t]
    p[t] <- 1 - exp(-beta_t[t] * I[t-1] / n_pop)
    
    new_cases[t]    ~ dbin(p[t],     S[t-1])
    S[t]            <- S[t-1] - new_cases[t]
    I[t]            <- I[t-1] + new_cases[t] - new_removals[t]
  }
})

# --- 3) compile & run MCMC ----------------------------
inits <- function() list(
  beta1 = runif(1,0.01,0.1),
  beta2 = runif(1,0.01,0.1),
  beta3 = runif(1,0.01,0.1),
  beta4 = runif(1,0.01,0.1),
  beta5 = runif(1,0.01,0.1),
  beta6 = runif(1,0.01,0.1)
)

mod     <- nimbleModel(multiTauCode,
                       data      = data_list,
                       constants = const_list,
                       inits     = inits())
cmod    <- compileNimble(mod)

conf    <- configureMCMC(mod, monitors = paste0("beta",1:6))
rmcmc   <- buildMCMC(conf)
crmcmc  <- compileNimble(rmcmc, project = mod)

set.seed(2025)
samples_multi <- runMCMC(crmcmc,
                         niter             = 15000,
                         nburnin           = 3000,
                         nchains           = 3,
                         samplesAsCodaMCMC = TRUE)

# --- 4) inspect results -------------------------------
print(summary(samples_multi), digits = 3)

mcmc_dens(samples_multi, pars = paste0("beta",1:6)) +
  ggtitle("Posterior densities of β1…β6 for 5 fixed changepoints")

mcmc_trace(samples_multi, pars = paste0("beta",1:6)) +
  ggtitle("Trace plots of β1…β6 for 5 fixed changepoints")

mcmc_acf(samples_multi, pars = paste0("beta",1:6)) + ggtitle("ACF: β")



## Question 9 ##



# 1) Pre‐process the Czech data (replace cz1 with your filtered dataset if needed)
cases   <- cz1$cases
lag_inf <- 7
removals <- c(rep(0, lag_inf), head(cases, -lag_inf))

# 2) Bundle into NIMBLE data/constants
pop_CR   <- N        # total Czech population you defined earlier
n_obs    <- length(cases)

data_list  <- list(
  new_cases    = cases,
  new_removals = removals
)
const_list <- list(
  n_obs = n_obs,
  n_pop = pop_CR,
  I0    = 1           # assume 1 initial infected
)

# 3) Define the fully‐Bayesian change‐point S→I model
cp_code <- nimbleCode({
  # 3.1) changepoint prior: uniform over days 1...n_obs
  for(t in 1:n_obs) p_tau[t] <- 1 / n_obs
  tau ~ dcat(p_tau[1:n_obs])
  
  # 3.2) piecewise transmission rates
  beta1 ~ dlnorm(0, 5)   # before or on tau
  beta2 ~ dlnorm(0, 5)   # after tau
  
  # 3.3) initial susceptibles
  S0 <- n_pop - I0
  
  # 3.4) time t = 1
  ind1_1    <- step(tau - 1 + 0.5)           # 1 if 1 ≤ tau, else 0
  beta_t[1] <- beta1 * ind1_1 + beta2 * (1 - ind1_1)
  p[1]       <- 1 - exp(-beta_t[1] * I0 / n_pop)
  new_cases[1] ~ dbin(p[1], S0)
  S[1]         <- S0 - new_cases[1]
  I[1]         <- I0 + new_cases[1] - new_removals[1]
  
  # 3.5) loop for t = 2 ... n_obs
  for(t in 2:n_obs) {
    ind1[t]    <- step(tau - t + 0.5)
    beta_t[t]  <- beta1 * ind1[t] + beta2 * (1 - ind1[t])
    p[t]       <- 1 - exp(-beta_t[t] * I[t-1] / n_pop)
    
    new_cases[t]   ~ dbin(p[t],    S[t-1])
    S[t]           <- S[t-1] - new_cases[t]
    I[t]           <- I[t-1] + new_cases[t] - new_removals[t]
  }
})

# 4) Compile the model
inits_fn <- function() list(
  beta1 = runif(1, 0.01, 0.1),
  beta2 = runif(1, 0.01, 0.1),
  tau   = sample(1:n_obs, 1)
)

nim_mod   <- nimbleModel(
  code      = cp_code,
  data      = data_list,
  constants = const_list,
  inits     = inits_fn()
)
c_nim_mod <- compileNimble(nim_mod)

# 5) Configure & build MCMC
mcmc_conf  <- configureMCMC(nim_mod, monitors = c("beta1","beta2","tau"))
mcmc_run   <- buildMCMC(mcmc_conf)
c_mcmc_run <- compileNimble(mcmc_run, project = nim_mod)

# 6) Run the MCMC
set.seed(2025)
samples_cp <- runMCMC(
  c_mcmc_run,
  niter             = 15000,
  nburnin           = 3000,
  nchains           = 3,
  samplesAsCodaMCMC = TRUE
)

# 7) Summaries & diagnostics
print(summary(samples_cp), digits = 3)

# Trace and density of β₁ (pre‐τ) vs. β₂ (post‐τ)
mcmc_trace(samples_cp, pars = c("beta1","beta2")) +
  ggtitle("Trace: β1 (pre-τ) & β2 (post-τ)")

mcmc_dens(samples_cp, pars = c("beta1","beta2")) +
  ggtitle("Density: β1 (pre-τ) & β2 (post-τ)")

mcmc_acf(samples_cp, pars = c("beta1","beta2")) + ggtitle("ACF: β")









