# -----------------------------------------------------------
# Omniscient manager control aka open-loop optimization
# Carl Walters and Chris Cahill 9 Oct 2022
# -----------------------------------------------------------
library(devtools)
library(TMB)
devtools::install_github("ChrisFishCahill/gg-qfc")
library(ggqfc)
library(tidyverse)
library(ggtext)
library(cowplot)
library(ggpmisc)
library(future)
library(furrr)
# -----------------------------------------------------------
# function to get recmult sequences:

get_recmult <- function(pbig, Rbig, sdr) {
  urand <- runif(n_year, 0, 1)
  Nrand <- rnorm(n_year, 0, 1)
  recmult <- rep(1, n_year)
  rlow <- (1 - pbig * Rbig) / (1 - pbig)
  if (rlow < 0) rlow <- 0
  for (t in 1:n_year) {
    recmult[t] <- rlow
    if (urand[t] < pbig) {
      recmult[t] <- Rbig
    }
    recmult[t] <- recmult[t] * exp(sdr * Nrand[t])
  }
  out <- tibble(
    year = 1:n_year,
    urand, Nrand, recmult
  )
  list(dat = out)
}

# testing recmult
# years <- 1:200
# n_year <- length(years)
# pbig <- 0.05
# Rbig <- 10
# sdr <- 0.4
# set.seed(1)
# sim <- get_recmult(n_year, pbig = 0.02, Rbig, sdr)
# sim$dat %>%
#   ggplot(aes(x = year, y = recmult)) +
#   geom_line() +
#   theme_qfc()

# -----------------------------------------------------------
# Set starting values
# leading parameters/values for simulation
years <- 1:2000
n_year <- length(years)
pbig <- 0.05
Rbig <- 10
sdr <- 0.4

ages <- 1:20 # slot 1 = recruits
cr <- 6
vbk <- .23
s <- .86
rinit <- 0.001
ro <- 1
uo <- 0.13
asl <- 0.5
ahm <- 6
upow <- 0.6
xinc <- 0.1

# simulate the om across these quantities
pbig <- 0.1
Rbig <- 0.99
sdr <- 0.3
ahv <- 5

# draw recruitment sequence
set.seed(1)
sim <- get_recmult(pbig, Rbig, sdr) 

tmb_data <- list(
  n_year = length(years),
  n_age = length(ages),
  vbk = vbk,
  s = s,
  cr = cr,
  rinit = rinit,
  ro = ro,
  uo = uo,
  asl = asl,
  ahv = ahv,
  ahm = ahm,
  upow = upow,
  ages = ages,
  recmult = sim$dat$recmult,
  obj_ctl = 0, # 0 = MAY, 1 = utility
  hcr = 2 # 0 = U(t), 1 = linear hcr, 2 = logistic hcr
)

if(tmb_data$hcr == 0){
  tmb_pars <- list(par = rep(0.5, length(years)))
}
if(tmb_data$hcr == 1){
  tmb_pars <- list(par = c(0.5, 0.5))
}
if(tmb_data$hcr == 2){
  tmb_pars <- list(par = rep(1, 3))
}
# compile and load the cpp
cppfile <- "src/om_hcr.cpp"
compile(cppfile)

dyn.load(dynlib("src/om_hcr"))
obj <- MakeADFun(tmb_data, tmb_pars,  silent = F, DLL = "om_hcr")
obj$fn()
obj$gr()

if(tmb_data$hcr == 0){
  lower = rep(0, length(years))
  upper = rep(1, length(years))
} 
if(tmb_data$hcr > 0){
  lower = rep(-Inf, length(tmb_pars$par))
  upper = rep(Inf, length(tmb_pars$par))
}

# run om simulation
opt <- nlminb(obj$par, obj$fn, obj$gr, 
              lower = lower,
              upper = upper
)
opt$par 
SD = sdreport( obj ) # standard errors








to_sim <- expand.grid(pbig = pbig, Rbig = Rbig, sdr = sdr, ahv = ahv, iter = iter)
to_sim <- to_sim %>% distinct()
glimpse(to_sim)

# set.seed(1)
out <- purrr::pmap(to_sim, run_om) # testing


future::plan(multisession)
system.time({
  out <- future_pmap(to_sim, run_om,
                     .options = furrr_options(seed = TRUE),
                     .progress = TRUE
  )
})
