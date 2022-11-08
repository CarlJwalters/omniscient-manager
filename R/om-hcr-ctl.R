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
library(splines)
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
rinit <- 0.6
ro <- 1
uo <- 0.13
asl <- 0.5
ahm <- 6
upow <- 0.6
ahv <- 5

# draw recruitment sequence
#set.seed(1)

pbig <- c(0.02)
Rbig <- c(7)
sdr <- c(0.25)
ahv <- c(5)

# Set up the B spline
vulb_seq <- seq(from=0, to=2, by=.001) 
B <- t(bs(vulb_seq, knots=2, degree=3, intercept = FALSE)) # creating the B-splines
num_basis <- nrow(B)
a = rnorm(num_basis)
Y_true <- as.vector(a%*%B)
plot(arm::invlogit(Y_true)~vulb_seq)

a = opt$par

# simulate recruitment sequence
# sim <- get_recmult(pbig, Rbig, sdr) 

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
  hcr = 3, # 0 = U(t), 1 = linear hcr, 2 = logistic hcr, 3 = spline
  B = as.matrix(t(B)), # need to transpose for TMB
  X = as.matrix(vulb_seq)
)

# set up the pars
if(tmb_data$hcr == 0){
  tmb_pars <- list(par = rep(0.5, length(years)))
}
if(tmb_data$hcr == 1){
  tmb_pars <- list(par = c(0.5, 0.5))
}
if(tmb_data$hcr == 2){
  tmb_pars <- list(par = jitter(rep(1, 3)))
}
if(tmb_data$hcr == 3){
  a <- rep(0, num_basis)
  # ao <- 0
  tmb_pars <- list(par = c(a))
}

if(tmb_data$hcr == 0){
  lower = rep(0, length(years))
  upper = rep(1, length(years))
} 
if(tmb_data$hcr > 0){
  lower = rep(-Inf, length(tmb_pars$par))
  upper = rep(Inf, length(tmb_pars$par))
}

#tmb_pars <- list(par = log(rep(0.1, length(a) + 1)))
# compile and load the cpp
cppfile <- "src/om_hcr.cpp"
compile(cppfile)

dyn.load(dynlib("src/om_hcr"))
obj <- MakeADFun(tmb_data, tmb_pars,  silent = F, DLL = "om_hcr")

obj$fn()
obj$gr()
obj$report()$`ut`

# run om simulation
opt <- nlminb(obj$par, obj$fn, obj$gr, 
              lower = lower, 
              upper = upper
              )
opt$objective

plot(obj$report()$`ut`~obj$report()$`vulb`)

# -67.91805 spline
#  -198.6367 logistic
# -198.63 line 


# re-run the optimization until convergence achieved
while (opt$convergence == 1) {
  tmb_pars <- list(par = opt$par)
  obj <- MakeADFun(tmb_data, tmb_pars, silent = T, DLL = "om_hcr")
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                lower = rep(0, length(years)),
                upper = rep(1, length(years))
  )
}
opt$objective
vulb = seq(from = 0, to = 10, by = 0.05)

Ut = opt$par[1] / (1 + exp(-opt$par[2]*(vulb-opt$par[3])))
lines(Ut~vulb, type = "l")
opt$par


# 201.26 for logistic
# 201.24 for linear
# 217.344
SD = sdreport( obj ) # standard errors
opt$par

# spline play
library(splines)
X <- seq(from=0, to=20, by=.1) # generating inputs
B <- t(bs(X, knots=3, degree=3, intercept = TRUE)) # creating the B-splines
num_data <- length(X); num_basis <- nrow(B)
a0 <- 0 # intercept
a <- rnorm(num_basis, 0, 1) # coefficients of B-splines
Y_true <- as.vector(a0*X + a%*%B) # generating the output
Y <- Y_true + rnorm(length(X),0,.2)
plot(Y~X)



Knots <- default.knots(vulb, K)

X.bs




opt <- nlminb(opt$par, obj$fn, obj$gr, 
              lower = lower, 
              upper = upper
             )
opt$convergence

plot(opt$par, type = "b")

plot(opt$par~obj$report(opt$par)$`vulb`)
dat <- data.frame(Ut = opt$par, vulb=obj$report(opt$par)$`vulb`, 
                  time = 1:n_year)

library(mgcv)
gam_setup = gam(dat$Ut ~ s(dat$vulb, bs = "cs"),
                fit=FALSE)

#Extrtact penelization matrices
S_vulb = gam_setup$smooth[[1]]$S[[1]]


devtools::install_github("AckerDWM/gg3D")
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
