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
ages <- 1:20
cr <- 6
vbk <- .23
s <- .86
rinit <- 0.01
ro <- 1.0
uo <- 0.13
asl <- 0.5
ahm <- 6
upow <- 0.6
ahv <- 5

pbig <- 0.05
Rbig <- 7
sdr <- 0.6
ahv <- 5

# simulate recruitment sequence
set.seed(13)
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
  objmode = 1, # 0 = MAY, 1 = utility
  hcrmode = 0, # 0 = U(t), 1 = linear hcr, 2 = logistic hcr, 3 = spline, 4 = rectilinear
  knots =  c(0,0.2,0.4,0.6,0.8,1.0,1.5,2,10)
)

# set up the pars
if (tmb_data$hcr == 0) {
  tmb_pars <- list(par = rep(0.1, length(years)))
} else if (tmb_data$hcr == 1) {
  tmb_pars <- list(par = c(0.1, 0.1))
} else if (tmb_data$hcr == 2) {
  tmb_pars <- list(par = rep(0.1, 3))
} else if (tmb_data$hcr == 3) {
  tmb_pars <- list(par = rep(0.1, length(tmb_data$knots)))
  #tmb_pars$par <- 0.2177*(tmb_data$knots-0.24)/(tmb_data$knots+1e-10)
  #tmb_pars$par <- ifelse(tmb_pars$par < 0, 0.02, tmb_pars$par)
} else if (tmb_data$hcr == 4) {
  tmb_pars <- list(par = c(0.02, 0.01, 0.1))
}

# set upper and lower bounds
if (tmb_data$hcr == 0) {
  lower <- rep(0, length(years))
  upper <- rep(1, length(years))
}
if (tmb_data$hcr > 0) {
  lower <- rep(-Inf, length(tmb_pars$par))
  upper <- rep(Inf, length(tmb_pars$par))
}
if(tmb_data$hcr == 3){
  lower <- rep(0, length(tmb_pars$par))
  upper <- rep(Inf, length(tmb_pars$par))
}
# compile and load the cpp
cppfile <- "src/om_hcr.cpp"
compile(cppfile)
dyn.load(dynlib("src/om_hcr"))
obj <- MakeADFun(tmb_data, tmb_pars, silent = F, DLL = "om_hcr")

# run om simulation
opt <- nlminb(obj$par, obj$fn, obj$gr, upper = upper, lower = lower)

while (opt$convergence == 1) {
  tmb_pars <- list(par = opt$par)
  obj <- MakeADFun(tmb_data, tmb_pars, silent = T, DLL = "om_hcr")
  opt <- nlminb(obj$par, obj$fn, obj$gr, upper = upper, lower = lower)
}

opt

sdreport(obj)

opt$objective

dat <- data.frame("Ut" = obj$report()$`ut`, 
                  "Vulb" = obj$report()$`vulb`, 
                  "hcr" = tmb_data$hcrmode
)

#------------------------------------------------------------------------------
# plot(obj$report(opt$par)$`ut` * obj$report(opt$par)$`vulb` ~
#   obj$report(opt$par)$`vulb`,
# type = "p", xlim = c(0, 10),
# xlab = "vulnerable biomass", ylab = "TAC"
# )

plot(obj$report()$`ut` ~ obj$report()$`vulb`,
  ylab = "Ut", xlab = "vulb", pch = 16, cex = 0.5
  #main = paste0("objective = ", round(-opt$objective, 2))
)
points(obj$report(opt$par)$`ut` ~ obj$report(opt$par)$`vulb`, col = "blue", 
       type = "p", cex = 0.5, pch = 16)
points(obj$report(opt$par)$`ut` ~ obj$report(opt$par)$`vulb`, col = "red", 
       type = "p", cex = 0.5, pch = 16)


plot(obj$report(opt$par)$`yield`, type = "b")
plot(obj$report(opt$par)$`ut` * obj$report(opt$par)$`vulb` ~ 
       obj$report(opt$par)$`vulb`, type = "p", xlab = "vulb", ylab = "TAC")
points(obj$report(opt$par)$`ut` * obj$report(opt$par)$`vulb` ~ obj$report(opt$par)$`vulb`, col = "blue", type = "b")

