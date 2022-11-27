# -----------------------------------------------------------
# Omniscient manager control aka open-loop optimization
# Cahill and Walters Nov 2022
# -----------------------------------------------------------
library(devtools)
library(TMB)
devtools::install_github("ChrisFishCahill/gg-qfc")
library(ggqfc)
library(tidyverse)
library(future)
library(furrr)
# -----------------------------------------------------------
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
  out <- dplyr::tibble(
    year = 1:n_year,
    urand, Nrand, recmult
  )
  list(dat = out)
}

get_fit <- function(hcrmode = NA, objmode = NA) {
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
    objmode = objmode, # 0 = MAY, 1 = utility
    hcrmode = hcrmode, # 0 = U(t), 1 = linear hcr, 2 = logistic hcr, 3 = spline, 4 = rectilinear
    knots = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2, 10)
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
  if (tmb_data$hcr == 3) {
    lower <- rep(0, length(tmb_pars$par))
    upper <- rep(Inf, length(tmb_pars$par))
  }
  # load the cpp
  dyn.load(TMB::dynlib("src/om_hcr"))
  #openmp(max=TRUE)
  obj <- TMB::MakeADFun(tmb_data, tmb_pars, silent = F, DLL = "om_hcr")

  # run om simulation
  opt <- nlminb(obj$par, obj$fn, obj$gr, upper = upper, lower = lower)

  while (opt$convergence == 1) {
    tmb_pars <- list(par = opt$par)
    obj <- MakeADFun(tmb_data, tmb_pars, silent = T, DLL = "om_hcr")
    opt <- nlminb(obj$par, obj$fn, obj$gr, upper = upper, lower = lower)
  }
  if (tmb_data$hcrmode == 0) {
    opt$SD <- NULL
  } else {
    opt$SD <- sdreport(obj)
  }
  dat <- data.frame(
    "Ut" = obj$report()$`ut`,
    "Vulb" = obj$report()$`vulb`,
    "hcr" = tmb_data$hcrmode,
    "obj" = -opt$objective,
    "convergence" = opt$convergence
  )
  opt$dat <- dat
  opt
}

# compile the cpp
cppfile <- "src/om_hcr.cpp"
compile(cppfile)

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
pbig <- 1.0 # 0.01, 0.1, 0.25, 0.5, 0.75, 1.0
Rbig <- 7
sdr <- 0.6

# testing recmult
# years <- 1:2000
# n_year <- length(years)
# pbig <- 0.05
# Rbig <- 10
# sdr <- 0.4
# set.seed(1)
# sim <- get_recmult(pbig = 0.02, Rbig, sdr)
# sim$dat %>%
#   ggplot(aes(x = year, y = recmult)) +
#   geom_line() +
#   theme_qfc()

# simulate recruitment sequence
set.seed(3)
sim <- get_recmult(pbig, Rbig, sdr)

#system.time({
#opt = get_fit(hcrmode = 0, objmode = 1)
#})
# estimate the hcr
system.time({
dat <- NULL
for (i in 0:4) {
  opt <- get_fit(hcrmode = i, objmode = 1)
  if (is.null(dat)) {
    dat <- opt$dat
  } else {
    dat <- rbind(dat, opt$dat)
  }
}
})
summary(warnings())

# to_fit = tibble(hcrmode = 0:4, objmode = 1L)

# set.seed(1)
# system.time({
# out <- purrr::pmap(to_fit, get_fit) # testing
# })
# 
# future::plan(multisession)
# system.time({
#   out <- future_pmap(to_fit, get_fit,
#                      .options = furrr_options(seed = TRUE),
#                      .progress = TRUE
#   )
# })

# extra code
# tmb_pars$par <- 0.2177*(tmb_data$knots-0.24)/(tmb_data$knots+1e-10)
# tmb_pars$par <- ifelse(tmb_pars$par < 0, 0.02, tmb_pars$par)
# plot(dat$Ut ~ dat$Vulb, xlim = c(0,5),
#     ylab = "Ut", xlab = "vulb", pch = 16, cex = 0.5,
#      col=as.factor(dat$hcr)
#     #main = paste0("objective = ", round(-opt$objective, 2))
# )

dat %>%
  filter(hcr == 1) %>%
  mutate(year = 1:length(years)) %>%
  ggplot(aes(x=years, y = Ut))+
  geom_line()

pd <- dat %>%
  mutate(obj = obj / max(obj)) %>%
  mutate(Utility = as.factor(case_when(
    hcr == "0" ~ paste0("OM = ", format(round(obj, 3), nsmall = 3)),
    hcr == "1" ~ paste0("linear = ", format(round(obj, 3), nsmall = 3)),
    hcr == "2" ~ paste0("logistic = ", format(round(obj, 3), nsmall = 3)),
    hcr == "3" ~ paste0("spline = ", format(round(obj, 3), nsmall = 3)),
    hcr == "4" ~ paste0("rectilinear = ", format(round(obj, 3), nsmall = 3))
  )))
my_levels <- unique(pd$Utility[rev(order(unlist(str_extract_all(pd$Utility, "\\(?[0-9,.]+\\)?"))))])
pd$Utility <- factor(pd$Utility, levels = my_levels)

p5 <-
  pd %>%
  ggplot(aes(x = Vulb, y = Ut, color = Utility)) +
  geom_point(size = 0.25) +
  scale_color_brewer(palette = "Paired") +
  ylab(expression(Exploitation ~ rate ~ U[t])) +
  xlab("Vulnerable biomass") +
  guides(color = guide_legend(reverse = TRUE)) +
  theme_qfc() +
  theme(
    legend.position = c(.8, .75),
    legend.title.align = 0.5
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_alpha(guide = "none") + 
  ggtitle(bquote(P[big]~`=`~ .(pbig)))
p5


pall <- cowplot::plot_grid(p, p1, p2, p3, p4, p5, nrow = 3, scale = 0.98)
ggsave("plots/pbigs.pdf", width = 8, height = 10)
