# ------------------------------------------------------------------------------
# Harvest strategies for fisheries with highly variable recruitment dynamics
#                  Cahill and Walters Fall 2022
# NOTE!!!

# This
# File
# Is 
# Deprecated
# see
# om_hcr.R instead
# ------------------------------------------------------------------------------
library(devtools)
library(TMB)
devtools::install_github("ChrisFishCahill/gg-qfc")
library(ggqfc)
library(tidyverse)
library(future)
library(furrr)

# ------------------------------------------------------------------------------
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

get_fit <- function(hcrmode = c(
                      "OM", "linear", "logistic", "spline", "rect",
                      "db_logistic", "exponential", "logit", 
                      "misery", "dfo"
                    ),
                    objmode = c("yield", "utility")) {
  hcrmode <- match.arg(hcrmode)
  objmode <- match.arg(objmode)
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
    upow = case_when(
      objmode == "yield" ~ 1,
      objmode == "utility" ~ upow
    ),
    ages = ages,
    recmult = sim_dat$dat$recmult,
    hcrmode = case_when(
      hcrmode == "OM" ~ 0,
      hcrmode == "linear" ~ 1,
      hcrmode == "logistic" ~ 2,
      hcrmode == "spline" ~ 3,
      hcrmode == "rect" ~ 4,
      hcrmode == "db_logistic" ~ 5,
      hcrmode == "exponential" ~ 6,
      hcrmode == "dfo" ~ 7, 
      hcrmode == "experimental" ~ 8, 
      hcrmode == "misery" ~ 9
    ),
    knots = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 10),
    dfopar = c(Umsy, Bmsy),
    vmult = exp(sd_survey * rnorm(length(years)) - 0.5 * (sd_survey)^2)
  )
  if (pbig > 0.4) {
    tmb_data$knots <- c(0, 1.0, 2.0, 5.0, 10)
  }
  if (hcrmode == "OM") {
    tmb_pars <- list(par = rep(0.1, length(years)))
  } else if (hcrmode == "linear") {
    tmb_pars <- list(par = c(0.1, 0.1))
  } else if (hcrmode == "logistic") {
    tmb_pars <- list(par = rep(0.1, 3))
  } else if (hcrmode == "spline") {
    tmb_pars <- list(par = rep(0.1, length(tmb_data$knots)))
  } else if (hcrmode == "rect") {
    tmb_pars <- list(par = c(0.02, 0.01, 0.1))
    if (objmode == "yield") {
      tmb_pars <- list(par = c(0.9, 0.01, 0.9))
    }
  } else if (hcrmode == "db_logistic") {
     tmb_pars <- list(par = c(0.2, rep(0.5, 4)))
     # tmb_pars <- list(par = c(0.3, 10, 0.7, 10, 0.6))
  } else if (hcrmode == "exponential") {
    tmb_pars <- list(par = rep(0.1, 3))
  } else if (hcrmode == "dfo") {
    tmb_pars <- list(par = rep(0.1, 3)) # just a filler for dfo policy
  } else if (hcrmode == "experimental"){
    tmb_pars <- list(par = rep(0, 3)) 
  } else if (hcrmode == "misery"){
    tmb_pars <- list(par = c(0.6, 1, 0.7))
  }
  if (hcrmode == "OM") {
    lower <- rep(0, length(years))
    upper <- rep(1, length(years))
  }
  if (hcrmode != "OM") {
    lower <- rep(-Inf, length(tmb_pars$par))
    upper <- rep(Inf, length(tmb_pars$par))
  }
  if (hcrmode == "spline") {
    lower <- rep(0, length(tmb_pars$par))
    upper <- rep(Inf, length(tmb_pars$par))
  }
  if (hcrmode == "db_logistic") {
    lower <- c(-Inf, 0, -Inf, 0, -Inf)
    upper <- c(1, Inf, Inf, Inf, Inf)
  }
  if (hcrmode == "misery") {
    lower <- c(-Inf, -Inf, -Inf)
    upper <- c(Inf, 100, Inf)
  }
  if (!"om_hcr" %in% names(getLoadedDLLs())) {
    dyn.load(TMB::dynlib("src/om_hcr"))
  }
  obj <- MakeADFun(tmb_data, tmb_pars, silent = F, DLL = "om_hcr")
  if (hcrmode != "dfo") {
    opt <- nlminb(obj$par, obj$fn, obj$gr, upper = upper, lower = lower)
    ctr <- 1
    if (opt$convergence == 1 && ctr < 5) {
      tmb_pars <- list(par = opt$par)
      obj <- MakeADFun(tmb_data, tmb_pars, silent = T, DLL = "om_hcr")
      opt <- nlminb(obj$par, obj$fn, obj$gr, upper = upper, lower = lower)
      ctr <- ctr + 1
    }
    pdHess <- sdreport(obj)$pdHess
    objective <- -opt$objective
    convergence <- opt$convergence
    pdHess <- ifelse(pdHess == TRUE, 0, 1)
  }
  dat <- dplyr::tibble(
    "Ut" = obj$report()$`ut`,
    "Vulb" = obj$report()$`vulb`,
    "Abar" = obj$report()$`abar`,
    "Wbar" = obj$report()$`wbar`,
    "hcr" = hcrmode,
    "obj" = ifelse(hcrmode != "dfo", objective, -obj$fn()),
    "convergence" = ifelse(hcrmode != "dfo", convergence, 0),
    "pdHess" = ifelse(hcrmode != "dfo", pdHess, 0),
    "criterion" = objmode
  )
  list(dat, opt$par)
}

#-------------------------------------------------------------------------------
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
asl <- 0.5 #0.5
ahm <- 6
upow <- 0.6 # note this only matters when objmode == "utility", else upow = 1
ahv <- 5
pbig <- 1
Rbig <- 9
sdr <- 0.6

sd_survey <- 1e-6

#-------------------------------------------------------------------------------
# now that we have starting parameters, run equilibrium analysis to calculate
# Umsy and Bmsy for DFO rule
mat <- 1 / (1 + exp(-asl * (ages - ahm)))
wt <- (1 - exp(-vbk * (ages)))^3
vul <- 1 / (1 + exp(-asl * (ages - ahv)))

Lo <- rep(NA, length(ages))
Lo[1] <- 1

for (a in 2:(length(ages) - 1)) {
  Lo[a] <- Lo[a - 1] * s
}

# plus group
Lo[length(ages)] <- Lo[length(ages) - 1] * s / (1 - s)

mwt <- mat * wt
sbro <- sum(Lo * mwt)
reca <- cr / sbro
recb <- (cr - 1) / (ro * sbro)
ln_ar <- log(reca)

Useq <- seq(from = 0.01, to = 1, length.out = 100)
Umsy <- msy <- 0

for (i in 1:length(Useq)) {
  Req <- Yeq <- sbrf <- ypr <- 0
  su <- 1
  for (a in 1:length(ages)) {
    if (a == length(ages)) {
      su <- su / (1 - s * (1 - Useq[i] * vul[a])) # plus group effect
    }
    sbrf <- sbrf + su * mwt[a]
    ypr <- ypr + su * Useq[i] * vul[a] * wt[a]
    su <- su * s * (1 - Useq[i] * vul[a])
  }
  Req <- (exp(ln_ar + 0.5 * sdr^2) * sbrf - 1.0) / (recb * sbrf) # Beverton-Holt prediction
  Yeq <- Req * ypr
  if (Yeq > msy) {
    msy <- Yeq
    Umsy <- Useq[i]
    Bmsy <- msy / Umsy
  }
}
Bmsy
Umsy

#-------------------------------------------------------------------------------
# compile the cpp
cppfile <- "src/om_hcr.cpp"
compile(cppfile)
dyn.load(TMB::dynlib("src/om_hcr"))

years <- 1:2000
n_year <- length(years)
set.seed(1)
pbig <- 0.1 # 0.01, 0.05, 0.1, 0.25, 0.5, 1
sim_dat <- get_recmult(pbig = pbig, Rbig, sdr)
opt <- get_fit(hcrmode = "misery", objmode = "yield")
opt[[1]]$obj[1]
opt[[2]]

unique(opt[[1]]$convergence)
unique(opt[[1]]$pdHess)

plot(opt[[1]]$Ut ~ opt[[1]]$Vulb,
  xlab = "vulnerable biomass", ylab = "ut"
)

plot(opt[[1]]$Ut ~ opt[[1]]$Wbar, main = unique(round(opt[[1]]$obj)))

my_dat <- data.frame(Ut = opt[[1]]$Ut, Vulb = opt[[1]]$Vulb, Wbar = opt[[1]]$Wbar)
library(plot3D)
png("plots/3dcarl.png",
 width = 10, height = 5, units = "in", res = 1200
)
par(mfrow=c(1,3))
scatter3D(my_dat$Vulb, my_dat$Wbar, my_dat$Ut, clab = "U(t)", 
          theta = 320, phi = 20, 
          zlab = "Exploitation rate", xlab = "Vulb", ylab = "Wbar",
          cex = 0.6, pch = 16, bty = "b2")

scatter3D(my_dat$Vulb, my_dat$Wbar, my_dat$Ut, clab = "U(t)", 
          theta = 200, phi = 20, 
          zlab = "Exploitation rate", xlab = "Vulb", ylab = "Wbar",
          cex = 0.6, pch = 16, bty = "b2")


scatter3D(my_dat$Vulb, my_dat$Wbar, my_dat$Ut, clab = "U(t)", 
          theta = 360, phi = 20, 
          zlab = "Exploitation rate", xlab = "Vulb", ylab = "Wbar",
          cex = 0.6, pch = 16, bty = "b2")

dev.off()
p <-
  my_dat %>%
  ggplot(aes(x = Vulb, y = Wbar)) +
  geom_point(shape = 21, aes(fill = Ut)) +
  scale_fill_gradient(low = "white", high = "black") +
  # scale_colour_gradient(low = "blue", high = "darkorange2") +
  theme_qfc()
p

p1 <-
  my_dat %>%
  ggplot(aes(x = Vulb, y = Ut)) +
  geom_point(color = "black") +
  theme_qfc()
p1

p2 <-
  my_dat %>%
  ggplot(aes(x = Wbar, y = Ut)) +
  geom_point(color = "black") +
  theme_qfc()
p2

bigplot <- cowplot::plot_grid(p, p1, p2, nrow = 3)
sim_dat <- get_recmult(pbig = pbig, Rbig, sdr)

which_rules <- c("linear", "db_logistic", "logistic")
system.time({
  dat <- NULL
  for (i in which_rules) {
    opt <- get_fit(hcrmode = i, objmode = "yield")
    if (is.null(dat)) {
      dat <- opt[[1]]
    } else {
      dat <- rbind(dat, opt[[1]])
    }
  }
})

unique(dat$convergence)
unique(dat$hcr[which(dat$convergence == 1)])
summary(warnings())

dat$Pbig <- pbig
pd <- dat %>%
  mutate(obj = obj / max(obj)) %>%
  mutate(Yield = paste0(hcr, " = ", format(round(obj, 3), nsmall = 3)))

my_levels <- unique(pd$Yield[rev(order(unlist(str_extract_all(pd$Yield, "\\(?[0-9,.]+\\)?"))))])
pd$Yield <- factor(pd$Yield, levels = my_levels)
p <-
  pd %>%
  filter(hcr != 0) %>%
  ggplot(aes(x = Vulb, y = Ut, color = Yield)) +
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
  ggtitle(bquote(P[big] ~ `=` ~ .(pd$Pbig)))
p
pall <- cowplot::plot_grid(p3, p4, nrow = 2)


pall <- cowplot::plot_grid(p, p1, p2, p3, p4, p5, nrow = 3, scale = 0.98)
ggsave("plots/upow.pdf", width = 8, height = 11)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Extra code:


# hcrmodes <- 0:4
# objmodes <- 1
# sdr <- 0.6
# pbigs <- c(0.01, 0.5, 1.0)
# sims <- 1
# seed <- 1
# to_fit <- expand.grid(hcrmodes, objmodes, pbigs, sims, seed)
# names(to_fit) <- c("hcrmode", "objmode", "pbig", "sim", "seed")
# to_fit <- to_fit %>% distinct()
# system.time({
#   out <- purrr::pmap_df(to_fit, get_fit)
# })
#
# future::plan(multisession)
# system.time({
#   out <- future_pmap_dfr(to_fit, get_fit,
#     .options = furrr_options(seed = TRUE),
#     .progress = TRUE
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
#
# dat %>%
#   filter(hcr == 0) %>%
#   mutate(year = 1:length(years)) %>%
#   ggplot(aes(x=years, y = Ut))+
#   geom_line()
# out <- out %>% group_by(Pbig) %>%
#   mutate(obj = obj / max(obj))
#
# plot_list <- list()
# for (i in unique(out$Pbig)) {
#   pd <- out %>%
#     filter(Pbig == i) %>%
#     mutate(Utility = as.factor(case_when(
#       hcr == "0" ~ paste0("OM = ", format(round(obj, 3), nsmall = 3)),
#       hcr == "1" ~ paste0("linear = ", format(round(obj, 3), nsmall = 3)),
#       hcr == "2" ~ paste0("logistic = ", format(round(obj, 3), nsmall = 3)),
#       hcr == "3" ~ paste0("spline = ", format(round(obj, 3), nsmall = 3)),
#       hcr == "4" ~ paste0("rectilinear = ", format(round(obj, 3), nsmall = 3))
#     )))
#   my_levels <- unique(pd$Utility[rev(order(unlist(str_extract_all(pd$Utility, "\\(?[0-9,.]+\\)?"))))])
#   pd$Utility <- factor(pd$Utility, levels = my_levels)
#   p <-
#     pd %>%
#     ggplot(aes(x = Vulb, y = Ut, color = Utility)) +
#     geom_point(size = 0.25) +
#     scale_color_brewer(palette = "Paired") +
#     ylab(expression(Exploitation ~ rate ~ U[t])) +
#     xlab("Vulnerable biomass") +
#     guides(color = guide_legend(reverse = TRUE)) +
#     theme_qfc() +
#     theme(
#       legend.position = c(.8, .75),
#       legend.title.align = 0.5
#     ) +
#     guides(colour = guide_legend(override.aes = list(size = 3))) +
#     scale_alpha(guide = "none") +
#     ggtitle(bquote(P[big] ~ `=` ~ .(pd$Pbig)))
#   plot_list[[which(unique(out$Pbig) == i)]] <- p
# }
#
# p1 <- plot_list[[1]]
# p2 <- plot_list[[2]]
# p3 <- plot_list[[3]]
# p4 <- plot_list[[4]]
# p5 <- plot_list[[5]]
# p6 <- plot_list[[6]]
#
