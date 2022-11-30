# -----------------------------------------------------------
# Omniscient manager control aka open-loop optimization
# Cahill and Walters Nov 2022
# TODO:  dfo precautionary rule
#        better plotting scheme
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
  out <- tibble(
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
    recmult = sim_dat$dat$recmult,
    objmode = objmode,
    hcrmode = hcrmode,
    knots = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 10)
  )
  if (pbig > 0.4) {
    tmb_data$knots <- c(0, 1.0, 2.0, 5.0, 10)
  }
  if (tmb_data$hcr == 0) {
    tmb_pars <- list(par = rep(0.1, length(years)))
  } else if (tmb_data$hcr == 1) {
    tmb_pars <- list(par = c(0.1, 0.1))
  } else if (tmb_data$hcr == 2) {
    tmb_pars <- list(par = rep(0.1, 3))
  } else if (tmb_data$hcr == 3) {
    tmb_pars <- list(par = rep(0.1, length(tmb_data$knots)))
    # tmb_pars$par <- 0.2177 * (tmb_data$knots - 0.24) / (tmb_data$knots + 1e-10)
    # tmb_pars$par <- ifelse(tmb_pars$par < 0, 0.02, tmb_pars$par)
  } else if (tmb_data$hcr == 4) {
    tmb_pars <- list(par = c(0.02, 0.01, 0.1))
  }
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
  if (!"om_hcr" %in% names(getLoadedDLLs())) {
    dyn.load(TMB::dynlib("src/om_hcr"))
  }
  obj <- MakeADFun(tmb_data, tmb_pars, silent = F, DLL = "om_hcr")
  opt <- nlminb(obj$par, obj$fn, obj$gr, upper = upper, lower = lower)
  ctr <- 1
  if (opt$convergence == 1 && ctr < 5) {
    tmb_pars <- list(par = opt$par)
    obj <- MakeADFun(tmb_data, tmb_pars, silent = T, DLL = "om_hcr")
    opt <- nlminb(obj$par, obj$fn, obj$gr, upper = upper, lower = lower)
    ctr <- ctr + 1
  }
  dat <- dplyr::tibble(
    "Ut" = obj$report()$`ut`,
    "Vulb" = obj$report()$`vulb`,
    "hcr" = tmb_data$hcrmode,
    "obj" = -opt$objective,
    "convergence" = opt$convergence
  )
  dat
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
pbig <- 1
Rbig <- 9
sdr <- 0.6

#-------------------------------------------------------------------------------
# equilibrium analysis to calculate Umsy and Umay for DFO rule
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
Umsy <- MSY <- 0

for (i in 1:length(Useq)) {
  Req <- Yeq <- sbrf <- ypr <- 0 
  su <- 1
  for(a in 1:length(ages)){
    if(a==length(ages)){su=su/(1-s*(1-Useq[i]*vul[a]))} # plus group effect
    sbrf = sbrf + su*mwt[a] 
    ypr = ypr + su*Useq[i]*vul[a]*wt[a]
    su = su * s * (1-Useq[i]*vul[a])
  }
  Req <- (exp(ln_ar + 0.5 * sdr^2) * sbrf - 1.0) / (recb * sbrf) # Beverton-Holt prediction
  Yeq <- Req * ypr
  if (Yeq > MSY) {
    MSY <- Yeq
    Umsy <- Useq[i]
  }
}
MSY
Umsy


#-------------------------------------------------------------------------------

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
# set.seed(333)
# sim <- get_recmult(pbig, Rbig, sdr)

# get_fit(hcrmode = 3, objmode = 1, pbig = 0.5, sim = 1, seed = 19)

years <- 1:2000
n_year <- length(years)
set.seed(1)
pbig <- 0.25 # 0.01, 0.05, 0.1, 0.25, 0.5, 1
sim_dat <- get_recmult(pbig = pbig, Rbig, sdr)

system.time({
  dat <- NULL
  for (i in 0:4) {
    opt <- get_fit(hcrmode = i, objmode = 1)
    if (is.null(dat)) {
      dat <- opt
    } else {
      dat <- rbind(dat, opt)
    }
  }
})
unique(dat$convergence)
summary(warnings())


dat$Pbig <- pbig
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
  ggtitle(bquote(P[big] ~ `=` ~ .(pd$Pbig)))
p5

pall <- cowplot::plot_grid(p, p1, p2, p3, p4, p5, nrow = 3, scale = 0.98)
ggsave("plots/pbigs.pdf", width = 8, height = 11)




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
