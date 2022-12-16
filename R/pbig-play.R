pbig <- 0.01 # play with this variable 

years <- 1:2000
n_year <- length(years)
Rbig <- 5
sdr <- 0.6
set.seed(1)
sim <- get_recmult(pbig = pbig, Rbig, sdr)
sim$dat %>%
  ggplot(aes(x = year, y = recmult)) +
  geom_line() +
  theme_qfc()
sdr <- 0.06
