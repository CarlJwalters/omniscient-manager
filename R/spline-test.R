# get a simple spline working in TMB
library(splines)

X <- seq(from=0, to=20, by=.01) # generating inputs
B <- t(bs(X, knots=1, degree=1, intercept = FALSE)) # creating the B-splines
num_data <- length(X); num_basis <- nrow(B)
a0 <- 0 # intercept
a <- rnorm(num_basis, 0, 1) # coefficients of B-splines
Y_true <- as.vector(a0*X + a%*%B) # generating the data
Y <- Y_true + rnorm(length(X),0,.2)
plot(Y~X)

attributes(B) <- attributes(B)["dim"] # convert to matrix (misery)

library(TMB)
cppfile <- "src/spline.cpp"
compile(cppfile)

tmb_data <- list(
  "n_data" = length(Y), 
  "Y" = Y, 
  "X" = as.matrix(X), 
  "B" = as.matrix(t(B)) # need to transpose for TMB
)

tmb_pars <- list("ao" = c(0), 
                 "a" = c(a), 
                 "ln_sig" = log(1)
)

dyn.load(dynlib("src/spline"))
obj <- MakeADFun(tmb_data, tmb_pars,  silent = F, DLL = "spline")

obj$fn()
obj$gr()



# run om simulation
opt <- nlminb(obj$par, obj$fn, obj$gr)
SD = sdreport( obj ) # standard errors
SD
opt <- nlminb(opt$par, obj$fn, obj$gr)
opt$convergence
library(TMBhelper)
TMBhelper::fit_tmb(obj, obj$fn, obj$gr, newtonsteps = 1)
