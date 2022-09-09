##--------------------------------------------
## Discrete version of a species distribution 
## model accounting for preferential sampling
##--------------------------------------------

require(TMB)

compile("models/1_ipp_discrete_toy_model.cpp")
dyn.load(dynlib("models/1_ipp_discrete_toy_model"))

## Create grid
grid_dim <- 25
n_cells <- grid_dim^2
loc_x <- expand.grid( "x"=1:grid_dim, "y"=1:grid_dim)
loc_x$cell <- 1:n_cells

## Simulate latent field of biomass (S_x)
source("r/function/sim_GF_Matern.R")
nu <- 1 # nu parameter of the matérn function
range_delta <- range_eta <- 10 # range parameter
SD_delta <- SD_eta <- 1 # marginal sd parameter
intercept_S <- runif(1,min = -1, max = 1)
delta_x <- sim_GF_Matern(loc_x, nu, range_delta, SD_delta^2)[[1]]$y # simulate random effect
S_x <- exp(intercept_S + delta_x) # biomass field

## Sample data preferentially
n_samp <- 300 # number of samples
intercept_l <- runif(1,min = -1, max = 1) # intercept of the point process
b <- runif(1,min = 0, max = 3) # preferential sampling parameter
eta_x <- sim_GF_Matern(loc_x, nu, range_eta, SD_eta^2)[[1]]$y # additionnal processes affecting data sampling
lambda <- exp(intercept_l + b*log(S_x) + eta_x) # intensity of sampling process
index_i <- sample(loc_x$cell, # Generate sampled locations
                  size=n_samp,
                  replace=T,
                  prob = lambda) # sampling depends on both biomass and additionnal processes
c_x <- rep(0,n_cells) # Number of data points in each cell
c_x <- do.call(c,lapply(1:n_cells, function(j){
  c_x[j] <- length(which(index_i == j))
}))


## Simulate observations (lognormal distribution)
SD_obs <- 0.5 # observations standard error
y_i <- do.call(c,lapply(1:length(index_i), function(j){
  exp_catch <- S_x[index_i[j]] # expected catch
  y_sci_i <- rlnorm(1,meanlog = log(exp_catch),sd = SD_obs)
  return(y_sci_i)
}))

# check: plot(log(y_i),S_x[index_i])

## Precision matrix
# Load from data
load(file="data/1_mesh.RData")
load(file="data/1_spde.RData")

# # # Build spde stuff from INLA functions (INLA library must be installed)
# library(INLA)
# mesh <- inla.mesh.create( loc_x[,c('x','y')] )
# spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]  # define sparse matrices for parametrisation of precision matrix

## Shape inputs (Data, Params, Map, Random) for TMB model
Options_vec <- c("ipp_lk"=0, # form of the ipp likelihood (0: raw form, 1: discretized form)
                 "SE"=1, # correction bias (0: no, 1: yes)
                 "RE"=1, # account for random effect (0: no, 1: yes)
                 "SP"=1) # account for sampling process (0:no, 1: yes)

Data <- list("Options_vec"=Options_vec,
             "n_cells" = n_cells,
             "c_x"=c_x,
             "y_i"=y_i,
             "index_i"=index_i-1,
             "spde"=spde)

Params <- list("intercept_S"=0,
               "deltainput_x"=rep(0,mesh$n),
               "intercept_l"=0,
               "b"=0,
               "etainput_x"=rep(0,mesh$n),
               "logtau"=c(0,0),
               "logkappa"=c(0,0),
               "logSD_obs"=log(1))

Map <- list()
if(Options_vec["RE"] == 1){ # estimate random effect
  Random <- c("deltainput_x","etainput_x")
}
if(Options_vec["RE"] == 0){ # do not estimate random effect
  Map[["deltainput_x"]] <- Map[["etainput_x"]] <- factor(rep(NA,mesh$n))
  Map[["logtau"]] <- Map[["logkappa"]] <- factor(rep(NA,2))
  Random=NULL
}

if(Options_vec["SP"] == 0){ # do not account for sampling process and preferential sampling
  Random <- c("deltainput_x")
  Map[["etainput_x"]] <- factor(rep(NA,mesh$n))
  Map[["logtau"]] <- Map[["logkappa"]] <- factor(c(1,NA))
  Map[["b"]] <- factor(NA)
  Map[["intercept_l"]] <- factor(NA)
}

## Fit TMB model
obj <- MakeADFun(data=Data,
                 parameters=Params,
                 random=Random,
                 map = Map,
                 DLL="1_ipp_discrete_toy_model")
runSymbolicAnalysis(obj)
obj$fn( obj$par )
system.time(opt <- nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr,
                           lower=-Inf, upper=+Inf,
                           control=list(trace=1, maxit=1000)))
Report <- obj$report()
rep <- sdreport(obj,bias.correct=TRUE)

# # check : 
# plot(Report$S_x,S_x)


# ## Plot matrices
# library(Matrix)
# # Hessian
# h <- obj$env$spHess(random=TRUE)
# image(h)
# 
# # Precision matrix (same as hessian here)
# kappa <- sqrt(8)/range_delta
# Q <- kappa ^ 4 * spde$M0 + 2 * kappa ^ 2 * spde$M1 + spde$M2
# # Q[which(Q < 1e-15)] <- 0
# Q@x[] <- 1
# n_per_row <- rowSums(Q)
# image(Q)
# hist(n_per_row)
# 
# # Cholesky matrix (I suppose...)
# L <- obj$env$L.created.by.newton
# image(L)

