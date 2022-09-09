##-------------------------------------------
## Spatio-temporal species distribution 
## model accounting for preferential sampling
##-------------------------------------------
# n.b. estimation includes a triangulation step to estimate the random effect
# (Cf. SPDE approach)
require(TMB)

compile("models/2_ipp_st_toy_model.cpp")
dyn.load(dynlib("models/2_ipp_st_toy_model"))
run_simulation <- F
# if False, load simulated latent field, data from RData files

if(run_simulation){
  
  ############################
  ## Simulation configurations
  ############################
  ## Create grid
  grid_dim <- 50
  n_cells <- grid_dim^2
  loc_x <- expand.grid( "x"=1:grid_dim, "y"=1:grid_dim)
  loc_x$cell <- 1:n_cells
  
  n_step <- 10 # number of time steps
  
  ## Latent field
  #--------------
  source("r/function/sim_GF_Matern.R")
  nu <- 1 # nu parameter of the matérn function
  range_delta <- range_eta <- grid_dim/2 # range parameter
  SD_delta <- SD_eta <- 1 # marginal sd parameter
  intercept_S <- rep(runif(1,min = -1, max = 1),n_step) # same intercept over the full period (but could be something else)
  rho_delta <- 0.5 # temporal auto-correlation parameter
  
  # Biomass field and random effect matrices
  S_x.t <- matrix(data = NA,nrow = n_cells,ncol = n_step)
  deltaAR1 <- matrix(data = NA,nrow = n_cells,ncol = n_step)
  
  # First time step random effect
  delta_x.0 <- sim_GF_Matern(loc_x, nu, range_delta, SD_delta^2)[[1]]$y # simulate random effect
  deltaAR1[,1] <- delta_x.0
  
  ## Sampling process
  #------------------
  n_samp <- rep(300,n_step) # number of samples per time step
  intercept_l <- rep(runif(1,min = -1, max = 1),n_step) # intercept of the point process
  b <- rep(runif(1,min = 0, max = 3),n_step) # preferential sampling parameter
  
  # Fishing intensity and random effect matrices
  lambda_x.t <- matrix(data = NA,nrow = n_cells,ncol = n_step) # fishing intensity (spatio-temporal)
  c_x.t <- matrix(data = NA,nrow = n_cells,ncol = n_step) # number of fishing points per cell and time step
  eta_x.t <- matrix(data = NA,nrow = n_cells,ncol = n_step) # additionnal processes affecting data sampling (spatial random effect for each time step)
  
  ## Data vectors
  #--------------
  t_i2 <- c() # time step of each observation
  index_i2 <- c() # grid cell of each observation
  # y_i: vector of all observations (simulated below)
  
  ##############
  ## Simulations
  ##############
  for(t in 1:n_step){
    
    print(t)
    
    ## Latent field
    #--------------
    if(t>=2){
      
      delta_x.t <- sim_GF_Matern(loc_x, nu, range_delta, SD_delta^2)[[1]]$y
      deltaAR1[, t] <- rho_delta * deltaAR1[, t - 1] + sqrt(1 - rho_delta^2) * delta_x.t
      
    }
    
    S_x.t[,t] = exp(intercept_S[t] + deltaAR1[, t])
    
    ## Sampling process
    #------------------
    # sampling depends on both biomass (b * log(S_x.t)) and additionnal processes (eta_x.t)
    eta_x.t[,t] <- sim_GF_Matern(loc_x, nu, range_eta, SD_eta^2)[[1]]$y
    lambda_x.t[,t] <- exp(intercept_l[t] + b[t]*log(S_x.t[,t]) + eta_x.t[,t]) # intensity of sampling process
    
    index_i <- sample(loc_x$cell, # samples' locations
                      size=n_samp[t],
                      replace=T,
                      prob = lambda_x.t[,t])
    index_i2 <- c(index_i2,
                  index_i)
    
    t_i <- rep(t,n_samp[t]) # samples' time step
    t_i2 <- c(t_i2,t_i)
    
    c_x.t[,t] <- do.call(c,lapply(1:n_cells, function(j){
      c_x.t[,t][j] <- length(which(index_i == j))
    }))
    
  }
  
  ## Observations (lognormal distribution)
  #--------------
  SD_obs <- 1 # observations' standard error
  y_i <- do.call(c,lapply(1:length(index_i2), function(j){
    exp_catch <- S_x.t[index_i2[j],t_i2[j]] # expected catch
    y_sci_i <- rlnorm(1,meanlog = log(exp_catch),sd = SD_obs)
    return(y_sci_i)
  }))
  
  # ## Plot and check
  # #----------------
  # # Values for scaling
  # par(mfrow = c(5, 3), mar = c(0, 0, 0.7, 0))
  # S_x.t.min <- min(as.vector(S_x.t))
  # S_x.t.max <- max(as.vector(S_x.t))
  # S_x.t.range <- S_x.t.max - S_x.t.min
  # 
  # c100 <- rev(terrain.colors(100))
  # 
  # for (j in 1:n_step){
  #   cols <- c100[1 + round(100 * (S_x.t[, j] - S_x.t.min)) / S_x.t.range ]
  #   plot(loc_x[,c("x","y")], col = cols, axes = FALSE, asp = 1, pch = 19, cex = 1.5,
  #        main = paste0("Time: ", j))
  # }
  # 
  # S_i <- do.call(c,lapply(1:length(index_i2), function(j){
  #   exp_catch <- S_x.t[index_i2[j],t_i2[j]] # expected catch
  #   return(exp_catch)
  # }))
  # plot(log(y_i),log(S_i))
  
  ###############
  ## Fitting step
  ###############
  ## Build interpolation mesh
  #--------------------------
  library(INLA)
  # domain boundary
  bound <- inla.nonconvex.hull(unique(as.matrix(loc_x[,c('x','y')])),convex=-0.05)
  bound2 <- inla.nonconvex.hull(unique(as.matrix(loc_x[,c('x','y')])),convex=-0.2)
  
  # mesh
  mesh <- inla.mesh.2d(
    boundary=list(bound,bound2),
    max.edge=c(grid_dim/10,grid_dim/5))
  
  # precision matrix components
  spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]  # define sparse matrices for parametrisation of precision matrix
  
  # link between data and mesh (A: projector matrix)
  A <- inla.spde.make.A(mesh,
                        loc=as.matrix(loc_x[index_i2,c("x","y")])) # sampled coordinates
  A <- as( A, "dgTMatrix")
  Aix_ij <- cbind(A@i,A@j)
  Aix_w <- A@x
  
  # link between prediction grid and mesh
  A_pred <- inla.spde.make.A(mesh,
                             loc=as.matrix(loc_x[,c("x","y")])) # sampled coordinates
  A_pred <- as( A_pred, "dgTMatrix")
  Aix_ij_pred <- cbind(A_pred@i,A_pred@j)
  Aix_w_pred <- A_pred@x
  
  ## Save simulated latent field and data
  list_simu.data <- list(n_cells=n_cells, # Domain definition
                         n_step=n_step,
                         grid_dim=grid_dim,
                         
                         range_delta=range_delta, # Latent field
                         range_eta=range_eta,
                         SD_delta=SD_delta,
                         SD_eta=SD_eta,
                         intercept_S=intercept_S,
                         rho_delta=rho_delta,
                         deltaAR1=deltaAR1,
                         S_x.t=S_x.t,
                         
                         n_samp=n_samp, # Sampling process
                         intercept_l=intercept_l,
                         b=b,
                         eta_x.t=eta_x.t,
                         lambda_x.t=lambda_x.t,
                         index_i2=index_i2,
                         t_i2=t_i2,
                         c_x.t=c_x.t,
                         
                         SD_obs=SD_obs, # Observation data
                         y_i=y_i,
                         
                         mesh=mesh, # Triangulation mesh and precision matrix
                         spde=spde,
                         Aix_ij=Aix_ij,
                         Aix_w=Aix_w,
                         Aix_ij_pred=Aix_ij_pred,
                         Aix_w_pred=Aix_w_pred)

  save(file = "data/2_simulation.RData", data = list_simu.data,version = 2)
  
  
}else{

  load("data/2_simulation.RData")
  n_cells <- list_simu.data$n_cells
  n_step <- list_simu.data$n_step
  index_i2 <- list_simu.data$index_i2
  t_i2 <- list_simu.data$t_i2
  y_i <- list_simu.data$y_i
  mesh <- list_simu.data$mesh
  spde <- list_simu.data$spde
  Aix_ij <- list_simu.data$Aix_ij
  Aix_w <- list_simu.data$Aix_w
  Aix_ij_pred <- list_simu.data$Aix_ij_pred
  Aix_w_pred <- list_simu.data$Aix_w_pred
  
}  
  
## Shape inputs (Data, Params, Map, Random) for TMB model
#--------------------------------------------------------
Options_vec <- c("SE"=0, # compute standard error and correction bias (0: no, 1: yes)
                 "RE"=1, # account for random effect (0: no, 1: yes)
                 "SP"=1) # account for sampling process (0:no, 1: yes)

Data <- list("Options_vec"=Options_vec,
             "n_cells"=n_cells, # number of cells in prediction grid
             "n_t"=n_step, 
             "n_nodes"=mesh$n, # number of mesh nodes
             "n_i"=length(y_i),
             "n_gf"=2, # number of random effect (here 2: delta and eta)
             "y_i"=y_i,
             "t_i"=t_i2-1,
             "Aix_ij"=Aix_ij,
             "Aix_w"=Aix_w,
             "Aix_ij_pred"=Aix_ij_pred,
             "Aix_w_pred"=Aix_w_pred,
             "spde"=spde)

Params <- list("intercept_S"=rep(0,n_step),
               "deltainput_x"=matrix(0,nrow=mesh$n,ncol=n_step),
               "intercept_l"=rep(0,n_step),
               "b"=rep(0,n_step),
               "etainput_x"=matrix(0,nrow=mesh$n,ncol=n_step),
               "logtau"=c(0,0),
               "logkappa"=c(0,0),
               "logSD_obs"=log(1),
               "rho_delta"=0)

Map <- list()
if(Options_vec["RE"] == 1){ # estimate random effect
  Random <- c("deltainput_x","etainput_x")
}

if(Options_vec["RE"] == 0){ # do not estimate random effect
  Map[["deltainput_x"]] <- Map[["etainput_x"]] <- factor(matrix(NA,nrow=mesh$n,ncol=n_step))
  Map[["logtau"]] <- Map[["logkappa"]] <- factor(rep(NA,2))
  Random=NULL
}

if(Options_vec["SP"] == 0){ # do not account for sampling process and preferential sampling
  Random <- c("deltainput_x")
  Map[["etainput_x"]] <- factor(matrix(NA,nrow=mesh$n,ncol=n_step))
  Map[["logtau"]] <- Map[["logkappa"]] <- factor(c(1,NA))
  Map[["b"]] <- factor(rep(NA,n_step))
  Map[["intercept_l"]] <- factor(rep(NA,n_step))
}

if( n_step == 1) Map[["rho_delta"]] <- factor(NA)


## Fit TMB model
#---------------
obj <- MakeADFun(data=Data,
                 parameters=Params,
                 random=Random,
                 map = Map,
                 DLL="2_ipp_st_toy_model",silent=TRUE)
runSymbolicAnalysis(obj)
obj$fn( obj$par )

lower_par <- rep(-50,length(obj$par))
upper_par <- rep(50,length(obj$par))
lower_par[which(names(obj$par) == "rho_delta")] <- -1
upper_par[which(names(obj$par) == "rho_delta")] <- 1


system.time(opt <- nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr,
                           lower=lower_par, upper=upper_par,
                           control=list(trace=0, maxit=1000)))
Report <- obj$report()
rep <- sdreport(obj,bias.correct=TRUE)

# # check :
# plot(log(Report$S_p[,1]),log(list_simu.data$S_x.t[,1]))

# ## Plot matrices
# library(Matrix)
# # Hessian
# h <- obj$env$spHess(random=TRUE)
# x11();image(h)
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
# # Cholesky matrix
# L <- obj$env$L.created.by.newton
# image(L)

