library(sybil)
source('~/GitHub/tnseq_modeling/falcon_tools.R')
source('~/GitHub/tnseq_modeling/sampling_tools.R')
library(gurobi)

# wt_sample_file <- "~/Documents/jensn lab/tnseq fitting/sampling_data/ko_0_samples.csv"
# wt_samples <- read_csv(wt_sample_file, col_names = TRUE)

# A_wt <- most_average_sample(wt_samples)

flux_diff <- function(w, A_wt, A){
  diff <- sum((A_wt - A)^2/w)
  return(diff)
}

# optimize_flux <- function(model, wt_flux, w){
#   init <- 
#   result <- optim()
# }

generate_parsimonious_fba_model <- function(model, minimization_idxs, fixed_idxs, fixed_vals){
  if (length(fixed_idxs) != length(fixed_vals)){
    return()
  }
  n <- length(model$varnames)
  
  # set Q matrix
  q_mtx <- matrix(data = 0, nrow = n, ncol = n)
  for (i in minimization_idxs){
    q_mtx[i,i] <- 1
  }
  model$Q <- q_mtx
  model$modelsense <- 'min'
  
  # fix reaction fluxes
  for (i in 1:length(fixed_idxs)){
    idx <- fixed_idxs[i]
    val <- fixed_vals[i]
    model$lb[idx] <- val
    model$ub[idx] <- val
  }
  
  return(model)
}

generate_moma_model <- function(model, wt_fluxes, idxs, weights){
  if (length(idxs) != length(weights)){return()}
  
  n <- length(model$varnames)
  
  q_mtx <- matrix(data = 0, nrow = n, ncol = n)
  model$modelsense <- 'min'
  
  for (i in 1:length(idxs)){
    idx <- idxs[i]
    weight <- weights[i]
    wt_flux <- wt_fluxes[idx]
    
    q_mtx[idx,idx] <- weight
    model$obj[idx] <- -2*weight*wt_flux
  }
  
  model$Q <- q_mtx
  return(model)
}

moma_optimization_wrapper <- function(model, wt_fluxes, idxs, weights, obj_idx){
  print('calculating ko fitness...')
  moma_model <- generate_moma_model(model, wt_fluxes, idxs, weights)
  ko_fitness <- matrix(data = 0, nrow = 1, ncol = length(idxs))
  for (i in 1:length(idxs)){
    # print(i)
    idx <- idxs[i]
    model <- moma_model
    model$lb[idx] <- 0
    model$ub[idx] <- 0
    results <- gurobi(model, params)
    if (!is.null(results$x)){
      ko_fitness[i] <- results$x[obj_idx]
    }
  }
  return(ko_fitness)
}

fitness_optimization <- function(model, wt_fluxes, idxs, weights, obj_idx, fit_obs){
  fit_pred <- moma_optimization_wrapper(model, wt_fluxes, gene_idxs, sample_weights, obj_idx)
  ss <- sum((fit_pred - fit_obs)^2)
  return(ss)
}

# make wt model 
params <- list(OutputFlag = 0)
model <- gurobi_read('GitHub/tnseq_modeling/data/ecoli_falcon.lp')
model$obj <- rep(0, length(model$obj))
vars <- model$varnames
biom_idx <- grep('Bio', vars)
model$obj[biom_idx] <- 1
gene_idxs <- grep('Ex_a_', vars)

# generate wt flux
result <- gurobi(model, params)
obj_flux <- result$objval
model$obj <- rep(0, length(model$obj))
parsim_fba_model <- generate_parsimonious_fba_model(model, gene_idxs, fixed_idxs = biom_idx, fixed_vals = obj_flux) # alter biomass
results <- gurobi(parsim_fba_model, params)
wt_fluxes <- results$x

# sample random weights
sample_weights <- sample((1:10)/10, length(gene_idxs), replace = TRUE)

# run all KOs

ko_fitness <- moma_optimization_wrapper(model, wt_fluxes, gene_idxs, sample_weights, biom_idx)

# moma_model <- generate_moma_model(model, wt_fluxes, gene_idxs, sample_weights)
# ko_fitness <- matrix(data = 0, nrow = 1, ncol = length(gene_idxs))
# for (i in 1:length(gene_idxs)){
#   idx <- gene_idxs[i]
#   model <- moma_model
#   model$lb[idx] <- 0
#   model$ub[idx] <- 0
#   results <- gurobi(model, params)
#   ko_fitness[i] <- results$x[biom_idx]
# }

# optimize 

result <- optim(par = rep(1, length(gene_idxs)), fitness_optimization, model = moma_model, wt_fluxes = wt_fluxes, idxs = gene_idxs, fit_obs = ko_fitness)