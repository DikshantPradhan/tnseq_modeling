

#' Fit energy requirements to data
#' @return Sybil model
fit_energy_cost <- function(model, energy_met_idx, ko_fitness, obs_fitness){
  new_e_cost <- function(e_old, df_obs, df_pred){
    e_new <- e_old*(1 - df_obs)/(1 - df_pred)
    return(e_new)
  }
  
  S <- model@S
  e <- S[,energy_met_idx]
}

