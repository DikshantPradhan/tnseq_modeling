
#' Reduce dimensionality of model; i.e. collapse pathways (R sets) into single reactions
#' @param model original Sybil model for organism metabolism
#' @param gene_sets known gene sets; defaults to empty
#' @param rxn_sets known rxn sets; defaults to empty
#' @return Sybil model expanded to include enzymatic activity
reduce_model_dimensionality <- function(model, coupling, flux){
  
  remove_reactions <- c()
  
  for (set in coupling){
    rxn_idxs <- vapply(set, function(x){get_set_idx(x)}, c(1))
    # calculate reduced reaction
    rxn_fluxes <- calculate_reaction_fluxes
    # add reaction to model
    model <- calculate_reduced_reaction
    # mark old reactions to be removed
    remove_reactions <- c(remove_reactions, rxn_idxs)
  }
  
  model <- rmReact(model = model, react = remove_reactions, rm_met = FALSE)
  
  return(model)
}

calculate_reaction_ratios <- function(rxn_idxs, flux){
  flux <- flux[,rxn_idxs]
  flux <- flux/flux[,1]
  
  for (i in 1:ncol(flux)){
    if (!all.equal(flux[,i])){print('inconsistent')}
  }
  
  return(flux[1,])
}

calculate_reduced_reactions <- function(rxn_idxs, ratios, react_id){
  S <- model@S
  
  rxns_S <- S[rxn_idxs,]*ratios
  new_react <- colSums(rxns_S)
  
  mets <- model@mets[which(new_react != 0)]
  coeffs <- new_react[which(new_react != 0)]
  
  model <- addReact(model, id = react_id, met = mets, Scoef = coeffs, reversible = any(model@rev[rxn_idxs]), lb = -1000, ub = 1000)
  
  return(model)
}