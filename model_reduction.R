
#' Reduce dimensionality of model; i.e. collapse pathways (R sets) into single reactions
#' @param model original Sybil model for organism metabolism
#' @param gene_sets known gene sets; defaults to empty
#' @param rxn_sets known rxn sets; defaults to empty
#' @return Sybil model expanded to include enzymatic activity
reduce_model_dimensionality <- function(model, coupling, flux){
  
  remove_reactions <- c()
  
  for (i in 1:length(coupling)){
    set <- coupling[[i]]
    if (length(set) == 1){next}
    rxn_idxs <- vapply(set, function(x){get_rxn_idx(model, x)}, c(1))
    # calculate reduced reaction
    rxn_ratios <- calculate_reaction_ratios(rxn_idxs, flux)
    # add reaction to model
    new_react_name <- paste("r_set_", i, sep = "")
    print(new_react_name)
    model <- calculate_reduced_reaction(model, rxn_idxs, rxn_ratios, new_react_name)
    # mark old reactions to be removed
    remove_reactions <- c(remove_reactions, rxn_idxs)
  }
  
  model <- rmReact(model = model, react = remove_reactions, rm_met = FALSE)
  
  return(model)
}

calculate_reaction_ratios <- function(rxn_idxs, flux){
  flux <- flux[,rxn_idxs]
  flux <- flux[which(flux[,1] != 0),]
  flux <- flux/flux[,1]
  
  for (i in 1:ncol(flux)){
    if (length(unique(flux[,i])) != 1){print('inconsistent'); print(rxn_idxs[i])}
  }
  
  return(flux[1,])
}

calculate_reduced_reaction <- function(model, rxn_idxs, ratios, react_id){
  S <- model@S
  
  rxns_S <- S[,rxn_idxs]*ratios
  new_react <- rowSums(rxns_S)
  
  met_idxs <- which(new_react != 0)
  mets <- model@met_id[met_idxs]
  coeffs <- new_react[met_idxs]
  
  print(paste(length(mets), length(coeffs)))
  print(paste(length(model@met_id), length(new_react)))
  
  rev <- any(model@react_rev[rxn_idxs])
  lb <- 0
  if (rev){lb <- -1000}
  
  model <- addReact(model = model, id = react_id, met = model@met_id, Scoef = new_react, reversible = rev, lb = lb, ub = 1000)
  
  return(model)
}