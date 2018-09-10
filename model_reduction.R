
#' Reduce dimensionality of model; i.e. collapse pathways (R sets) into single reactions
#' @param model original Sybil model for organism metabolism
#' @param gene_sets known gene sets; defaults to empty
#' @param rxn_sets known rxn sets; defaults to empty
#' @return Sybil model expanded to include enzymatic activity
reduce_model_dimensionality <- function(model, coupling, flux){
  
  remove_reactions <- c()
  pathway_ratios <- c()
  
  for (i in 1:length(coupling)){
    set <- coupling[[i]]
    if (length(set) == 1){next}
    rxn_idxs <- vapply(set, function(x){get_rxn_idx(model, x)}, c(1))
    # calculate reduced reaction
    rxn_ratios <- calculate_reaction_ratios(rxn_idxs, flux)
    # add reaction to model
    new_react_name <- paste("pathway_", i, sep = "")
    print(i); print(rxn_ratios)
    model <- calculate_reduced_reaction(model, rxn_idxs, rxn_ratios, new_react_name)
    # mark old reactions to be removed
    remove_reactions <- c(remove_reactions, rxn_idxs)
    pathway_ratios[i] <- list(rxn_ratios)
  }
  
  model <- rmReact(model = model, react = remove_reactions, rm_met = FALSE)
  
  
  list(model = model, pathway_ratios = pathway_ratios)
}

calculate_reaction_ratios <- function(rxn_idxs, flux){
  flux <- flux[,rxn_idxs]
  flux <- flux[which(flux[,1] != 0),]
  flux <- flux/flux[,1]
  
  for (i in 1:ncol(flux)){
    # if (length(unique(as.vector(flux[,i]))) > 1){print('inconsistent'); print(rxn_idxs[i])}
  }
  
  return(flux[1,])
}

calculate_reduced_reaction <- function(model, rxn_idxs, ratios, react_id){
  S <- model@S
  
  rxns_S <- S[,rxn_idxs]
  
  for (i in 1:length(ratios)){
    rxns_S[,i] <- rxns_S[,i]*ratios[i]
  }
  
  new_react <- rowSums(rxns_S)
  
  met_idxs <- which(new_react != 0)
  mets <- model@met_id[met_idxs]
  coeffs <- new_react[met_idxs]
  
  # print(paste(length(mets), length(coeffs)))
  # print(paste(length(model@met_id), length(new_react)))
  
  # rev <- any(model@react_rev[rxn_idxs])
  lbs <- model@lowbnd[rxn_idxs]*ratios
  lb <- max(lbs)
  ubs <- model@uppbnd[rxn_idxs]*ratios
  ub <- min(ubs)
  
  model <- addReact(model = model, id = react_id, met = model@met_id, Scoef = new_react, lb = lb, ub = ub,
                    gprAssoc = concatenate_gprs(rxn_idxs, model@gpr))
  
  return(model)
}

concatenate_gprs <- function(idxs, gpr){
  gprs <- c()
  for (i in idxs){
    new_gpr <- gpr[i]
    if (nchar(new_gpr) == 0){next}
    if (grepl('or', new_gpr)){
      new_gpr <- paste('(', new_gpr, ')', sep = '')
    }
    gprs <- c(gprs, new_gpr)
  }
  
  new_gprs <- paste(gprs, collapse = ' and ')
  return(new_gprs)
}

#' Calculate the fluxes of the reactions in the original model based on fluxes in the reduced model
#' @param flux vector of flux values from reduced model
#' @param sets sets of reactions which were used to reduce the original model
#' @param ratios flux ratios of reactions in each set
#' @param react_id_split string used in naming the reduced pathways which will be used to identify the associated sets
#' @return vector of fluxes for each of the original, un-reduced reactions
extrapolate_from_reduced_flux <- function(flux, sets, ratios, react_id_split = 'pathway_'){
  if ((nrow(flux) != length(sets)) & (length(sets) != length(ratios))){
    print('error in data')
    return()
  }
  
  extrapolate_rxn_flux <- function(set, ratio, flux){
    # convert ratio to matrix and multiply by flux
    output_flux <- matrix(data = ratio, nrow = length(ratio), ncol = 1)*flux
    rownames(output_flux) <- set
    return(output_flux)
  }
  
  rxns <- unlist(sets)
  output_flux <- matrix(data = 0, nrow = length(rxns), ncol = 1)
  rownames(output_flux) <- rxns
  
  pathway_names <- names(flux)
  for (i in 1:length(flux)){
    if (grepl(react_id_split) %in% pathway_names[i]){ # reduced pathway
      set <- sets[[strsplit(pathway_names[i], split = react_id_split)[[1]][2]]]
      new_flux <- extrapolate_rxn_flux(set, ratio[[i]], flux[i])
      
      for (j in 1:nrow(new_flux)){
        output_flux[rownames(new_flux)[j]] <- new_flux[j,1]
      }
      
    }
    else { # single reaction
      output_flux[pathway_names[i]] <- flux[i]
    }
  }
}

#' Calculate the fluxes of the reactions in the reduced model based on fluxes in the original model
#' @param flux vector of flux values from reduced model
#' @param sets sets of reactions which were used to reduce the original model
#' @param ratios flux ratios of reactions in each set
#' @param react_id_split string used in naming the reduced pathways which will be used to identify the associated sets
#' @return vector of fluxes for each of the original, un-reduced reactions
extrapolate_from_reduced_flux <- function(flux, sets, ratios = NULL, react_id_split = 'pathway_'){
  
  output_flux <- matrix(data = 0, nrow = length(sets), ncol = 1)
  path_names <- c()
  # rownames(output_flux) <- 
  
  for (i in 1:length(sets)){
    rxn <- sets[[i]][1]
    
    if (length(sets[[i]]) == 1){
      path_names <- c(path_names, rxn)
    }
    else {
      path_names <- c(path_names, paste(react_id_split, i, sep = ''))
    }
    
    output_flux[i] <- flux[rxn]
  }
  
  rownames(output_flux) <- path_names
  return(output_flux)
}