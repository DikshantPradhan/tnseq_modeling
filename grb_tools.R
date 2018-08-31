## MODELS

#' description
#' @return GRB model
#' @seealso 
#' @export
#' @examples
#' 
GRB_ecoli_model <- function(){
  data(Ec_core);
  model=Ec_core;
  model <- changeBounds(model, 11, lb = 0)
  # model <- changeBounds(model, 13, lb = 0, ub = 0)
  model <- rmReact(model = model, react = 13)
  for (i in findExchReact(model)@react_pos){
    model <- changeBounds(model, i, lb = -1000, ub = 1000)
    # if (model@lowbnd[i] == 0){
    #   model <- changeBounds(model, i, lb = -1000)
    # }
  }
  
  ecoli <- as_GRBmodel(model)
  ecoli$show_output(FALSE)
  return(ecoli)
}

#' description
#' @return GRB falcon model
#' @seealso 
#' @export
#' @examples
#' 
GRB_ecoli_falcon_model <- function(){
  sybil_ecoli <- get_ecoli_model()
  ecoli_falcon_model <- GRB_generate_falcon_model(sybil_ecoli)
  return(ecoli_falcon_model)
}

#' Gurobi S. cerevisiae model
#' @param model filename of RData file for sybil model; defaults to 'yeast_model.RData'
#' @return GRB model
#' @seealso 
#' @export
#' @examples
#' 
GRB_yeast_model <- function(model = 'yeast_model.RData'){
  #maranas_model_exch_add_biom_rm/maranas_model_lipid_exch
  load(model)
  
  # setwd("~/GitHub/PathwayMining/")
  
  yeast <- as_GRBmodel(yeast_model)
  yeast$show_output(FALSE)
  
  return(yeast)
}

#' Gurobi S. cerevisiae falcon model
#' @param model filename of RData file for sybil model; defaults to 'yeast_model.RData'
#' @return GRB falcon model
#' @seealso 
#' @export
#' @examples
#' 
GRB_yeast_falcon_model <- function(model = 'yeast_model.RData'){
  
  load(model)
  
  sybil_yeast <- yeast_model
  yeast_falcon_model <- GRB_generate_falcon_model(sybil_yeast)
  return(yeast_falcon_model)
}

#' Gurobi S. mutans model
#' @param model filename of RData file for sybil model; defaults to 'mutans_model.RData'
#' @return GRB model
#' @seealso 
#' @export
#' @examples
#' 
GRB_mutans_model <- function(model = 'mutans_model.RData'){
  
  load(model)
  
  mutans <- as_GRBmodel(mutans)
  mutans$show_output(FALSE)
  
  return(mutans)
}

#' Gurobi S. mutans falcon model
#' @param model filename of RData file for sybil model; defaults to 'mutans_model.RData'
#' @return GRB falcon model
#' @seealso 
#' @export
#' @examples
#' 
GRB_mutans_falcon_model <- function(model = 'mutans_model.RData'){
  load(model)
  
  sybil_mutans <- mutans
  mutans_falcon_model <- GRB_generate_falcon_model(sybil_mutans)
  return(mutans_falcon_model)
}

#' Gurobi P. aeruginosa model
#' @param model filename of RData file for sybil model; defaults to 'pao_model.RData'
#' @return GRB model
#' @seealso 
#' @export
#' @examples
#' 
GRB_pao_model <- function(model = 'pao_model.RData'){
  load(model)
  
  pao <- as_GRBmodel(pao_model)
  pao$show_output(FALSE)
  
  return(pao)
}

#' Gurobi P. aeruginosa falcon model
#' @param model filename of RData file for sybil model; defaults to 'pao_model.RData'
#' @return GRB falcon model
#' @seealso 
#' @export
#' @examples
#' 
GRB_pao_falcon_model <- function(model = 'pao_model.RData'){
  load(model)
  
  pao_falcon <- GRB_generate_falcon_model(pao_model)
  pao_falcon$show_output(FALSE)
  
  return(pao_falcon)
}

## FUNCTIONS

#' Generate Gurobi falcon model. First generates model in sybil then converts the model to Gurobi in order to add additional
#' constraints on flux
#' @import grb
#' @param sybil_model base model to build falcon model from
#' @param falcon_model boolean value indicating whether or not the base model has enzyme activity already. If TRUE, then the function will simply add contstraints to the falcon model
#' @param r0_gene_set list of genes passed in to the falcon-model building method (generate_falcon_model)
#' @param r0_rxn_set_list list of reaction sets passed in to the falcon-model building method (generate_falcon_model)
#' @return Gurobi model based on metabolic model with enzymatic acctivity
#' @seealso 
#' @export
#' @examples
#' 
GRB_generate_falcon_model <- function(sybil_model, falcon_model = FALSE, r0_gene_set = c(), r0_rxn_set_list = c()){
  
  sybil_falcon_model <- sybil_model
  
  if (!falcon_model){ # if user has not passed in a falcon model already
    sybil_falcon_model <- generate_falcon_model(sybil_model, r0_gene_set, r0_rxn_set_list)
  }
  
  grb_falcon_model <- as_GRBmodel(sybil_falcon_model)
  grb_falcon_model$show_output(FALSE)
  
  ## ADD NECESSARY CONSTRAINTS TO MODEL
  vars <- grb_falcon_model$get_names()$VarName
  
  split_fwd_rxns <- vars[grep('fwd', vars)]
  split_rev_rxns <- vars[grep('rev', vars)]
  
  split_rxns <- sapply(split_fwd_rxns, function(x) strsplit(x, ' ')[[1]][1])
  split_rxns <- unique(split_rxns)
  
  if (length(split_fwd_rxns) != length(split_rev_rxns)){
    print('fed rev matchup error')
    return()
  }
  
  for (rxn in split_rxns){
    #a_rxn <- paste('a', rxn, sep = '_')
    I <- paste('I', rxn, sep = '_')
    .Call("GRB_addvar", grb_falcon_model$exptr, 0L, integer(0), numeric(0), 1.0, 0.0, 1.0, 'B', I, PACKAGE = 'grb')
  }
  
  .Call("GRB_updatemodel", grb_falcon_model$exptr, PACKAGE = 'grb')
  vars <- grb_falcon_model$get_names()$VarName
  n <- grb_falcon_model$get_sizes()$NumVars
  
  # add bounds on split conversion reactions (gene -> activity_[rxn])
  for (i in 1:length(split_fwd_rxns)){
    fwd <- split_fwd_rxns[i]
    rev <- split_rev_rxns[i]
    
    # get bounds (care about fwd_ub & rev_lb)
    fwd_ub <- grb_falcon_model$getattr("UB")[[fwd]]
    fwd_lb <- grb_falcon_model$getattr("LB")[[fwd]]
    rev_ub <- grb_falcon_model$getattr("UB")[[rev]]
    rev_lb <- grb_falcon_model$getattr("LB")[[rev]]
    
    a_rxn <- strsplit(fwd, ' ')[[1]][1]
    I <- paste('I', a_rxn, sep = '_')
    i_idx <- which(vars == I)
    fwd_idx <- which(vars == fwd)
    rev_idx <- which(vars == rev)
    
    #print(paste(fwd, fwd_ub, fwd_lb, rev, rev_ub, rev_lb, I))
    fwd_vec <- c(fwd_ub, -1)
    rev_vec <- c(-1*rev_ub, -1)
    fwd_idxs <- c(i_idx, which(vars == fwd))
    rev_idxs <- c(i_idx, which(vars == rev))
    
    fwd_name <- paste(fwd, 'I', sep = ' ')
    rev_name <- paste(rev, 'I', sep = ' ')
    
    #.Call("GRB_updatemodel", grb_falcon_model$exptr) ##
    .Call("GRB_addconstr", grb_falcon_model$exptr, 2L, as.integer(fwd_idxs-1), fwd_vec, ">=", 0.0,
          fwd_name, PACKAGE = 'grb')
    #.Call("GRB_addconstr", grb_falcon_model$exptr, 2L, as.integer(rev_idxs-1), rev_vec, "<=", (-1*rev_lb),
    #  rev_name)
    .Call("GRB_addconstr", grb_falcon_model$exptr, 2L, as.integer(rev_idxs-1), rev_vec, ">=", (-1*rev_ub),
          rev_name, PACKAGE = 'grb')
    .Call("GRB_updatemodel", grb_falcon_model$exptr, PACKAGE = 'grb')
    #grb_falcon_model$addconstr(paste(fwd, '*', I, sep = ''), sense="<=", rhs= fwd_ub, name = paste(a_rxn, 'fwd', sep = '_')) # bound on fwd conversion
    #grb_falcon_model$addconstr(paste(rev, '*(1 - ', I, ')',  sep = ''), sense=">=", rhs= rev_lb, name = paste(a_rxn, 'rev', sep = '_')) # bound on rev conversion
    #.Call("GRB_updatemodel", grb_falcon_model$exptr)
  }
  
  # for each reaction w fwd and rev components:
  #   add constraint so that only one can run at a time
  #   binary vars I_fwd, I_rev <- {0, 1} and I_fwd + I_rev = 1
  
  .Call("GRB_updatemodel", grb_falcon_model$exptr, PACKAGE = 'grb')
  return(grb_falcon_model)
}

#' get the position of a reaction along its axis on the S matrix
#' @param model Gurobi model
#' @param rxn reaction id to get the index of
#' @return integer indicating the index of the reaction along the S matrix
#' @seealso 
#' @export
#' @examples
#' 
GRB_get_rxn_idx <- function(model, rxn){
  vars <- model$get_names()$VarName
  return(return(which(vars == rxn)))
}

#' wrapper function for flux coupling
#' @param i integer indicating index of the reaction to suppress
#' @param vars vector of reaction ids
#' @param model_og Gurobi model
#' @param reaction_indexes integers indicating reactions to run flux coupling over; defaults to 1:length(vars)
#' @param compare_mtx boolean to indicate whether or not to couple known sets together if any individual reactions are coupled; FALSE by default
#' @param init_coupling_mtx boolean matrix of known coupling between reactions; empty by default
#' @param file_output filename to output coupling vector to
#' @return list of integers which indicate positions of TRUE values in couling matrix
#' @seealso 
#' @examples
#' 
GRB_flux_coupling_raptor_wrapper <- function(i, vars, model_og, reaction_indexes = 1:length(vars), compare_mtx = FALSE, init_coupling_mtx = c(), file_output = NULL){
  print(paste('suppression index:', i))
  
  model <- model_og$copy()
  
  # block i
  model$setattr("UB", setNames(0, vars[i]))
  model$setattr("LB", setNames(0, vars[i]))
  
  coupling_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes, compare_mtx = compare_mtx, known_set_mtx = init_coupling_mtx, stored_obs = 4000)$coupled
  
  coupling_idxs <- which(coupling_mtx)
  if (!is.null(file_output)){
    write(paste(c(i,coupling_idxs), collapse = ','),file = file_output,append=TRUE)
  }
  return(coupling_idxs)
}

#' Function for PACT method
#' @param model_og model upon which to run PACT
#' @param suppression_idxs list of integers indicating indexes of reactions to iteratively block and calculate coupling for. DEFAULT is -1, which indicates that all reactions should be suppressed
#' @param reaction_indexes list of integers indicating indexes of reaction to check coupling for. DEFAULT is empty list, which indicates that all reactions in the model should be considered
#' @param compare_known_init_sets boolean indicating whether or not to calculate R0 sets in order to run further optimizations during flux coupling and PACT. DEFAULT: FALSE
#' @param optimize_suppr boolean indicating whether or not to optimize the reactions to suppress based on R0 sets. DEFAULT: FALSE
#' @param optimize_rxns boolean indicating whether or not to optimize the reactions couple. DEFAULT: FALSE
#' @param cores integer indicating the number of cores to run on, if intending to parallelize
#' @param avoid_idxs list of integers indicating indexes to specifically avoid suppressing during PACT method
#' @param file_output filename to output coupling vector to
#' @return a coupling list: list of list of integers which each indicate the couplings induced by each suppression. The index position of each list indicates the index of the reaction which was suppressed and the integers in eacch list indicate the positions of TRUE values in the coupling array
#' @seealso 
#' @export
#' @examples
#' 
GRB_generate_set_lists_cluster <- function(model_og, suppression_idxs = -1, reaction_indexes = c(),
                                           compare_known_init_sets = FALSE, optimize_suppr = FALSE, optimize_rxns = FALSE, cores = 1, avoid_idxs = c(), file_output = NULL){
  
  n <- model_og$get_sizes()$NumVars
  vars <- model_og$get_names()$VarName
  
  if (suppression_idxs[1] == -1){
    if (length(reaction_indexes) > 0){
      suppression_idxs = reaction_indexes
    }
    else {
      suppression_idxs = 1:n
    }
  }
  
  # dim: rxns_row, rxns_col, deletions
  model <- model_og$copy()
  
  suppr_vector <- Matrix(data = FALSE, nrow = 1, ncol = n, sparse = TRUE)
  suppr_vector[suppression_idxs] <- TRUE
  
  #init_coupling_mtx <- c()
  if (compare_known_init_sets){
    init_coupling_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes)$coupled
    init_coupling_mtx <- fill_coupling_matrix(init_coupling_mtx)
  }
  
  if (compare_known_init_sets & optimize_suppr){
    i <- 1
    while (i <= n){
      if (suppr_vector[i]){ # if tagged to be suppressed
        set_idx <- which(init_coupling_mtx[,i])[1] # which is first reaction (row) i is coupled to
        if (!is.na(set_idx)){
          rxn_idxs <- which(init_coupling_mtx[set_idx,]) # other reactions in set
          # only suppress first reaction in set since, theoretically, suppressing any should have the same effect
          suppr_vector[rxn_idxs] <- FALSE
          suppr_vector[rxn_idxs[1]] <- TRUE
        }
        else {
          suppr_vector[i] <- FALSE
        }
        
      }
      i <- i+1
    }
    
    if (optimize_rxns){
      reaction_indexes <- which(suppr_vector[1,])
    }
  }
  
  if (length(avoid_idxs) > 0){
    suppr_vector[avoid_idxs] <- FALSE
  }
  
  print(paste("# of suppressions:", length(which(suppr_vector[1,])), sep = " "))
  coupling <- mclapply(which(suppr_vector[1,]), function(x) GRB_flux_coupling_raptor_wrapper(x, vars, model_og, reaction_indexes = reaction_indexes, compare_mtx = compare_known_init_sets, init_coupling_mtx = init_coupling_mtx, file_output = file_output), mc.cores = cores)
  
  return(coupling)
}

#' Convert the output of GRB_generate_set_lists_cluster to a coupling matrix. G sets or R sets are fully coupled to each other in G* or R* sets
#' @param coupling_list list of list of integers which each indicate the couplings induced by each suppression in PACT
#' @param n_react number of reactions
#' @param vars names of the reactions; length should equal n_react
#' @param init_sets known sets to indicate coupling for in the matrix
#' @return boolean matrix indicating coupling between reactions
#' @seealso 
#' @export
#' @examples
#' 
full_coupling_matrix_from_coupling_vector_list <- function(coupling_list, n_react, vars, init_sets = NULL){
  
  coupling_matrix <- Matrix(data = FALSE, nrow = n_react, ncol = n_react)
  rownames(coupling_matrix) <- vars
  colnames(coupling_matrix) <- vars
  if (!is.null(init_sets)){
    coupling_matrix <- fill_coupling_matrix_from_sets(coupling_matrix, init_sets)
  }
  
  for (i in 1:length(coupling_list)){
    if (is.null(coupling_list[[i]])){next}
    coupling_matrix[coupling_list[[i]]] <- TRUE
  }
  
  coupling_matrix <- fill_coupling_matrix(coupling_matrix)
  
  return(coupling_matrix)
}

#' Convert the output of GRB_generate_set_lists_cluster to a coupling matrix
#' @param coupling_list list of list of integers which each indicate the couplings induced by each suppression in PACT
#' @param n_react number of reactions
#' @param vars names of the reactions; length should equal n_react
#' @param init_sets known sets to indicate coupling for in the matrix
#' @return boolean matrix indicating coupling between reactions
#' @seealso 
#' @export
#' @examples
#' 
coupling_matrix_from_coupling_vector_list <- function(coupling_list, n_react, vars, init_sets){
  
  coupling_matrix <- Matrix(data = FALSE, nrow = n_react, ncol = n_react)
  rownames(coupling_matrix) <- vars
  colnames(coupling_matrix) <- vars
  for (i in 1:length(coupling_list)){
    if (is.null(coupling_list[[i]])){next}
    print(i)
    intermediate_mtx <- Matrix(data = FALSE, nrow = n_react, ncol = n_react)
    rownames(intermediate_mtx) <- vars
    colnames(intermediate_mtx) <- vars
    intermediate_mtx[coupling_list[[i]]] <- TRUE
    intermediate_mtx <- fill_coupling_matrix_from_sets(intermediate_mtx, init_sets)
    coupling_matrix[which(intermediate_mtx)] <- TRUE
  }
  
  return(coupling_matrix)
}

#' identify which G sets and R sets within the same G* or R* set never couple to each other
#' @param full_coupling_mtx boolean matrix indicating fully coupled G* or R* sets
#' @param coupling_mtx  boolean matrix indicating only observed couplings
#' @param n_react number of reactions
#' @return boolean matrix with TRUE values indicating which reactions in the same G* or R* set never couple
#' @seealso 
#' @export
#' @examples
#' 
identify_intermediate_uncoupled <- function(full_coupling_mtx, coupling_mtx, n_react){
  uncoupled_matrix <- Matrix(data = FALSE, nrow = n_react, ncol = n_react)
  rownames(uncoupled_matrix) <- rownames(coupling_mtx)
  colnames(uncoupled_matrix) <- colnames(coupling_mtx)
  
  uncoupled <- which(full_coupling_mtx & !coupling_mtx)
  uncoupled_matrix[uncoupled] <- TRUE
  
  return(uncoupled_matrix)
}

#' flux-balance optimization function for Gurobi
#' @param model_og Gurobi model
#' @param obj integer indicating the index of the objective reaction
#' @param suppress list of integers indicating which reactions to suppress before optimizing
#' @param max boolean indicating whether to maximize or minimize the objective. TRUE indicates maximization, FALSE indicates minimization. DEFAULT: TRUE
#' @return
#' @seealso 
#' @export
#' @examples
#' 
GRB_maximize <- function(model_og, obj, suppress = c(), max = TRUE){ # suppress is characters
  model <- model_og$copy()
  
  n <- model$get_sizes()$NumVars
  vars <- model$get_names()$VarName
  # clear obj
  model$setattr("Obj", setNames(rep(0.0, times = n), vars))
  
  # set suppressions
  if (length(suppress) > 0){
    suppr_idxs <- which(vars %in% suppress)
    model$setattr("UB", setNames(rep(0.0, times = length(suppr_idxs)), vars[suppr_idxs]))
    model$setattr("LB", setNames(rep(0.0, times = length(suppr_idxs)), vars[suppr_idxs]))
  }
  
  # set obj
  model$setattr("Obj", setNames(1.0, vars[obj]))
  model$set_model_sense(maximize=TRUE)
  if (!max){
    model$set_model_sense(minimize=TRUE)
  }
  model$optimize()
  sol <- model$get_solution()
  obj_max <- sol$ObjVal
  return(obj_max)
}
