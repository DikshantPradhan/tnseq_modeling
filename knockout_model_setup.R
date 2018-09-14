#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
objective <- args[2]
knockout_idx <- as.numeric(args[3])
output_file <- args[4]
obj_flux <- as.numeric(args[5])

library(gurobi)
# print(input_file)
model <- gurobi_read(input_file)
obj_idx <- which(model$varnames == objective)
model$obj <- rep(0, length(model$varnames))
model$obj[obj_idx] <- 1

# flux setting
# model$lb[obj_idx] <- obj_flux
# model$ub[obj_idx] <- obj_flux

#knockout
if (knockout_idx != 0){
  model$lb[knockout_idx] <- 0
  model$ub[knockout_idx] <- 0
}

gurobi_write(model, output_file)