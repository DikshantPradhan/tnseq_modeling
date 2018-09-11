# Functions related to sampling
library(sybilcycleFreeFlux)
library(plotrix)

warmup = 1500

ACHR_sampler <- function(model, W=3000, nPnts=5000, steps=1, Floor = FALSE, matrix = TRUE){
  sample = ACHR(model,W,nPoints=nPnts,stepsPerPoint=steps)
  sample = t(sample$Points)
  colnames(sample) <- model@react_id
  sample_df <- as.data.frame(sample)
  if (Floor){
    for (i in 1:nrow(sample_df)){
      for (j in 1:ncol(sample_df)){
        if (abs(sample_df[i,j]) < 1.0e-9){
          sample_df[i,j] <- 0
        }
      }
    }
  }
  if (matrix){
    sample_df <- as.matrix(sample_df)
  }

  return(sample_df)
}

rescale_sample <- function(sample, model, rxn_idx = 0){ # rxn_idx is idx of blocked reaction, 0 if nothing is blocked

  opt <- fluxVar(model, percentage = 99)
  model_fva <- opt@lp_obj
  fva_min <- model_fva[1:95]
  fva_max <- model_fva[96:190]

  rxn_list = c(1:ncol(sample))
  if (rxn_idx != 0){
    rxn_list <- rxn_list[-rxn_idx]
  }

  for(i in rxn_list){ # rescale to flux variability range
    rescaled <- rescale(c(fva_min[i], sample[,i], fva_max[i]), c(-1,1))
    sample[,i] <- rescaled[-c(1, length(rescaled))]
  }

  return(sample)
}

suppressed_model <- function(model, rxn_idx){
  # model@lowbnd[rxn_idx] <- 0
  # model@uppbnd[rxn_idx] <- 0
  if (rxn_idx == 0){
    return(model)
  }

  for (i in rxn_idx){
    model <- changeBounds(model, rxn_idx, lb = 0, ub = 0)
  }

  return(model)
}

maxDiff_dist <- function(sample, model){

  rescaled <- rescale_sample(sample, model = model)
  maxDiff <- array(0, dim = c(nrow(sample)-1, ncol(sample)))

  for (i in 1:ncol(rescaled)){
    sorted <- sort(rescaled[,i])
    max = 0
    for (j in 1:(length(sorted)-1)){
      temp_max <- sorted[j+1] - sorted[j]
      maxDiff[j,i] <- temp_max
    }
    # maxDiff <- c(maxDiff, max)
  }

  colnames(maxDiff) <- colnames(sample)
  return(maxDiff)
}

color_binning <- function(x, n_bins){
  x[is.nan(x)] <- 0
  min <- min(x)
  max <- max(x)
  
  bins <- .bincode(x, breaks = seq(from = min-1, to = max+1, length.out = n_bins+1))
  color_range <- rev(bluered(n_bins))
  colors <- vapply(bins, function(x){color_range[x]}, c('col'))
  
  return(colors)
}

sample_diff <- function(wt_samples, ko_samples, sample_dim = 2){
  wt_mean <- apply(wt_samples, sample_dim, mean)
  wt_sd <- apply(wt_samples, sample_dim, sd)
  ko_mean <- apply(ko_samples, sample_dim, mean)
  ko_sd <- apply(ko_samples, sample_dim, sd)
  
  names <- colnames(wt_samples)
  d_mean <- (ko_mean - wt_mean)
  d_sd <- (ko_sd - wt_sd)
  
  sd_colors <- color_binning(d_sd, n_bins = 10)
  
  ko_change <- data.frame(names = names, d_mean = d_mean, d_sd = d_sd, sd_col = sd_colors)
  
  ggplot(data=ko_change, aes(x = names, y=d_mean)) +
    geom_bar(stat="identity", fill = sd_colors) #+ scale_color_manual(values=sd_colors) + theme_minimal()
}