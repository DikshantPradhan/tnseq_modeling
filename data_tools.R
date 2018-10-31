## data parsing for falcon model
library(readr)
library(readxl)

parse_brenda_kcat <- function(filename, organism){
  full_data <- read_tsv(filename, col_names = FALSE)
  organisms_list <- full_data$X3
  idxs <- grep(organism, organisms_list)
  kcats <- list(EC = full_data$X1[idxs], enzyme = full_data$X2[idxs], kcat = full_data$X4[idxs])
}

parse_swiss_prot <- function(filename){
  full_data <- read_csv(filename, col_names = TRUE)
  data <- list(gene = full_data$kegg3, EC = full_data$kegg4)
  return(data)
}

combine_kcat_prot <- function(kcat, prot, organism){
  len <- length(kcat$kcat)
  
  genes <- c()
  for (i in 1:len){
    ec <- strsplit(kcat$EC[i],split = 'EC')[[1]][2]
    print(ec)
    gene_idx <- grep(ec, prot$EC)
    print(gene_idx)
    if (length(gene_idx) == 0){next}
    gene <- prot$gene[gene_idx]
    # print(gene)
    # gene <- tolower(paste(strsplit(gene, '_')[[1]], collapse = ''))
    genes[i] <- list(gene) # how to deal with multiple genes???
  }
  
  kcat$genes <- genes
  return(kcat)
}

combine_EC_kcat <- function(EC_rxn_filename, kcat_filename, phylogeny = c()){
  EC_rxn <- read_excel(EC_rxn_filename)
  n <- length(EC_rxn$Abbreviation)
  
  # clear unknowns
  clear_str <- c('-', 'no', 'No')
  for (str in clear_str){
    idxs <- which(grepl(str,EC_rxn$`EC Number`))
    EC_rxn$KCAT[idxs] <- 0
  }
  idxs <- which(is.na(EC_rxn$`EC Number`))
  EC_rxn$KCAT[idxs] <- 0
  
  for (i in 1:n){
    if (!is.na(EC_rxn$KCAT[i])){next}
    EC <- EC_rxn$`EC Number`[i] #strsplit(EC_rxn$`EC Number`[i], split = 'EC')[[1]][2]
    EC <- paste(c(EC, '$'), collapse = '')
    
    flag <- FALSE
    phylo_idx <- 1
    
    while (!flag){
      print(phylogeny[phylo_idx])
      kcat <- 0
      if (phylo_idx > length(phylogeny)){
        EC_rxn$KCAT[i] <- kcat
        flag <- TRUE
        next
      }
      kcat_data <- parse_brenda_kcat(kcat_filename, phylogeny[phylo_idx])
      
      idx <- grep(EC, kcat_data$EC)
      if (length(idx) == 1){
        kcat <- kcat_data$kcat[idx]
        EC_rxn$KCAT[i] <- kcat
        flag <- TRUE
      }
      if (length(idx) > 1){
        kcat <- max(kcat_data$kcat[idx])
        EC_rxn$KCAT[i] <- kcat
        flag <- TRUE
      }
      if (length(idx) == 0){
        phylo_idx <- phylo_idx + 1
      }
    }
    
  }
  
  kcat_data <- list(rxn = EC_rxn$Abbreviation, kcat = EC_rxn$KCAT)
  
  return(kcat_data)
}

kcat_from_rxn <- function(kcat_data, rxn){
  idx <- grep(rxn, kcat_data$rxn)
  if (length(idx) < 1){return(1)}
  if (length(idx) > 1){
    idx <- idx[1]
  }
  kcat <- kcat_data$kcat[idx]
  if (length(kcat) < 1){return(1)}
  if (kcat == 0){return(1)}
  return(kcat)
}

phylogeny <- c('mutans', 'streptococcus', 'streptococcaceae', 'lactobacillales', 'bacilli', 'firmicutes', 'bacteria', 'prokaryotes')

# data <- parse_brenda_kcat('~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', 'mutans')
# strep_kcat <- parse_brenda_kcat('~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', 'streptococcus')
# mutans_kcat <- parse_brenda_kcat('~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', 'bacteria')
# prot <- parse_swiss_prot('~/databases/table.txt')
# data <- combine_kcat_prot(kcat, prot, 'ECDH10')

# mutans_kcat <- combine_EC_kcat("~/Documents/jensn lab/gecko/iSMUv01_rxn_EC.xlsx", '~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', phylogeny)
mutans_falcon <- generate_falcon_model(mutans, kcat_data = mutans_kcat, enzyme_constraint = TRUE)
# for (gene in Ec_core@allGenes){
#   num <- strsplit(gene, split = 'b')[[1]][2]
#   print(grep(num, prot$`Gene names`))
# }
