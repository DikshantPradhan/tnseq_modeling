## data parsing for falcon model
library(readr)
library(readxl)
library(dplyr)

parse_brenda_kcat <- function(filename, organism){
  full_data <- read_tsv(filename, col_names = FALSE)
  organisms_list <- full_data$X3
  idxs <- grep(organism, organisms_list)
  kcats <- data.frame(EC = full_data$X1[idxs], enzyme = full_data$X2[idxs], kcat = full_data$X4[idxs])
}

parse_swiss_prot <- function(filename){
  full_data <- read_csv(filename, col_names = TRUE)
  data <- data.frame(gene = full_data$kegg3, EC = full_data$kegg4)
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
  
  source <- c()
  
  for (i in 1:n){
    if (!is.na(EC_rxn$KCAT[i])){next}
    EC <- EC_rxn$`EC Number`[i] #strsplit(EC_rxn$`EC Number`[i], split = 'EC')[[1]][2]
    EC <- paste(c(EC, '$'), collapse = '')
    
    flag <- FALSE
    phylo_idx <- 1
    
    while (!flag){
      new_source <- phylogeny[phylo_idx]
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
  
  kcat_data <- data.frame(rxn = EC_rxn$Abbreviation, kcat = EC_rxn$KCAT, source = source)
  
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

concatenate_kcat_data <- function(kcat_filename, organisms, concatenate_by = 'enzyme', remove = c('EC')){
  if (length(organisms) < 2){
    print('single organism provided')
    return()
  }
  
  data <- parse_brenda_kcat(kcat_filename, organisms[1])
  data <- data[, !(names(data) %in% remove)]
  
  for (i in 2:length(organisms)){
    new_data <- parse_brenda_kcat(kcat_filename, organisms[i])
    new_data <- new_data[, !(names(new_data) %in% remove)]
    data <- left_join(data, new_data, by = concatenate_by)
  }
  
  names(data) <- c(concatenate_by, organisms)
  
  return(data)
}

# phylogeny <- c('mutans', 'streptococcus', 'streptococcaceae', 'lactobacillales', 'bacilli', 'firmicutes', 'bacteria', 'prokaryotes')

streptococcus_list <- c('streptococcus pneumoniae', 'streptococcus mutans', 'streptococcus pyogenes', 'streptococcus oralis',
                        'streptococcus agalactiae', 'streptococcus equi', 'streptococcus thermophilus', 'streptococcus gordonii',
                        'streptococcus suis')

pneumoniae_data <- parse_brenda_kcat('~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', 'streptococcus pneumoniae')# streptococcus pneumoniae
mutans_data <- parse_brenda_kcat('~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', 'streptococcus mutans')# streptococcus mutans
# pyogenes_data <- parse_brenda_kcat('~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', 'streptococcus pyogenes')# streptococcus pyogenes
# oralis_data <- parse_brenda_kcat('~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', 'streptococcus oralis')# streptococcus oralis
# agalactiae_data <- parse_brenda_kcat('~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', 'streptococcus agalactiae')# streptococcus agalactiae
# equi_data <- parse_brenda_kcat('~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', 'streptococcus equi')# streptococcus equi
# thermophilus_data <- parse_brenda_kcat('~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', 'streptococcus thermophilus')# streptococcus thermophilus
# # streptococcus sp.
# gordonii_data <- parse_brenda_kcat('~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', 'streptococcus gordonii')# streptococcus gordonii
# suis_data <- parse_brenda_kcat('~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', 'streptococcus suis')# streptococcus suis
# # streptococcus dysgalactiae subsp. equisimilis
# 
# EC <- unique(c(pneumoniae_data$EC, mutans_data$EC, pyogenes_data$EC, oralis_data$EC, agalactiae_data$EC, equi_data$EC, thermophilus_data$EC,
#         gordonii_data$EC, suis_data$EC))

strep_data <- concatenate_kcat_data('~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', 
                                    organisms = streptococcus_list)

streptomyces_list <- c('streptomyces lividans', 'streptomyces exfoliatus', 'streptomyces ambofaciens', 'streptomyces avermitilis',
                       'streptomyces coelicolor', 'streptomyces griseus', 'streptomyces pristinaespiralis', 'streptomyces viridifaciens',
                       'streptomyces hygroscopicus', 'streptomyces wedmorensis', 'streptomyces glaucescens', 'streptomyces antibioticus',
                       'streptomyces griseorubens', 'streptomyces griseosporeus', 'streptomyces lavendulae', 'streptomyces rochei',
                       'streptomyces toyocaensis', 'streptomyces nogalater', 'streptomyces clavuligerus', 'streptomyces vinaceus', 
                       'streptomyces sahachiroi', 'streptomyces nodosus', 'streptomyces tubercidicus', 'streptomyces globisporus',
                       'streptomyces coeruleorubidus', 'streptomyces exfoliatus', 'streptomyces toxytricini', 'streptomyces rugosporus',
                       'streptomyces noursei', 'streptomyces cinnamonensis', 'streptomyces albus', 'streptomyces seoulensis', 
                       'streptomyces venezuelae', 'streptomyces fradiae', 'streptomyces niveus', 'streptomyces sahachiroi',
                       'streptomyces carzinostaticus', 'streptomyces davaonensis', 'streptomyces collinus', 'streptomyces mobaraensis',
                       'streptomyces diastaticus', 'streptomyces kanamyceticus', 'streptomyces niveus', 'streptomyces galilaeus', 
                       'streptomyces scabiei', 'streptomyces anulatus')
streptomyces_list <- unique(streptomyces_list)

# streptomyces lividans
# streptomyces exfoliatus
# streptomyces ambofaciens
# streptomyces avermitilis
# streptomyces coelicolor
# streptomyces griseus
# streptomyces pristinaespiralis
# streptomyces viridifaciens
# streptomyces sp.
# streptomyces hygroscopicus
# streptomyces wedmorensis
# streptomyces glaucescens
# streptomyces antibioticus
# streptomyces griseorubens
# streptomyces griseosporeus
# streptomyces lavendulae
# streptomyces rochei
# streptomyces toyocaensis
# streptomyces nogalater
# streptomyces clavuligerus
# streptomyces vinaceus
# streptomyces sahachiroi
# streptomyces nodosus
# streptomyces tubercidicus
# streptomyces globisporus
# streptomyces coeruleorubidus
# streptomyces exfoliatus
# streptomyces toxytricini
# streptomyces rugosporus
# streptomyces noursei
# streptomyces cinnamonensis
# streptomyces albus
# streptomyces seoulensis
# streptomyces venezuelae
# streptomyces fradiae
# streptomyces niveus
# streptomyces sahachiroi
# streptomyces carzinostaticus
# streptomyces davaonensis
# streptomyces collinus
# streptomyces mobaraensis
# streptomyces diastaticus
# streptococcus dysgalactiae subsp. equisimilis
# streptomyces kanamyceticus
# streptomyces niveus
# streptomyces galilaeus
# streptomyces scabiei
# streptomyces roseochromogenus subsp. oscitans
# streptomyces anulatus
streptomyces_data <- concatenate_kcat_data('~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', 
                                           organisms = streptomyces_list)

num_entries <- c()
for (i in 1:20){
  num_entries[i] <- length(which(!is.na(streptomyces_data[i,])))
}

species <- which(!is.na(streptomyces_data[9,]))
species <- species[-c(1)]

reduced_streptomyces_data <- streptomyces_data[3:9,species]
d <- dist(t(reduced_streptomyces_data), method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 
plot(fit) # display dendogram

# mutans_kcat <- combine_EC_kcat("~/Documents/jensn lab/gecko/iSMUv01_rxn_EC.xlsx", '~/Documents/jensn lab/gecko/brenda_max/max_KCAT.txt', phylogeny)
# mutans_falcon <- generate_falcon_model(mutans, kcat_data = mutans_kcat, enzyme_constraint = TRUE)
