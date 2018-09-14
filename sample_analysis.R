library(readr)
library(ggplot2)
library(gplots)
library(data.tree)
library(pryr)

# ko_0_samples <- read_csv("~/Documents/jensn lab/tnseq fitting/sampling_data/ko_0_samples.csv", col_names = TRUE)
# ko_5_samples <- read_csv("~/Documents/jensn lab/tnseq fitting/sampling_data/ko_5_samples.csv", col_names = FALSE)
# ko_18_samples <- read_csv("~/Documents/jensn lab/tnseq fitting/sampling_data/ko_18_samples.csv", col_names = FALSE)
# ko_19_samples <- read_csv("~/Documents/jensn lab/tnseq fitting/sampling_data/ko_19_samples.csv", col_names = FALSE)
# ko_33_samples <- read_csv("~/Documents/jensn lab/tnseq fitting/sampling_data/ko_33_samples.csv", col_names = FALSE)
# 
# ko_0_samples <- as.matrix(ko_0_samples)
# ko_5_samples <- as.matrix(ko_5_samples)
# ko_18_samples <- as.matrix(ko_18_samples)
# ko_19_samples <- as.matrix(ko_19_samples)
# ko_33_samples <- as.matrix(ko_33_samples)
# sample_diff(ko_0_samples, ko_18_samples)
# sample_diff(ko_0_samples, ko_19_samples)
# sample_diff(ko_0_samples, ko_33_samples)
ko_0_samples <- read_csv("~/Documents/jensn lab/tnseq fitting/sampling_data/ko_0_samples.csv", col_names = TRUE)
ko_0_samples <- as.matrix(ko_0_samples)
ko_2_samples <- read_csv("~/Documents/jensn lab/tnseq fitting/sampling_data/ko_2_samples.csv", col_names = TRUE)
ko_2_samples <- as.matrix(ko_2_samples)
ko_5_samples <- read_csv("~/Documents/jensn lab/tnseq fitting/sampling_data/ko_5_samples.csv", col_names = TRUE)
ko_5_samples <- as.matrix(ko_5_samples)

rxn_1 <- 'bio00001'
rxn_2 <- 'R00209'
rxn_3 <- 'R00248'
rxn_4 <- 'R04439'
rxn_5 <- 'R04440'

par(mfrow=c(3,5))
hist(ko_0_samples[,grep(rxn_1, colnames(ko_0_samples))])
hist(ko_0_samples[,grep(rxn_2, colnames(ko_0_samples))])
hist(ko_0_samples[,grep(rxn_3, colnames(ko_0_samples))])
hist(ko_0_samples[,grep(rxn_4, colnames(ko_0_samples))])
hist(ko_0_samples[,grep(rxn_5, colnames(ko_0_samples))])
hist(ko_2_samples[,grep(rxn_1, colnames(ko_2_samples))])
hist(ko_2_samples[,grep(rxn_2, colnames(ko_2_samples))])
hist(ko_2_samples[,grep(rxn_3, colnames(ko_2_samples))])
hist(ko_2_samples[,grep(rxn_4, colnames(ko_2_samples))])
hist(ko_2_samples[,grep(rxn_5, colnames(ko_2_samples))])
hist(ko_5_samples[,grep(rxn_1, colnames(ko_5_samples))])
hist(ko_5_samples[,grep(rxn_2, colnames(ko_5_samples))])
hist(ko_5_samples[,grep(rxn_3, colnames(ko_5_samples))])
hist(ko_5_samples[,grep(rxn_4, colnames(ko_5_samples))])
hist(ko_5_samples[,grep(rxn_5, colnames(ko_5_samples))])

ko_0_change <- sample_diff(wt_samples, ko_0_samples)
ko_2_change <- sample_diff(ko_0_samples, ko_2_samples)
ko_5_change <- sample_diff(ko_0_samples, ko_5_samples)

plot_sample_diff(ko_0_change)
plot_sample_diff(ko_2_change)
plot_sample_diff(ko_5_change)
