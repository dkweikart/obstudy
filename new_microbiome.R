setwd("~/PSU/Obesity Study/Microbiome_code")
library(phyloseq)

OTU.df <- data.frame(read.table(file="otus_new.txt", header=TRUE, row.names=1, sep="\t"))
tax.df <- data.frame(read.table(file="PhylTaxTable.txt", header=TRUE, row.names=1, sep="\t"))
sample <- read.csv(file="Diet.csv", header=TRUE, row.names=1, sep=",")


tax.df$Species <- make.names (ifelse(tax.df$Genus==tax.df$Species,tax.df$Genus,paste(tax.df$Genus, tax.df$Species, sep="_")),  
                              unique=FALSE)
tax.df$Species <- make.names (ifelse(grepl("_",tax.df$Species,fixed=TRUE ), tax.df$Species, paste(tax.df$Species, "unclassified", sep = "_")), 
                              unique = FALSE)

##Create phyloseq object##
OTU.p <- otu_table(OTU.df, taxa_are_rows =TRUE)
TAX.p <- tax_table(as.matrix(tax.df))
sample.p <- sample_data(sample)
physeq <- phyloseq(OTU.p,TAX.p,sample.p)

newPhy = subset_samples(physeq, sample_names(physeq)!="DC006169")

###ANCOM
#load libraries
library(exactRankTests)
library(nlme)
library(dplyr)
library(ggplot2)
library(compositions)
library(phyloseq)
library(data.table)
library(tidyverse)
library(magrittr)

# ANCOM 2.1 from: https://github.com/FrederickHuangLin/ANCOM/blob/master/scripts/ancom_v2.1.R
# Place in working directory
source("ancom_v2.1.R")

##Create filter to 5% prevalence/ greater than 8 cages
physeq_5 <- filter_taxa(newPhy, function(x){sum(x > 0) > 8}, prune = TRUE)

# Compute prevalence of each feature, store as data.frame
prevdf5 <- apply(X = otu_table(physeq_5),
                  MARGIN = ifelse(taxa_are_rows(physeq_5), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf5 <- data.frame(Prevalence = prevdf5,
                       TotalAbundance = taxa_sums(physeq_5),
                       tax_table(physeq_5))

prevalenceThreshold <- 0.05 * nsamples(physeq_5)
keepTaxa2 <- rownames(prevdf5)[(prevdf5$Prevalence >= prevalenceThreshold)]
physeq_5_5 <- prune_taxa(keepTaxa2, physeq_5)

#setup variables for ANCOM
sample_var <- c("SeqNum")
lib_cut = 1000
neg_lb=FALSE   
main_var <- c("Diet")
treats <- c("HIHI","HINO","LOLO","LONO","NOHI","NOLO","NONO")

for (fac in treats){
  #set pair for comparison, change here and output file name
  physeq_test <- subset_samples(physeq_5_5, Diet %in% c("HF", fac))
  
  feature_table <- as.data.frame(otu_table(physeq_test))
  meta_data <- as.data.table(data.frame(sample_data(physeq_test)), keep.rownames = "SeqNum")
  
  # Data Pre-Processing
  processed_data <- feature_table_pre_process(feature_table, meta_data, sample_var, group_var=NULL, 
                                              out_cut = 0.05, zero_cut = 0.95, lib_cut, neg_lb)
  
  #Run ANCOM - could not include random formula for inoculum or else it would crash
  #Didn't seem to be a memory issue, but not sure
  ANCOM_results <- ANCOM(processed_data$feature_table, processed_data$meta_data, 
                         struc_zero = NULL, main_var, p_adj_method = "BH", 
                         alpha = 0.05, adj_formula = NULL, rand_formula = NULL)
  
  #Process results
  results_5_5 <- data.table(ANCOM_results$out)
  results_0.6 <- results_5_5[detected_0.6==TRUE,]
  result_dat <- ANCOM_results$fig$data
  
  theresults <- filter(result_dat, taxa_id %in% results_0.6$taxa_id)
  theresults <- cbind(theresults, results_0.6[,detected_0.9:detected_0.6])
  theresults <- dplyr::rename(theresults, meanCLR=x)
  theresults <- dplyr::rename(theresults, W_score=y)
  
  spec_nam <- as.data.table(as.data.frame(tax_table(physeq_test)), keep.rownames = "OTU")
  theresults <- cbind(taxa_id=theresults$taxa_id, Species=spec_nam[spec_nam$OTU %in% theresults$taxa_id, Species], 
                      theresults[,c(2:8)])
  
  fn <- paste("ANCOM_", fac, ".sig", sep = "")
  fwrite(theresults, file = fn, sep = "\t") 
}


###OTU ANALYSIS DESeq2
library(DESeq2)
library(phyloseq)

#import files and create phyloseq object (see this in ANCOM)
#Filter sequences

##Create filter to 5% prevalence/ greater than 8 cages
physeq_5 <- filter_taxa(newPhy, function(x){sum(x > 0) > 8}, prune = TRUE)


# Compute prevalence of each feature, store as data.frame
prevdf5 <- apply(X = otu_table(physeq_5),
                  MARGIN = ifelse(taxa_are_rows(physeq_5), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf5 <- data.frame(Prevalence = prevdf5,
                       TotalAbundance = taxa_sums(physeq_5),
                       tax_table(physeq_5))

prevalenceThreshold <- 0.05 * nsamples(physeq_5)
keepTaxa2 <- rownames(prevdf5)[(prevdf5$Prevalence >= prevalenceThreshold)]
physeq_5_5 <- prune_taxa(keepTaxa2, physeq_5)

#DESeq2 Analysis###

#convert to factors
sample$Diet <- relevel(as.factor(sample$Diet), ref = "HF")

#Run analysis looking at differences in Diet
OTU5_5.matrix <- as.data.frame(otu_table(physeq_5_5))
DE_Diet_data5_5 <- DESeqDataSetFromMatrix(countData = OTU5_5.matrix, colData = sample[-52,], design = ~Diet)
DE_Diet5_5 <- DESeq(DE_Diet_data5_5)

alpha_cut <- 0.05
#Compare each treatment to water

treats <- factor(factor(c("HIHI","HINO","LOLO","LONO","NOHI","NOLO","NONO")))

for (fac in  treats){
  res5_5 <- results(DE_Diet5_5, contrast = c("Diet", fac, "HF"), cooksCutoff = FALSE)
  sigtab5_5 <- res5_5[which(res5_5$padj < alpha_cut), ]
  sigtab5_5 = cbind(as(sigtab5_5, "data.frame"), as(tax_table(physeq_5_5)[rownames(sigtab5_5), ], "matrix"))
  fn_DESeq <- paste("DESEQ_", fac, ".csv", sep = "")
  write.csv(sigtab5_5,file = fn_DESeq)}


####ANCOM at Genus and Phylum level
physeq_phylum <- tax_glom(physeq_5_5, taxrank = "Phylum")
physeq_gen <- tax_glom(physeq_5_5, taxrank = "Genus")

#load libraries
library(exactRankTests)
library(nlme)
library(dplyr)
library(ggplot2)
library(compositions)
library(phyloseq)
library(data.table)
library(tidyverse)
library(magrittr)

# ANCOM 2.1 from: https://github.com/FrederickHuangLin/ANCOM/blob/master/scripts/ancom_v2.1.R
# Place in working directory
source("ancom_v2.1.R")

# Compute prevalence of each feature, store as data.frame
prevdf5 <- apply(X = otu_table(physeq_phylum),
                 MARGIN = ifelse(taxa_are_rows(physeq_phylum), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf5 <- data.frame(Prevalence = prevdf5,
                      TotalAbundance = taxa_sums(physeq_phylum),
                      tax_table(physeq_phylum))

prevalenceThreshold <- 0.05 * nsamples(physeq_phylum)
keepTaxa2 <- rownames(prevdf5)[(prevdf5$Prevalence >= prevalenceThreshold)]
physeq_5_5 <- prune_taxa(keepTaxa2, physeq_phylum)

#setup variables for ANCOM
sample_var <- c("SeqNum")
lib_cut = 1000
neg_lb=FALSE   
main_var <- c("Diet")
treats <- c("HIHI","HINO","LOLO","LONO","NOHI","NOLO","NONO")

for (fac in treats){
  #set pair for comparison, change here and output file name
  physeq_test <- subset_samples(physeq_5_5, Diet %in% c("HF", fac))
  
  feature_table <- as.data.frame(otu_table(physeq_test))
  meta_data <- as.data.table(data.frame(sample_data(physeq_test)), keep.rownames = "SeqNum")
  
  # Data Pre-Processing
  processed_data <- feature_table_pre_process(feature_table, meta_data, sample_var, group_var=NULL, 
                                              out_cut = 0.05, zero_cut = 0.95, lib_cut, neg_lb)
  
  #Run ANCOM - could not include random formula for inoculum or else it would crash
  #Didn't seem to be a memory issue, but not sure
  ANCOM_results <- ANCOM(processed_data$feature_table, processed_data$meta_data, 
                         struc_zero = NULL, main_var, p_adj_method = "BH", 
                         alpha = 0.05, adj_formula = NULL, rand_formula = NULL)
  
  #Process results
  results_5_5 <- data.table(ANCOM_results$out)
  results_0.6 <- results_5_5[detected_0.6==TRUE,]
  result_dat <- ANCOM_results$fig$data
  
  theresults <- filter(result_dat, taxa_id %in% results_0.6$taxa_id)
  theresults <- cbind(theresults, results_0.6[,detected_0.9:detected_0.6])
  theresults <- dplyr::rename(theresults, meanCLR=x)
  theresults <- dplyr::rename(theresults, W_score=y)
  
  spec_nam <- as.data.table(as.data.frame(tax_table(physeq_test)), keep.rownames = "OTU")
  theresults <- cbind(taxa_id=theresults$taxa_id, Species=spec_nam[spec_nam$OTU %in% theresults$taxa_id, Species], 
                      theresults[,c(2:8)])
  
  fn <- paste("ANCOMphy_", fac, ".sig", sep = "")
  fwrite(theresults, file = fn, sep = "\t") 
}

#GENUS
# Compute prevalence of each feature, store as data.frame
prevdf5 <- apply(X = otu_table(physeq_gen),
                 MARGIN = ifelse(taxa_are_rows(physeq_gen), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf5 <- data.frame(Prevalence = prevdf5,
                      TotalAbundance = taxa_sums(physeq_gen),
                      tax_table(physeq_gen))

prevalenceThreshold <- 0.05 * nsamples(physeq_gen)
keepTaxa2 <- rownames(prevdf5)[(prevdf5$Prevalence >= prevalenceThreshold)]
physeq_5_5 <- prune_taxa(keepTaxa2, physeq_gen)

#setup variables for ANCOM
sample_var <- c("SeqNum")
lib_cut = 1000
neg_lb=FALSE   
main_var <- c("Diet")
treats <- c("HIHI","HINO","LOLO","LONO","NOHI","NOLO","NONO")

for (fac in treats){
  #set pair for comparison, change here and output file name
  physeq_test <- subset_samples(physeq_5_5, Diet %in% c("HF", fac))
  
  feature_table <- as.data.frame(otu_table(physeq_test))
  meta_data <- as.data.table(data.frame(sample_data(physeq_test)), keep.rownames = "SeqNum")
  
  # Data Pre-Processing
  processed_data <- feature_table_pre_process(feature_table, meta_data, sample_var, group_var=NULL, 
                                              out_cut = 0.05, zero_cut = 0.95, lib_cut, neg_lb)
  
  #Run ANCOM - could not include random formula for inoculum or else it would crash
  #Didn't seem to be a memory issue, but not sure
  ANCOM_results <- ANCOM(processed_data$feature_table, processed_data$meta_data, 
                         struc_zero = NULL, main_var, p_adj_method = "BH", 
                         alpha = 0.05, adj_formula = NULL, rand_formula = NULL)
  
  #Process results
  results_5_5 <- data.table(ANCOM_results$out)
  results_0.6 <- results_5_5[detected_0.6==TRUE,]
  result_dat <- ANCOM_results$fig$data
  
  theresults <- filter(result_dat, taxa_id %in% results_0.6$taxa_id)
  theresults <- cbind(theresults, results_0.6[,detected_0.9:detected_0.6])
  theresults <- dplyr::rename(theresults, meanCLR=x)
  theresults <- dplyr::rename(theresults, W_score=y)
  
  spec_nam <- as.data.table(as.data.frame(tax_table(physeq_test)), keep.rownames = "OTU")
  theresults <- cbind(taxa_id=theresults$taxa_id, Species=spec_nam[spec_nam$OTU %in% theresults$taxa_id, Species], 
                      theresults[,c(2:8)])
  
  fn <- paste("ANCOMgen_", fac, ".sig", sep = "")
  fwrite(theresults, file = fn, sep = "\t") 
}

###DESeq2 at phylum and genus level
physeq_phylum <- tax_glom(physeq_5_5, taxrank = "Phylum")
physeq_gen <- tax_glom(physeq_5_5, taxrank = "Genus")

#DESeq2 Analysis###
#convert to factors
sample$Diet <- relevel(as.factor(sample$Diet), ref = "HF")

#Run analysis looking at differences in Diet
OTU5_5.matrix <- as.data.frame(otu_table(physeq_phylum))
DE_Diet_data5_5 <- DESeqDataSetFromMatrix(countData = OTU5_5.matrix, colData = sample[-52,], design = ~Diet)
DE_Diet5_5 <- DESeq(DE_Diet_data5_5)

alpha_cut <- 0.05
#Compare each treatment to water

treats <- factor(factor(c("HIHI","HINO","LOLO","LONO","NOHI","NOLO","NONO")))

for (fac in  treats){
  res5_5 <- results(DE_Diet5_5, contrast = c("Diet", fac, "HF"), cooksCutoff = FALSE)
  sigtab5_5 <- res5_5[which(res5_5$padj < alpha_cut), ]
  sigtab5_5 = cbind(as(sigtab5_5, "data.frame"), as(tax_table(physeq_5_5)[rownames(sigtab5_5), ], "matrix"))
  fn_DESeq <- paste("DESEQphy_", fac, ".csv", sep = "")
  write.csv(sigtab5_5,file = fn_DESeq)}

#genus
#DESeq2 Analysis#
#convert to factors
sample$Diet <- relevel(as.factor(sample$Diet), ref = "HF")

#Run analysis looking at differences in Diet
OTU5_5.matrix <- as.data.frame(otu_table(physeq_gen))
DE_Diet_data5_5 <- DESeqDataSetFromMatrix(countData = OTU5_5.matrix, colData = sample[-52,], design = ~Diet)
DE_Diet5_5 <- DESeq(DE_Diet_data5_5)

alpha_cut <- 0.05
#Compare each treatment to water

treats <- factor(factor(c("HIHI","HINO","LOLO","LONO","NOHI","NOLO","NONO")))

for (fac in  treats){
  res5_5 <- results(DE_Diet5_5, contrast = c("Diet", fac, "HF"), cooksCutoff = FALSE)
  sigtab5_5 <- res5_5[which(res5_5$padj < alpha_cut), ]
  sigtab5_5 = cbind(as(sigtab5_5, "data.frame"), as(tax_table(physeq_5_5)[rownames(sigtab5_5), ], "matrix"))
  fn_DESeq <- paste("DESEQgen_", fac, ".csv", sep = "")
  write.csv(sigtab5_5,file = fn_DESeq)}











####HEATMAP
#RELATIVE ABUNDACE Species
library(phyloseq)
library(data.table)
library(dplyr)
library(ggplot2)

#Create phyloseq object (see above)
##Create filter to 5% prevalence/ 8 samples
physeq_5 <- filter_taxa(newPhy, function(x){sum(x > 0) > 8}, prune = TRUE)


# Compute prevalence of each feature, store as data.frame
prevdf5 <- apply(X = otu_table(physeq_5),
                 MARGIN = ifelse(taxa_are_rows(physeq_5), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf5 <- data.frame(Prevalence = prevdf5,
                      TotalAbundance = taxa_sums(physeq_5),
                      tax_table(physeq_5))

prevalenceThreshold <- 0.05 * nsamples(physeq_5)
keepTaxa2 <- rownames(prevdf5)[(prevdf5$Prevalence >= prevalenceThreshold)]
physeq_5_5 <- prune_taxa(keepTaxa2, physeq_5)

#Get the average relative abundance for each of the diet conditions for each member 
#of the taxonomic level, in this case phylum - can be changed to other
physeq_5_5 <- tax_glom(physeq_5_5, taxrank = "Species")
physeq_5_5 <- transform_sample_counts(physeq_5_5, function(x) x/sum(x))

physeq_HF <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "HF")))
physeq_NONO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "NONO")))
physeq_LONO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "LONO")))
physeq_HINO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "HINO")))
physeq_NOLO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "NOLO")))
physeq_NOHI <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "NOHI")))
physeq_LOLO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "LOLO")))
physeq_HIHI <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "HIHI")))


#update Phylum here if using other taxonomic level
Diet.pct <- as.data.table(data.frame(OTU=rownames(tax_table(physeq_5_5)), Phylum=tax_table(physeq_5_5)[,"Species"], 
                                     HF = round(rowMeans(physeq_HF)*100,3), 
                                     HIHI = round(rowMeans(physeq_HIHI)*100,3), HINO = round(rowMeans(physeq_HINO)*100,3), 
                                     LOLO = round(rowMeans(physeq_LOLO)*100,3), LONO = round(rowMeans(physeq_LONO)*100,3),
                                     NOHI = round(rowMeans(physeq_NOHI)*100,3), NOLO = round(rowMeans(physeq_NOLO)*100,3),
                                     NONO = round(rowMeans(physeq_NONO)*100,3)))

fwrite(Diet.pct, file="percent_abund_species.txt", sep = "\t")

###Percent Abundance Genus
#Create phyloseq object (see above)
##Create filter to have minimum 5 counts (0.001%) and 5% prevalence
physeq_5 <- filter_taxa(newPhy, function(x){sum(x > 0) > 8}, prune = TRUE)


# Compute prevalence of each feature, store as data.frame
prevdf5 <- apply(X = otu_table(physeq_5),
                 MARGIN = ifelse(taxa_are_rows(physeq_5), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf5 <- data.frame(Prevalence = prevdf5,
                      TotalAbundance = taxa_sums(physeq_5),
                      tax_table(physeq_5))

prevalenceThreshold <- 0.05 * nsamples(physeq_5)
keepTaxa2 <- rownames(prevdf5)[(prevdf5$Prevalence >= prevalenceThreshold)]
physeq_5_5 <- prune_taxa(keepTaxa2, physeq_5)

#Get the average relative abundance for each of the diet conditions for each member 
#of the taxonomic level, in this case phylum - can be changed to other
physeq_5_5 <- tax_glom(physeq_5_5, taxrank = "Genus")
physeq_5_5 <- transform_sample_counts(physeq_5_5, function(x) x/sum(x))

physeq_HF <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "HF")))
physeq_NONO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "NONO")))
physeq_LONO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "LONO")))
physeq_HINO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "HINO")))
physeq_NOLO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "NOLO")))
physeq_NOHI <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "NOHI")))
physeq_LOLO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "LOLO")))
physeq_HIHI <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "HIHI")))


#update Phylum here if using other taxonomic level
Diet.pct <- as.data.table(data.frame(OTU=rownames(tax_table(physeq_5_5)), Phylum=tax_table(physeq_5_5)[,"Genus"], 
                                     HF = round(rowMeans(physeq_HF)*100,3), 
                                     HIHI = round(rowMeans(physeq_HIHI)*100,3), HINO = round(rowMeans(physeq_HINO)*100,3), 
                                     LOLO = round(rowMeans(physeq_LOLO)*100,3), LONO = round(rowMeans(physeq_LONO)*100,3),
                                     NOHI = round(rowMeans(physeq_NOHI)*100,3), NOLO = round(rowMeans(physeq_NOLO)*100,3),
                                     NONO = round(rowMeans(physeq_NONO)*100,3)))

fwrite(Diet.pct, file="percent_abundance_genus.txt", sep = "\t")

###Percent Abundance Phylum
#Create phyloseq object (see above)
##Create filter to have minimum 5 counts (0.00%) and 5% prevalence
physeq_5 <- filter_taxa(newPhy, function(x){sum(x > 0) > 8}, prune = TRUE)


# Compute prevalence of each feature, store as data.frame
prevdf5 <- apply(X = otu_table(physeq_5),
                 MARGIN = ifelse(taxa_are_rows(physeq_5), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf5 <- data.frame(Prevalence = prevdf5,
                      TotalAbundance = taxa_sums(physeq_5),
                      tax_table(physeq_5))

prevalenceThreshold <- 0.05 * nsamples(physeq_5)
keepTaxa2 <- rownames(prevdf5)[(prevdf5$Prevalence >= prevalenceThreshold)]
physeq_5_5 <- prune_taxa(keepTaxa2, physeq_5)

#Get the average relative abundance for each of the diet conditions for each member 
#of the taxonomic level, in this case phylum - can be changed to other
physeq_5_5 <- tax_glom(physeq_5_5, taxrank = "Phylum")
physeq_5_5 <- transform_sample_counts(physeq_5_5, function(x) x/sum(x))

physeq_HF <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "HF")))
physeq_NONO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "NONO")))
physeq_LONO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "LONO")))
physeq_HINO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "HINO")))
physeq_NOLO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "NOLO")))
physeq_NOHI <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "NOHI")))
physeq_LOLO <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "LOLO")))
physeq_HIHI <- data.frame (otu_table(subset_samples(physeq_5_5, Diet == "HIHI")))


#update Phylum here if using other taxonomic level
Diet.pct <- as.data.table(data.frame(OTU=rownames(tax_table(physeq_5_5)), Phylum=tax_table(physeq_5_5)[,"Phylum"], 
                                     HF = round(rowMeans(physeq_HF)*100,3), 
                                     HIHI = round(rowMeans(physeq_HIHI)*100,3), HINO = round(rowMeans(physeq_HINO)*100,3), 
                                     LOLO = round(rowMeans(physeq_LOLO)*100,3), LONO = round(rowMeans(physeq_LONO)*100,3),
                                     NOHI = round(rowMeans(physeq_NOHI)*100,3), NOLO = round(rowMeans(physeq_NOLO)*100,3),
                                     NONO = round(rowMeans(physeq_NONO)*100,3)))

fwrite(Diet.pct, file="percent_abundance_phylum.txt", sep = "\t")



####HEATMAP 2
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)

#Import data
da_tab <- fread("Diff_abund.csv")
rel_phy <- fread("percent_abundance_phylum.txt")
rel_gen <- fread("percent_abundance_genus.txt")
rel_sp <- fread("percent_abund_species.txt")

#Get all of the significant relative abundances
rel_phy_sig <- rel_phy[Phylum %in% da_tab$Name,-1]
gen_names <- da_tab[Level=="Genus", "Name"]
rel_gen_sig <- rel_gen[Genus %in% gen_names$Name,-1]
sp_names <- da_tab[Level=="Species", "OTU"]
rel_sp_sig <- rel_sp[OTU %in% sp_names$OTU,]

#Reform columns and sorting to match between levels and then join
rel_sp_sig$Name <- paste(rel_sp_sig$OTU, rel_sp_sig$Species, sep = "_")
rel_sp_sig <- rel_sp_sig[,-c(1:2)]
rel_sp_sig <- cbind("Name"=rel_sp_sig$Name, rel_sp_sig[,c(1:8)])

rel_phy_sig <- dplyr::rename(rel_phy_sig, Name=Phylum)
rel_gen_sig <- dplyr::rename(rel_gen_sig, Name=Genus)

rel_phy_sig <- setorder(rel_phy_sig, Name)
rel_gen_sig <- setorder(rel_gen_sig, Name)

rel_all <- rbind(rel_gen_sig, rel_phy_sig, rel_sp_sig)

#Calculate normalized values for differences between treatments and water
rel_nums <- rel_all[,-c(1:2)]
HF <- rel_all$HF
diff_nums <- rel_nums - HF
norm_nums <- diff_nums/apply(abs(diff_nums), 1, max)

#Prepare objects for plotting
plot_vals <- cbind(Name=rel_all$Name, Level = da_tab$Level, norm_nums)
plot_sigs <- cbind(Name=rel_all$Name, da_tab[,c(3:10)])
colnames(plot_sigs) <- colnames(plot_vals)

plot_vals$Level <- factor(plot_vals$Level, levels = c("Phylum", "Genus", "Species"))
plot_vals <- plot_vals[order(plot_vals$Level, plot_vals$Name, decreasing = TRUE)]
ordernames <- plot_vals$Name

plot_vals <- data.table::melt(plot_vals, variable.name = "Treatment", id.vars = c("Name", "Level"), value.name = "Rel_Change")
plot_sigs <- data.table::melt(plot_sigs, variable.name = "Treatment", id.vars = c("Name", "Level"), value.name = "Sig")

plot_vals$Level <- factor(plot_vals$Level, levels = c("Phylum", "Genus", "Species"))
plot_sigs$Level <- factor(plot_sigs$Level, levels = c("Phylum", "Genus", "Species"))
treat_ord <- c("HIHI","HINO","LOLO","LONO", "NOHI","NOLO","NONO")
plot_vals$Treatment <- factor(plot_vals$Treatment, levels = treat_ord)
plot_sigs$Treatment <- factor(plot_sigs$Treatment, levels = treat_ord)

plot_vals <- plot_vals[order(plot_vals$Level, plot_vals$Name, plot_vals$Treatment, decreasing = TRUE)]
plot_sigs <- plot_sigs[order(plot_sigs$Level, plot_sigs$Name, plot_sigs$Treatment, decreasing = TRUE)]

plot_vals$Name <- factor(plot_vals$Name, levels = ordernames)
plot_sigs$Name <- factor(plot_sigs$Name, levels = ordernames)

plot_sigs <- plot_sigs[is.na(plot_sigs$Sig)==FALSE,]

sig_dat <- data.frame(Treatment = plot_sigs$Treatment, Name = plot_sigs$Name)
Lab_txt <- c("HIHI","HINO","LOLO","LONO", "NOHI","NOLO","NONO")

tiff("differential abundance.tiff", units="in", width=7, height=6, res=300)
ggplot(plot_vals, aes(x = Treatment, y = Name)) +
  geom_tile(aes(fill = Rel_Change)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  scale_x_discrete(label=Lab_txt) + 
  labs(y = "Taxa") +
  geom_text(data = sig_dat, label = plot_sigs$Sig)
dev.off()
