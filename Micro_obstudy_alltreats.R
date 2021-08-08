#load in files
library(phyloseq)
OTU.df <- data.frame(read.table(file="otus_new", header=TRUE, row.names=1, sep="\t"))
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

#Remove Mock
newPhy = subset_samples(physeq, sample_names(physeq)!="DC006169")

###ANCOM all treatments compared 
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

##Create filter to have minimum 51 counts (0.001%) and 5% prevalence
physeq_51 <- filter_taxa(newPhy, function(x) sum(x) >=51, prune=TRUE)

# Compute prevalence of each feature, store as data.frame
prevdf51 <- apply(X = otu_table(physeq_51),
                  MARGIN = ifelse(taxa_are_rows(physeq_51), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf51 <- data.frame(Prevalence = prevdf51,
                       TotalAbundance = taxa_sums(physeq_51),
                       tax_table(physeq_51))

prevalenceThreshold <- 0.05 * nsamples(physeq_51)
keepTaxa2 <- rownames(prevdf51)[(prevdf51$Prevalence >= prevalenceThreshold)]
physeq_51_5 <- prune_taxa(keepTaxa2, physeq_51)

#setup variables for ANCOM
sample_var <- c("SeqNum")
lib_cut = 1000
neg_lb=FALSE   
main_var <- c("Diet")
treats <- c("NO/NO", "LO/NO", "HI/NO", "NO/LO", "NO/HI", "LO/LO", 
            "HI/HI")

for (fac in treats){
  #set pair for comparison, change here and output file name
  physeq_test <- subset_samples(physeq_51_5, Diet %in% c("HF", fac))
  
  feature_table <- as.data.frame(otu_table(physeq_test))
  meta_data <- as.data.table(data.frame(sample_data(physeq_test)), keep.rownames = "SeqNum")
  
  # Data Pre-Processing
  processed_data <- feature_table_pre_process(feature_table, meta_data, sample_var, group_var=NULL, 
                                              out_cut = 0.05, zero_cut = 0.95, lib_cut, neg_lb)
  
  #Run ANCOM 
  ANCOM_results <- ANCOM(processed_data$feature_table, processed_data$meta_data, 
                         struc_zero = NULL, main_var, p_adj_method = "BH", 
                         alpha = 0.05, adj_formula = NULL, rand_formula = NULL)
  
  #Process results
  results_51_5 <- data.table(ANCOM_results$out)
  results_0.6 <- results_51_5[detected_0.6==TRUE,]
  result_dat <- ANCOM_results$fig$data
  
  theresults <- filter(result_dat, taxa_id %in% results_0.6$taxa_id)
  theresults <- cbind(theresults, results_0.6[,detected_0.9:detected_0.6])
  theresults <- dplyr::rename(theresults, meanCLR=x)
  theresults <- dplyr::rename(theresults, W_score=y)
  
  spec_nam <- as.data.table(as.data.frame(tax_table(physeq_test)), keep.rownames = "OTU")
  theresults <- cbind(taxa_id=theresults$taxa_id, Species=spec_nam[spec_nam$OTU %in% theresults$taxa_id, Species], 
                      theresults[,c(2:8)])
  
  fn <- paste("ANCOM_", fac, ".sig", sep = "")
  fwrite(theresults, file = "fn", sep = "\t") 
}

###DESEQ2 
#All treatments 
library(DESeq2)
library(phyloseq)

##Create filter to have minimum 51 counts (0.001%) and 5% prevalence
physeq_51 <- filter_taxa(physeq, function(x) sum(x) >=51, prune=TRUE)

# Compute prevalence of each feature, store as data.frame
prevdf51 <- apply(X = otu_table(physeq_51),
                  MARGIN = ifelse(taxa_are_rows(physeq_51), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf51 <- data.frame(Prevalence = prevdf51,
                       TotalAbundance = taxa_sums(physeq_51),
                       tax_table(physeq_51))

prevalenceThreshold <- 0.05 * nsamples(physeq_51)
keepTaxa2 <- rownames(prevdf51)[(prevdf51$Prevalence >= prevalenceThreshold)]
physeq_51_5 <- prune_taxa(keepTaxa2, physeq_51)

#DESeq2 Analysis###

#convert to factors
sample$Diet <- relevel(as.factor(sample$Diet), ref = "HF")

#Run analysis looking at differences in Diet
OTU51_5.matrix <- as.data.frame(otu_table(physeq_51_5))
DE_Diet_data51_5 <- DESeqDataSetFromMatrix(countData = OTU51_5.matrix, colData = sample, design = ~Diet)
DE_Diet51_5 <- DESeq(DE_Diet_data51_5)

alpha_cut <- 0.05
#Compare each treatment to water

treats <- factor(factor(c("NO/NO", "LO/NO", "HI/NO", "NO/LO", "NO/HI", "LO/LO", 
                          "HI/HI")))

for (fac in  treats){
  res51_5 <- results(DE_Diet51_5, contrast = c("Diet", fac, "HF"), cooksCutoff = FALSE)
  sigtab51_5 <- res51_5[which(res51_5$padj < alpha_cut), ]
  sigtab51_5 = cbind(as(sigtab51_5, "data.frame"), as(tax_table(physeq_51_5)[rownames(sigtab51_5), ], "matrix"))
  fn_DESeq <- paste("DESEQ_", fac, ".csv", sep = "")
  write.csv(sigtab51_5,file = "fn_DESeq") 
}

####PERMANOVA
library(phyloseq)
library(microbiome)
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
# Pick relative abundances (compositional) and sample metadata
pseq.rel <- microbiome::transform(newPhy, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

#Visualize
p <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "Diet", size = 3)
print(p)

p <- plot_landscape(pseq.rel, method = "PCoA", distance = "bray", col = "Diet", size = 3)
print(p)

f <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "Ferm", size = 3)
print(f)

r <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "Roast", size = 3)
print(r)

NMDS <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "Sex", size = 3)
print(NMDS)

pCoA <- plot_landscape(pseq.rel, method = "PCoA", distance = "bray", col = "Sex", size = 3)
print(pCoA)

# samples x species as input
library(vegan)
permanova <- adonis(t(otu) ~ Diet,
                    data = meta, permutations=99, method = "bray")

# P-value
print(as.data.frame(permanova$aov.tab)["Diet", "Pr(>F)"])

permanova.ferm <- adonis(t(otu) ~ Ferm,
                         data = meta, permutations=99, method = "bray")
print(as.data.frame(permanova.ferm$aov.tab)["Ferm", "Pr(>F)"])

permanova.roast <- adonis(t(otu) ~ Roast,
                          data = meta, permutations=99, method = "bray")
print(as.data.frame(permanova.roast$aov.tab)["Roast", "Pr(>F)"])

permanova.sex <- adonis(t(otu) ~ Sex,
                        data = meta, permutations=99, method = "bray")
print(as.data.frame(permanova.sex$aov.tab)["Sex", "Pr(>F)"])

#Check homogeniety
# Note the assumption of similar multivariate spread among the groups
# ie. analogous to variance homogeneity
# Here the groups have signif. different spreads and
# permanova result may be potentially explained by that.
dist <- vegdist(t(otu))
anova(betadisper(dist, meta$Diet))

anova(betadisper(dist, meta$Ferm))

anova(betadisper(dist, meta$Roast))

anova(betadisper(dist, meta$Sex))

permutest(betadisper(dist, meta$Diet), pairwise = TRUE)
permutest(betadisper(dist, meta$Ferm), pairwise = TRUE)
permutest(betadisper(dist, meta$Roast), pairwise = TRUE)
permutest(betadisper(dist, meta$Sex), pairwise = TRUE)

###Alpha Diversity
library(phyloseq)
library(DivNet)
library(breakaway)
library(lme4)
library(dplyr)
library(magrittr)
library(picante)
library(data.table)
library(ggsignif)
library(ggplot2)
library(multcomp)
library(tibble)
library(pairwiseAdonis)
library(DescTools)
library(mctoolsr)

#Create phyloseq object
physeq_gen <- tax_glom(physeq, taxrank = "Genus")
otu_gen <- as.data.table(otu_table(physeq_gen), keep.rownames = "OTU")
tax_gen <- as.data.table(as.data.frame(tax_table(physeq_gen)), keep.rownames = "OTU")

fwrite(otu_gen, file= "otus_new_gen", col.name = TRUE, sep = "\t")
fwrite(tax_gen, file= "tax_new_gen.txt", col.name = FALSE, sep = "\t")

OTU.matrix <- data.matrix(read.table(file="otus_new_gen", header=TRUE, row.names=1, sep="\t"))
tax.matrix <- as.matrix(read.table(file="tax_new_gen.txt", header=TRUE, row.names=1, sep="\t"))
sample <- read.csv(file="Diet.csv", header=TRUE, row.names=1, sep=",")

#Create phyloseq object
OTU.p <- otu_table(OTU.matrix, taxa_are_rows =TRUE)
TAX.p <- tax_table(tax.matrix)
sample.p <- sample_data(sample)
tree.p <- read_tree(treefile="FF2.tree")
physeq_gen <- phyloseq(OTU.p,TAX.p,sample.p, tree.p)


##Run DivNet diversity tests##
a_div <- physeq_gen %>% divnet(ncores = 6)
estimates <- a_div$shannon %>% summary %$% estimate
ses <- a_div$shannon %>% summary %$% error
sample2 <- sample
sample2$Diet <- relevel(as.factor(sample2$Diet), ref="HF")

dv <- DivNet::divnet(physeq_gen, X = NULL)



#calculate significance accounting for repeated measures
set.seed(3488)
sample2 <- sample
a_div_test <- betta_random(chats = estimates,
                           ses = ses,
                           X = model.matrix(~Diet, data = sample2),
                           groups=sample2$Sample)

#correct p-values via fdr
adjusted.shannon <- as.data.frame(a_div_test$table)
adjusted.shannon <- dplyr::rename(adjusted.shannon, pvalues='p-values')
adjusted.shannon$qvalues<- round(p.adjust(adjusted.shannon$pvalues, "BH"),3)
adjusted.shannon$Sig <- NA
adjusted.shannon$Sig[adjusted.shannon$qvalues<0.05] <- "*"
adjusted.shannon$Sig[adjusted.shannon$qvalues<0.01] <- "**"
adjusted.shannon$Sig[adjusted.shannon$qvalues<0.001] <- "***"

write ("\nGenus Level DivNet Treatment Shannon Test \n", file="ANOVA_results.txt", append=TRUE)
capture.output(adjusted.shannon, file="ANOVA_results_Shannon.txt", append=TRUE)
adjusted.shannon
sample2$Shannon<-estimates

#Repeat procedure for simpson diversity
estimates2 <- a_div$simpson %>% summary %$% estimate
ses2 <- a_div$simpson %>% summary %$% error
a_div_test2 <- betta_random(chats = estimates2,
                            ses = ses2,
                            X = model.matrix(~Diet, data = sample2),
                            groups=sample$Sample)

adjusted.simpson <- as.data.frame(a_div_test2$table)
adjusted.simpson <- dplyr::rename(adjusted.simpson, pvalues='p-values')
adjusted.simpson$qvalues<- round(p.adjust(adjusted.simpson$pvalues, "BH"),3)
adjusted.simpson$Sig <- NA
adjusted.simpson$Sig[adjusted.simpson$qvalues<0.05] <- "*"
adjusted.simpson$Sig[adjusted.simpson$qvalues<0.01] <- "**"
adjusted.simpson$Sig[adjusted.simpson$qvalues<0.001] <- "***"

write ("\nGenus Level DivNet Treatment Simmpson Test \n", file="ANOVA_results.txt", append=TRUE)
capture.output(adjusted.simpson, file="ANOVA_results_Simpson.txt", append=TRUE)
adjusted.simpson
sample2$Simpson <- estimates2

##Faith's Phylogenetic Diversity##
set.seed(3488)
physeq_gen.rar <- rarefy_even_depth(physeq_gen, replace = FALSE)
OTU_gen.rar <- as.data.table(otu_table(physeq_gen.rar), keep.rownames = "OTU")
OTU_gen.rar.t <- data.table::transpose(OTU_gen.rar, keep.names = "Group", make.names = "OTU")
Faith_d.rar <- picante::pd(OTU_gen.rar.t, tree.p, include.root = FALSE)
sample2$Faith <- Faith_d.rar$PD

#Fit model with random effects
fit.faith <- lmerTest::lmer(Faith ~ Diet + (1 | Sample), data = sample2)
faith_test <- glht(fit.faith, linfct = mcp(Diet = "Dunnett"))

write("Genus Level lmer Treatment and Faith's Diversity Test \n", file= "ANOVA_results.txt", append=TRUE)
capture.output(summary(faith_test, test=adjusted("BH")), file="ANOVA_results_Faith.txt", append=TRUE)

#Create object to aid plotting
faith_sum <- summary(faith_test, test=adjusted("BH"))
adjusted.faith <- data.frame(qvalues=faith_sum$test$pvalues) 
adjusted.faith$Sig<- NA
adjusted.faith$Sig[adjusted.faith$qvalues<0.05] <- "*"
adjusted.faith$Sig[adjusted.faith$qvalues<0.01] <- "**"
adjusted.faith$Sig[adjusted.faith$qvalues<0.001] <- "***"

#use inverse simpson for actual plotting
inv_simp <- 1/estimates2
sample2$InvSimp <- inv_simp

#Create labels for plots and data frames
Treatments_abrev <- c("HF","NO/NO", "LO/NO", "HI/NO", "NO/LO", "NO/HI",
                      "LO/LO", "HI/HI", "MOCK")
Treatments_full <- c("HF","NO/NO", "LO/NO", "HI/NO", "NO/LO", "NO/HI",
                     "LO/LO", "HI/HI", "MOCK")
Treatments_test_full <- Treatments_full[-1]
Treatments_test_abrev <- Treatments_abrev[-1]

#Create objects for plotting significance of Shannon values
adjusted.shannon$max_shannon <- tapply(sample2$Shannon, sample2$Diet, max)
adjusted.shannon$coord <- adjusted.shannon$max_shannon +0.25
Shannon_sig <- adjusted.shannon[-1,]
rownames(Shannon_sig) <- Treatments_test_full
Shannon_sig <- Shannon_sig[is.na(Shannon_sig$Sig)==FALSE, ]
Shan_label <- data.frame(Treatment = rownames(Shannon_sig), Shannon = Shannon_sig$coord)

#Create objects for plotting significance of Simpson values
adjusted.simpson$max_simpson <- tapply(sample2$Simpson, sample2$Diet, max)
adjusted.simpson$coord <- adjusted.simpson$max_simpson +0.1
Simpson_sig <- adjusted.simpson[-1,]
rownames(Simpson_sig) <- Treatments_test_full
Simpson_sig <- Simpson_sig[is.na(Simpson_sig$Sig)==FALSE, ]
Simp_label <- data.frame(Treatment = rownames(Simpson_sig), Simpson = Simpson_sig$coord)

#Create objects for plotting Inverse Simpson
adjusted.invsimp <- data.frame(qvalues = adjusted.simpson$qvalues, Sig = adjusted.simpson$Sig, row.names = Treatments_full)
adjusted.invsimp$max_invsimpson <- tapply(sample2$InvSimp, sample2$Diet, max)
adjusted.invsimp$coord <- adjusted.invsimp$max_invsimpson + 1
InvSimpson_sig <- adjusted.invsimp[-1,]
InvSimpson_sig <- InvSimpson_sig[is.na(InvSimpson_sig$Sig)==FALSE, ]
InvSimp_label <- data.frame(Treatment = rownames(InvSimpson_sig), invSimp = InvSimpson_sig$coord)

#Create objects for plotting Faith's Diversity
sample3 <- data.frame(Treatment = as.character(sample2$Diet), Faith = sample2$Faith, row.names = rownames(sample2))
sample3 <- sample3[sample3$Treatment != "HF",]
adjusted.faith$max_faith <- tapply(sample3$Faith, sample3$Treatment, max)
adjusted.faith$coord <- adjusted.faith$max_faith +1
Faith_sig <- adjusted.faith
rownames(Faith_sig) <- Treatments_test_full
Faith_sig <- Faith_sig[is.na(Faith_sig$Sig)==FALSE, ]
Faith_label <- data.frame(Treatment = rownames(Faith_sig), Faith = Faith_sig$coord)

plot_sample <- data.frame(Group=rownames(sample2), Treatment=sample2$Diet, 
                          Shannon=sample2$Shannon, Simpson=sample2$Simpson, invSimp=sample2$InvSimp, Faith=sample2$Faith)


tiff("Genus_Shannon_Diet.tiff", units="in", width=5, height=5, res=300)
ggplot(plot_sample, aes(Treatment, Shannon)) + 
  geom_boxplot() + 
  ggtitle("Genus-level Shannon Diversity") + 
  theme (plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(label=Treatments_abrev) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text (data = Shan_label, label = Shannon_sig$Sig)
dev.off()

tiff("Genus_Simpson_Diet.tiff", units="in", width=5, height=5, res=300)
ggplot(plot_sample, aes(Treatment, Simpson)) + 
  geom_boxplot() +
  ggtitle("Genus-level Simpson Diversity") +
  theme (plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(label=Treatments_abrev) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text (data = Simp_label, label = Simpson_sig$Sig)
dev.off()

tiff("Genus_Inverse_Simpson_Diet.tiff", units="in", width=5, height=5, res=300)
ggplot(plot_sample, aes(Treatment, invSimp)) + 
  geom_boxplot() +
  ggtitle("Genus-level Inverse Simpson Diversity") +
  theme (plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(label=Treatments_abrev) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text (data = InvSimp_label, label = InvSimpson_sig$Sig)
ylab("Inverse Simpson")
dev.off()

tiff("Genus_Faith_Diet.tiff", units="in", width=5, height=5, res=300)
ggplot(plot_sample, aes(Treatment, Faith)) + 
  geom_boxplot() +
  ggtitle("Genus-level Faith Diversity") +
  theme (plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(label=Treatments_abrev) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text (data = Faith_label, label = Faith_sig$Sig)
dev.off()

####Beta-diversity Analysis####
sample.a <- as.data.frame(sample2)

##Bray-Curtis Distance via DivNet##

#Create distance object and perform ordination
b_div_bray <- as.dist(a_div$`bray-curtis`)
b_div_bray_pcoa <- ordinate(physeq_gen, method="PCoA", distance=b_div_bray) 


tiff("DivNet_Genus_Bray_Diet.tiff", units="in", width=5, height=5, res=300)
phyloseq::plot_ordination(physeq_gen, b_div_bray_pcoa, type ="samples", color = "Diet") +
  ggtitle("Bray-Curtis Dissimilarity") +
  theme (plot.title = element_text(hjust = 0.5))
dev.off()

#Test for overall significance
write("Genus Level DivNet Diet and Bray-Curtis Diversity Test", 
      file= "PERMANOVA_results_overall.txt", append=TRUE)
set.seed(3488)
capture.output(adonis(formula = b_div_bray ~ Diet, data = sample.a), 
               file="PERMANOVA_results_overall.txt", append=TRUE)
#error in `colnames<-`(`*tmp*`, value = colnames(lhs)) :attempt to set 'colnames' on an object with less than two dimensions

#Test for pair-wise significance
set.seed(3488)
adjusted.bray <- calc_pairwise_permanovas(b_div_bray, sample.a, "Diet")
adjusted.bray <- adjusted.bray[adjusted.bray$X1=="HF",-c(5:6)]
adjusted.bray$qval <- p.adjust(adjusted.bray$pval, method = "BH")
adjusted.bray$Sig <- NA
adjusted.bray$Sig[adjusted.bray$qval<0.05] <- "*"
adjusted.bray$Sig[adjusted.bray$qval<0.01] <- "**"
adjusted.bray$Sig[adjusted.bray$qval<0.001] <- "***"

write("Genus Level DivNet Treatment and Bray-Curtis Diversity Test", file= "PERMANOVA_results.txt", append=TRUE)
capture.output(adjusted.bray, file="PERMANOVA_results.txt", append=TRUE)

##Aitchison Distance Analysis##

#CLR transformation of counts
physeq_gen_clr <- microbiome::transform(physeq_gen, "clr")
otu.clr <- t(as(otu_table(physeq_gen_clr), "matrix"))
#calculate aitchison distance (Euclidan distance of clr transformed data)
clr_ait <- dist(otu.clr, method='euc')
b_div_ait_pcoa <- ordinate(physeq_gen, method="PCoA", distance=clr_ait)


tiff("Genus_Aitchison_Diet.tiff", units="in", width=5, height=5, res=300)
phyloseq::plot_ordination(physeq_gen, b_div_ait_pcoa, type ="samples", color = "Diet") +
  ggtitle("Aitchison Distance") +
  theme (plot.title = element_text(hjust = 0.5))
dev.off()

#Test for overall significance
write("\nGenus Level Diet and Aitchison Diversity Test", 
      file= "PERMANOVA_results_overall.txt", append=TRUE)
set.seed(3488)
capture.output(adonis(formula = clr_ait ~ Diet, data = sample.a), 
               file="PERMANOVA_results_overall.txt", append=TRUE)

#Test for pair-wise significance
set.seed(3488)
adjusted.ait <- calc_pairwise_permanovas(clr_ait, sample.a, "Diet")
adjusted.ait <- adjusted.ait[adjusted.ait$X1=="HF",-c(5:6)]
adjusted.ait$qval <- p.adjust(adjusted.ait$pval, method = "BH")
adjusted.ait$Sig <- NA
adjusted.ait$Sig[adjusted.ait$qval<0.05] <- "*"
adjusted.ait$Sig[adjusted.ait$qval<0.01] <- "**"
adjusted.ait$Sig[adjusted.ait$qval<0.001] <- "***"

write("\nGenus Level Diet and Aitchison Diversity Test", file= "PERMANOVA_results.txt", append=TRUE)
capture.output(adjusted.ait, file="PERMANOVA_results.txt", append=TRUE)

##Unifrac analysis##

dist_gen_rar_uniw <- phyloseq::UniFrac(physeq_gen.rar, weighted = TRUE)
gen_rar_uniw_pcoa <- ordinate(physeq_gen.rar, method="PCoA", distance=dist_gen_rar_uniw)

tiff("Genus_UnifracW_Diet.tiff", units="in", width=5, height=5, res=300)
phyloseq::plot_ordination(physeq_gen.rar, gen_rar_uniw_pcoa, type ="samples", color = "Diet") +
  ggtitle("Weighted UniFrac") +
  theme (plot.title = element_text(hjust = 0.5))
dev.off()

#Test for overall significance
write("\nGenus Level Diet and Weighted Unifrac Diversity Test", 
      file= "PERMANOVA_results_overall.txt", append=TRUE)
set.seed(3488)
capture.output(adonis(formula = dist_gen_rar_uniw ~ Treatment, data = sample.a), 
               file="PERMANOVA_results_overall.txt", append=TRUE)

#Test for pair-wise significance
set.seed(3488)
adjusted.uniw <- calc_pairwise_permanovas(dist_gen_rar_uniw, sample.a, "Treatment")
adjusted.uniw <- adjusted.uniw[adjusted.uniw$X1=="Water",-c(5:6)]
adjusted.uniw$qval <- p.adjust(adjusted.uniw$pval, method = "BH")
adjusted.uniw$Sig <- NA
adjusted.uniw$Sig[adjusted.uniw$qval<0.05] <- "*"
adjusted.uniw$Sig[adjusted.uniw$qval<0.01] <- "**"
adjusted.uniw$Sig[adjusted.uniw$qval<0.001] <- "***"

write("\nGenus Level Treatment and Weighted Unifrac Diversity Test", file= "PERMANOVA_results.txt", append=TRUE)
capture.output(adjusted.uniw, file="PERMANOVA_results.txt", append=TRUE)


