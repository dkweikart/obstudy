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


###RELATIVE ABUNDACE Species
library(phyloseq)
library(data.table)
library(dplyr)
library(ggplot2)

#Create phyloseq object (see above)
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

#Get the average relative abundance for each of the diet conditions for each member 
#of the taxonomic level, in this case phylum - can be changed to other
physeq_51_5 <- tax_glom(physeq_51_5, taxrank = "Species")
physeq_51_5 <- transform_sample_counts(physeq_51_5, function(x) x/sum(x))

physeq_HF <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "HF")))
physeq_NONO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "NO/NO")))
physeq_LONO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "LO/NO")))
physeq_HINO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "HI/NO")))
physeq_NOLO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "NO/LO")))
physeq_NOHI <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "NO/HI")))
physeq_LOLO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "LO/LO")))
physeq_HIHI <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "HI/HI")))


#update Phylum here if using other taxonomic level
Diet.pct <- as.data.table(data.frame(OTU=rownames(tax_table(physeq_51_5)), Phylum=tax_table(physeq_51_5)[,"Species"], 
                                     HF = round(rowMeans(physeq_HF)*100,3), 
                                     HIHI = round(rowMeans(physeq_HIHI)*100,3), HINO = round(rowMeans(physeq_HINO)*100,3), 
                                     LOLO = round(rowMeans(physeq_LOLO)*100,3), LONO = round(rowMeans(physeq_LONO)*100,3),
                                     NOHI = round(rowMeans(physeq_NOHI)*100,3), NOLO = round(rowMeans(physeq_NOLO)*100,3),
                                     NONO = round(rowMeans(physeq_NONO)*100,3)))

fwrite(Diet.pct, file="Perc_abund_Spec.txt", sep = "\t")

###Percent Abundance Genus
#Create phyloseq object (see above)
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

#Get the average relative abundance for each of the diet conditions for each member 
#of the taxonomic level, in this case phylum - can be changed to other
physeq_51_5 <- tax_glom(physeq_51_5, taxrank = "Genus")
physeq_51_5 <- transform_sample_counts(physeq_51_5, function(x) x/sum(x))

physeq_HF <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "HF")))
physeq_NONO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "NO/NO")))
physeq_LONO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "LO/NO")))
physeq_HINO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "HI/NO")))
physeq_NOLO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "NO/LO")))
physeq_NOHI <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "NO/HI")))
physeq_LOLO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "LO/LO")))
physeq_HIHI <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "HI/HI")))


#update Phylum here if using other taxonomic level
Diet.pct <- as.data.table(data.frame(OTU=rownames(tax_table(physeq_51_5)), Phylum=tax_table(physeq_51_5)[,"Genus"], 
                                     HF = round(rowMeans(physeq_HF)*100,3), 
                                     HIHI = round(rowMeans(physeq_HIHI)*100,3), HINO = round(rowMeans(physeq_HINO)*100,3), 
                                     LOLO = round(rowMeans(physeq_LOLO)*100,3), LONO = round(rowMeans(physeq_LONO)*100,3),
                                     NOHI = round(rowMeans(physeq_NOHI)*100,3), NOLO = round(rowMeans(physeq_NOLO)*100,3),
                                     NONO = round(rowMeans(physeq_NONO)*100,3)))

fwrite(Diet.pct, file="Perc_abund_Gen.txt", sep = "\t")

###Percent Abundance Phylum
#Create phyloseq object (see above)
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

#Get the average relative abundance for each of the diet conditions for each member 
#of the taxonomic level, in this case phylum - can be changed to other
physeq_51_5 <- tax_glom(physeq_51_5, taxrank = "Phylum")
physeq_51_5 <- transform_sample_counts(physeq_51_5, function(x) x/sum(x))

physeq_HF <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "HF")))
physeq_NONO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "NO/NO")))
physeq_LONO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "LO/NO")))
physeq_HINO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "HI/NO")))
physeq_NOLO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "NO/LO")))
physeq_NOHI <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "NO/HI")))
physeq_LOLO <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "LO/LO")))
physeq_HIHI <- data.frame (otu_table(subset_samples(physeq_51_5, Diet == "HI/HI")))


#update Phylum here if using other taxonomic level
Diet.pct <- as.data.table(data.frame(OTU=rownames(tax_table(physeq_51_5)), Phylum=tax_table(physeq_51_5)[,"Phylum"], 
                                     HF = round(rowMeans(physeq_HF)*100,3), 
                                     HIHI = round(rowMeans(physeq_HIHI)*100,3), HINO = round(rowMeans(physeq_HINO)*100,3), 
                                     LOLO = round(rowMeans(physeq_LOLO)*100,3), LONO = round(rowMeans(physeq_LONO)*100,3),
                                     NOHI = round(rowMeans(physeq_NOHI)*100,3), NOLO = round(rowMeans(physeq_NOLO)*100,3),
                                     NONO = round(rowMeans(physeq_NONO)*100,3)))

fwrite(Diet.pct, file="Perc_abund_Phy.txt", sep = "\t")



####HEATMAP 2
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)

#Import data
da_tab <- fread("Diff_abund.csv")
rel_phy <- fread("Perc_abund_Phy.txt")
rel_gen <- fread("Perc_abund_Gen.txt")
rel_sp <- fread("Perc_abund_Spec.txt")

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
treat_ord <- c("NONO","LONO","HINO","NOLO","NOHI", "LOLO","HIHI" )
plot_vals$Treatment <- factor(plot_vals$Treatment, levels = treat_ord)
plot_sigs$Treatment <- factor(plot_sigs$Treatment, levels = treat_ord)

plot_vals <- plot_vals[order(plot_vals$Level, plot_vals$Name, plot_vals$Treatment, decreasing = TRUE)]
plot_sigs <- plot_sigs[order(plot_sigs$Level, plot_sigs$Name, plot_sigs$Treatment, decreasing = TRUE)]

plot_vals$Name <- factor(plot_vals$Name, levels = ordernames)
plot_sigs$Name <- factor(plot_sigs$Name, levels = ordernames)

plot_sigs <- plot_sigs[is.na(plot_sigs$Sig)==FALSE,]

sig_dat <- data.frame(Treatment = plot_sigs$Treatment, Name = plot_sigs$Name)
Lab_txt <- c("NONO","LONO","HINO","NOLO","NOHI", "LOLO","HIHI")

tiff("ANCOM_Diff_abundance.tiff", units="in", width=7, height=6, res=300)
ggplot(plot_vals, aes(x = Treatment, y = Name)) +
  geom_tile(aes(fill = Rel_Change)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  scale_x_discrete(label=Lab_txt) + 
  labs(y = "Taxa") +
  geom_text(data = sig_dat, label = plot_sigs$Sig)
dev.off()

