setwd("~/PennState/Obesity Study/Microbiome_code")

# Load libraries
library(microbiome)
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names

#Load the data
taxa <- read_excel("taxa.xlsx") # taxon x sample matrix
otus <- read_excel("otus.xlsx", sheet = "t.otus") # lineage x taxon matrix
metadata <- read_excel("DaphneSampleKey.xlsx") # metadata x sample data matrix

###Make Phyloseq Object
otus <- otus %>%
  tibble::column_to_rownames("otu") 

taxa <- taxa %>% 
  tibble::column_to_rownames("otu")

metadata <- metadata %>% 
  tibble::column_to_rownames("File") 

otus <- as.matrix(otus)
taxa <- as.matrix(taxa)

OTU = otu_table(otus, taxa_are_rows = TRUE)
TAX = tax_table(taxa)
samples = sample_data(metadata)

phy.o <- phyloseq(OTU, TAX, samples)
phy.o


sample_names(phy.o)
rank_names(phy.o)
sample_variables(phy.o)

total = median(sample_sums(phy.o))
standf = function(x, t=total) round(t * (x / sum(x)))
phy.o = transform_sample_counts(phy.o, standf)

###No MOck
newPhy.o = subset_samples(phy.o, sample_names(phy.o)!="DC006169")

sample_names(newPhy.o)
rank_names(newPhy.o)
sample_variables(newPhy.o)

total = median(sample_sums(newPhy.o))
standf = function(x, t=total) round(t * (x / sum(x)))
newPhy.o = transform_sample_counts(newPhy.o, standf)

# Pick relative abundances (compositional) and sample metadata
pseq.rel <- microbiome::transform(newPhy.o, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

#Visualize
p <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "Diet", size = 3)
print(p)

# samples x species as input
library(vegan)
permanova <- adonis(t(otu) ~ Diet,
                    data = meta, permutations=99, method = "bray")

# P-value
print(as.data.frame(permanova$aov.tab)["Diet", "Pr(>F)"])

#Check homogeniety
# Note the assumption of similar multivariate spread among the groups
# ie. analogous to variance homogeneity
# Here the groups have signif. different spreads and
# permanova result may be potentially explained by that.
dist <- vegdist(t(otu))
anova(betadisper(dist, meta$Diet))

permutest(betadisper(dist, meta$Diet), pairwise = TRUE)


#Investigate top factors 
coef <- coefficients(permanova)["Diet1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
