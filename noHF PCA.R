###PCA no HF group
library("FactoMineR")
library("factoextra")

#Load File 
noHF <- read.csv('noHFD.csv', header=T)
View(noHF)


#Subset into male and female
treat.mal <- noHF[1:35,3:28]
treat.fem <- noHF[36:70,3:28]
head(treat.mal)
head(treat.fem)

#Load library
library(ggfortify)
pca_res <- prcomp(treat.mal[,-1], scale. = TRUE)

#Color indiv. by group in males
fviz_pca_ind(pca_res, label="none", habillage=treat.mal$Treament, 
             geom="point", addEllipses = TRUE)

#FEMALES
pca_res.f <- prcomp(treat.fem[,-1], scale. = TRUE)
fviz_pca_ind(pca_res.f, label="none", habillage=treat.fem$Treament, 
             geom="point", addEllipses = TRUE)
