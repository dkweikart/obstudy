setwd("~/PSU/Obesity Study/New_Figs")
###PCA
# Load
library("FactoMineR")
library("factoextra")

# CVA on treatment-sex interaction
# use new dataset that has a ST column

df <- read.table("CVA_obstudy_1.csv", sep=',', header =TRUE)
head(df)
str(df)
df <- df[,1:30]

# set factors ... columns 1:5 and 32:33
df$Pool <- as.factor(df$Pool)
df$Sex <- as.factor(df$Sex)
df$Treament <- as.factor(df$Treament)
df$ST <- as.factor(df$ST)
df$ST <- as.factor(df$ST_No)
df$Ferm <- as.factor(df$Ferm)
df$Roast <- as.factor(df$Roast)

dim(df)
str(df)

rownames(df) <- df$ST_No
head(df)

df <- df[,-5]
head(df)

library(candisc)

# model for Sex-Treatment effect

df.lm <- lm(as.matrix(df[, -c(1:4, 31, 32)]) ~ ST, data = df)
summary(manova(df.lm))  # sign. ST effects

# CVA
df.cva <- candisc(df.lm)
df.cva

plot(df.cva, which = 1, points.1d = TRUE)
plot(df.cva, which = 1:2, 
     conf = 0.95, var.col = 'grey',
     ellipse = TRUE, ellipse.prob = 0.68)


## PCA with ST coloring
ST.pca <- PCA(df[,-c(1,31,32)], scale.unit = TRUE,
              quali.sup = c(1,2,3))

# fviz_pca doesn't appear to allow you the customization needed
# so I changed to standard ggplot2
###Like the colors

fviz_pca_ind(ST.pca, axes = c(1,2), 
             geom.ind = c("point"), # show points only (but not "text")
             pointsize = 2, 
             col.ind = df$Sex, # color by groups
             # label = ,
             habillage = df$Sex,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.level = 0.95, 
             legend.title = "Groups",
             mean.point = FALSE,
             repel = TRUE) +
  theme_minimal() +
  scale_shape_manual(values = rep(c(15, 16, 17, 18, 19, 20, 21, 22), 2))

# ggplot 2
library(ggplot2)
library(ggthemes)

# prepare data set
# extract coordinates of individuals and add ST column

ST.pca.df <- as.data.frame(ST.pca$ind$coord[,1:2])
ST.pca.df$ST <- c(rep("M_HF",5),
                  rep("M_NO/NO",5),
                  rep("M_LO/NO",5),
                  rep("M_HI/NO",5),
                  rep("M_NO/LO",5),
                  rep("M_NO/HI",5),
                  rep("M_LO/LO",5),
                  rep("M_HI/HI",5),
                  rep("F_HF",5),
                  rep("F_NO/NO",5),
                  rep("F_LO/NO",5),
                  rep("F_HI/NO",5),
                  rep("F_NO/LO",5),
                  rep("F_NO/HI",5),
                  rep("F_LO/LO",5),
                  rep("F_HI/HI",5))
ST.pca.df$Sex <- c(rep("Male", 40),
                   rep("Female", 40))
ST.pca.df$Treatment <- c(rep("HF",5),
                         rep("NO/NO",5),
                         rep("LO/NO",5),
                         rep("HI/NO",5),
                         rep("NO/LO",5),
                         rep("NO/HI",5),
                         rep("LO/LO",5),
                         rep("HI/HI",5),
                         rep("HF",5),
                         rep("NO/NO",5),
                         rep("LO/NO",5),
                         rep("HI/NO",5),
                         rep("NO/LO",5),
                         rep("NO/HI",5),
                         rep("LO/LO",5),
                         rep("HI/HI",5))

ST.pca.df
str(ST.pca.df)

ST.plot <- ggplot(ST.pca.df) +
  geom_point(aes(x = Dim.1, y = Dim.2, 
                 shape = Treatment,
                 color = Sex), size = 2) + 
  scale_shape_manual(values = c(16, 0, 2, 4, 5, 6, 8, 7)) +
  stat_ellipse(geom = 'polygon',
               aes(x = Dim.1, y = Dim.2,
                   fill = Sex,
                   color = Sex),
               alpha = 0.25, 
               type = 'norm') +
  labs(x = "PC 1, 30%", y = 'PC 2, 22%') +
  theme_classic() +
  theme(legend.position = 'top',
        legend.box = 'vertical')
ST.plot

# ggpubr
library(ggpubr)

ggscatter(ST.pca.df, x = "Dim.1", y = "Dim.2",
          color = "Sex", palette = c("#00AFBB", "#E7B800"),
          shape = "Treatment", 
          ellipse = TRUE, ellipse.type = 'norm',
          ellipse.alpha = 0.25, 
          mean.point = FALSE,
          xlab = 'PC 1, 30%', ylab = 'PC 2, 22%')
