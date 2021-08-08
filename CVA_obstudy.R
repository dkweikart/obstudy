##CVA with Sex and Treatment combined
library(candisc)

#Load file
all <- read.csv('CVA_obstudy.csv', header=T)
View(all)

#Remove last two columns
all <- all[,1:30]

#Set variables as factors 
all$Pool <- as.factor(all$Pool)
all$Sex <- as.factor(all$Sex)
all$Treament <- as.factor(all$Treament)
all$SexTreat <- as.factor(all$SexTreat)

#Create a linear model
all.mod <- lm(as.matrix(all[,-c(1:4)]) ~ SexTreat, data = all)

# MANOVA on the lm, as CVA uses a MANOVA model
summary(manova(all.mod), test = 'Wilks')  

#Perform canonoical discriminant analysis and linear model
all.cva <- candisc(all.mod)
# look at CVA summary
all.cva # provides brief summary, incl. Bartlett's test => 2 sign. dims
summary(all.cva)  # provides different results - product means and std coeffs
all.cva$scores
plot(all.cva$eigenvalues)

# plot CVA plot with overlaid confidence ellipses
col1 <- c('#fee090','#d73027','green','#4575b4','#fc8d59','#91bfdb')
pch1 <- c(15,16,17,8,6,3)

# confidence ellipses built into candisc package
plot(all.cva, which = 1:2, conf = 0.95, col = col1,
     pch = pch1, ellipse=FALSE, prefix = "Can", suffix=TRUE)

plot(all.cva, which = 1:2, conf = 0.95, col = col1,
     pch = pch1, ellipse=TRUE, ellipse.prob = 0.68, fill.alpha=0.1,
     prefix = "Can", suffix=TRUE)

#1 dimensional
plot(all.cva, which = 1, conf = 0.95, col = col1,
     pch = pch1, ellipse=FALSE, prefix = "Can", suffix=TRUE)

plot(all.cva, which = 1, conf = 0.95, col = col1,
     pch = pch1, ellipse=TRUE, ellipse.prob = 0.68, fill.alpha=0.1,
     prefix = "Can", suffix=TRUE)



