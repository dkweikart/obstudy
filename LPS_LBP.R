###Load packages###
library(dplyr)
library(ggpubr)
library(agricolae)

####Inflammation Markers####
###Assign Data to Object####
lipo <- read.csv('LPS_LBP.csv', header = TRUE)
View(lipo)

#Subset into male and female
lipo.mal <- lipo[1:40,]
lipo.fem <- lipo[1:40,]

###Summary###
#By Treatment for LPS in males
group_by(lipo.mal, Treament) %>%
  summarise(
    count = n(),
    mean = mean(LPS, na.rm = TRUE),
    sd = sd(LPS, na.rm = TRUE)
  )

#By Treatment for LBP in males
group_by(lipo.mal, Treament) %>%
  summarise(
    count = n(),
    mean = mean(LBP, na.rm = TRUE),
    sd = sd(LBP, na.rm = TRUE)
  )

#By Treatment for LPS in females 
group_by(lipo.fem, Treament) %>%
  summarise(
    count = n(),
    mean = mean(LPS, na.rm = TRUE),
    sd = sd(LPS, na.rm = TRUE)
  )

#By Treatment for LBP in females
group_by(lipo.fem, Treament) %>%
  summarise(
    count = n(),
    mean = mean(LBP, na.rm = TRUE),
    sd = sd(LBP, na.rm = TRUE)
  )


###Visualize data with boxplot###
boxplot(LPS ~ Treament, data = lipo.mal,
        xlab = "Treatment", ylab = "Concentration",
        frame = FALSE)

boxplot(LBP ~ Treament, data = lipo.mal,
        xlab = "Treatment", ylab = "Concentration",
        frame = FALSE)


boxplot(LPS ~ Treament, data = lipo.fem,
        xlab = "Treatment", ylab = "Concentration",
        frame = FALSE)

boxplot(LBP ~ Treament, data = lipo.fem,
        xlab = "Treatment", ylab = "Concentration",
        frame = FALSE)



###Compute ANOVA###
#LPS in males by treatment
mal.LPS.aov <- aov(LPS ~ Treament, data = lipo.mal)
#LBP in males by treatment
mal.LBP.aov <- aov(LBP ~ Treament, data = lipo.mal)
#LPS in females by treatment
fem.LPS.aov <- aov(LPS ~ Treament, data = lipo.fem)
#LBP in females by treatment
fem.LBP.aov <- aov(LBP ~ Treament, data = lipo.fem)



# Summary of the analysis#
summary(mal.LPS.aov)
summary(mal.LBP.aov)
summary(fem.LPS.aov)
summary(fem.LBP.aov)

####Tukey Post Hoc Test####
TukeyHSD(mal.LPS.aov)
TukeyHSD(mal.LBP.aov)
TukeyHSD(fem.LPS.aov)
TukeyHSD(fem.LBP.aov)



#Group based on significance (assigns letters based on sig)#
HSD.test(mal.LPS.aov, "Treament", console=TRUE)
HSD.test(mal.LBP.aov, "Treament", console=TRUE)
HSD.test(fem.LPS.aov, "Treament", console=TRUE)
HSD.test(fem.LBP.aov, "Treament", console=TRUE)


