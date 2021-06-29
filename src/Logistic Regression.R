###################################################
## Project: Applied Bio Statistics Semester Project
## Script purpose: Logistic Regression Modelling
## Date: May 2021
## Author: Kamran Nejad-Sattary
## 
## Dataset:  gold_target1.dat
## Source: N.R. Sahoo and H.S. Pandalai (1999). "Integration of Sparse Geologic
## Information in Gold Targeting Using Logistic Regression Analysis in the 
## Hutti-Maski Schist Belt, Raichur, Karnataka, India - A Case Study,"
## Natural Resources Research, Vol. 8, #3, pp. 233-250
## 
## Description: Logistic Regression of Presence/Absent of gold deposit
## as a function of water/chemical factors: As, Sb, and presence/absence
## of lineament.
## 
## Variables/Columns
## As level     1-8
## Sb level     10-16
## Lineament Proximity  24   1=Present, 0 if absent  (0.5km)  
## Gold deposit proximity  32   1=Present, 0=absent  (0.5km)  
###################################################
rm(list=ls())

# set session variables
library(ResourceSelection)
library(DescTools)
library(stargazer)
set.seed(42)

my_dir = "/media/kamran/hdd/Dropbox/EPFL/Semester 2/Applied Biostats/AppBioStats/"
setwd(my_dir)

output_dir = './output/'
data_url = "http://users.stat.ufl.edu/~winner/data/gold_target1.dat"

# read in data
master = read.table(url(data_url))
colnames(master) <- c("as_lvl", "sb_lvl", "lineament_proxim", "gold_proxim")

# Check for Multicollinearity
cor(master)
cor.test(master$as_lvl, master$lineament_proxim)

# EDA
boxplot_eda = function(){
  
  pdf(paste0(output_dir, "eda_boxplots.pdf"),
      width = 5,
      height = 7)
  
  par(mfrow=c(1,2), mar = c(8, 4, 2, 2))
  boxplot(as_lvl ~ gold_proxim, data = master, ylab ="Arsenic (As) Level (ppm)", xlab="", xaxt='n')
  text(x = 1:length(unique(master$gold_proxim)), y = par("usr")[3]-7.5,
       labels = c('No Proximity','Gold Proximity\n(0.5 km)'),
       cex = 1.2, xpd = NA, srt = 90)
  
  boxplot(sb_lvl ~ gold_proxim, data = master, ylab ="Antimony (Sb) Level (ppm)", xlab="", xaxt='n')
  text(x = 1:length(unique(master$gold_proxim)), y = par("usr")[3]-3.3,
       labels = c('No Proximity','Gold Proximity\n(0.5 km)'),
       cex = 1.2, xpd = NA, srt = 90)
  
  dev.off()
}
boxplot_eda()

barplot_eda = function(){
  
  pdf(paste0(output_dir, "eda_barplots.pdf"),
      width = 8,
      height = 5)
  
  par(mfrow=c(1,1), mar = c(5, 7, 0, 2))
  
  counts <- table(master$lineament_proxim, factor(master$gold_proxim, levels=c(0,1), labels=c('No Proximity','Gold Proximity\n(0.5 km)')))
  barplot(counts, main="", xlab="Count of Respective Category", col=c("#8fba8d","#c0b3c9"), beside=T, horiz=TRUE, las=1)

  legend("center",
         legend = c('No Proximity','Lineament Proximity (0.5 km)'),
         pch = 15,
         col = c("#8fba8d","#c0b3c9"))
  dev.off()
  
}
barplot_eda()

quantile(master[master$gold_proxim==0,]$as_lvl)
quantile(master[master$gold_proxim==0,]$sb_lvl)

quantile(master[master$gold_proxim==1,]$as_lvl)
quantile(master[master$gold_proxim==1,]$sb_lvl)

# Logistic Regression
master$lineament_proxim <- factor(master$lineament_proxim)
mylogit <- glm(gold_proxim ~ as_lvl + sb_lvl + lineament_proxim, data = master, family = "binomial", method="glm.fit")

stargazer(mylogit, font.size='small', no.space=TRUE, single.row=TRUE, report=('vc*s'))
summary(mylogit)

pdf(paste0(output_dir, "modelfit.pdf"),
    width = 6,
    height = 6)
par(mfrow=(c(2,2)))
plot(mylogit, which=c(1,2,3,5))
dev.off()

preds = predict(mylogit, newdata = master, type = "response")
print(paste0("Accuracy:", sum(as.numeric(preds > 0.5) == master$gold_proxim) / nrow(master)))
print(paste0("Misclassified Cases: ", sum(as.numeric(preds > 0.5) != master$gold_proxim), "/", nrow(master)))

# we cannot reject the null hypothesis that the model fits the data well
hl <- ResourceSelection::hoslem.test(master$gold_proxim, fitted(mylogit))
hl

# Logistic R2
DescTools::PseudoR2(mylogit, which='Nagelkerke')

drop1(mylogit, test='Chisq')

# compute predictive odds
exp(coef(mylogit))

lineament_proxim_data <- with(master, data.frame(as_lvl = mean(as_lvl), sb_lvl = mean(sb_lvl), lineament_proxim = factor(0:1)))
predict(mylogit, newdata = lineament_proxim_data, type = "response")



