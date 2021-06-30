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

### set session variables
library(ResourceSelection)
library(DescTools)
library(stargazer)
library(ggplot2)
library(car)
library(tidyverse)

set.seed(42)

my_dir = "/media/kamran/hdd/Dropbox/EPFL/Semester 2/Applied Biostats/AppBioStats/"
setwd(my_dir)

output_dir = './output/'
data_url = "http://users.stat.ufl.edu/~winner/data/gold_target1.dat"

### read in data
master = read.table(url(data_url))
colnames(master) <- c("as_lvl", "sb_lvl", "lineament_proxim", "gold_proxim")

### Check for Multicollinearity (point-biserial), vif below under the model part
cor(master)
cor.test(master$as_lvl, master$lineament_proxim)


### EDA


# Univariate numerical analysis
# produce summary statistics without disaggregating
summary(master)

length(master[master$gold_proxim==1,]$as_lvl)
length(master[master$lineament_proxim==1,]$as_lvl)

# produce summary statistics with disaggregating by gold proximity
summary(master[master$gold_proxim==0,])
summary(master[master$gold_proxim==1,])

# Numerical analysis of each class combination
quantile(master[master$gold_proxim==0,]$as_lvl)
quantile(master[master$gold_proxim==0,]$sb_lvl)

quantile(master[master$gold_proxim==1,]$as_lvl)
quantile(master[master$gold_proxim==1,]$sb_lvl)

# Distributions of elements by gold proximity
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

# Counts of lineament proximity by gold proximity
barplot_eda = function(){
  
  pdf(paste0(output_dir, "eda_barplots.pdf"),
      width = 8,
      height = 5)
  
  par(mfrow=c(1,1), mar = c(5, 7, 0, 2))
  
  counts <- table(master$lineament_proxim, factor(master$gold_proxim, levels=c(0,1), labels=c('No Proximity','Gold Proximity\n(0.5 km)')))
  barplot(counts, main="", xlab="Count of Respective Category", col=c("#c0b3c9", "#8fba8d"), beside=T, horiz=TRUE, las=1)

  legend("center",
         legend = c('No Proximity','Lineament Proximity (0.5 km)'),
         pch = 15,
         col = c("#c0b3c9","#8fba8d"))
  dev.off()
  
}
barplot_eda()

# Bivariate scatter plot of arsenic and antimony levels
scatterplot_eda = function(){
    
  pdf(paste0(output_dir, "eda_scatterplots.pdf"),
      width = 8,
      height = 5)
  p = ggplot(master, aes(x=as_lvl, y=sb_lvl, color=factor(gold_proxim))) +
        geom_point(size=2) + 
        geom_smooth(method='lm', se=F) + 
        xlab("Arsenic (As) Level (ppm)") + 
        ylab("Antimony (Sb) Level (ppm)") +
        theme_bw() + 
        labs(color="Gold Proximity (0.5 km)") +
        scale_color_manual(values = c("#77b874","#b099bf"), labels = c("Yes", "No"))
  print(p)
  dev.off()
}
scatterplot_eda()

### Model Fitting

# Logistic Regression
master$lineament_proxim <- factor(master$lineament_proxim)
mylogit <- glm(gold_proxim ~ as_lvl + sb_lvl + lineament_proxim, data = master, family = "binomial", method="glm.fit")

# check model VIF
car::vif(mylogit)

stargazer(mylogit, font.size='small', no.space=TRUE, single.row=TRUE, report=('vc*s'))
summary(mylogit)

# Save diagnostic plots
pdf(paste0(output_dir, "modelfit.pdf"),
    width = 6,
    height = 6)
par(mfrow=(c(2,2)))
plot(mylogit, which=c(1,2,3,5))
dev.off()

preds = predict(mylogit, newdata = master, type = "response")
print(paste0("Accuracy:", sum(as.numeric(preds > 0.5) == master$gold_proxim) / nrow(master)))
print(paste0("Misclassified Cases: ", sum(as.numeric(preds > 0.5) != master$gold_proxim), "/", nrow(master)))

# logit vs predictors (linear assumption plot)
# method taken from source: http://www.sthda.com/english/articles/36-classification-methods-essentials/148-logistic-regression-assumptions-and-diagnostics-in-r/
logit_eda = function(){
  
  logit.data = master
  logit.data$`Arsenic Level` = logit.data$as_lvl
  logit.data$`Antimony Level` = logit.data$sb_lvl
  logit.data[,c('gold_proxim','lineament_proxim','as_lvl','sb_lvl')] = NULL
  model.probs <- predict(mylogit, type = "response")
  predicted.classes <- ifelse(model.probs > 0.5, "pos", "neg")
  logit.data$logit = log(model.probs/(1-model.probs))
  
  logit.plot.data = logit.data %>% gather(key = "predictors", value = "predictor.value", -logit)
  
  pdf(paste0(output_dir, "eda_logits.pdf"),
      width = 8,
      height = 5)
  
  p = ggplot(logit.plot.data, aes(logit, predictor.value))+
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(method = "loess") + 
    theme_bw() + 
    facet_wrap(~predictors, scales = "free_y") +
    xlab("Outcome in Logit Scale") + ylab("Predictor Value")
  print(p)
  
  dev.off()
}
logit_eda()

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



