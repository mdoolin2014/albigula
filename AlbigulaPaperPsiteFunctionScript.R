#Code to process data from M.L. Doolin, S.B. Weinstein, M.D. Dearing —— 
#Neotoma albigula parasite-microbiome study, Journal of Parasitology
#Submitted for publication February 2022. 
#Parasitism metrics, Dry matter digestibility and Oxalate analyses. 

##R version 4.0.2 (2020-06-22)
##Platform: x86_64-apple-darwin17.0 (64-bit)
##Running under: macOS Mojave 10.14.3


############## Parasite prevalence and intensity #####

##### Confidence intervals for population-level prevalence #####
library(coda)
library(rjags)
library(prevalence) #parasite prevalence package
library(lme4)

#Bayesian confidence interval for prevalence, using above prevalence package
?propCI

#e.g., May animals with no infection
propCI(4, 15, method = "all", level= 0.95, sortby= "level")


#Upload the csv's I can use with all of the parasite infection classes
prev <- read.csv("Psite Prev.csv") #All infection classes
prev
prev$Parasite <- factor(prev$Parasite, levels=c("Coinfection", "Coccidian", "Pinworm", "No Infection"))

#Subset df to remove coinfections
prev1 <- prev[prev$Parasite != "Coinfection",]
prev1$Parasite = factor(prev1$Parasite, levels=c("Coccidian", "Pinworm", "No Infection"))

#Subset that df to remove coccidians, and make it only pinworm and uninfected.
prev2 <- prev1[which(prev1$Parasite != "Coccidian"),]
prev2$Parasite = factor(prev2$Parasite, levels=c("Pinworm","No Infection"))


##### Making dataframes and plots to look at prevalence and egg intensity. #####

#Defining where to draw CI info from. Or could just put these variables in the calls for the plot.
b=prev$HiCI
a=prev$LoCI

p1 <- ggplot(data=prev, aes(x=Month, y=Prevalence, group=Parasite, shape=Parasite, colour=Month)) +
  geom_point(data=prev, aes(size=10), position=position_dodge(0.6)) + 
  scale_shape_manual(values=c(15, 18, 19, 17)) +
  theme_bw() + scale_colour_manual(values=c("#000000", "#999999")) +
  scale_y_continuous(limits=c(0,.8)) + labs(title="Parasite Prevalence", x="Parasite", y="Prevalence") + 
  theme(plot.title=element_text(hjust=0.5,size=20,face="bold")) + 
  guides(fill=guide_legend(title="Parasite Infection")) +
  geom_errorbar(data=prev, aes(ymin=a, ymax=b, width=0.25), position=position_dodge(0.6))
p1


#Now, make a data frame with the intensities for nematodes eggs per gram.
#Redoing this because wtf was I doing here? These numbers don't make any sense.
NemaInt <- data.frame("ID" = c(577, 580, 594, 599, 554, 557, 559, 586), 
                      "AvgIntens" = c(0.104, 0.549, 0.959, 0.41, 0.754, 0.136, 0.397, 0.49))
NemaInt$Parasite <- "Nematode"
NemaInt$Parasite <- as.factor(NemaInt$Parasite)

#Plot the nematode egg intensities. 
ggplot(data=NemaInt, aes(Parasite, AvgIntens))+
  geom_boxplot(data=NemaInt, aes(x=Parasite, y=AvgIntens, fill=Parasite)) +
  theme_bw() + scale_fill_manual(values=c("dodgerblue2")) +
  scale_y_continuous(limits=c(0,1)) + labs(title="Fecal Egg Intensity",
                                           y="Eggs Per Gram Feces") + theme(plot.title=element_text(hjust=0.5,size=20,face="bold")) +
  guides(fill=guide_legend(title="Parasite Infection")) +
  geom_jitter(shape=19, size=5, position=position_jitter(0.1))



##### Is there a difference in body weight based on parasitism? #####

#plotting parasite infection by weight over all animals in both months
ggplot(alldat1, aes(psite, rat_wt)) + 
  geom_boxplot(data=alldat1, aes(x=psite, y=rat_wt, fill=psite)) +
  scale_y_continuous(limits=c(90,250)) + labs(title="Rat weight by psite infection",
                                              x="Parasite infection", y="Rat weight (g)") + 
  theme(plot.title=element_text(hjust=0.5,size=20,face="bold")) +
  guides(fill=guide_legend(title="Parasite Infection")) +
  geom_jitter(shape=21, position=position_jitter(0.1))


#Choose the best model for fit, looks like there are differences based on month,
# but is there a difference based on parasitism?
null=lm(rat_wt~1, data=alldat1)
full=lm(rat_wt~psite + month + sex + psite*month, data=dat)
step(full, data=dat, direction="backward")




############## Functional metrics #####
##### Is there a difference in October oxalate intake based on psite? - ANCOVA ####

#Tutorial from: https://www.datanovia.com/en/lessons/ancova-in-r/#:~:text=The%20Analysis%20of%20Covariance%20(ANCOVA,two%20or%20more%20independent%20groups.
#Using their code, but putting in my own data.

library(tidyverse)
library(ggpubr)
library(rstatix)  #This is how they're doing the stats, instead of base stats package
library(broom)

#Make a scatterplot to check linearity of the groups within the grouping variable, assess 
# by eye and also use the equations of the lines to assess. 
ggscatter(oxdata1, x = "rat_wt", y = "oxalate_intake",
          color = "parasite", add = "reg.line")+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = parasite))

#Look for interaction btwn grouping variable and independent variable
oxdata1 %>% anova_test(oxalate_intake ~ rat_wt*parasite)

#Make a model to look at the residuals
model <- lm(oxalate_intake ~ rat_wt + parasite, data = oxdata1)
#Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted) # Remove details
head(model.metrics, 3)

#Assess normality with shapiro-wilkes test
shapiro_test(model.metrics$.resid)
#Not significant, so can proceed with parametrics

#Also want to check homogeneity of variances, with Levene's test
model.metrics %>% levene_test(.resid ~ parasite)
#Also not signif bc, so we can proceed.

#Also need to check for outliers, with standardized residuals greater than 3.
model.metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()
##There are none.

#After those tests, we can now move forward in computing the ANCOVA
res.aov <- oxdata1 %>% anova_test(oxalate_intake ~ rat_wt + parasite)
get_anova_table(res.aov)

#Post-hoc testing with bonferroni multiple testing correction
install.packages("emmeans")
library(emmeans) #This is for estimated marginal means (i.e. least square means)
pwc <- oxdata1 %>% 
  emmeans_test(
    oxalate_intake ~ parasite, covariate = rat_wt,
    p.adjust.method = "bonferroni")

#Look at the means and s.e., and other descriptive stuff.
get_emmeans(pwc)

#Making a plot to see what is sig diff, if anything. 
pwc <- pwc %>% add_xy_position(x = "parasite", fun = "mean_se")
ggline(get_emmeans(pwc), x = "parasite", y = "emmean") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
  stat_pvalue_manual(pwc, hide.ns = TRUE, tip.length = FALSE) +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))

#Visualizing the data to see if oxalate intake varies by CAC versus OX diet. 
ggplot(data=octdat1, aes(x=trial, y=oxalate_intake, fill=parasite))+
  geom_boxplot(outlier.shape=NA) +
  theme_bw(base_size = 16) + 
  geom_point(pch=21, alpha= 0.8, position=position_jitterdodge(0.2), size=5) 




##### Is there a difference in October DMD based on psite? #####
t.test(dmd~as.factor(parasite), data=octdat1) 


