
#Set working directory

setwd("~/Desktop/MAC/Environomics_FSP/Analysis/Proof_of_concept_paper/JULY22/R")

#Load packages

library(lme4)
library(performance)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(ggplot2)
library(tidyr)

#------------------------------------
#Evenness of sequencing depth
#------------------------------------
  #Experiment 3 and 4
    #lmer

#The within-pool co-efficient of variation (CV) was calculated for each sample. 
#The co-efficient of variation (CV) is calculated as the standard deviation of within-pool read proportions divided by the mean of the within-pool read proportions. 
#A CV of 0 means samples were sequenced exactly evenly. 

#Load dataset

evenness_cv_exp34<-read.csv("evenness_cv_exp34.csv",header=TRUE)
head(evenness_cv_exp34)

#We performed a linear mixed model in R v4.2.3 with the CV as the response variable, with capture or shotgun (“Sequencing”) and library preparation input and size (“Treatment”) as fixed effects, and pool (“Batch”) as a random effect.

cv.mod<-lmer(CV ~ Treatment * Sequencing + (1|Batch), data=evenness_cv_exp34,REML=TRUE)
summary(cv.mod)

#Model estimated means and 95% confidence intervals were calculated:

emm<- emmeans(cv.mod, ~Treatment*Sequencing, lmer.df="satterthwaite")
emm

#Further simplify the model by removing the interaction between “Treatment” and “Sequencing”

cv.mod2<-lmer(CV ~ Treatment + Sequencing + (1|Batch), data=evenness_cv_exp34,REML=TRUE)
anova(cv.mod, cv.mod2)

#The interaction is non-signficiant. Further simplify the model by removing the “Sequencing” term:

cv.mod3<-lmer(CV ~ Treatment + (1|Batch), data=evenness_cv_exp34,REML=TRUE)
anova(cv.mod2, cv.mod3)

#“Sequencing” is not a significant factor. Further simplify the model by removing “Treatment”:

cv.mod4<-lmer(CV ~ (1|Batch), data=evenness_cv_exp34,REML=TRUE)
anova(cv.mod3, cv.mod4)

#”Treatment” is not a significant factor. The random factor “Batch” contributes marginally to variance. 

VarCorr(cv.mod4)
summary(cv.mod4)

#Model residuals were approximately normally distributed according to the Q-Q plot, however, the Shapiro-Wilk test indicated significant deviations from normality:

hist(residuals(cv.mod4))
qqnorm(residuals(cv.mod4))
shapiro.test(residuals(cv.mod4))

#As such, we performed a glmm with a beta family distribution and logit link transformation (below).

#------------------------------------

  #Experiment 3 and 4
   #glmm

#We performed a generalised linear mixed model in R v4.2.3 with the CV as the response variable, with capture or shotgun (“Sequencing”) and library preparation input and size (“Treatment”) as fixed effects, and pool (“Batch”) as a random effect:

cv.mod5<-glmmTMB(CV ~ Treatment * Sequencing + (1|Batch), family=beta_family(link="logit"), data=evenness_cv_exp34)

#The interaction term between Treatment and Sequencing was tested:

cv.mod6<-glmmTMB(CV ~ Treatment + Sequencing + (1|Batch), family=beta_family(link="logit"), data=evenness_cv_exp34)
anova(cv.mod5, cv.mod6)

#The interaction between Treatment and Sequencing was non-significant and was removed from the model. Next, the Treatment effect was tested:

cv.mod7<-glmmTMB(CV ~ Sequencing + (1|Batch), family=beta_family(link="logit"), data=evenness_cv_exp34)
anova(cv.mod6, cv.mod7)

#Treatment was a non-significant term and could be removed from the model. Next, the Sequencing effect was tested:

cv.mod8<-glmmTMB(CV ~ (1|Batch), family=beta_family(link="logit"), data=evenness_cv_exp34)
anova(cv.mod7, cv.mod8)

#Sequencing was a non-significant term and could be removed from the model.

summary(cv.mod8)

#Model residuals remain approximately normally distributed according to the Q-Q plot, however, the Shapiro-Wilk test still indicated significant deviations from normality:

hist(residuals(cv.mod8))
qqnorm(residuals(cv.mod8))
shapiro.test(residuals(cv.mod8))

#------------------------------------
  #Evenness
    #Experiment 5

evenness_cv_exp5<-read.csv("evenness_cv_exp5.csv",header=TRUE)

#The within-pool co-efficient of variation (CV) was calculated for each sample. We performed a linear model in R v4.2.3 with the CV as the response variable, with capture conditions (“Treatment”) as a fixed effect. Model residuals were normally distributed. Pool could not be included as a random effect because there was only one observation per pool:

cv.mod<-lm(CV ~ Treatment, data= evenness_cv_exp5)
summary(cv.mod)
hist(residuals(cv.mod))
qqnorm(residuals(cv.mod))
shapiro.test(residuals(cv.mod))

#Model estimated means and 95% confidence intervals were calculated:

emm<- emmeans(cv.mod, ~Treatment, lmer.df="satterthwaite")
emm

#------------------------------------
#Enrichment (Aim1)
#------------------------------------
#Target recovery
#lmer

exp1<-read.csv("exp1_w_year.csv", header=TRUE)
head(exp1)

#We performed a linear mixed model with “Treatment” (shotgun or capture) as a fixed effect, “Batch” (pool) as a random effect, and “Sample nested in Batch” as a random effect. “Sample nested in Batch” was removed as an effect because it contributed 0 to variance. Note that phylogenetic effects were not included i.e., tribe and genus because this experiment only covered one tribe and two genera.

exp1.mod<-lmer(Recovery ~ Treatment * Year + (1 | Batch) + (1| Batch:Sample), data=exp1)

VarCorr(exp1.mod)

#Remove Batch:Sample term as it contributes 0 to variance.

exp1.mod2<-lmer(Recovery ~ Treatment * Year + (1 | Batch), data=exp1)
anova(exp1.mod,exp1.mod2)

#Test removing Batch:

exp1.mod3<-lm(Recovery ~ Treatment * Year, data=exp1)
anova(exp1.mod2,exp1.mod3)

#Batch is also a non-significant term and can be removed from the model. 
#Test removing the interaction between Treatment and Year

exp1.mod4<-lm(Recovery ~ Treatment + Year, data=exp1)
anova(exp1.mod3,exp1.mod4)

#The interaction between treatment and year is non-significant and can be removed from the model. 
#Test removing Year:

exp1.mod5<-lm(Recovery ~ Treatment, data=exp1)
anova(exp1.mod4,exp1.mod5)

#Year is non-significant and can be removed from the model. 

AIC(exp1.mod5)
summary(exp1.mod5)

#Model estimated means and 95% confidence intervals were calculated:

emm<- emmeans(exp1.mod5, ~Treatment, lmer.df="satterthwaite")
emm

#Model residuals were approximately normally distributed according to the histogram and Q-Q plot, however, the Shapiro-Wilk test indicated significant deviations from normality:

hist(residuals(exp1.mod5))
qqnorm(residuals(exp1.mod5))
shapiro.test(residuals(exp1.mod5))

#------------------------------------

#Target recovery
#glmm

#Change target recovery to a proportion

exp1$Recovery_prop <- exp1$Recovery / 100

#We performed a glmm with a beta family distribution and logit link transformation in R v4.2.3 with the proportion of target recovered as the response variable (“Recovery”, with capture or shotgun (“Sequencing”) and library preparation input and size (“Treatment”) as fixed effects, and pool (“Batch”) as a random effect:

exp1.mod6<-glmmTMB(Recovery_prop ~ Treatment * Year + (1 | Batch) + (1| Batch:Sample), family=beta_family(link="logit"), data=exp1)

#Test removing Sample nested within Batch:

exp1.mod7<-glmmTMB(Recovery_prop ~ Treatment * Year + (1 | Batch), family=beta_family(link="logit"), data=exp1)
anova(exp1.mod6, exp1.mod7)

#Sample nested within Batch contributes significantly to variance and cannot be removed from the model. Test removing the Batch effect:

exp1.mod8<-glmmTMB(Recovery_prop ~ Treatment * Year + (1| Batch:Sample), family=beta_family(link="logit"), data=exp1)
anova(exp1.mod6, exp1.mod8)

#Batch is non-significant and can be removed from the model.
#Test removing the interaction term between treatment and year:

exp1.mod9<-glmmTMB(Recovery_prop ~ Treatment + Year + (1| Batch:Sample), family=beta_family(link="logit"), data=exp1)
anova(exp1.mod8, exp1.mod9)

#The interaction term between treatment and year is non-significant and can be removed from the model. 
#Test removing Year:

exp1.mod10<-glmmTMB(Recovery_prop ~ Treatment + (1| Batch:Sample), family=beta_family(link="logit"), data=exp1)
anova(exp1.mod9, exp1.mod10)

#Year is insignificant and can be removed from the model.

summary(exp1.mod10)
AIC(exp1.mod10)
emm<- emmeans(exp1.mod10, ~Treatment, type="response")
emm

#Model residuals remain approximately normally distributed according to the Q-Q plot, however, the Shapiro-Wilk test still indicated significant deviations from normality:

hist(residuals(exp1.mod10))
qqnorm(residuals(exp1.mod10))
shapiro.test(residuals(exp1.mod10))

#------------------------------------

#Target Recovery
#Non-parametric

#The absolute difference in the percentage of reads on target between the shotgun and capture datasets was also calculated for each sample (n = 32):

wide <- pivot_wider(exp1, id_cols = c(Sample,Batch),names_from = Treatment, values_from = Recovery)
wide$Diff <-wide$shotgun - wide$capture

#A non-parametric Wilcoxon signed rank test with continuity correction was performed in R v4.2.3: 

wilcox.test(wide$shotgun, wide$capture, exact=FALSE)

#------------------------------------

#Chloroplast recovery
#lmer

#We also performed a linear mixed model with chloroplast recovery (“Chlp”) as the response variable, “Treatment” (shotgun or capture) as a fixed effect, “Batch” (pool) as a random effect, and “Sample nested in Batch” as a random effect. “Sample nested in Batch” was removed as an effect because it contributed 0 to variance. 

exp1.mod11<-lmer(Chlp ~Treatment * Age + (1|Batch) + (1|Batch:Sample), data=exp1)

#Test removing Sample nested within Batch:

exp1.mod12<-lmer(Chlp ~Treatment * Age + (1|Batch), data=exp1)
anova(exp1.mod11, exp1.mod12)

#Sample nested within Batch is a significant term and cannot be removed from the model. Next test removing Batch:

exp1.mod13<-lmer(Chlp ~Treatment * Age + (1|Batch:Sample), data=exp1)
anova(exp1.mod11, exp1.mod13)

#Batch can be removed from the model.
#Test removing the interacton term between treatment and year. 

exp1.mod14<-lmer(Chlp ~Treatment + Age + (1|Batch:Sample), data=exp1)
anova(exp1.mod13, exp1.mod14)

#The interaction between Treatment and Age is significant and cannot be removed from the model.

summary(exp1.mod14)
AIC(exp1.mod14)

#Model estimated means and 95% confidence intervals were calculated:

emm<- emmeans(exp1.mod14, ~Treatment, lmer.df="satterthwaite")
emm

#Model residuals indicated significant deviations from normality:

hist(residuals(exp1.mod14))
qqnorm(residuals(exp1.mod14))
shapiro.test(residuals(exp1.mod14))

#------------------------------------

#Chloroplast recovery
#glmm

exp1$Chlp_prop <- exp1$Chlp / 100

#We performed a glmm with a beta family distribution and logit link transformation in R v4.2.3 with the proportion of chloroplast recovered as the response variable (“Chlp_prop”, as the response variable, “Treatment” (shotgun or capture) as a fixed effect, “Batch” (pool) as a random effect, and “Sample nested in Batch” as a random effect:

exp1.mod15<-glmmTMB(Chlp_prop ~Treatment * Age + (1|Batch) + (1|Batch:Sample), family=beta_family(link="logit"),  data=exp1)

#Test removing Sample nested in Batch from the model:

exp1.mod16<-glmmTMB(Chlp_prop ~Treatment * Age + (1|Batch), family=beta_family(link="logit"),  data=exp1)
anova(exp1.mod15, exp1.mod16)

#Sample nested in Batch cannot be removed from the model. Test removing Batch:

exp1.mod17<-glmmTMB(Chlp_prop ~Treatment * Age + (1|Batch:Sample), family=beta_family(link="logit"),  data=exp1)
anova(exp1.mod15, exp1.mod17)

#Batch can be removed from the model.
#Test removing the interaction between Treatment and Age. 

exp1.mod18<-glmmTMB(Chlp_prop ~Treatment + Age + (1|Batch:Sample), family=beta_family(link="logit"),  data=exp1)
anova(exp1.mod17, exp1.mod18)

#The interaction between treatment and age is non-significant and can be removed from the model. 
#Test removing Age. 

exp1.mod19<-glmmTMB(Chlp_prop ~Treatment + (1|Batch:Sample), family=beta_family(link="logit"),  data=exp1)
anova(exp1.mod18, exp1.mod19)

#Age is significant and cannot be removed from the model. 

AIC(exp1.mod18)
summary(exp1.mod18)
emm<- emmeans(exp1.mod18, ~Treatment, type="response")
emm

#Model residuals still indicate significant deviations from normality:

hist(residuals(exp1.mod18))
qqnorm(residuals(exp1.mod18))
shapiro.test(residuals(exp1.mod18))

#------------------------------------

#Chloroplast recovery
#Non-parametric

#A non-parametric Wilcoxon signed rank test with continuity correction was also performed in R v4.2.3: 

wide <- pivot_wider(exp1, id_cols = c(Sample,Batch),names_from = Treatment, values_from = Chlp)
wilcox.test(wide$shotgun, wide$capture, exact=FALSE)

#------------------------------------
#Between replicates (Aim 2)
#------------------------------------
#Target Recovery
#Non-parametric

exp2<-read.csv("exp2_w_year.csv",header=TRUE)
head(exp2)

#A non-parametric Wilcoxon signed rank test with continuity correction was performed in R v4.2.3 to compare target recovery between replicates within a sample: 

wide <- pivot_wider(exp2, id_cols = c(Sample),names_from = Replicate, values_from = Recovery)
wide
wilcox.test(wide$'1', wide$'2', exact=FALSE)

#------------------------------------

#Target Recovery
#lmer

#To partition the variance components, we performed a linear mixed model with the percentage of target base-pairs recovered (“Recovery”) as the response variable with “Sample”, “Day” (day 1 or day 2) and “Batch” (pool) nested within “Day” as random effects:

exp2.mod <- lmer (Recovery ~ Age + (1|Tribe) + (1|Tribe:Genus) + (1 | Tribe:Genus:Sample) + (1 | Day) + (1 | Day:Batch), data = exp2, REML = TRUE)
VarCorr(exp2.mod)
check_singularity(exp2.mod)
isSingular(exp2.mod)

#Test removing Day from the model:

exp2.mod2 <- lmer (Recovery ~ Age + (1|Tribe) + (1|Tribe:Genus) + (1 | Tribe:Genus:Sample) + (1 | Day:Batch), data = exp2, REML = TRUE)
isSingular(exp2.mod2)
anova(exp2.mod, exp2.mod2)

#Day can be removed from the model. Test removing Tribe from the model:

exp2.mod3 <- lmer (Recovery ~ Age + (1|Tribe:Genus) + (1 | Tribe:Genus:Sample) + (1 | Day:Batch), data = exp2, REML = TRUE)
anova(exp2.mod2, exp2.mod3)

#Tribe is non-significant and can be removed from the model. Test removing Genus nested in Tribe from the model:

exp2.mod4 <- lmer (Recovery ~  Age + (1 | Tribe:Genus:Sample) + (1 | Day:Batch), data = exp2, REML = TRUE)
anova(exp2.mod3, exp2.mod4)

#Genus nested in Tribe significant and cannot be removed from the model. Test removing Sample nested within Genus and Tribe:

exp2.mod5 <- lmer (Recovery ~ Age + (1|Tribe:Genus) + (1 | Day:Batch), data = exp2, REML = TRUE)
anova(exp2.mod3, exp2.mod5)

#Sample nested within Genus and Tribe is significant and cannot be removed from the model. Test removing Batch nested in Day:

exp2.mod6 <- lmer (Recovery ~ Age + (1|Tribe:Genus) + (1 | Tribe:Genus:Sample), data = exp2, REML = TRUE)
anova(exp2.mod3, exp2.mod6)

#Batch nested within Day is significant and cannot be removed from the model.
#Test removing Age:

exp2.mod7 <- lmer (Recovery ~ (1|Tribe:Genus) + (1 | Tribe:Genus:Sample) + (1 | Day:Batch), data = exp2, REML = TRUE)
anova(exp2.mod3, exp2.mod7)

#Age is insignificant and can be removed from the model. 

summary(exp2.mod7)

#The intraclass correlation was calculated using the variance components of the final model:

vc <- as.data.frame(VarCorr(exp2.mod7))
vc
sigma_sample  <- vc$vcov[vc$grp == "Tribe:Genus:Sample"]
sigma_genus <- vc$vcov[vc$grp == "Tribe:Genus"]
sigma_batch   <- vc$vcov[vc$grp == "Day:Batch"]
sigma_resid   <- vc$vcov[vc$grp == "Residual"]
total_var <- sigma_genus + sigma_sample + sigma_batch + sigma_resid
ICC <- (sigma_genus + sigma_sample) / total_var
ICC
AIC(exp2.mod7)

#Model residuals are approximately normally distributed according to the histogram and Q-Q plot, however, the Shapiro-Wilk test still indicated significant deviations from normality:

hist(residuals(exp2.mod7))
qqnorm(residuals(exp2.mod7))
shapiro.test(residuals(exp2.mod7))

#------------------------------------

#Target Recovery
#glmm

#We performed a generalised linear mixed model with the proportion of target base-pairs recovered (“Recovery_prop”) as the response variable with “Sample”, “Day” (day 1 or day 2) and “Batch” (pool) nested within “Day” as random effects. A beta family distribution was used with a logit link transformation:

#Convert percentage target recovery to a proportion between 0 and 1:

exp2$Recovery_prop <- exp2$Recovery / 100
exp2.mod8 <- glmmTMB (Recovery_prop ~ Age + (1|Tribe) + (1|Tribe:Genus) + (1 | Tribe:Genus:Sample) + (1 | Day) + (1 | Day:Batch), family=beta_family(link="logit"), data = exp2)

#Test removing the random effect Tribe:

exp2.mod9 <- glmmTMB (Recovery_prop ~ Age + (1|Tribe:Genus) + (1 | Tribe:Genus:Sample) + (1 | Day) + (1 | Day:Batch), family=beta_family(link="logit"), data = exp2)
anova(exp2.mod8, exp2.mod9)

#Tribe is non-significant and can be removed from the model. Test removing Genus nested within Tribe:

exp2.mod10 <- glmmTMB (Recovery_prop ~ Age + (1 | Tribe:Genus:Sample) + (1 | Day) + (1 | Day:Batch), family=beta_family(link="logit"), data = exp2)
anova(exp2.mod9, exp2.mod10)

#Genus nested in tribe is significant and cannot be removed from the model. Test removing Sample nested within Genus nested within Tribe:

exp2.mod11 <- glmmTMB (Recovery_prop ~ Age + (1|Tribe:Genus) + (1 | Day) + (1 | Day:Batch), family=beta_family(link="logit"), data = exp2)
anova(exp2.mod9, exp2.mod11)

#Sample nested in Genus nested in tribe is significant and cannot be removed from the model. Test removing Day:

exp2.mod12 <- glmmTMB (Recovery_prop ~ Age + (1|Tribe:Genus) + (1 | Tribe:Genus:Sample) + (1 | Day:Batch), family=beta_family(link="logit"), data = exp2)
anova(exp2.mod9, exp2.mod12)

#Day is non-significant and can be removed from the model. Test removing Batch nested within Day:

exp2.mod13 <- glmmTMB (Recovery_prop ~ Age + (1|Tribe:Genus) + (1 | Tribe:Genus:Sample), family=beta_family(link="logit"), data = exp2)
anova(exp2.mod12, exp2.mod13)

#Batch nested within Day is significant and cannot be removed from the model. 
#Test removing Age:

exp2.mod14 <- glmmTMB (Recovery_prop ~  (1|Tribe:Genus) + (1 | Tribe:Genus:Sample) + (1 | Day:Batch), family=beta_family(link="logit"), data = exp2)
anova(exp2.mod12, exp2.mod14)

#Age is non-significant and can be removed from the model. 

summary(exp2.mod14)
AIC(exp2.mod14)

#Model residuals are remain approximately normally distributed according to the histogram and Q-Q plot, however, the Shapiro-Wilk test still indicated significant deviations from normality:

hist(residuals(exp2.mod14))
qqnorm(residuals(exp2.mod14))
shapiro.test(residuals(exp2.mod14))

#------------------------------------

#Chloroplast Recovery
#Non-parametric

#A non-parametric Wilcoxon signed rank test with continuity correction was performed in R v4.2.3 to compare chloroplast recovery between replicates within a sample: 

wide <- pivot_wider(exp2, id_cols = c(Sample),names_from = Replicate, values_from = Chlp)
wide
wilcox.test(wide$'1', wide$'2', exact=FALSE)

#------------------------------------

#Chloroplast Recovery
#lmer 

exp2.mod15 <- lmer(Chlp ~ Age + (1|Tribe) + (1|Tribe:Genus) + (1 | Tribe:Genus:Sample) +(1|Day)+ (1 | Day:Batch), data = exp2,REML = TRUE)

#Test removing Tribe from the model:

exp2.mod16 <- lmer(Chlp ~  Age + (1|Tribe:Genus) + (1 | Tribe:Genus:Sample) +(1|Day)+ (1 | Day:Batch), data = exp2,REML = TRUE)
anova(exp2.mod15, exp2.mod16)

#Tribe is non-significant and can be removed from the model. Test removing Genus nested in Tribe:

exp2.mod17 <- lmer(Chlp ~   Age + (1 | Tribe:Genus:Sample) +(1|Day)+ (1 | Day:Batch), data = exp2,REML = TRUE)
anova(exp2.mod16, exp2.mod17)

#Genus nested in tribe is significant and cannot be removed from the model. Test removing Sample nested in Genus nested in Tribe:

exp2.mod18 <- lmer(Chlp ~  Age + (1|Tribe:Genus)  +(1|Day)+ (1 | Day:Batch), data = exp2,REML = TRUE)
anova(exp2.mod16, exp2.mod18)

#Sample nested in Genus nested in tribe is significant and cannot be removed from the model. Test removing Day:

exp2.mod19 <- lmer(Chlp ~  Age + (1|Tribe:Genus) + (1 | Tribe:Genus:Sample) + (1 | Day:Batch), data = exp2,REML = TRUE)
anova(exp2.mod16, exp2.mod19)

#Day is non-significant and can be removed from the model. Test removing Batch nested in Day:

exp2.mod20 <- lmer(Chlp ~  Age + (1|Tribe:Genus) + (1 | Tribe:Genus:Sample), data = exp2,REML = TRUE)
anova(exp2.mod19, exp2.mod20)

#Batch nested in Day is non-significant and can be removed from the model.
#Test removing Age:

exp2.mod21 <- lmer(Chlp ~   (1|Tribe:Genus) + (1 | Tribe:Genus:Sample), data = exp2,REML = TRUE)
anova(exp2.mod20, exp2.mod21)

#Age is non-significant and can be removed from the model.

summary(exp2.mod21)
AIC(exp2.mod21)

#The intraclass correlation was calculated using the variance components of the final model:

VarCorr(exp2.mod21)
vc <- as.data.frame(VarCorr(exp2.mod21))
vc
sigma_sample  <- vc$vcov[vc$grp == "Tribe:Genus:Sample"]
sigma_genus  <- vc$vcov[vc$grp == "Tribe:Genus"]
total_var <- sigma_genus + sigma_sample + sigma_resid 
ICC <- (sigma_genus + sigma_sample) / total_var
ICC

#Model residuals are approximately normally distributed according to the histogram and Q-Q plot, however, the Shapiro-Wilk test still indicated significant deviations from normality:

hist(residuals(exp2.mod21))
qqnorm(residuals(exp2.mod21))
shapiro.test(residuals(exp2.mod21))

#------------------------------------

#Chloroplast Recovery
#glmm 

#We performed a generalised linear mixed model with the proportion of off-target chloroplast base-pairs recovered (“Chlp_prop”) as the response variable with “Sample”, “Day” (day 1 or day 2) and “Batch” (pool) nested within “Day” as random effects:

exp2$Chlp_prop <- exp2$Chlp / 100
exp2.mod22 <- glmmTMB(Chlp_prop ~ Age + (1|Tribe) + (1|Tribe:Genus) + (1 | Tribe:Genus:Sample) +(1|Day)+ (1 | Day:Batch), family=beta_family(link="logit"), data = exp2, subset = Chlp_prop != 0)

#Test removing Tribe:

exp2.mod23 <- glmmTMB(Chlp_prop ~ Age + (1|Tribe:Genus) + (1 | Tribe:Genus:Sample) +(1|Day)+ (1 | Day:Batch), family=beta_family(link="logit"), data = exp2, subset = Chlp_prop != 0)
anova(exp2.mod22, exp2.mod23)

#Tribe is non-significant and can be removed from the model. Test removing Genus nested in Tribe:

exp2.mod24 <- glmmTMB(Chlp_prop ~  Age + (1 | Tribe:Genus:Sample) +(1|Day)+ (1 | Day:Batch), family=beta_family(link="logit"), data = exp2, subset = Chlp_prop != 0)
anova(exp2.mod23, exp2.mod24)

#Genus nested in Tribe is non-significant and can be removed from the model. Test removing Sample nested in Genus nested in Tribe:

exp2.mod25 <- glmmTMB(Chlp_prop ~ Age + (1|Day)+ (1 | Day:Batch), family=beta_family(link="logit"), data = exp2, subset = Chlp_prop != 0)
anova(exp2.mod24, exp2.mod25)

#Sample nested in genus nested in Tribe is non-significant and can be removed from the model. Test removing Day:

exp2.mod26 <- glmmTMB(Chlp_prop ~ Age  + (1 | Day:Batch), family=beta_family(link="logit"), data = exp2, subset = Chlp_prop != 0)
anova(exp2.mod25, exp2.mod26)

#Day is non-significant and can be removed from the model. Test removing Batch within Day:

exp2.mod27 <- glmmTMB(Chlp_prop ~ Age, family=beta_family(link="logit"), data = exp2, subset = Chlp_prop != 0)
anova(exp2.mod26, exp2.mod27)

#Batch within Day is non-significant and can be removed from the model. 

summary(exp2.mod27)
AIC(exp2.mod27)

#Model residuals seem to be even less normally distributed according to the histogram and Q-Q plot, and the Shapiro-Wilk test still indicated significant deviations from normality:

hist(residuals(exp2.mod27))
qqnorm(residuals(exp2.mod27))
shapiro.test(residuals(exp2.mod27))

#------------------------------------
#Between treatments (Aim 3)
#------------------------------------

exp3<-read.csv("exp3_w_year.csv", header=TRUE)
head(exp3)

#Target recovery
#lmer

#We performed a linear mixed model with the percentage target base-pairs recovered (“Ast353”) as the response variable with “Size” (mini or standard) as a fixed factor, “Tribe”, “Genus within Tribe” and “Sample nested in Genus nested in Tribe” as a random effects, and “Batch” (pool) as a random effect: 

exp3.mod<-lmer(Ast353 ~ Size * Age + (1|Tribe) + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp3, REML=TRUE)

#Test removing Tribe from the model:

exp3.mod2<-lmer(Ast353 ~ Size * Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp3, REML=TRUE)
anova(exp3.mod, exp3.mod2)

#Tribe does not contribute to variance so is removed from the model.Genus within Tribe is then removed from the model:

exp3.mod3<-lmer(Ast353 ~ Size * Age + (1|Tribe:Genus:Sample) + (1|Batch), data=exp3, REML=TRUE)
anova(exp3.mod2, exp3.mod3)

#Genus also does not significantly contribute to variance. The Sample term is simplified:

exp3.mod4<-lmer(Ast353 ~ Size * Age + (1|Sample) + (1|Batch), data=exp3, REML=TRUE)
anova(exp3.mod3, exp3.mod4)

#Test removing the factor Sample:

exp3.mod5<-lmer(Ast353 ~ Size * Age  + (1|Batch), data=exp3, REML=TRUE)
anova(exp3.mod4, exp3.mod5)

#Sample cannot be removed from the model. Test removing Batch:

exp3.mod6<-lmer(Ast353 ~ Size * Age + (1|Sample), data=exp3, REML=TRUE)
anova(exp3.mod4, exp3.mod6)

#Batch cannot be removed from the model.
#Test removing the interaction term between size and age:

exp3.mod7<-lmer(Ast353 ~ Size + Age + (1|Sample) + (1|Batch), data=exp3, REML=TRUE)
anova(exp3.mod4, exp3.mod7)

#The interaction between size and age is non-significant and can be removed from the model. 
#Test removing Age:

exp3.mod8<-lmer(Ast353 ~ Size + (1|Sample) + (1|Batch), data=exp3, REML=TRUE)
anova(exp3.mod7, exp3.mod8)

#Age is significant and cannot be removed from the model.

summary(exp3.mod7)


#Library reaction size does not have a significant effect on target recovery.
#The model estimated marginal means were calculated for each library size treatment:

emm<- emmeans(exp3.mod7, ~Age | Size, at = list(Age = seq(min(exp3$Age), max(exp3$Age), length=100)), lmer.df="satterthwaite")
emm
emm_df<-as.data.frame(emm)
emm_df

ggplot() +
  geom_point(data = exp3,
             aes(x = Age, y = Ast353, colour = Size),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Age, ymin = lower.CL, ymax = upper.CL, fill = Size),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Age, y = emmean, colour = Size),
            linewidth = 1) +
  theme_classic()
      
#The variance components were estimated:

vc <- as.data.frame(VarCorr(exp3.mod7))
vc$percent <- vc$vcov / sum(vc$vcov) * 100
vc
AIC(exp3.mod7)

#Model residuals are normally distributed according to the histogram, Q-Q plot, and the Shapiro-Wilk test:

hist(residuals(exp3.mod7))
qqnorm(residuals(exp3.mod7))
shapiro.test(residuals(exp3.mod7))

#------------------------------------

#Target recovery
#glmm

exp3$Ast353_prop <- exp3$Ast353 / 100

#We performed a generalised linear mixed model with the proportion of target base-pairs recovered (“Ast353_prop”) as the response variable with “Size” (mini or standard) as a fixed factor, “Tribe”, “Genus within Tribe” and “Sample nested in Genus nested in Tribe” as a random effects, and “Batch” (pool) as a random effect: 

exp3.mod9<-glmmTMB(Ast353_prop ~ Size * Age + (1|Tribe) + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp3)

#Test removing Tribe from the model:

exp3.mod10<-glmmTMB(Ast353_prop ~ Size * Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp3)
anova(exp3.mod9, exp3.mod10)

#Tribe is non-significant and can be removed from the model. Test removing Genus nested within Tribe:

exp3.mod11<-glmmTMB(Ast353_prop ~ Size * Age + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp3)
anova(exp3.mod10, exp3.mod11)

#Genus nested in Tribe is also non-significant and can be removed from the model. Test removing Sample nested within Genus nested within Tribe:

exp3.mod12<-glmmTMB(Ast353_prop ~ Size * Age  + (1|Batch), family=beta_family(link="logit"), data=exp3)
anova(exp3.mod11, exp3.mod12)

#Sample nested within Genus nested in Tribe is significant and cannot be removed from the model. Test removing Batch:

exp3.mod13<-glmmTMB(Ast353_prop ~ Size * Age + (1|Tribe:Genus:Sample), family=beta_family(link="logit"), data=exp3)
anova(exp3.mod11, exp3.mod13)

#Batch is significant and cannot be removed from the model. 
#Test removing the interaction between size and age:

exp3.mod14<-glmmTMB(Ast353_prop ~ Size + Age + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp3)
anova(exp3.mod11, exp3.mod14)

#The interacton between size and age is insignificant and can be removed from the model. 
#Test removing age:

exp3.mod15<-glmmTMB(Ast353_prop ~ Size  + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp3)
anova(exp3.mod14, exp3.mod15)

#Age is significant and cannot be removed from the model. 

summary(exp3.mod14)

#Library size does not significantly influence target recovery.

emm <- emmeans(exp3.mod14, ~Size + Age, type="response")
emm
pairs(emm, adjust="tukey")
AIC(exp3.mod14)

emm<- emmeans(exp3.mod14, ~Age | Size, at = list(Age = seq(min(exp3$Age), max(exp3$Age), length=100)), type="response")
emm
emm_df<-as.data.frame(emm)
emm_df

ggplot() +
  geom_point(data = exp3,
             aes(x = Age, y = Ast353_prop, colour = Size),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Age, ymin = asymp.LCL, ymax = asymp.UCL, fill = Size),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Age, y = response, colour = Size),
            linewidth = 1) +
  theme_classic()

#Model residuals are now non-normally distributed according to the Shapiro-Wilk test, though the histogram and Q-Q plot show residuals are approximately normally distributed:

hist(residuals(exp3.mod14))
qqnorm(residuals(exp3.mod14))
shapiro.test(residuals(exp3.mod14))

#------------------------------------

#Chloroplast recovery
#lmer

#We performed a linear mixed model with the percentage of chloroplast base pairs recovered (“Chlp”) as the response variable with “Size” (mini or standard) as a fixed factor, “Tribe”, “Genus within Tribe” and “Sample nested in Genus nested in Tribe” as a random effects, and “Batch” (pool) as a random effect: 

exp3.mod16<-lmer(Chlp ~ Size * Age + (1|Tribe) + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp3, REML=TRUE)

#Test removing Tribe as a factor:

exp3.mod17<-lmer(Chlp ~ Size * Age  + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp3, REML=TRUE)
anova(exp3.mod16, exp3.mod17)

#Tribe does not significantly contribute to the model. Next, remove Genus within Tribe:

exp3.mod18<-lmer(Chlp ~ Size * Age  +  (1|Tribe:Genus:Sample) + (1|Batch), data=exp3, REML=TRUE)
anova(exp3.mod17, exp3.mod18)

#Genus within Tribe cannot be removed. Next, test removing Sample from the model:

exp3.mod19<-lmer(Chlp ~ Size * Age  + (1|Tribe:Genus)  + (1|Batch), data=exp3, REML=TRUE)
anova(exp3.mod17, exp3.mod19)

#Sample cannot be removed. Next, test removing Batch from the model:

exp3.mod20<-lmer(Chlp ~ Size * Age  + (1|Tribe:Genus) + (1|Tribe:Genus:Sample), data=exp3, REML=TRUE)
anova(exp3.mod17, exp3.mod20)

#Batch cannot be removed.
#Test removing the interaction between size and age:

exp3.mod21<-lmer(Chlp ~ Size + Age  + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp3, REML=TRUE)
anova(exp3.mod17, exp3.mod21)

#The interaction between size and age is significant and cannot be removed.

summary(exp3.mod17)
vc <- as.data.frame(VarCorr(exp3.mod17))
vc$percent <- vc$vcov / sum(vc$vcov) * 100
vc
emm<- emmeans(exp3.mod17, ~Size * Age, lmer.df="satterthwaite")
emm
pairs(emm, adjust="none")
AIC(exp3.mod17)

emm<- emmeans(exp3.mod17, ~Age | Size, at = list(Age = seq(min(exp3$Age), max(exp3$Age), length=100)), lmer.df="satterthwaite")
emm
emm_df<-as.data.frame(emm)
emm_df

ggplot() +
  geom_point(data = exp3,
             aes(x = Age, y = Chlp, colour = Size),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Age, ymin = lower.CL, ymax = upper.CL, fill = Size),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Age, y = emmean, colour = Size),
            linewidth = 1) +
  theme_classic()

#Model residuals are approximately normally distributed according to the histogram, Q-Q plot, but the Shapiro-Wilk test indicates that statistically they are non-normal:

hist(residuals(exp3.mod17))
qqnorm(residuals(exp3.mod17))
shapiro.test(residuals(exp3.mod17))

#------------------------------------

#Chloroplast recovery    
#glmm

#We performed a linear mixed model with the proportion of chloroplast base pairs recovered (“Chlp_prop”) as the response variable with “Size” (mini or standard) as a fixed factor, “Tribe”, “Genus within Tribe” and “Sample nested in Genus nested in Tribe” as a random effects, and “Batch” (pool) as a random effect: 

exp3$Chlp_prop <- exp3$Chlp / 100
exp3.mod22<-glmmTMB(Chlp_prop ~ Size * Age + (1|Tribe) + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp3)

#Test removing Tribe as a factor:

exp3.mod23<-glmmTMB(Chlp_prop ~ Size * Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp3)
anova(exp3.mod22, exp3.mod23)

#Tribe is non-significant and can be removed from the model. Test removing Genus nested in Tribe:

exp3.mod24<-glmmTMB(Chlp_prop ~ Size * Age + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp3)
anova(exp3.mod23, exp3.mod24)

#Genus nested in Tribe is significant and cannot be removed from the model. Test removing Sample nested in Genus nested in Tribe:

exp3.mod25<-glmmTMB(Chlp_prop ~ Size * Age + (1|Tribe:Genus) + (1|Batch), family=beta_family(link="logit"), data=exp3)
anova(exp3.mod23, exp3.mod25)

#Sampled nested in Genus nested in Tribe is significant and cannot be removed from the model. Test removing Batch:

exp3.mod26<-glmmTMB(Chlp_prop ~ Size * Age  + (1|Tribe:Genus) + (1|Tribe:Genus:Sample), family=beta_family(link="logit"), data=exp3)
anova(exp3.mod23, exp3.mod26)

#Batch is significant and cannot be removed from the model. 
#Test removing the interaction between size and age:

exp3.mod27<-glmmTMB(Chlp_prop ~ Size + Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp3)
anova(exp3.mod23, exp3.mod27)

#The interaction betwen size and age is significant and can't be removed from the model.

summary(exp3.mod23)
emm <- emmeans(exp3.mod23, ~Size, type="response")
emm
pairs(emm, adjust="none")
AIC(exp3.mod23)

emm<- emmeans(exp3.mod23, ~Age | Size, at = list(Age = seq(min(exp3$Age), max(exp3$Age), length=100)), type="response")
emm
emm_df<-as.data.frame(emm)
emm_df

ggplot() +
  geom_point(data = exp3,
             aes(x = Age, y = Chlp_prop, colour = Size),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Age, ymin = asymp.LCL, ymax = asymp.UCL, fill = Size),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Age, y = response, colour = Size),
            linewidth = 1) +
  theme_classic()

#Model residuals remain approximately normally distributed according to the histogram, Q-Q plot, but the Shapiro-Wilk test indicates that statistically they are non-normal:

hist(residuals(exp3.mod23))
qqnorm(residuals(exp3.mod23))
shapiro.test(residuals(exp3.mod23))

#------------------------------------
#Between treatments (Aim 4)
#------------------------------------

#Target recovery
#lmer

exp4<-read.csv("exp4_w_year.csv", header=TRUE)
head(exp4)

#We performed a linear mixed model with the percentage target base-pairs recovered (“Ast353”) as the response variable with “Input” (2 ng, 5 ng, 20 ng, 40 ng) as a fixed, categorical factor, “Tribe”, “Genus within Tribe”, “Sample nested in tribe nested in Genus” and “Batch” (pool) as random effects:

exp4.mod<-lmer(Ast353 ~ Input_cat * Age + (1|Tribe) + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp4, REML=TRUE)

#Test removing Tribe as a factor:

exp4.mod2<-lmer(Ast353 ~ Input_cat * Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp4, REML=TRUE)
anova(exp4.mod, exp4.mod2)

#Tribe can be removed as it contributes nothing to variance. Now test removing Genus within Tribe:

exp4.mod3<-lmer(Ast353 ~ Input_cat * Age  + (1|Tribe:Genus:Sample) + (1|Batch), data=exp4, REML=TRUE)
anova(exp4.mod2, exp4.mod3)

#Genus within Tribe can be removed. Simplify the sample term:

exp4.mod4<-lmer(Ast353 ~ Input_cat * Age  + (1|Sample) + (1|Batch), data=exp4, REML=TRUE)
anova(exp4.mod3, exp4.mod4)

#Now test removing Sample from the model:

exp4.mod5<-lmer(Ast353 ~ Input_cat * Age  + (1|Batch), data=exp4, REML=TRUE)
anova(exp4.mod4, exp4.mod5)

#Sample cannot be removed from the model. Now try removing Batch:

exp4.mod6<-lmer(Ast353 ~ Input_cat * Age  + (1|Sample) , data=exp4, REML=TRUE)
anova(exp4.mod4, exp4.mod6)

#Batch cannot be removed from the model. 
#Test removing the interaction term between Input_cat and age:

exp4.mod7<-lmer(Ast353 ~ Input_cat + Age  + (1|Sample) + (1|Batch), data=exp4, REML=TRUE)
anova(exp4.mod4, exp4.mod7)

#The interaction term between input_cat and age is not significant and can be removed.
#Test removing age:

exp4.mod8<-lmer(Ast353 ~ Input_cat  + (1|Sample) + (1|Batch), data=exp4, REML=TRUE)
anova(exp4.mod7, exp4.mod8)

#Age is significant and cannot be removed from the model. 
#Test removing Input_cat:

exp4.mod9<-lmer(Ast353 ~  Age  + (1|Sample) + (1|Batch), data=exp4, REML=TRUE)
anova(exp4.mod7, exp4.mod9)

#Input_cat is significant and cannot be removed from the model.

summary(exp4.mod7)
AIC(exp4.mod7)

#Model residuals are approximately normally distributed according to the histogram and Q-Q plot, but the Shapiro-Wilk test indicates that statistically they are non-normal:

hist(residuals(exp4.mod7))
qqnorm(residuals(exp4.mod7))
shapiro.test(residuals(exp4.mod7))

#Calculate the variance components:

vc <- as.data.frame(VarCorr(exp4.mod7))
vc$percent <- vc$vcov / sum(vc$vcov) * 100
vc

#The model estimated marginal means were calculated for each library input treatment:

emm<-emmeans(exp4.mod7, ~Input_cat, lmer.df="satterthwaite")
emm

#Pairwise contrasts were performed between input treatment groups with a Tukey correction:

pairs(emm, adjust="tukey")

#Plot the model:

emm<- emmeans(exp4.mod7, ~Age + Input_cat, at = list(Age = seq(min(exp3$Age), max(exp3$Age), length=100)), lmer.df="satterthwaite")
emm
emm_df<-as.data.frame(emm)
emm_df

ggplot() +
  geom_point(data = exp4,
             aes(x = Age, y = Ast353, colour = Input_cat),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Age, ymin = lower.CL, ymax = upper.CL, fill = Input_cat),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Age, y = emmean, colour = Input_cat),
            linewidth = 1) +
  theme_classic()

#We also performed a linear mixed model with the percentage target base-pairs recovered (“Ast353”) as the response variable with “Input” (2, 5, 20, 40) as a fixed, continuous variable, and “Sample” and “Batch” (pool) as random effects:

exp4$Input_cont <- as.numeric(exp4$Input_cont)
exp4.mod10<-lmer(Ast353 ~ Input_cont * Age + (1|Sample)+ (1|Batch),data=exp4, REML=TRUE)

#Sample and batch were presumed to be significant random factors as in the previous analysis.
#Test removing the interaction between Input_cont and age:

exp4.mod11<-lmer(Ast353 ~ Input_cont + Age + (1|Sample)+ (1|Batch),data=exp4, REML=TRUE)
anova(exp4.mod10, exp4.mod11)

#The interaction between Input_cont and age is non-significant and can be removed from model.
#Test removing age:

exp4.mod12<-lmer(Ast353 ~ Input_cont + (1|Sample)+ (1|Batch),data=exp4, REML=TRUE)
anova(exp4.mod11, exp4.mod12)

#Age is significant and cannot be removed from the model. 
#Test removing Input_cont:

exp4.mod13<-lmer(Ast353 ~ Age + (1|Sample)+ (1|Batch),data=exp4, REML=TRUE)
anova(exp4.mod11, exp4.mod13)

#Input_cont is non-significant and can be removed from the model.

summary(exp4.mod13)

emm<- emmeans(exp4.mod13, ~Age, at = list(Age = seq(min(exp4$Age), max(exp4$Age), length=100)), lmer.df="satterthwaite")
emm
emm_df<-as.data.frame(emm)
emm_df

ggplot() +
  geom_point(data = exp4,
             aes(x = Age, y = Ast353),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Age, ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Age, y = emmean),
            linewidth = 1) +
  theme_classic()

#Model residuals are approximately normally distributed according to the histogram and Q-Q plot, but the Shapiro-Wilk test indicates that statistically they are non-normal:

hist(residuals(exp4.mod13))
qqnorm(residuals(exp4.mod13))
shapiro.test(residuals(exp4.mod13))

#------------------------------------

#Target recovery
#glmm

#We performed a generalised linear mixed model with the proportion target base-pairs recovered (“Ast353_prop”) as the response variable with “Input” (2 ng, 5 ng, 20 ng, 40 ng) as a fixed, categorical factor, “Tribe”, “Genus within Tribe”, “Sample nested in tribe nested in Genus” and “Batch” (pool) as random effects:

exp4$Ast353_prop <- exp4$Ast353 / 100
exp4.mod14<-glmmTMB(Ast353_prop ~ Input_cat * Age + (1|Tribe) + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp4)

#Test removing Tribe:

exp4.mod15<-glmmTMB(Ast353_prop ~ Input_cat * Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp4)
anova(exp4.mod14, exp4.mod15)

#Tribe is non-significant and can be removed from the model. Test removing Genus nested within Tribe:

exp4.mod16<-glmmTMB(Ast353_prop ~ Input_cat * Age + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp4)
anova(exp4.mod15, exp4.mod16)

#Genus nested within Tribe is non-significant and can be removed from the model. Test removing Sample nested in Genus nested within Tribe:

exp4.mod17<-glmmTMB(Ast353_prop ~ Input_cat * Age + (1|Batch), family=beta_family(link="logit"), data=exp4)
anova(exp4.mod16, exp4.mod17)

#Sample nested in genus nested within Tribe is significant and cannot be removed from the model. Test removing Batch:

exp4.mod18<-glmmTMB(Ast353_prop ~ Input_cat * Age + (1|Tribe:Genus:Sample), family=beta_family(link="logit"), data=exp4)
anova(exp4.mod16, exp4.mod18)

#Batch is significant and cannot be removed from the model. 
#Test removing the interaction between Input_cat and age:

exp4.mod19<-glmmTMB(Ast353_prop ~ Input_cat + Age + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp4)
anova(exp4.mod16, exp4.mod19)

#The interaction between Input_cat and age is not significant and can be removed. 
#Test removing Age:

exp4.mod20<-glmmTMB(Ast353_prop ~ Input_cat + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp4)
anova(exp4.mod19, exp4.mod20)

#Age is significant and cannot be removed from the model.

summary(exp4.mod19)
AIC(exp4.mod19)
emm <- emmeans(exp4.mod19, ~Input_cat, type = "response")
emm
pairs(emm, adjust = "tukey")

#Plot model
emm<- emmeans(exp4.mod19, ~Age | Input_cat, at = list(Age = seq(min(exp4$Age), max(exp4$Age), length=100)), type="response")
emm
emm_df<-as.data.frame(emm)
emm_df

ggplot() +
  geom_point(data = exp4,
             aes(x = Age, y = Ast353_prop, colour = Input_cat),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Age, ymin = asymp.LCL, ymax = asymp.UCL, fill = Input_cat),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Age, y = response, colour = Input_cat),
            linewidth = 1) +
  theme_classic()

#Model residuals remain approximately normally distributed according to the histogram and Q-Q plot, but the Shapiro-Wilk test indicates that statistically they are still non-normal:

hist(residuals(exp4.mod19))
qqnorm(residuals(exp4.mod19))
shapiro.test(residuals(exp4.mod19))

#------------------------------------

#Chloroplast recovery
#lmer

#We performed a linear mixed model with the percentage chloroplast base-pairs recovered (“Chlp”) as the response variable with “Input” (2 ng, 5 ng, 20 ng, 40 ng) as a fixed, categorical factor, “Tribe”, “Genus within Tribe”, “Sample nested in tribe nested in Genus” and “Batch” (pool) as random effects:

exp4.mod21<-lmer(Chlp ~ Input_cat * Age  + (1|Tribe) +(1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp4, REML=TRUE)

#Test removing Tribe from the model:

exp4.mod22<-lmer(Chlp ~ Input_cat * Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp4, REML=TRUE)
anova(exp4.mod21, exp4.mod22)

#Tribe can be removed from the model. Next, test removing Genus within Tribe:

exp4.mod23<-lmer(Chlp ~ Input_cat * Age + (1|Tribe:Genus:Sample) + (1|Batch), data=exp4, REML=TRUE)
anova(exp4.mod22, exp4.mod23)

#Genus within Tribe cannot be removed from the model. Next test removing Sample from the model:

exp4.mod24<-lmer(Chlp ~ Input_cat * Age + (1|Tribe:Genus) + (1|Batch), data=exp4, REML=TRUE)
anova(exp4.mod22, exp4.mod24)

#Sample cannot be removed from the model. Next test removing Batch from the model:

exp4.mod25<-lmer(Chlp ~ Input_cat * Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample), data=exp4, REML=TRUE)
anova(exp4.mod22, exp4.mod25)

#Batch cannot be removed from the model. 
#Test removing the interaction term between Input_cat and age:

exp4.mod26<-lmer(Chlp ~ Input_cat + Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp4, REML=TRUE)
anova(exp4.mod22, exp4.mod26)

#The interaction term between input_cat and age is significant and cannot be removed from the model.

summary(exp4.mod22)

#The model estimated marginal means were calculated for each library input treatment:

emm<- emmeans(exp4.mod22, ~Input_cat, lmer.df="satterthwaite")
emm

#Pairwise contrasts were performed between input treatment groups with a Tukey correction:

pairs(emm, adjust="tukey")

#Calculate the variance components:

vc <- as.data.frame(VarCorr(exp4.mod22))
vc$percent <- vc$vcov / sum(vc$vcov) * 100
vc
AIC(exp4.mod22)

#Plot model
emm<- emmeans(exp4.mod22, ~Age | Input_cat, at = list(Age = seq(min(exp4$Age), max(exp4$Age), length=100)), lmer.df="satterthwaite")
emm
emm_df<-as.data.frame(emm)
emm_df

ggplot() +
  geom_point(data = exp4,
             aes(x = Age, y = Chlp, colour = Input_cat),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Age, ymin = lower.CL, ymax = upper.CL, fill = Input_cat),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Age, y = emmean, colour = Input_cat),
            linewidth = 1) +
  theme_classic()

#Model residuals are approximately normally distributed according to the histogram and Q-Q plot, but the Shapiro-Wilk test indicates that statistically they are non-normal:

hist(residuals(exp4.mod22))
qqnorm(residuals(exp4.mod22))
shapiro.test(residuals(exp4.mod22))

#------------------------------------

#Chloroplast recovery
#glmm

#We performed a generalised linear mixed model with the proportion chloroplast base-pairs recovered (“Chlp_prop”) as the response variable with “Input” (2 ng, 5 ng, 20 ng, 40 ng) as a fixed, categorical factor, “Tribe”, “Genus within Tribe”, “Sample nested in tribe nested in Genus” and “Batch” (pool) as random effects:

exp4$Chlp_prop <- exp4$Chlp / 100
exp4.mod27<-glmmTMB(Chlp_prop ~ Input_cat * Age  + (1|Tribe) +(1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp4)

#Test removing Tribe as a factor:

exp4.mod28<-glmmTMB(Chlp_prop ~ Input_cat * Age +(1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp4)
anova(exp4.mod27, exp4.mod28)

#Tribe is non-significant and can be removed as a factor. Test removing Genus nested Tribe:

exp4.mod29<-glmmTMB(Chlp_prop ~ Input_cat * Age + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp4)
anova(exp4.mod28, exp4.mod29)

#Genus nested in Tribe is significant and cannot be removed as a factor. Test removing Sample nested in Genus nested Tribe:

exp4.mod30<-glmmTMB(Chlp_prop ~ Input_cat * Age +(1|Tribe:Genus)+ (1|Batch), family=beta_family(link="logit"), data=exp4)
anova(exp4.mod28, exp4.mod30)

#Sample nested in Genus nested in Tribe is significant and cannot be removed as a factor. Test removing Batch:

exp4.mod31<-glmmTMB(Chlp_prop ~ Input_cat  * Age+(1|Tribe:Genus) + (1|Tribe:Genus:Sample), family=beta_family(link="logit"), data=exp4)
anova(exp4.mod28, exp4.mod31)

#Batch significant and cannot be removed as a factor. 
#Test removing the interacton between Input_cat and age:

exp4.mod32<-glmmTMB(Chlp_prop ~ Input_cat + Age +(1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp4)
anova(exp4.mod28, exp4.mod32)

#The interaction between input_cat and age is significant and cannot be removed from the model.

summary(exp4.mod28)
emm <- emmeans(exp4.mod28, ~Input_cat, type="response")
emm
AIC(exp4.mod28)

#Plot model
emm<- emmeans(exp4.mod28, ~Age | Input_cat, at = list(Age = seq(min(exp4$Age), max(exp4$Age), length=100)), type="response")
emm
emm_df<-as.data.frame(emm)
emm_df

ggplot() +
  geom_point(data = exp4,
             aes(x = Age, y = Chlp_prop, colour = Input_cat),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Age, ymin = asymp.LCL, ymax = asymp.UCL, fill = Input_cat),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Age, y = response, colour = Input_cat),
            linewidth = 1) +
  theme_classic()

#Model residuals remain approximately normally distributed according to the histogram and Q-Q plot, but the Shapiro-Wilk test indicates that statistically they are still non-normal:

hist(residuals(exp4.mod28))
qqnorm(residuals(exp4.mod28))
shapiro.test(residuals(exp4.mod28))

#------------------------------------
#Between treatments (Aim 5)
#------------------------------------

exp5<-read.csv("exp5_w_year.csv", header=TRUE)
head(exp5)

#Target recovery
#lmer

#We performed a linear mixed model with the percentage target base-pairs recovered (“Ast353”) as the response variable with “Treatment” (16rxn, 32rxn, half_bait) as a fixed, categorical factor, and “Sample” and “Pool” (batch) as random effects:

exp5.mod<-lmer(Ast353 ~ Treatment * Age + (1|Tribe) + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Pool), data=exp5, REML=TRUE)

#Remove Tribe as a random factor:

exp5.mod2<-lmer(Ast353 ~ Treatment* Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Pool), data=exp5, REML=TRUE)
anova(exp5.mod, exp5.mod2)

#Tribe is nonsignificant as a factor. Remove Genus nested within Tribe. 

exp5.mod3<-lmer(Ast353 ~ Treatment * Age + (1|Tribe:Genus:Sample) + (1|Pool), data=exp5, REML=TRUE)
anova(exp5.mod2, exp5.mod3)

#Genus is a nonsignificant factor. Simplify the terms:

exp5.mod4<-lmer(Ast353 ~ Treatment * Age + (1|Sample) + (1|Pool), data=exp5, REML=TRUE)
anova(exp5.mod3, exp5.mod4)

#Model 3 = Model 4. Remove Sample:

exp5.mod5<-lmer(Ast353 ~ Treatment * Age+ (1|Pool), data=exp5, REML=TRUE)
anova(exp5.mod4, exp5.mod5)

#Sample is a significant factor and can’t be removed. Remove Pool as a factor:

exp5.mod6<-lmer(Ast353 ~ Treatment * Age + (1|Pool), data=exp5, REML=TRUE)
anova(exp5.mod4, exp5.mod6)

#Pool is a significant factor and can’t be removed. 
#Remove the interaction between treatment and age:

exp5.mod7<-lmer(Ast353 ~ Treatment + Age + (1|Sample) + (1|Pool), data=exp5, REML=TRUE)
anova(exp5.mod4, exp5.mod7)

#The interaction between treatment and age is significant and cannot be removed.

summary(exp5.mod4)
AIC(exp5.mod4)
vc <- as.data.frame(VarCorr(exp5.mod4))
vc$percent <- vc$vcov / sum(vc$vcov) * 100
vc
emm<-emmeans(exp5.mod4, ~Treatment, lmer.df="satterthwaite")
pairs(emm, adjust="tukey")

#Plot model
emm<- emmeans(exp5.mod4, ~Age | Treatment, at = list(Age = seq(min(exp5$Age), max(exp5$Age), length=100)), lmer.df="satterthwaite")
emm
emm_df<-as.data.frame(emm)
emm_df

ggplot() +
  geom_point(data = exp5,
             aes(x = Age, y = Ast353, colour = Treatment),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Age, ymin = lower.CL, ymax = upper.CL, fill = Treatment),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Age, y = emmean, colour = Treatment),
            linewidth = 1) +
  theme_classic()

#Model residuals are approximately normally distributed according to the histogram and Q-Q plot, but the Shapiro-Wilk test indicates that statistically they are non-normal:

qqnorm(residuals(exp5.mod4))
hist(residuals(exp5.mod4))
shapiro.test(residuals(exp5.mod4))

#------------------------------------

#Target recovery
#glmm

#We performed a generalised linear mixed model with the proportion of target base-pairs recovered (“Ast353_prop”) as the response variable with “Treatment” (16rxn, 32rxn, half_bait) as a fixed, categorical factor, and “Sample” and “Pool” (batch) as random effects:

exp5$Ast353_prop <- exp5$Ast353/100
exp5.mod8<-glmmTMB(Ast353_prop ~ Treatment * Age + (1|Tribe) + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Pool), family=beta_family(link="logit"), data=exp5, subset = Ast353_prop != 0)

#Remove Tribe as a factor:

exp5.mod9<-glmmTMB(Ast353_prop ~ Treatment * Age+ (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Pool), family=beta_family(link="logit"), data=exp5, subset = Ast353_prop != 0)
anova(exp5.mod8, exp5.mod9)

#Tribe is non-significant and can be removed as a factor. Remove Genus nested within Tribe:

exp5.mod10<-glmmTMB(Ast353_prop ~ Treatment * Age+ (1|Tribe:Genus:Sample) + (1|Pool), family=beta_family(link="logit"), data=exp5, subset = Ast353_prop != 0)
anova(exp5.mod9, exp5.mod10)

#Genus is non-significant and can be removed as a factor. Remove Sample nested in Genus and Tribe:

exp5.mod11<-glmmTMB(Ast353_prop ~ Treatment * Age+ (1|Pool), family=beta_family(link="logit"), data=exp5, subset = Ast353_prop != 0)
anova(exp5.mod10, exp5.mod11)

#Sample is significant and cannot be removed as a factor. Remove Pool:

exp5.mod12<-glmmTMB(Ast353_prop ~ Treatment * Age + (1|Tribe:Genus:Sample), family=beta_family(link="logit"), data=exp5, subset = Ast353_prop != 0)
anova(exp5.mod10, exp5.mod12)

#Pool is significant and cannot be removed as a factor. 
#Remove the interaction between treatment and age:

exp5.mod13<-glmmTMB(Ast353_prop ~ Treatment + Age+ (1|Tribe:Genus:Sample) + (1|Pool), family=beta_family(link="logit"), data=exp5, subset = Ast353_prop != 0)
anova(exp5.mod10, exp5.mod13)

#The interaction between treatment and age is non-significant and can be removed. 
#Remove age:

exp5.mod14<-glmmTMB(Ast353_prop ~ Treatment + (1|Tribe:Genus:Sample) + (1|Pool), family=beta_family(link="logit"), data=exp5, subset = Ast353_prop != 0)
anova(exp5.mod13, exp5.mod14)

#Age is significant and cannot be removed from the model.
#Remove treatment:

exp5.mod15<-glmmTMB(Ast353_prop ~ Age+ (1|Tribe:Genus:Sample) + (1|Pool), family=beta_family(link="logit"), data=exp5, subset = Ast353_prop != 0)
anova(exp5.mod13, exp5.mod15)

#Treatment is non-significant and can be removed from the model.

summary(exp5.mod13)
emm <- emmeans(exp5.mod13, ~Treatment, type="response")
emm
pairs(emm, adjust="tukey")

#Plot model
emm<- emmeans(exp5.mod13, ~Age | Treatment, at = list(Age = seq(min(exp5$Age), max(exp5$Age), length=100)), type="response")
emm
emm_df<-as.data.frame(emm)
emm_df

ggplot() +
  geom_point(data = exp5,
             aes(x = Age, y = Ast353_prop, colour = Treatment),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Age, ymin = asymp.LCL, ymax = asymp.UCL, fill = Treatment),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Age, y = response, colour = Treatment),
            linewidth = 1) +
  theme_classic()

#Model residuals remain approximately normally distributed according to the histogram and Q-Q plot, but the Shapiro-Wilk test indicates that statistically they are still non-normal:

hist(residuals(exp5.mod13))
qqnorm(residuals(exp5.mod13))
shapiro.test(residuals(exp5.mod13))

AIC(exp5.mod13)

#------------------------------------

#Chloroplast recovery
#lmer

#We performed a linear mixed model with the percentage of chloroplast base-pairs recovered (“Chlp”) as the response variable with “Treatment” (16rxn, 32rxn, half_bait) as a fixed, categorical factor, and “Sample” and “Pool” (batch) as random effects:

exp5.mod16<-lmer(Chlp ~ Treatment * Age + (1|Tribe) + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Pool), data=exp5, REML=TRUE)

#Test removing Tribe from the model:

exp5.mod17<-lmer(Chlp ~ Treatment  * Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Pool), data=exp5, REML=TRUE)
anova(exp5.mod16, exp5.mod17)

#Tribe is non-significant and can be removed from the model. Test removing Genus nested in Tribe:

exp5.mod18<-lmer(Chlp ~ Treatment * Age +  (1|Tribe:Genus:Sample) + (1|Pool), data=exp5, REML=TRUE)
anova(exp5.mod17, exp5.mod18)

#Genus nested in Tribe is significant and cannot be removed from the model. Test removing Sample nested in Genus nested in Tribe:

exp5.mod19<-lmer(Chlp ~ Treatment  * Age + (1|Tribe:Genus) + (1|Pool), data=exp5, REML=TRUE)
anova(exp5.mod17, exp5.mod19)

#Sample nested in Genus nested in Tribe is significant and cannot be removed from the model. Test removing Pool:

exp5.mod20<-lmer(Chlp ~ Treatment  * Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample), data=exp5, REML=TRUE)
anova(exp5.mod17, exp5.mod20)

#Pool is significant and cannot be removed from the model.
#Test removing the interaction between treatment and age:

exp5.mod21<-lmer(Chlp ~ Treatment  + Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Pool), data=exp5, REML=TRUE)
anova(exp5.mod17, exp5.mod21)

#The interaction between treatment and age is not significant and can be removed from the model.
#Remove Age:

exp5.mod22<-lmer(Chlp ~ Treatment + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Pool), data=exp5, REML=TRUE)
anova(exp5.mod21, exp5.mod22)

#Age is significant and cannot be removed from the model.
#Remove Treatment:

exp5.mod23<-lmer(Chlp ~ Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Pool), data=exp5, REML=TRUE)
anova(exp5.mod21, exp5.mod23)

#Treatment is not significant and can be removed from the model. 

summary(exp5.mod21)
AIC(exp5.mod21)

#Plot model
emm<- emmeans(exp5.mod21, ~Age | Treatment, at = list(Age = seq(min(exp5$Age), max(exp5$Age), length=100)), lmer.df="satterthwaite")
emm
emm_df<-as.data.frame(emm)
emm_df

ggplot() +
  geom_point(data = exp5,
             aes(x = Age, y = Chlp, colour = Treatment),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Age, ymin = lower.CL, ymax = upper.CL, fill = Treatment),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Age, y = emmean, colour = Treatment),
            linewidth = 1) +
  theme_classic()

#Model residuals are approximately normally distributed according to the histogram and Q-Q plot, but the Shapiro-Wilk test indicates that statistically they are non-normal:

hist(residuals(exp5.mod21))
qqnorm(residuals(exp5.mod21))
shapiro.test(residuals(exp5.mod21))

#------------------------------------

#Chloroplast recovery
#glmm

#We performed a generalised linear mixed model with the proportion of chloroplast base-pairs recovered (“Chlp_prop”) as the response variable with “Treatment” (16rxn, 32rxn, half_bait) as a fixed, categorical factor, and “Sample” and “Pool” (batch) as random effects:

exp5_nz_n1 <- exp5[exp5$Chlp > 0 & exp5$Chlp <= 1, ]
exp5_nz_n1$Chlp_prop <- exp5_nz_n1$Chlp / 100
exp5.mod24<-glmmTMB(Chlp_prop ~ Treatment * Age + (1|Tribe) + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Pool), family=beta_family(link="logit"), data=exp5_nz_n1)

#Note this model is rank-deficient.

#Test removing Tribe as a factor:

exp5.mod25<-glmmTMB(Chlp_prop ~ Treatment* Age + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Pool), family=beta_family(link="logit"), data=exp5_nz_n1)
anova(exp5.mod24, exp5.mod25)

#Tribe is non-significant and can be removed as a factor. Test removing Genus within Tribe:

exp5.mod26<-glmmTMB(Chlp_prop ~ Treatment* Age  + (1|Tribe:Genus:Sample) + (1|Pool), family=beta_family(link="logit"), data=exp5_nz_n1)
anova(exp5.mod25, exp5.mod26)

#Genus within Tribe is non-significant and can be removed as a factor. Test removing Sample within Genus within Tribe:

exp5.mod27<-glmmTMB(Chlp_prop ~ Treatment* Age  + (1|Pool), family=beta_family(link="logit"), data=exp5_nz_n1)
anova(exp5.mod26, exp5.mod27)

#Sample within Genus within Tribe is non-significant and can be removed as a factor. Test removing Pool:

exp5.mod28<-glmmTMB(Chlp_prop ~ Treatment* Age, family=beta_family(link="logit"), data=exp5_nz_n1)
anova(exp5.mod27, exp5.mod28)

#Pool is non-significant and cannot be removed as a factor.
#Remove the interaction between treatment and age:

exp5.mod29<-glmmTMB(Chlp_prop ~ Treatment+ Age, family=beta_family(link="logit"), data=exp5_nz_n1)
anova(exp5.mod28, exp5.mod29)

#The interaction between treatment and age is non-signficant and can be removed. 
#Remove Age:

exp5.mod30<-glmmTMB(Chlp_prop ~ Treatment, family=beta_family(link="logit"), data=exp5_nz_n1)
anova(exp5.mod29, exp5.mod30)

#Age is significant and cannot be removed from the model. 
#Remove Treatment:

exp5.mod31<-glmmTMB(Chlp_prop ~  Age, family=beta_family(link="logit"), data=exp5_nz_n1)
anova(exp5.mod29, exp5.mod31)

#Treatment is nonsignificant and can be removed from the model. 

summary(exp5.mod29)
AIC(exp5.mod29)
emm <- emmeans(exp5.mod29, ~Treatment, type="response")
emm
pairs(emm, adjust="tukey")

#Plot model
emm<- emmeans(exp5.mod29, ~Age | Treatment, at = list(Age = seq(min(exp5$Age), max(exp5$Age), length=100)), type="response")
emm
emm_df<-as.data.frame(emm)
emm_df

ggplot() +
  geom_point(data = exp5_nz_n1,
             aes(x = Age, y = Chlp_prop, colour = Treatment),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Age, ymin = asymp.LCL, ymax = asymp.UCL, fill = Treatment),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Age, y = response, colour = Treatment),
            linewidth = 1) +
  theme_classic()

#Model residuals remain approximately normally distributed according to the histogram and Q-Q plot, but the Shapiro-Wilk test indicates that statistically they are still non-normal:

hist(residuals(exp5.mod29))
qqnorm(residuals(exp5.mod29))
shapiro.test(residuals(exp5.mod29))

#------------------------------------
#Target recovery vs chloroplast recovery
#------------------------------------

#lmer

exp2345<-read.csv("exp2_3_4_5_w_year.csv", header=TRUE)
head(exp2345)

#We performed a linear mixed model with the percentage chloroplast base-pairs recovered (“Chlp”) as the response variable with the percentage of target base-pairs recovered (“Ast353”) as a continuous predictor, “Experiment” (2, 3_4, and 5) as a fixed, categorical factor, and “Sample” and “Pool” (batch) as random effects:

exp2345.mod<-lmer(Chlp ~ Ast353 * Age * Experiment + Treatment + (1|Tribe) + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp2345, REML=TRUE)

#Test removing Tribe as a factor:

exp2345.mod2<-lmer(Chlp ~ Ast353 * Age* Experiment + Treatment + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp2345, REML=TRUE)
anova(exp2345.mod, exp2345.mod2)

#Tribe is non-significant and can be removed from the model. Test removing Genus within Tribe:

exp2345.mod3<-lmer(Chlp ~ Ast353 * Age* Experiment + Treatment + (1|Tribe:Genus:Sample) + (1|Batch), data=exp2345, REML=TRUE)
anova(exp2345.mod2, exp2345.mod3)

#Genus within Tribe is significant and cannot be removed from the model. Test removing Sample within Genus within Tribe:

exp2345.mod4<-lmer(Chlp ~ Ast353 * Age* Experiment + Treatment + (1|Tribe:Genus)  + (1|Batch), data=exp2345, REML=TRUE)
anova(exp2345.mod2, exp2345.mod4)

#Sample within Genus within Tribe is significant and cannot be removed from the model. Test removing Batch:

exp2345.mod5<-lmer(Chlp ~ Ast353 * Age* Experiment + Treatment + (1|Tribe:Genus) + (1|Tribe:Genus:Sample), data=exp2345, REML=TRUE)
anova(exp2345.mod2, exp2345.mod5)

#Batch is significant and cannot be removed from the model. Test the interaction between target recovery/Age and experiment:

exp2345.mod6<-lmer(Chlp ~ Ast353 * Age + Experiment + Treatment + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp2345, REML=TRUE)
anova(exp2345.mod2, exp2345.mod6)

#The interaction between target recovery/Age and experiment is significant and cannot be removed from the model. Test removing Treatment:

exp2345.mod7<-lmer(Chlp ~ Ast353 * Age* Experiment + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp2345, REML=TRUE)
anova(exp2345.mod2, exp2345.mod7)

#Treatment is non-significant and can be removed from the model. 
#Test removing the interaction between target recovery/experiment and age:

exp2345.mod8<-lmer(Chlp ~ Ast353 + Age* Experiment + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), data=exp2345, REML=TRUE)
anova(exp2345.mod7, exp2345.mod8)

#The interaction between target recovery/experiment and age is significant and cannot be removed from the model.


summary(exp2345.mod7)
vc <- as.data.frame(VarCorr(exp2345.mod7))
vc$percent <- vc$vcov / sum(vc$vcov) * 100
vc
emm<-emmeans(exp2345.mod7, ~Ast353 | Experiment | Age, at = list (Ast353 = seq (min(exp2345$Ast353), max(exp2345$Ast353), length=100)), lmer.df="satterthwaite")
emm_df<-as.data.frame(emm)

ggplot() +
  geom_point(data = exp2345,
             aes(x = Ast353, y = Chlp, colour = Experiment),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Ast353, ymin = lower.CL, ymax = upper.CL, fill = Experiment),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Ast353, y = emmean, colour = Experiment),
            linewidth = 1) +
  theme_classic()

AIC(exp2345.mod7)

#Model residuals are approximately normally distributed according to the histogram and Q-Q plot, but the Shapiro-Wilk test indicates that statistically they are non-normal:

hist(residuals(exp2345.mod7))
qqnorm(residuals(exp2345.mod7))
shapiro.test(residuals(exp2345.mod7))

#------------------------------------

#glmm

#We performed a linear mixed model with the proportion of chloroplast base-pairs (“Chlp_prop”) as the response variable with the proportion of target base-pairs recovered (“Ast353_prop”) as a continuous predictor, “Experiment” (2, 3_4, and 5) as a fixed, categorical factor, and “Sample” and “Pool” (batch) as random effects:

exp2345$Chlp_prop <- exp2345$Chlp / 100
exp2345$Ast353_prop <- exp2345$Ast353 / 100
exp2345_nz_n1 <- exp2345[exp2345$Chlp_prop > 0 & exp2345$Chlp_prop <= 1, ]
exp2345_n0_n1 <- exp2345_nz_n1[exp2345_nz_n1$Ast353_prop > 0 & exp2345_nz_n1$Ast353_prop <= 1, ]
head(exp2345_n0_n1)

exp2345.mod9<-glmmTMB(Chlp_prop ~ Ast353_prop *Age * Experiment + Treatment + (1|Tribe) + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp2345_n0_n1)

#Test removing batch from the model:

exp2345.mod10<-glmmTMB(Chlp_prop ~ Ast353_prop *Age * Experiment + Treatment + (1|Tribe) + (1|Tribe:Genus) + (1|Tribe:Genus:Sample), family=beta_family(link="logit"), data=exp2345_n0_n1)
anova(exp2345.mod9,exp2345.mod10)

#Batch is significant and cannot be removed from the model.
#Test removing Sample:

exp2345.mod11<-glmmTMB(Chlp_prop ~ Ast353_prop *Age * Experiment + Treatment + (1|Tribe) + (1|Tribe:Genus) + (1|Batch), family=beta_family(link="logit"), data=exp2345_n0_n1)
anova(exp2345.mod9,exp2345.mod11)

#Sample is significant and cannot be removed from the model. 
#Test removing Genus:

exp2345.mod12<-glmmTMB(Chlp_prop ~ Ast353_prop *Age * Experiment + Treatment + (1|Tribe)  + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp2345_n0_n1)
anova(exp2345.mod9,exp2345.mod12)

#Genus is significant and cannot be removed from the model. 
#Test removing Tribe:

exp2345.mod13<-glmmTMB(Chlp_prop ~ Ast353_prop *Age * Experiment + Treatment  + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp2345_n0_n1)
anova(exp2345.mod9,exp2345.mod13)

#Tribe is insignificant and can be removed from the model. 
#Test removing Treatment:

exp2345.mod14<-glmmTMB(Chlp_prop ~ Ast353_prop *Age * Experiment + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp2345_n0_n1)
anova(exp2345.mod13,exp2345.mod14)

#Treatment is not significant and can be removed from the model. 
#Test the interaction between target recovery/age and experiment:

exp2345.mod15<-glmmTMB(Chlp_prop ~ Ast353_prop *Age + Experiment + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp2345_n0_n1)
anova(exp2345.mod14,exp2345.mod15)

#The interaction between target recovery/age and experiment is significant and cannot be removed from the model. 
#Test the interaction between target recovery/experiment and age:

exp2345.mod16<-glmmTMB(Chlp_prop ~ Ast353_prop + Age * Experiment + (1|Tribe:Genus) + (1|Tribe:Genus:Sample) + (1|Batch), family=beta_family(link="logit"), data=exp2345_n0_n1)
anova(exp2345.mod14,exp2345.mod16)

#The interaction between target recovery/experiment and age is significant and cannot be removed from the model.

summary(exp2345.mod14)
AIC(exp2345.mod14)
emm <- emmeans(exp2345.mod14, ~Ast353_prop | Experiment | Age, type="response")
emm
emm<-emmeans(exp2345.mod14, ~Ast353_prop | Experiment | Age, type="response", at = list (Ast353_prop = seq (min(exp2345_n0_n1$Ast353_prop), max(exp2345_n0_n1$Ast353_prop), length=100)))
emm_df<-as.data.frame(emm)

ggplot() +
  geom_point(data = exp2345_n0_n1,
             aes(x = Ast353_prop, y = Chlp_prop, colour = Experiment),
             alpha = 0.4) +
  geom_ribbon(data = emm_df,
              aes(x = Ast353_prop, ymin = asymp.LCL, ymax = asymp.UCL, fill = Experiment),
              alpha = 0.2) +
  geom_line(data = emm_df,
            aes(x = Ast353_prop, y = response, colour = Experiment),
            linewidth = 1) +
  theme_classic()

#Model residuals remain approximately normally distributed according to the histogram and Q-Q plot, but the Shapiro-Wilk test indicates that statistically they are still non-normal:

hist(residuals(exp2345.mod14))
qqnorm(residuals(exp2345.mod14))
shapiro.test(residuals(exp2345.mod14))

