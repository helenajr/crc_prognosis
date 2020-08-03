#Load data
samples <- 
  read.delim("~/Desktop/Xena/tcga_target_gtex_colorectal.tsv", 
             header=TRUE)
View(samples)
#Load required packages
library(survival)
library(survminer)
library(dplyr)
library(forcats)
glimpse(samples)

#Summarise - number of samples for each sample type
count(samples, X_sample_type)

#Boxplot to compare normal and cancer distributions
normal <- filter(samples, X_sample_type == "Normal Tissue")
summary(normal$T)
count(normal, T >0)
cancer <- filter(samples, X_sample_type == "Primary Tumor")
summary(cancer$T)
count(cancer, T >0)
count(cancer, T >4.411)
boxplot(normal$T, cancer$T,
        at = c(1,2),
        names = c("Normal", "Cancer"),
        xlab = "Tissue type",
        ylab = "T expression (log2(norm_count+1))")


#Dichotamise T expression
cancer <- cancer %>% 
  mutate(T_group = ifelse(T > 4.4110, "> normal range","normal range"))
cancer$T_group <- factor(cancer$T_group)
str(cancer)

#Fit survival data using Kap-Meier method
surv_object <- Surv(time = cancer$OS.time, event = cancer$OS)
surv_object

#Fit Kaplan Meier curves
fit1 <- survfit(surv_object ~ T_group, data = cancer)
summary(fit1)

#Display curve
ggsurvplot(fit1, data = cancer, pval = TRUE, xlab = "Time(days)")

#Find n in groups and number of deleted obs
fit1

#CoxPH model T group univariable
cancer$T_group <- relevel(cancer$T_group, ref = "normal range")
fit.coxph <- coxph(surv_object ~ T_group, data = cancer)
ggforest(fit.coxph, data = cancer)

####Expected count measure####
#Boxplot to compare normal and cancer distributions
summary(normal$ENSG00000164458.9)
count(normal, ENSG00000164458.9 >0)
summary(cancer$ENSG00000164458.9)
count(cancer, ENSG00000164458.9 >0)
count(cancer, ENSG00000164458.9 >4.411)
boxplot(normal$ENSG00000164458.9, cancer$ENSG00000164458.9,
        at = c(1,2),
        names = c("Normal", "Cancer"),
        xlab = "Tissue type",
        ylab = "T expression (log2(norm_count+1))")

#Dichotamise T expression
cancer <- cancer %>% 
  mutate(T_group_exp = ifelse(ENSG00000164458.9 > 3.917, "> normal range","normal range"))
cancer$T_group_exp <- factor(cancer$T_group_exp)
str(cancer)

#Fit survival data using Kap-Meier method
surv_object <- Surv(time = cancer$OS.time, event = cancer$OS)
surv_object

#Fit Kaplan Meier curves
fit1 <- survfit(surv_object ~ T_group, data = cancer)
summary(fit1)

#Display curve
ggsurvplot(fit1, data = cancer, pval = TRUE, xlab = "Time(days)")

#Find n in groups and number of deleted obs
fit1

#CoxPH model T group univariable
cancer$T_group <- relevel(cancer$T_group, ref = "normal range")
fit.coxph <- coxph(surv_object ~ T_group, data = cancer)
ggforest(fit.coxph, data = cancer)