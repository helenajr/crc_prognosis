#Load data
samples <- 
  read.delim("~/Desktop/crc_prognosis/tcga_target_gtex_colorectal.tsv", 
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
count(cancer, T >4.411)
boxplot(normal$T, cancer$T,
        at = c(1,2),
        names = c("Normal", "Cancer"),
        xlab = "Tissue type",
        ylab = "T expression (log2(norm_count+1))")

####Kaplan Meier plots####
#Dichotamise T expression based on normal range
cancer <- cancer %>% 
  mutate(TBXT_group = ifelse(T > 4.4110, " > Normal range"," Normal range"))
cancer$TBXT_group <- factor(cancer$TBXT_group)
str(cancer)

#Fit survival data using Kap-Meier method
surv_object <- Surv(time = cancer$OS.time, event = cancer$OS)
surv_object

#Fit Kaplan Meier curves
fit1 <- survfit(surv_object ~ TBXT_group, data = cancer)
summary(fit1)

#Display curve
ggsurvplot(fit1, data = cancer, pval = TRUE, xlab = "Time(days)")

#Find n in groups and number of deleted obs
fit1

####Dichotamise T expression based on arbitary cut-off####
cancer <- cancer %>% 
  mutate(TBXT_group = ifelse(T > 5.7, " High"," Low"))
cancer$TBXT_group <- factor(cancer$TBXT_group)
str(cancer)

#Fit survival data using Kap-Meier method
surv_object <- Surv(time = cancer$OS.time, event = cancer$OS)
surv_object

#Fit Kaplan Meier curves
fit1 <- survfit(surv_object ~ TBXT_group, data = cancer)
summary(fit1)

#Display curve
ggsurvplot(fit1, data = cancer, pval = TRUE, xlab = "Time(days)")

#Find n in groups and number of deleted obs
fit1

