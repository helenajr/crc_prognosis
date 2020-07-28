#Load data
Texpression <- read.delim("~/Desktop/Xena/2020-21-07_coadread.tsv")
View(Texpression)
#Look at sample type
count(Texpression, sample_type)
#new data frame filtered for only prim tumour
prim_tumor <- filter(Texpression, sample_type == "Primary Tumor")

#Load required packages
library(survival)
library(survminer)
library(dplyr)
library(forcats)
glimpse(Texpression)
glimpse(gtex)

#Summarise - number of samples in TCGA coadread
length(prim_tumor$T)
length(na.omit(prim_tumor$T))
summary(prim_tumor$T)

#Specify factor levels and convert other values to NA
Texpression$gender <- factor(Texpression$gender, c("FEMALE", "MALE"))
str(Texpression$gender)
Texpression <- Texpression %>% 
  mutate(stage_group = fct_collapse(Texpression$pathologic_stage, 
                                    Stage1 = c("Stage I", "Stage IA"),
                                    Stage2 = c("Stage II", "Stage IIA",
                                               "Stage IIB", "Stage IIC"),
                                    Stage3 = c("Stage III", "Stage IIIA",
                                               "Stage IIIB", "Stage IIIC"),
                                    Stage4 = c("Stage IV", "Stage IVA", "Stage IVB")))

Texpression %>% select(pathologic_stage, stage_group)


#Dichotamise T expression
prim_tumor <- prim_tumor %>% 
  mutate(T_group = ifelse(T > 2.525, "> normal range","normal range"))
prim_tumor$T_group <- factor(prim_tumor$T_group)

#Fit survival data using Kap-Meier method
surv_object <- Surv(time = prim_tumor$OS.time, event = prim_tumor$OS)
surv_object

#Fit Kaplan Meier curves
fit1 <- survfit(surv_object ~ T_group, data = prim_tumor)
summary(fit1)

#Display curve
ggsurvplot(fit1, data = prim_tumor, pval = TRUE, xlab = "Time(days)")

#Find n in groups and number of deleted obs
fit1

#CoxPH model T group univariable
Texpression$T_group <- relevel(Texpression$T_group, ref = "normal range")
fit.coxph <- coxph(surv_object ~ T_group, data = Texpression)
ggforest(fit.coxph, data = Texpression)

#Dichotamise other variables
Texpression <- Texpression %>% 
  mutate(age_group = ifelse(age_at_initial_pathologic_diagnosis > 68, "old","young"))
Texpression$age_group <- factor(Texpression$age_group)
Texpression <- Texpression %>% 
  mutate(di_stage = ifelse(stage_group == "Stage4", "late", "early"))
Texpression$di_stage <- factor(Texpression$di_stage)

#CoxPH model multivariable
fit.coxph <- coxph(surv_object ~ T_group + gender + age_group + di_stage, data = Texpression)
ggforest(fit.coxph, data = Texpression)

#% Texpression >0
(sum(na.omit(Texpression$T>0)))/(length(na.omit(Texpression$T)))*100
(sum(na.omit(gtex_colon$T>0)))/(length(na.omit(gtex_colon$T)))*100

#chi squared analysis age
highT_old <- sum(na.omit(Texpression$T_group == "> normal range" 
            &(Texpression$age_group =="old")))
highT_old

highT_young <- sum(na.omit(Texpression$T_group == "> normal range" 
                           &(Texpression$age_group =="young")))
highT_young

lowT_old <- sum(na.omit(Texpression$T_group == "normal range" 
                         &(Texpression$age_group =="old")))
lowT_old

lowT_young <- sum(na.omit(Texpression$T_group == "normal range" 
                        &(Texpression$age_group =="young")))
lowT_young

age_obs_table <- matrix(c(233, 184, 7, 6), nrow = 2, 
                        ncol =2, byrow = T)
rownames(age_obs_table) <- c('Low T', 'High T')
colnames(age_obs_table) <- c('Young', 'Old')
age_obs_table

chisq.test(age_obs_table)

#chi squared gender
highT_male <- sum(na.omit(Texpression$T_group == "> normal range" 
                         &(Texpression$gender =="MALE")))
highT_male

highT_female <- sum(na.omit(Texpression$T_group == "> normal range" 
                           &(Texpression$gender =="FEMALE")))
highT_female

lowT_male <- sum(na.omit(Texpression$T_group == "normal range" 
                        &(Texpression$gender =="MALE")))
lowT_male

lowT_female <- sum(na.omit(Texpression$T_group == "normal range" 
                          &(Texpression$gender =="FEMALE")))
lowT_female

gender_obs_table <- matrix(c(221, 196, 10, 3), nrow = 2, 
                        ncol =2, byrow = T)
rownames(gender_obs_table) <- c('Low T', 'High T')
colnames(gender_obs_table) <- c('Male', 'Female')
gender_obs_table

chisq.test(gender_obs_table)

#Chi squared stage
highT_early <- sum(na.omit(Texpression$T_group == "> normal range" 
                          &(Texpression$di_stage =="early")))
highT_early

highT_late <- sum(na.omit(Texpression$T_group == "> normal range" 
                            &(Texpression$di_stage =="late")))
highT_late

lowT_early <- sum(na.omit(Texpression$T_group == "normal range" 
                         &(Texpression$di_stage =="early")))
lowT_early

lowT_late <- sum(na.omit(Texpression$T_group == "normal range" 
                           &(Texpression$di_stage =="late")))
lowT_late

stage_obs_table <- matrix(c(364, 57, 9, 4), nrow = 2, 
                           ncol =2, byrow = T)
rownames(stage_obs_table) <- c('Low T', 'High T')
colnames(stage_obs_table) <- c('Early', 'Late')
stage_obs_table

chisq.test(stage_obs_table)
#Warning message appears when expected counts in a category are fewer
#than 5. Still ok to use - see help.


