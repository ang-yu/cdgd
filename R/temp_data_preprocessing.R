
rm(list=ls(all=TRUE))

library("caret")
library("nnet")
library("gbm")
library("xtable")

options(scipen=999)

data <- readRDS("/Users/Ang/Desktop/Research/Counterfactual covariances/Data/data_cleaned_ge.rds")
colSums(is.na(data))/nrow(data)

data$gender <- as.numeric(data$gender)-1

data$parental_presence <- as.numeric(data$parental_presence)-1

data$n_sib <- as.numeric(data$n_sib)-1

data$urban <- as.numeric(data$urban)-1

data$edu_exp <- as.numeric(data$edu_exp)-2
data$edu_exp[data$edu_exp==-1] <- 0

data$age <- as.numeric(data$age)-1

data$friend_edu_exp <- as.numeric(data$friend_edu_exp)-2
data$friend_edu_exp[data$friend_edu_exp==-1] <- 0

data$sig_other_exp1 <- NA
data$sig_other_exp1[data$sig_other_expec==1 & !is.na(data$sig_other_expec)] <- 1
data$sig_other_exp1[data$sig_other_expec!=1 & !is.na(data$sig_other_expec)] <- 0
data$sig_other_exp2 <- NA
data$sig_other_exp2[data$sig_other_expec==2 & !is.na(data$sig_other_expec)] <- 1
data$sig_other_exp2[data$sig_other_expec!=2 & !is.na(data$sig_other_expec)] <- 0
data$sig_other_exp3 <- NA
data$sig_other_exp3[data$sig_other_expec==3 & !is.na(data$sig_other_expec)] <- 1
data$sig_other_exp3[data$sig_other_expec!=3 & !is.na(data$sig_other_expec)] <- 0
data$sig_other_exp4 <- NA
data$sig_other_exp4[data$sig_other_expec==4 & !is.na(data$sig_other_expec)] <- 1
data$sig_other_exp4[data$sig_other_expec!=4 & !is.na(data$sig_other_expec)] <- 0

data$foreign_lang <- as.numeric(data$foreign_lang)-1

data$SMSA1 <- NA
data$SMSA1[data$SMSA==0 & !is.na(data$SMSA)] <- 1
data$SMSA1[data$SMSA!=0 & !is.na(data$SMSA)] <- 0
data$SMSA2 <- NA
data$SMSA2[data$SMSA==1 & !is.na(data$SMSA)] <- 1
data$SMSA2[data$SMSA!=1 & !is.na(data$SMSA)] <- 0
data$SMSA3 <- NA
data$SMSA3[data$SMSA==2 & !is.na(data$SMSA)] <- 1
data$SMSA3[data$SMSA!=2 & !is.na(data$SMSA)] <- 0
data$SMSA4 <- NA
data$SMSA4[data$SMSA==3 & !is.na(data$SMSA)] <- 1
data$SMSA4[data$SMSA!=3 & !is.na(data$SMSA)] <- 0

data$mother_seperate <- as.numeric(data$mother_seperate)-1

data$school_satis1 <- NA
data$school_satis1[(data$school_satisfaction==1 | data$school_satisfaction==2) & !is.na(data$school_satisfaction)] <- 1
data$school_satis1[data$school_satisfaction!=1 & data$school_satisfaction!=2 & !is.na(data$school_satisfaction)] <- 0
data$school_satis2 <- NA
data$school_satis2[data$school_satisfaction==3 & !is.na(data$school_satisfaction)] <- 1
data$school_satis2[data$school_satisfaction!=3 & !is.na(data$school_satisfaction)] <- 0
data$school_satis3 <- NA
data$school_satis3[data$school_satisfaction==4 & !is.na(data$school_satisfaction)] <- 1
data$school_satis3[data$school_satisfaction!=4 & !is.na(data$school_satisfaction)] <- 0

data$fm_foreign_born <- NA
data$fm_foreign_born[(data$f_foreign_born==2 | data$m_foreign_born==2) & !is.na(data$f_foreign_born) & !is.na(data$m_foreign_born)] <- 1
data$fm_foreign_born[(data$f_foreign_born==1 | data$m_foreign_born==1) & !is.na(data$f_foreign_born) & !is.na(data$m_foreign_born)] <- 0

data$region1 <- NA
data$region1[data$region==1 & !is.na(data$region)] <- 1
data$region1[data$region!=1 & !is.na(data$region)] <- 0
data$region2 <- NA
data$region2[data$region==2 & !is.na(data$region)] <- 1
data$region2[data$region!=2 & !is.na(data$region)] <- 0
data$region3 <- NA
data$region3[data$region==3 & !is.na(data$region)] <- 1
data$region3[data$region!=3 & !is.na(data$region)] <- 0
data$region4 <- NA
data$region4[data$region==4 & !is.na(data$region)] <- 1
data$region4[data$region!=4 & !is.na(data$region)] <- 0

data$m_work <- as.numeric(data$m_work)-1

data$race1 <- NA
data$race1[data$race=="Other" & !is.na(data$race)] <- 1
data$race1[data$race!="Other" & !is.na(data$race)] <- 0
data$race2 <- NA
data$race2[data$race=="Black" & !is.na(data$race)] <- 1
data$race2[data$race!="Black" & !is.na(data$race)] <- 0
data$race3 <- NA
data$race3[data$race=="Hispanic" & !is.na(data$race)] <- 1
data$race3[data$race!="Hispanic" & !is.na(data$race)] <- 0

data$completion <- as.numeric(data$completion)-1

data <- na.omit(data)

data$pincome_1 <- NA
data$pincome_1[data$parental_income_rank>=0.75 & !is.na(data$parental_income_rank)] <- 1
data$pincome_1[data$parental_income_rank<0.75 & !is.na(data$parental_income_rank)] <- 0

data$pincome_2 <- NA
data$pincome_2[data$parental_income_rank<=0.25 & !is.na(data$parental_income_rank)] <- 1
data$pincome_2[data$parental_income_rank>0.25 & !is.na(data$parental_income_rank)] <- 0

data <- data[data$pincome_1==1 | data$pincome_2==1, ]

#table(data$completion[data$parental_income_log<8])   # very low parental income is associated with very low probability of college completion
# note that at very high parental income, the distribution of treatment is much more balanced
#data <- data[data$parental_income_log>=8,]  # sample size down by 245


### Unconditional decomposition
# The Q variable should be listed first to ensure consistency with conditional decomposition
result_parametric <- cdgd0_parametric(Y="adult_income_rank",D="completion",G="pincome_1",
                    X=c("AFQT","gender","medu","parental_presence",
                        "n_sib","urban","edu_exp","age","friend_edu_exp","rotter_score","rosenberg_irt_score",
                        "sig_other_exp1","sig_other_exp2","sig_other_exp3","sig_other_exp4","foreign_lang",
                        "SMSA1","SMSA2","SMSA3","SMSA4","mother_seperate",
                        "school_satis1","school_satis2","school_satis3","fm_foreign_born",
                        "region1","region2","region3","region4","m_work","race1","race2","race3"),
                    data=data)

set.seed(1)
result_gbm <- cdgd0(Y="adult_income_rank",D="completion",G="pincome_1",
      X=c("AFQT","gender","medu","parental_presence",
          "n_sib","urban","edu_exp","age","friend_edu_exp","rotter_score","rosenberg_irt_score",
          "sig_other_exp1","sig_other_exp2","sig_other_exp3","sig_other_exp4","foreign_lang",
          "SMSA1","SMSA2","SMSA3","SMSA4","mother_seperate",
          "school_satis1","school_satis2","school_satis3","fm_foreign_born",
          "region1","region2","region3","region4","m_work","race1","race2","race3"),
      data=data,algorithm="gbm")

set.seed(1)
result_nnet <- cdgd0(Y="adult_income_rank",D="completion",G="pincome_1",
                    X=c("AFQT","gender","medu","parental_presence",
                        "n_sib","urban","edu_exp","age","friend_edu_exp","rotter_score","rosenberg_irt_score",
                        "sig_other_exp1","sig_other_exp2","sig_other_exp3","sig_other_exp4","foreign_lang",
                        "SMSA1","SMSA2","SMSA3","SMSA4","mother_seperate",
                        "school_satis1","school_satis2","school_satis3","fm_foreign_born",
                        "region1","region2","region3","region4","m_work","race1","race2","race3"),
                    data=data,algorithm="nnet")

set.seed(1)
result_ranger <- cdgd0(Y="adult_income_rank",D="completion",G="pincome_1",
                     X=c("AFQT","gender","medu","parental_presence",
                         "n_sib","urban","edu_exp","age","friend_edu_exp","rotter_score","rosenberg_irt_score",
                         "sig_other_exp1","sig_other_exp2","sig_other_exp3","sig_other_exp4","foreign_lang",
                         "SMSA1","SMSA2","SMSA3","SMSA4","mother_seperate",
                         "school_satis1","school_satis2","school_satis3","fm_foreign_born",
                         "region1","region2","region3","region4","m_work","race1","race2","race3"),
                     data=data,algorithm="ranger")


### Conditional decomposition
cond_result_parametric <- cdgd1_parametric(Y="adult_income_rank",D="completion",G="pincome_1",Q="AFQT",
                         X=c("gender","medu","parental_presence",
                             "n_sib","urban","edu_exp","age","friend_edu_exp","rotter_score","rosenberg_irt_score",
                             "sig_other_exp1","sig_other_exp2","sig_other_exp3","sig_other_exp4","foreign_lang",
                             "SMSA1","SMSA2","SMSA3","SMSA4","mother_seperate",
                             "school_satis1","school_satis2","school_satis3","fm_foreign_born",
                             "region1","region2","region3","region4","m_work","race1","race2","race3"),
                         data=data)

set.seed(1)
cond_result_gbm <- cdgd1(Y="adult_income_rank",D="completion",G="pincome_1",Q="AFQT",
                    X=c("gender","medu","parental_presence",
                        "n_sib","urban","edu_exp","age","friend_edu_exp","rotter_score","rosenberg_irt_score",
                        "sig_other_exp1","sig_other_exp2","sig_other_exp3","sig_other_exp4","foreign_lang",
                        "SMSA1","SMSA2","SMSA3","SMSA4","mother_seperate",
                        "school_satis1","school_satis2","school_satis3","fm_foreign_born",
                        "region1","region2","region3","region4","m_work","race1","race2","race3"),
                    data=data,algorithm="gbm")

set.seed(1)
cond_result_nnet <- cdgd1(Y="adult_income_rank",D="completion",G="pincome_1",Q="AFQT",
                     X=c("gender","medu","parental_presence",
                         "n_sib","urban","edu_exp","age","friend_edu_exp","rotter_score","rosenberg_irt_score",
                         "sig_other_exp1","sig_other_exp2","sig_other_exp3","sig_other_exp4","foreign_lang",
                         "SMSA1","SMSA2","SMSA3","SMSA4","mother_seperate",
                         "school_satis1","school_satis2","school_satis3","fm_foreign_born",
                         "region1","region2","region3","region4","m_work","race1","race2","race3"),
                     data=data,algorithm="nnet")

set.seed(1)
cond_result_ranger <- cdgd1(Y="adult_income_rank",D="completion",G="pincome_1",Q="AFQT",
                       X=c("gender","medu","parental_presence",
                           "n_sib","urban","edu_exp","age","friend_edu_exp","rotter_score","rosenberg_irt_score",
                           "sig_other_exp1","sig_other_exp2","sig_other_exp3","sig_other_exp4","foreign_lang",
                           "SMSA1","SMSA2","SMSA3","SMSA4","mother_seperate",
                           "school_satis1","school_satis2","school_satis3","fm_foreign_born",
                           "region1","region2","region3","region4","m_work","race1","race2","race3"),
                       data=data,algorithm="ranger")







