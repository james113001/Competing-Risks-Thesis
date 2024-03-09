library(cmprsk)
library(survival)
library(skimr)
library(dplyr)
library(Survout)
library(survminer)
library(rpivotTable)
library(gtsummary)
library(knitr)
library(ezfun)
library(randomForestSRC)
library(boot)

# data(cancer, package="survival")
# data("peakVO2", package = "randomForestSRC")
# colnames(veteran)
# unique(pbc$status)
# skim(ovarian)
# 
# dim(survival::colon)
# dim(wihs)
# 

dim(survival::udca)

pbc2 <- na.omit(pbc)
skim(udca)
table(udca$stage)
skim(melanoma$status)
table(pbc$status)
#
# 
# clinicalpbc<- pbc %>% filter(!is.na(trt)) 
# 
# clinicalpbc$status<- factor(clinicalpbc$status, levels=c(0,1,2), labels=c("censored", "transplant", "death"))
# clinicalpbc$edema<- factor(clinicalpbc$edema, levels=c(0,0.5,1), labels=c("no edema", "untreated or successfully treated", "edema despite therapy"))
# clinicalpbc$trt<- factor(clinicalpbc$trt, levels=c(1,2), labels=c("D-penicillmain", "placebo"))
# skim(clinicalpbc$trt)
# 
# attach(clinicalpbc)
# table(status, trt)
# 
# tapply(time, list(trt, status), mean) #time given in days
# round(tapply(time, list(trt, status), mean), digits=2)
# 
# rpivotTable(data=clinicalpbc, rows="trt", cols="status", vals="time",
#             aggregatorName = "Average",rendererName = "Table")
# 
# 
# 
# crrpbc <- pbc %>% filter(!is.na(treatment))
# covariates<- crrpbc %>% 
#   select(-c(id, time, status)) %>%
#   mutate(sex = factor2ind(sex,"f"))
# 
# skim(covariates)
# attach(crrpbc)
# mod1<-crr(time, status, covariates[, c("trt","age","sex","stage")], cencode = 2)
# summary(mod1)
# mvcrrres(mod1) %>% 
#   kable()
# 
# fit <- 
#   coxph(
#     Surv(time, ifelse(status == 2, 1, 0)) ~ trt + age + sex + stage, 
#     data = covariates
#   )
# 
# tbl_regression(fit, exp = TRUE)
# 
# ggsurvplot(survfit(fit, data = crrpbc), color = "#2E9FDF",
#            ggtheme = theme_minimal())
# 
# skim(crrpbc)


## Random Forest
data(cancer, package = "survival")
melanoma
skim(melanoma)

## Analysis 1
## modified Gray's weighted log-rank splitting
## (equivalent to cause=c(1,1) and splitrule="logrankCR")
o1 <- rfsrc(Surv(time, status) ~ ., melanoma, ntree = 1000)

## Analysis 2
## log-rank cause-1 (death) specific splitting and targeted VIMP
o2 <- rfsrc(Surv(time, status) ~ ., melanoma, 
            splitrule = "logrank", cause = c(1,0), importance = TRUE)

## Analysis 3
## log-rank cause-2 (transplant) specific splitting and targeted VIMP
o3 <- rfsrc(Surv(time, status) ~ ., melanoma, 
            splitrule = "logrank", cause = c(0,1), importance = TRUE)

## extract VIMP from the log-rank forests: event-specific
## extract minimal depth from the Gray log-rank forest: non-event specific
vimpOut <- data.frame(md = max.subtree(o1)$order[, 1],
                      vimp.death = 100 * o2$importance[ ,1],
                      vimp.transplant = 100 * o3$importance[ ,2])

print(vimpOut[order(vimpOut$md), ], digits = 2)


pdf("melanoma.pdf", width = 8, height = 8)
par(cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,6,1,1), mgp = c(4, 1, 0))
plot.competing.risk(o1)
dev.off()

o1$predicted
plot(o1)

get.mv.error(o1)
get.auc(y, prob)

jk.o1 <- subsample(o1)
pdf("VIMPsur.pdf", width = 15, height = 20)
par(oma = c(0.5, 10, 0.5, 0.5))
par(cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,17,1,1), mgp = c(4, 1, 0))
plot(jk.o1, xlab = "Variable Importance (x 100)", cex = 1.2)
dev.off()
