library(boot)
library(skimr)
library(table1)
library(survival)
library(dplyr)
library(cmprsk)
library(randomForestSRC)
library(riskRegression)
library(prodlim)
library(pec)
library(party)
library(BART)
library(feather)
library(rsample) 




#Import dataset
data(wihs)
cmpdata <- wihs
skim(cmpdata)
recoded <- cmpdata
##0 = Alive, 1 = HAART Initiate, 2= death
table(recoded$status)
recoded$black <- as.factor(recoded$black)
recoded$idu <- as.factor(recoded$idu)

#table 1
cmpdata.labeled <- recoded
table(recoded$status)
cmpdata.labeled$status <- 
  factor(recoded$status, 
         levels=c(0,1,2),
         labels=c("Censored",
                  "HAART Initiate", # Reference
                  "Death"))
cmpdata.labeled$idu <- 
  factor(cmpdata$idu, 
         levels=c(0,1),
         labels=c("No history", # Reference
                  "History"))
cmpdata.labeled$black <- 
  factor(cmpdata$black, 
         levels=c(0,1),
         labels=c("Not black", # Reference
                  "Black"))
table1(~ time + idu + ageatfda + cd4nadir + black | status, data=cmpdata.labeled)

recoded <- recoded %>% mutate_at(c("ageatfda", "cd4nadir"), 
                                 ~(scales::rescale(.) %>% as.vector))
skim(recoded)
write_feather(recoded, "wihs.feather")


#Cox PH (Cause-specific Hazards model)
print("COX PROPORTIONAL HAZARDS")
cmpdatacen <- recoded %>%  #censoring competing risks (Death)
  mutate(status = replace(status, status != 1, 0))

#fit model
coxmodel<- coxph(Surv(time, status == 1) ~ idu + ageatfda + cd4nadir + black, 
                 data = cmpdatacen, x = TRUE)
cox.zph(coxmodel)

#C-index
c<-(cindex(coxmodel, formula = Surv(time, status) ~
               idu + ageatfda + cd4nadir + black,
             data=cmpdatacen, splitMethod="boot632", B=1000,verbose = FALSE, keep.index = TRUE,
           keep.matrix = TRUE, keep.pvalues = TRUE))
print(c)
cvector<-unlist(c$BootstrapCrossValCindexMat)
print(sd(cvector))
brier<-pec(coxmodel, data = cmpdatacen, splitMethod="boot632", B=1000,
           verbose = FALSE)
print(brier)
ls<-brier$Boot632Err[2]
print(sd(unlist(ls), na.rm = TRUE))
print(quantile(unlist(ls), c(.025, .975), na.rm = TRUE))
print("END OF COXPH")
cat("\n\n\n\n") #for creating new lines in printed output


saveRDS(cvector, file = "cwihs.rds")
saveRDS(cvector2, file = "cwihs2.rds")
saveRDS(cvector3, file = "cwihs3.rds")
save.image()
t<-t.test(cvector, cvector3, paired = TRUE, alternative = "two.sided")
diff<-(cvector-cvector3)
mean<-mean(diff)
variance<-var(diff)
tstat<-mean/sqrt(variance*(1/1000 + 999))
p<-pt(tstat, 999, lower.tail = FALSE)
p
p.adjust(p, method = "fdr", n = length(p))

#Fine & Gray model (Subdistribution model)
print("FINE AND GRAY MODEL")

#visualise CIF
CIF <- cuminc(ftime   = cmpdata.labeled$time,  # failure time variable (years)
              fstatus = cmpdata.labeled$status,  # variable with distinct codes for different causes of failure
              #group   = cmpdata.labeled$black,  # estimates will calculated within groups
              ## strata  = ,  # Tests will be stratified on this variable.
              rho     = 0, # Power of the weight function used in the tests.
              cencode = 2, # value of fstatus variable which indicates the failure time is censored.
              ## subset = ,
              ## na.action = na.omit
)
plot(CIF, color=1:6)
title("WIHS")
#c-index     must standardize td-c-index date across all methods
fgr <- FGR(Hist(time,status)~idu + ageatfda + cd4nadir + black,
           data=recoded, cause = 1) ##Fine-Gray
c2<-(cindex(fgr, data = recoded, splitMethod="boot632", B=1000,
             verbose = FALSE, keep.index = TRUE,
            keep.matrix = TRUE, keep.pvalues = TRUE))
cvector2<-unlist(c2$BootstrapCrossValCindexMat)
print(sd(cvector2))
brier<-pec(fgr, data = recoded, splitMethod="boot632", cause = 1, B=1000,
           verbose = FALSE)
print(brier)
ls<-brier$Boot632Err[2]
print(sd(unlist(ls), na.rm = TRUE))
print(quantile(unlist(ls), c(.025, .975), na.rm = TRUE))
print("END OF FINE & GRAY")
cat("\n\n\n\n") #for creating new lines in printed output


##Random Survival Forest (Ishwaran, 2013)
## Analysis 1
## modified Gray's weighted log-rank splitting
## (equivalent to cause=c(1,1) and splitrule="logrankCR")
print("RANDOM SURVIVAL FOREST")
fitform <- Surv(time, status) ~ idu + ageatfda + cd4nadir + black
o1 <- rfsrc(fitform, recoded, ntree = 1000)

##C-index
c_index_test <- cindex(o1, fitform, data=recoded, splitMethod="boot632",
                        cause = 1, 
                        B=1000, verbose = FALSE, keep.index = TRUE,
                       keep.matrix = TRUE, keep.pvalues = TRUE)
print(c_index_test)
cvector3<-unlist(c_index_test$BootstrapCrossValCindexMat)
print(sd(cvector3))
brier<-pec(o1, data = recoded, splitMethod="boot632", cause = 1, B=1000,
           verbose = FALSE)
print(brier)
ls<-brier$Boot632Err[2]
print(sd(unlist(ls), na.rm = TRUE))
print(quantile(unlist(ls), c(.025, .975), na.rm = TRUE))

print("END OF RANDOM SURVIVAL FOREST")
cat("\n\n\n\n")




#####BART (Bayesian additive regression trees)
print("BAYESIAN ADDITIVE REGRESSION TREES")

set.seed(99)


split <- initial_split(recoded, prop = .7)
train <- training(split)
times <- pmax(1, ceiling(train$time)) ## years
delta <- train$status
test  <- testing(split)


train <- data.matrix(train[,-c(1,2)])
test <- data.matrix(test[,-c(1,2)])

pre <- crisk.pre.bart(x.train=train, times=times, delta=delta, x.test=test)

post <- crisk.bart(x.train=train, times=times, delta=delta, x.test=test)


Cindex=function(risk, times, delta=NULL)
{   
  N=length(risk)
  if(N!=length(times))
    stop('risk and times must be the same length')
  if(length(delta)==0) delta=rep(1, N)
  else if(N!=length(delta))
    stop('risk and delta must be the same length')
  
  l=0
  k=0
  for(i in 1:N) {
    h=which((times[i]==times & delta[i]>delta) |
              (times[i]<times & delta[i]>0))
    if(length(h)>0) {
      l=l+sum(risk[i]<risk[h])
      k=k+length(h)
    }
  }
  return(l/k)
} 


ttrain <- post$tx.train
times<-ttrain[,1]
ttrain[,1]
ttest <- post$tx.test
times<-ttest[,1]
print("Bootstrap c-index")
print(boot(data=c(post$prob.test.mean, times), statistic=Cindex, R=1000))
print("END OF BART")
cat("\n\n\n\n")