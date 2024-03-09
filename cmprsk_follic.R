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
data(follic)
cmpdata <- follic
skim(cmpdata)
recoded <- cmpdata
##0 = Alive, 1 = Relapse, 2= Death
table(recoded$status)
recoded$clinstg <- as.factor(recoded$clinstg)
recoded$ch <- as.factor(recoded$ch)

#table 1
cmpdata.labeled <- recoded
cmpdata.labeled$status <- 
  factor(recoded$status, 
         levels=c(0,1,2),
         labels=c("Alive",
                  "Relapse", # Reference
                  "Death"))
table1(~ age + hgb + clinstg + ch + rt + time | status, data=cmpdata.labeled)
recoded <- recoded %>% mutate_at(c("age", "hgb"), 
                                        ~(scales::rescale(.) %>% as.vector))
write_feather(recoded, "follic.feather")


#Cox PH (Cause-specific Hazards model)
print("COX PROPORTIONAL HAZARDS")
cmpdatacen <- recoded %>%  #censoring competing risks (Death)
  mutate(status = replace(status, status != 1, 0))
skim(cmpdatacen)
#fit model
coxmodel<- coxph(Surv(time, status) ~ age + hgb + clinstg + ch,
                 data = cmpdatacen, x = TRUE)
cox.zph(coxmodel)

#C-index
c<-(cindex(coxmodel, formula = Surv(time, status) ~
               age + hgb + clinstg + ch,
             data=cmpdatacen, splitMethod="boot632", B=1000,verbose = FALSE, keep.index = TRUE,
             keep.matrix = TRUE, keep.pvalues = TRUE))
print(c)
cvector<-unlist(c$BootstrapCrossValCindexMat)
print(sd(cvector))

brier<-pec(coxmodel, data = cmpdatacen, splitMethod="boot632", cause = 1, B=1000,
           verbose = FALSE)
print(brier)
ls<-brier$Boot632Err[2]
print(sd(unlist(ls), na.rm = TRUE))
print(quantile(unlist(ls), c(.025, .975), na.rm = TRUE))
print("END OF COXPH")
cat("\n\n\n\n") #for creating new lines in printed output



#Fine & Gray model (Subdistribution model)
print("FINE AND GRAY MODEL")
#visualise CIF
CIF <- cuminc(ftime   = cmpdata.labeled$time,  # failure time variable (years)
              fstatus = cmpdata.labeled$status,  # variable with distinct codes for different causes of failure
              # group   = cmpdata.labeled$clinstg,  # estimates will calculated within groups
              ## strata  = ,  # Tests will be stratified on this variable.
              rho     = 0, # Power of the weight function used in the tests.
              cencode = 1, # value of fstatus variable which indicates the failure time is censored.
              ## subset = ,
              ## na.action = na.omit
)
plot(CIF,color=1:6)
title("Follic")
#c-index     
fgr <- FGR(Hist(time,status)~age + hgb + clinstg + ch,
           data=recoded, cause = 1) ##Fine-Gray
c2<-(cindex(fgr, data = recoded, splitMethod="boot632", B=1000,
             verbose = FALSE, keep.index = TRUE,
             keep.matrix = TRUE, keep.pvalues = TRUE))
print(c2)
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
fitform <- Surv(time, status) ~ age + hgb + clinstg + ch
o1 <- rfsrc(fitform, recoded, ntree = 1000)


##C-index
c3<-(cindex(o1, fitform, data=recoded, splitMethod="boot632",
                        cause = 1, 
                        B=1000, verbose = FALSE, keep.index = TRUE,
             keep.matrix = TRUE, keep.pvalues = TRUE))
print(c3)
cvector3<-unlist(c3$BootstrapCrossValCindexMat)
print(sd(cvector3))

brier<-pec(o1, data = recoded, splitMethod="boot632", cause = 1, B=1000,
           verbose = FALSE)
print(brier)
ls<-brier$Boot632Err[2]
print(sd(unlist(ls), na.rm = TRUE))
print(quantile(unlist(ls), c(.025, .975), na.rm = TRUE))
print("END OF RANDOM SURVIVAL FOREST")
cat("\n\n\n\n")


saveRDS(cvector, file = "cfol.rds")
saveRDS(cvector2, file = "cfol2.rds")
saveRDS(cvector3, file = "cfol3.rds")
save.image()
t<-t.test(cvector, cvector3, paired = TRUE, alternative = "two.sided")
diff<-(cvector-cvector3)
mean<-mean(diff)
variance<-var(diff)
tstat<-mean/sqrt(variance*(1/1000 + 999))
p<-pt(tstat, 999, lower.tail = FALSE)
p
p.adjust(p, method = "fdr", n = length(p))



#####BART (Bayesian additive regression trees)
print("BAYESIAN ADDITIVE REGRESSION TREES")
set.seed(99)

split <- initial_split(recoded, prop = .7)
train <- training(split)
times <- pmax(1, ceiling(train$time)) ## years
delta <- train$status
test  <- testing(split)


train <- data.matrix(train[,-c(5:7)])
test <- data.matrix(test[,-c(5:7)])

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
print("Bootstrap c-index")
print(boot(data=c(post$prob.train.mean, times), statistic=Cindex,
           R=1000))
print("END OF BART")
cat("\n\n\n\n")



