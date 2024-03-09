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
library(survsim)



#Import dataset
cmpdata <- melanoma
skim(cmpdata)
recoded <- cmpdata
##recode so that 0 = Alive, 1 = Melanoma death, 2= Non-melanoma death
recoded$status[recoded$status==2] <- 0
recoded$status[recoded$status>1] <- 2
recoded$sex <- as.factor(recoded$sex)
recoded$ulcer <- as.factor(recoded$ulcer)
table(recoded$status)

#table 1
cmpdata.labeled <- recoded
table(recoded$status)
cmpdata.labeled$status <- 
  factor(recoded$status, 
         levels=c(0,1,2),
         labels=c("Alive",
                  "Melanoma death", # Reference
                  "Non-melanoma death"))
cmpdata.labeled$ulcer <- 
  factor(recoded$ulcer, 
         levels=c(0,1),
         labels=c("Ulcer-free", # Reference
                  "Ulcer"))
cmpdata.labeled$time <- recoded$time / 365

table1(~ time + sex + age + year + thickness + ulcer | status, data=cmpdata.labeled)

recoded <- recoded %>% mutate_at(c("age", "year", "thickness"), 
                                        ~(scales::rescale(.) %>% as.vector))
write_feather(recoded, "melanoma.feather")




#Cox PH (Cause-specific Hazards model)
print("COX PROPORTIONAL HAZARDS")
cmpdatacen <- recoded %>%  #censoring competing risks (non-melanoma death)
  mutate(status = replace(status, status != 1, 0))
# table(cmpdatacen$status)
# table(recoded$status)
#fit model
coxmodel<- coxph(Surv(time, status == 1) ~ sex + age + thickness + ulcer, 
      data = cmpdatacen, x = TRUE)
cox.zph(coxmodel)

cf=calPlot(coxmodel,data=cmpdatacen,legend.legend="Cause-specific Cox regression")

# modelsum <- summary(coxmodel)
# print(mean(modelsum$residuals^2))
#C-index
c_melcox<- cindex(coxmodel, formula = Surv(time, status) ~
               sex + age + thickness + ulcer,
             data=cmpdatacen, splitMethod="boot632", B=1000, confInt = TRUE,verbose = FALSE, keep.index = TRUE,
             keep.matrix = TRUE, keep.pvalues = TRUE)
print(c_melcox)
cvector<-unlist(c_melcox$BootstrapCrossValCindexMat)
print(sd(cvector))

brier<-pec(coxmodel, data = cmpdatacen, splitMethod="boot632", B=1000,
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
              #group   = cmpdata.labeled$ulcer,  # estimates will calculated within groups
              ## strata  = ,  # Tests will be stratified on this variable.
              rho     = 0, # Power of the weight function used in the tests.
              cencode = 2, # value of fstatus variable which indicates the failure time is censored.
              ## subset = ,
              ## na.action = na.omit
)
plot(CIF,color=1:6)
title("Melanoma")
#c-index     
fgr <- FGR(Hist(time,status)~sex + age + thickness + ulcer,
         data=recoded, cause = 1) ##Fine-Gray
c<-(cindex(fgr, data = recoded, splitMethod="boot632", cause = 1, B=1000,
            verbose = FALSE, keep.index = TRUE,
           keep.matrix = TRUE, keep.pvalues = TRUE))
print(c)
cvector2<-unlist(c$BootstrapCrossValCindexMat)
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
fitform <- Surv(time, status) ~ sex + age + thickness + ulcer
o1 <- rfsrc(fitform, recoded, ntree = 1000)

##C-index
c2<-(cindex(o1, fitform, data=recoded, splitMethod="boot632",
                        cause = 1, 
                        B= 1000, verbose = FALSE, keep.index = TRUE,
            keep.matrix = TRUE, keep.pvalues = TRUE))
print(c2)
cvector3<-unlist(c2$BootstrapCrossValCindexMat)
print(sd(cvector3))
brier<-pec(o1, fitform, data = recoded, splitMethod="boot632", B=1000,
           verbose = TRUE)
print(o1)

print(brier)
ls<-brier$BootCvErr[2]
print(sd(unlist(ls), na.rm = TRUE))
print(quantile(unlist(ls), c(.025, .975), na.rm = TRUE))
print("END OF RANDOM SURVIVAL FOREST")
cat("\n\n\n\n")


saveRDS(cvector, file = "cmel.rds")
saveRDS(cvector2, file = "cmel2.rds")
saveRDS(cvector3, file = "cmel3.rds")
save.image()
t<-t.test(cvector2, cvector3, paired = TRUE, alternative = "two.sided")
t
diff<-(cvector2-cvector3)
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
times <- pmax(1, ceiling(train$time/365)) ## years
delta <- train$status
test  <- testing(split)
timestest <- pmax(1, ceiling(test$time/365)) ## years
deltatest <- test$status


testdf <- test[,-c(1,2,5)]


train <- data.matrix(train[,-c(1,2,5)])
test <- data.matrix(test[,-c(1,2,5)])

pre <- crisk.pre.bart(x.train=train, times=times, delta=delta, x.test=test)

post <- crisk.bart(x.train=train, times=times, delta=delta, x.test=test)

post$cif.test.mean

# Output predicted probabilities on grid from BART, where post is the BART object that was returned
cif.bart.pred <- matrix(post$cif.test.mean, nrow=length(test),
                        ncol=length(timestest),byrow=TRUE)
# select timepoint and timeindex to use for Brier score or C index, as well as the timepoint from the prediction model

for(tindx in (1:4)) {
  
  timept <- tindx*26
  timeindx <- max(which(timestest <=timept))
  
  # run cindex function, passing bart predictions as matrix, selecting the timepoint of interest, but also passing it the survival times (times.test) and event indicator (delta.test)
    
  cindex.bart <- cindex(object=as.matrix(cif.bart.pred[,timeindx],length(test),1),
                        formula=Surv(timestest,deltatest)~1,
                        data=testdf,
                        cens.model="marginal",
                        eval.times=timestest[timeindx],
                        cause=1)
  cstat.bart[tindx] <- cindex.bart$AppCindex$matrix
  
}

#x1=Score(as.list(post$prob.train.mean),formula=fitform,data=recoded,
 #        cause =1,split.method="bootcv",B=100)
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


