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

#####Simulate Data
### A cohort with 50 subjects, with a maximum follow-up time of 100 days and two
### covariates, following Bernoulli distributions, and a corresponding beta of
### 0.1698695 and 0.0007010932 for each event for the first covariate and a
### corresponding beta of 0.3735146 and 0.5591244 for each event for the
### second covariate. Notice that the time to censorship is assumed to follow a
### log-normal distribution.
sim.data <- crisk.sim(n=1500, foltime=10, 
                      dist.ev=c("weibull","weibull"),
                      anc.ev=c(1.2, 1.12),
                      beta0.ev=c(1.80342, 1.35374),
                      dist.cens="unif",
                      anc.cens=10,
                      beta0.cens=0,
                      z=list(c("unif", 0.2,2), c("unif", 0.1, 1.6)),
                      beta=list(c(0.1698695,0.0007010932),
                                c(0.3735146,0.5591244),
                                c(0.0012239,0.8361621),
                                c(0.4139212,0.4362511),
                                c(0.3749190,0.2748191)),
                      x=list(c("bern", 0.381), 
                             c("bern", 0.564),
                             c("normal", 0.564, 0.167),
                             c("normal", 0.456, 0.02),
                             c("normal", 0.183, 0.06)),
                      nsit=2)
colnames(sim.data) <- c("nid", "status", "time", "status2","start","stop",
                        "z", "x_binary1","x_binary2", "x_cont1", "x_cont2", "x_cont3")
sim.data <- sim.data %>%
  select(-c(nid,status2,start,stop,z))
sim.data$x_binary1 <- as.factor(sim.data$x_binary1)
sim.data$x_binary2 <- as.factor(sim.data$x_binary2)

sim.data$status[is.na(sim.data$status)] <- 0
skim(sim.data)
table(sim.data$status)
write_feather(sim.data, "simdata.feather")

sim.data<-read_feather("simdata.feather")

#Import dataset
cmpdata <- sim.data
recoded <- cmpdata
table(recoded$status)
skim(recoded)
#table 1
cmpdata.labeled <- recoded
cmpdata.labeled$status <- 
  factor(recoded$status, 
         levels=c(0,1,2),
         labels=c("Alive",
                  "Focal Event", # Reference
                  "Competing Event"))
table1(~ . | status, data=cmpdata.labeled)



#Cox PH (Cause-specific Hazards model)
print("COX PROPORTIONAL HAZARDS")
cmpdatacen <- recoded %>%  #censoring competing risks (non-melanoma death)
  mutate(status = replace(status, status != 1, 0))
# table(cmpdatacen$status)
# table(recoded$status)
#fit model
coxmodel<- coxph(Surv(time, status == 1) ~ x_binary1 + x_binary2 + x_cont1 + x_cont2 + x_cont3, 
                 data = cmpdatacen, x = TRUE)
cox.zph(coxmodel)

# modelsum <- summary(coxmodel)
# print(mean(modelsum$residuals^2))
#C-index
c_simcox<- cindex(coxmodel, formula = Surv(time, status) ~.,
                  data=cmpdatacen, splitMethod="boot632", B=1000, confInt = TRUE,verbose = FALSE, 
                  keep.index = TRUE,
                  keep.matrix = TRUE, keep.pvalues = TRUE)
print(c_simcox)
cvector<-unlist(c_simcox$BootstrapCrossValCindexMat)
print(sd(cvector))

brier<-pec(coxmodel, data = recoded, splitMethod="boot632", cause = 1, B=1000,
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
plot(CIF, color=1:6)
title("Sim")
#c-index     
fgr <- FGR(Hist(time,status)~ x_binary1 + x_binary2 + x_cont1 + x_cont2 + x_cont3,
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
fitform <- Surv(time, status) ~ x_binary1 + x_binary2 + x_cont1 + x_cont2 + x_cont3
o1 <- rfsrc(fitform, recoded, ntree = 1000)


##C-index
c3<-(cindex(o1, fitform, data=recoded, splitMethod="boot632",
             cause = 1, 
             B= 1000, verbose = FALSE, keep.index = TRUE,
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




#####BART (Bayesian additive regression trees)
print("BAYESIAN ADDITIVE REGRESSION TREES")

set.seed(99)

split <- initial_split(scaled.cmpdata, prop = .7)
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
print("Bootstrap c-index")
print(boot(data=c(post$prob.train.mean, times), statistic=Cindex,
           R=1000))
print("END OF BART")
cat("\n\n\n\n")
