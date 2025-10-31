Functions are provided to estimate VUS and cause-specific ExCs for overall and cause-specific evaluations of competing risk models. 

Input:

To estimate both VUS and ExC, users need to provide a data frame by specifying "mydata" that consists of five columns in the following order:

mydata[,1]=observation times

mydata[,2]=cause indicators, where 0 indicates right censoring, 1 for event 1 and 2 for event 2.

mydata[,3]=CIF estimates of a cause-1 event obtained from the model to be evaluated. The CIF should be estimated at "predicttime".

mydata[,4]=CIF estimates of a cause-2 event obtained from the model to be evaluated. The CIF should be estimated at "predicttime".

mydata[,5]=estimates of survival function obtained from the model to be evaluated and estimated at "predicttime".

Both VUS and ExC also require to provide a "predicttime", and the models will be evaluated at this time point since both are time-dependent model evaluation metrics.

Cause-specific ExC additionally requires users to specify the cause of interest by specifying "cause=".

Output:

Two functions are included in VUS_new.R. VUS(mydata, predicttime) will produce a list that contains 

$est: VUS estimate

$asd: estimate of asymptotic standard deviation

$infun: estimate of influence function.

To save time and if asymptotic standard deviation is not of interest, a simple version vus.sim(mydata, predicttime) is provided that only estimates VUS.

To estimate cause-specific ExC, call ExCl_pdi(mydata, predictime, cause). The output is the estimate of cause-specific of ExC for the cause of interest.

Example based on melanoma data without cross fitting:

source("VUS_new.R")

source("AUC_est.R")

library(riskRegression)

library(survival)

library(prodlim)

data(Melanoma)

original=Melanoma

#obsdata$sex1=ifelse(as.numeric(obsdata$sex)==1, 0, 1)

#obsdata$ulcer1=ifelse(as.numeric(obsdata$ulcer)==1, 0, 1)

folds=createFolds(1:nrow(original),k=2)

obsdata=original[folds[[1]], ] # train

test=original[folds[[2]], ]

predicttime=2000

obsdata1=obsdata

obsdata1$status=ifelse(obsdata$status==1, 1, 0)

obsdata2=obsdata

obsdata2$status=ifelse(obsdata$status==2, 1, 0)

\# model fitting

csc1=CSC(Hist(time, status)~sex+age+logthick+ulcer, data=obsdata, cause=1)

csc2=CSC(Hist(time, status)~sex+age+logthick+ulcer, data=obsdata, cause=2)

cscrisk1=predictRisk(csc1, test[,c("sex","age","logthick","ulcer")], times=predicttime, cause=1)

cscrisk2=predictRisk(csc2, test[,c("sex","age","logthick","ulcer")], times=predicttime, cause=2)

cscsurv=1-cscrisk1-cscrisk2

\# data frame construction

cscdata=data.frame(cbind(test$time, test$status, cscrisk1, cscrisk2, cscsurv))

VUS(mydata=cscdata, predicttime=predicttime) # estimate VUS and its standard error

ExCl_pdi(mydata=cscdata, predicttime, cause=1) # estimate ExC for cause-1 event

ExCl_pdi(mydata=cscdata, predicttime, cause=2) # estiamte ExC for cause-2 event
