# based on the true VUS definition from Li & Fine 2008
VUS=function(mydata,predicttime){
  # mydata supposed to be a dataframe
  # mydata[,1] supposed to be observed time
  # mydata[,2] supposed to be cause indicator
  # mydata[, 3:5] supposed to be CIF1, CIF2, S at predicttime for each subject
  
  colnames(mydata)=c("obstime","causesta","cif1","cif2","s")
  mydata$trueclass=ifelse((mydata$obstime<=predicttime & mydata$causesta==1), 1, ifelse((mydata$obstime<=predicttime & mydata$causesta==2),2,ifelse((mydata$obstime>predicttime),3,NA))) 
  # for those censored before predicttime, trueclass is NA
  cpmdata=mydata[order(mydata$obstime),] # order by obseved time
  cpmdataobs<-na.omit(cpmdata) # remove those censored before predicttime
  n.nmiss=length(cpmdataobs$trueclass) # number of subjects not missing
  
  obsclass1=which(cpmdataobs$trueclass==1)                      ## T<t0 and event=1 and T<C
  obsclass2=which(cpmdataobs$trueclass==2)                      ## T<t0 and event=2 and T<C
  obsclass3=which(cpmdataobs$trueclass==3)                      ## T>t0
  
  nclass1=length(obsclass1)
  nclass2=length(obsclass2)
  nclass3=length(obsclass3)
  
  if(nclass1==0 | nclass2==0 | nclass3==0){
    return(0)
  }
  
  Tevent1=cpmdataobs$obstime[obsclass1] # event time for event=1
  Tevent2=cpmdataobs$obstime[obsclass2] # event time for event=2
  Tevent3=rep(predicttime, nclass3) # censored at predicttime
  Teventobs=c(Tevent1, Tevent2, Tevent3)
  obsevent.npts=length(Teventobs)
  
  # estimate survival function of censoring
  dim.time = length(cpmdata$obstime)
  censoring.index=which(cpmdata$causesta==0)
  censoring.survival=NULL # KM estimator of censoring
  censoring.time=unique(sort(cpmdata$obstime[censoring.index])) # sorted unique time of censoring
  censoring.npts=length(censoring.time)
  
  timepool=matrix(rep(cpmdata$obstime, each=censoring.npts),censoring.npts, dim.time) # each row is a copy of obstime, censoring.npts copies/rows
  statuspool=matrix(rep(cpmdata$causesta, each=censoring.npts),censoring.npts, dim.time) # each row is a copy of causesta, censoring.npts copies/rows
  censor.Y=apply(ifelse(timepool>=censoring.time, 1,0), 1, sum) # at risk set
  censor.d=apply(ifelse(timepool==censoring.time & statuspool==0, 1,0),1,sum) # event
  censoring.survival=cumprod(1-censor.d/censor.Y)
  censoring.time=c(0, censoring.time)
  censoring.survival=c(1,censoring.survival)
  
  censortimepool=matrix(rep(censoring.time, each=obsevent.npts), obsevent.npts, length(censoring.time))
  # each row is a copy of censoring.time, obsevent.npts number of rows
  
  max_index = function(uniquetime.event_time) 
  {
    vec.len=length(uniquetime.event_time)
    uniquetime = uniquetime.event_time[1:vec.len-1]
    event_time = uniquetime.event_time[vec.len]
    if(event_time<=min(uniquetime))
      return(1)
    else
      return(max(which(uniquetime<event_time)))
  }
  # for n event times, from the first n-1 points, find the last one that is <= the nth time
  max_index_noevent = function(uniquetime.event_time) 
  {
    vec.len=length(uniquetime.event_time)
    uniquetime = uniquetime.event_time[1:vec.len-1]
    event_time = uniquetime.event_time[vec.len]
    if(sum(uniquetime>=event_time)==0)
      return(vec.len-1)
    else
      return(min(which(uniquetime>=event_time)))
    
  }
  # for n event times, from the first n-1 times, find the first one that is >= the nth time.
  
  eventatcensordist.index = apply(cbind(censortimepool, Teventobs), 1, max_index)
  probcen.event=censoring.survival[eventatcensordist.index]
  
  probcen.event1=probcen.event[1:nclass1]
  probcen.event2=probcen.event[(nclass1+1):(nclass1+nclass2)]
  noeventatcensordist.index = apply(cbind(censortimepool, Teventobs), 1, max_index_noevent)
  probcen.noevent=censoring.survival[noeventatcensordist.index]
  probcen.event3=probcen.noevent[(nclass1+nclass2+1):(nclass1+nclass2+nclass3)]
  probcensorweig=kronecker(probcen.event1%*%t(probcen.event2),probcen.event3)
  
  class1comp=cpmdataobs[obsclass1,]
  class2comp=cpmdataobs[obsclass2,]
  class3comp=cpmdataobs[obsclass3,]
  
  class1comp$dist1=(class1comp$cif1-1)^2+(class1comp$cif2)^2+(class1comp$s)^2
  class1comp$dist2=(class1comp$cif1)^2+(class1comp$cif2-1)^2+(class1comp$s)^2
  class1comp$dist3=(class1comp$cif1)^2+(class1comp$cif2)^2+(class1comp$s-1)^2
  
  class2comp$dist1=(class2comp$cif1-1)^2+(class2comp$cif2)^2+(class2comp$s)^2
  class2comp$dist2=(class2comp$cif1)^2+(class2comp$cif2-1)^2+(class2comp$s)^2
  class2comp$dist3=(class2comp$cif1)^2+(class2comp$cif2)^2+(class2comp$s-1)^2
  
  class3comp$dist1=(class3comp$cif1-1)^2+(class3comp$cif2)^2+(class3comp$s)^2
  class3comp$dist2=(class3comp$cif1)^2+(class3comp$cif2-1)^2+(class3comp$s)^2
  class3comp$dist3=(class3comp$cif1)^2+(class3comp$cif2)^2+(class3comp$s-1)^2
  
  Identcla3=rep(1, nclass3)
  Identcla12=matrix(1, nclass1, nclass2)
  class1m1=kronecker(matrix(rep(class1comp$dist1,nclass2), nclass1, nclass2), Identcla3)
  class1m2=kronecker(matrix(rep(class1comp$dist2,nclass2), nclass1, nclass2), Identcla3)
  class1m3=kronecker(matrix(rep(class1comp$dist3,nclass2), nclass1, nclass2), Identcla3)
  # inner: each column is a copy of y[obsclass1], nclass2 copies.
  # stack Identcla3 copies of the inner matrix
  class2m1=kronecker(matrix(rep(class2comp$dist1,each=nclass1), nclass1, nclass2), Identcla3)
  class2m2=kronecker(matrix(rep(class2comp$dist2,each=nclass1), nclass1, nclass2), Identcla3)
  class2m3=kronecker(matrix(rep(class2comp$dist3,each=nclass1), nclass1, nclass2), Identcla3)
  # inner: each row is a copy of y[obsclass2], nclass2 copies
  # stack Identcla3 copies of inner matrix
  class3m1=kronecker(Identcla12,class3comp$dist1)
  class3m2=kronecker(Identcla12,class3comp$dist2)
  class3m3=kronecker(Identcla12,class3comp$dist3)
  # each matrix is nclass1*nclass2 with all elements equal to y[obsclass3][i]
  # stack all matrices
  m1=class1m1+class2m2+class3m3
  m2=class1m1+class2m3+class3m2
  m3=class1m2+class2m1+class3m3
  m4=class1m2+class2m3+class3m1
  m5=class1m3+class2m1+class3m2
  m6=class1m3+class2m2+class3m1
  
  
  Num.IPCW=sum(as.numeric(m1<=m2 & m1<=m3 & m1<=m4 & m1<=m5 & m1<=m6)/probcensorweig)
  
  event1n=rep(1, nclass1)
  event2n=rep(1, nclass2)
  event3n=rep(1, nclass3)
  Demon.IPCW=sum(kronecker(event1n%*%t(event2n),event3n)/probcensorweig)
  
  weightedHUM=Num.IPCW/Demon.IPCW
  # return(weightedHUM)
  
  ##################################################################################
  # Asymptotic Error Estimation
  ##################################################################################
  
  class1=which(cpmdata$trueclass==1)                      ## T<t0 and event=1 and T<C
  
  class2=which(cpmdata$trueclass==2)                      ## T<t0 and event=2 and T<C
  
  class3=which(cpmdata$trueclass==3)                      ## T>t0 
  
  
  cpmdata$Tobs=cpmdata$obstime  
  
  cpmdata$Tobs[class3]=predicttime
  
  
  censortimepooltot=matrix(rep(censoring.time, each=dim.time), dim.time, length(censoring.time))
  
  
  
  Tobsatcensordist.index = apply(cbind(censortimepooltot, cpmdata$Tobs), 1, max_index)   
  
  probobscen.event=censoring.survival[Tobsatcensordist.index] 
  
  noTobsatcensordist.index = apply(cbind(censortimepooltot, cpmdata$Tobs), 1, max_index_noevent)  
  
  noevent.index=noTobsatcensordist.index[class3] 
  
  probobscen.event[class3]=censoring.survival[noevent.index] 
  
  
  
  
  #############################################################################################################################################################
  ###calculate matingale matrix for survival function of censoring G()                                                                                      ##
  
  #############################################################################################################################################################
  
  times=cpmdata$obstime # already sorted at the beginning
  
  status=cpmdata$causesta
  
  
  computeMC=function(times, status){
    indCens=as.numeric(status==0) # indicator of censoring
    n=length(times) # number of observations
    
    
    ##calculate survival function for censoring cases
    hatSurvCens=1-cumsum(indCens)/n
    
    hatSurv=1-(1:n)/n # survival function of X=min(T,C)
    
    temp1=cbind(times, indCens, hatSurv, hatSurvCens)
    temp1=rbind(c(0, 0, 1,1), temp1)
    colnames(temp1)=c("T", "Ind4Censor", "hatSurv", "hatSurvCens")
    
    
    ##calculate hazard rate for censoring cases
    lambdaC=(temp1[-1, "Ind4Censor"])/(n:1)
    temp1=cbind(temp1, c(0, lambdaC))
    colnames(temp1)[ncol(temp1)]="lambdaC"
    
    
    ##calculate cumulative hazard rate for censoring cases
    LambdaC=cumsum(lambdaC)
    temp1=cbind(temp1, c(0, LambdaC))
    colnames(temp1)[ncol(temp1)]="LambdaC"
    
    temp2=temp1[-1,]
    
    
    
    ##calculate martingale matrix
    # temp2: T, ind4censor, hatsurv, hatsurvcens, lambdac, LambdaC, ordered by T
    hatMC=matrix(NA, n,n)
    for (i in 1:n){
      hatMC[,i]=temp2[i, 2]*as.numeric(temp2[i, 1]<=temp2[, 1])- c(temp2[0:i, 6], rep(temp2[i, 6], (n-i))) # N(t)-Lambda(t) at all time points for subject i
    }
    
    
    ##calculate d(martingale matrix) 
    
    dhatMC=rbind(hatMC[1,], hatMC[-1, ]-hatMC[-nrow(hatMC),]) # dMc(t)
    
    # column: subject. Row: time points
    
    
    ##calculate d(martingale matrix)/survival prob
    
    dhatbySurv=function(v){
      n=length(v)
      v/c(1, 1-(1:(n-1))/n)
    }
    
    dhatdivSurv=apply(dhatMC, 2, dhatbySurv)
    
    
    ##calculate intergration of d(martingale matrix)/survival prob
    
    MatInt0TCidhatMC=apply(dhatdivSurv, 2, cumsum)
    colnames(MatInt0TCidhatMC)=paste("M_{C_", 1: length(times),"}", sep="")
    rownames(MatInt0TCidhatMC)=times
    return(MatInt0TCidhatMC)
  }
  
  MatInt0TcidhatMCksurEff =computeMC(times=cpmdata$obstime,status=cpmdata$causesta)
  
  #############################################################################################################################################################
  ###develop influence function of VUS                                                                                                                      ##
  
  #############################################################################################################################################################
  cpmdata$dist1=(cpmdata$cif1-1)^2+(cpmdata$cif2)^2+(cpmdata$s)^2
  cpmdata$dist2=(cpmdata$cif1)^2+(cpmdata$cif2-1)^2+(cpmdata$s)^2
  cpmdata$dist3=(cpmdata$cif1)^2+(cpmdata$cif2)^2+(cpmdata$s-1)^2
  
  cpmdata$weight=1/probobscen.event
  
  # event1=cpmdata[class1,]
  # event2=cpmdata[class2,]
  # eventfree=cpmdata[class3,]
  
  
  n23=nclass2*nclass3
  n13=nclass1*nclass3
  n12=nclass1*nclass2
  
  ##########
  # A_{ijk}
  ##########
  Aijk=function(i,cpmdata){
    event1=cpmdata[class1,]
    event2=cpmdata[class2,]
    eventfree=cpmdata[class3,]
    i1=rep(event1$dist1[i], n23)
    i2=rep(event1$dist2[i], n23)
    i3=rep(event1$dist3[i], n23)
    j1=rep(event2$dist1, each=nclass3)
    j2=rep(event2$dist2, each=nclass3)
    j3=rep(event2$dist3, each=nclass3)
    k1=rep(eventfree$dist1, nclass2)
    k2=rep(eventfree$dist2, nclass2)
    k3=rep(eventfree$dist3, nclass2)
    d1=i1+j2+k3
    d2=i1+j3+k2
    d3=i2+j1+k3
    d4=i2+j3+k1
    d5=i3+j1+k2
    d6=i3+j2+k1
    wi=rep(event1$weight[i], n23)
    wj=rep(event2$weight, each=nclass3)
    wk=rep(eventfree$weight, nclass2)
    return(as.numeric(d1<=d2 & d1<=d3 & d1<=d4 & d1<=d5 & d1<=d6)*wi*wj*wk)
    # should return a vector: Aijk with fixed i
  }
  
  Bijk=function(i,cpmdata){
    event1=cpmdata[class1,]
    event2=cpmdata[class2,]
    eventfree=cpmdata[class3,]
    wi=rep(event1$weight[i], n23)
    wj=rep(event2$weight, each=nclass3)
    wk=rep(eventfree$weight, nclass2)
    return(wi*wj*wk)
  }
  
  MatAtijk=matrix(NA, n23, nclass1)
  MatBtijk=matrix(NA, n23, nclass1)
  
  for(i in 1:nclass1){
    MatAtijk[,i]=Aijk(i,cpmdata)
    MatBtijk[,i]=Bijk(i,cpmdata)
  }
  
  
  ##########
  # A_{jik}
  ##########
  Ajik=function(i,cpmdata){
    event1=cpmdata[class1,]
    event2=cpmdata[class2,]
    eventfree=cpmdata[class3,]
    j1=rep(event1$dist1, each=nclass3)
    j2=rep(event1$dist2, each=nclass3)
    j3=rep(event1$dist3, each=nclass3)
    i1=rep(event2$dist1[i], n13)
    i2=rep(event2$dist2[i], n13)
    i3=rep(event2$dist3[i], n13)
    k1=rep(eventfree$dist1, nclass1)
    k2=rep(eventfree$dist2, nclass1)
    k3=rep(eventfree$dist3, nclass1)
    d1=j1+i2+k3
    d2=j1+i3+k2
    d3=j2+i1+k3
    d4=j2+i3+k1
    d5=j3+i1+k2
    d6=j3+i2+k1
    wj=rep(event1$weight, each=nclass3)
    wi=rep(event2$weight[i], n13)
    wk=rep(eventfree$weight, nclass1)
    return(as.numeric(d1<=d2 & d1<=d3 & d1<=d4 & d1<=d5 & d1<=d6)*wj*wi*wk)
    # should return a vector: Aijk with fixed i
  }
  
  Bjik=function(i,cpmdata){
    event1=cpmdata[class1,]
    event2=cpmdata[class2,]
    eventfree=cpmdata[class3,]
    wj=rep(event1$weight, each=nclass3)
    wi=rep(event2$weight[i], n13)
    wk=rep(eventfree$weight, nclass1)
    return(wi*wj*wk)
  }
  
  MatAtjik=matrix(NA, n13, nclass2)
  MatBtjik=matrix(NA, n13, nclass2)
  
  for(i in 1:nclass2){
    MatAtjik[,i]=Ajik(i,cpmdata)
    MatBtjik[,i]=Bjik(i,cpmdata)
  }
  
  ##########
  # A_{jki}
  ##########
  Ajki=function(i,cpmdata){
    event1=cpmdata[class1,]
    event2=cpmdata[class2,]
    eventfree=cpmdata[class3,]
    j1=rep(event1$dist1, each=nclass2)
    j2=rep(event1$dist2, each=nclass2)
    j3=rep(event1$dist3, each=nclass2)
    k1=rep(event2$dist1, nclass1)
    k2=rep(event2$dist2, nclass1)
    k3=rep(event2$dist3, nclass1)
    i1=rep(eventfree$dist1[i], n12)
    i2=rep(eventfree$dist2[i], n12)
    i3=rep(eventfree$dist3[i], n12)
    d1=j1+k2+i3
    d2=j1+k3+i2
    d3=j2+k1+i3
    d4=j2+k3+i1
    d5=j3+k1+i2
    d6=j3+k2+i1
    wj=rep(event1$weight, each=nclass2)
    wk=rep(event2$weight, nclass1)
    wi=rep(eventfree$weight[i], n12)
    return(as.numeric(d1<=d2 & d1<=d3 & d1<=d4 & d1<=d5 & d1<=d6)*wj*wk*wi)
    # should return a vector: Aijk with fixed i
  }
  
  Bjki=function(i,cpmdata){
    event1=cpmdata[class1,]
    event2=cpmdata[class2,]
    eventfree=cpmdata[class3,]
    wj=rep(event1$weight, each=nclass2)
    wk=rep(event2$weight, nclass1)
    wi=rep(eventfree$weight[i], n12)
    return(wj*wk*wi)
  }
  
  MatAtjki=matrix(NA, n12, nclass3)
  MatBtjki=matrix(NA, n12, nclass3)
  
  for(i in 1:nclass3){
    MatAtjki[,i]=Ajki(i,cpmdata)
    MatBtjki[,i]=Bjki(i,cpmdata)
  }
  
  ####################################
  # A_{jks}(dMc)
  ####################################
  # use MatAtjki as MatAtjks
  martingale1=MatInt0TcidhatMCksurEff[class1,]
  martingale2=MatInt0TcidhatMCksurEff[class2,]
  martingale3=MatInt0TcidhatMCksurEff[class3,]
  if(is.vector(martingale1)){
    martingale1=matrix(martingale1, nrow=1, ncol=100, byrow=T)
  }
  if(is.vector(martingale2)){
    martingale2=matrix(martingale2, nrow=1, ncol=100, byrow=T)
  }
  if(is.vector(martingale3)){
    martingale3=matrix(martingale3, nrow=1, ncol=100, byrow=T)
  }
  
  nobs=nrow(cpmdata)

  AdM=rep(0, nobs)
  BdM=rep(0, nobs)
  
  for(i in 1:nobs){
    # print(i)
    xj=kronecker(rep(martingale1[,i], each=nclass2), t(rep(1, nclass3)))
    xk=kronecker(rep(martingale2[,i], nclass1), t(rep(1, nclass3)))
    xs=kronecker(rep(1, n12), t(martingale3[,i]))
    AdM[i]=sum(MatAtjki*(xj+xk+xs))/(nobs*nobs*nobs)
    BdM[i]=sum(MatBtjki*(xj+xk+xs))/(nobs*nobs*nobs)
  }
  
  ########################################
  # Final Step for Influence Function Est
  ########################################
  aveMatAt=(sum(MatAtijk))/(nobs*nobs*nobs) # estimator for numerator in VUS.hat, Q in proof
  aveMatBt=(sum(MatBtijk))/(nobs*nobs*nobs) # estimator for denominator in VUS.hat, Q0in proof
  
  cpmdata$ijks=0
  cpmdata$ijks[class1]=(colSums(MatAtijk)/(nobs*nobs)-aveMatAt-aveMatAt/aveMatBt*(colSums(MatBtijk)/(nobs*nobs)-aveMatBt))/aveMatBt
  
  cpmdata$jiks=0
  cpmdata$jiks[class2]=(colSums(MatAtjik)/(nobs*nobs)-aveMatAt-aveMatAt/aveMatBt*(colSums(MatBtjik)/(nobs*nobs)-aveMatBt))/aveMatBt
  
  cpmdata$jkis=0
  cpmdata$jkis[class3]=(colSums(MatAtjki)/(nobs*nobs)-aveMatAt-aveMatAt/aveMatBt*(colSums(MatBtjki)/(nobs*nobs)-aveMatBt))/aveMatBt
  
  cpmdata$jksi=(AdM-aveMatAt/aveMatBt*BdM)
  
  cpmdata$infun=cpmdata$ijks+cpmdata$jiks+cpmdata$jkis+cpmdata$jksi
  
  return(list(est=weightedHUM, asd=sd(cpmdata$infun)/sqrt(nobs), infun=cpmdata$infun))
}

vus.sim=function(mydata, predicttime){
  # mydata supposed to be a dataframe
  # mydata[,1] supposed to be observed time
  # mydata[,2] supposed to be cause indicator
  # mydata[, 3:5] supposed to be CIF1, CIF2, S at predicttime for each subject
  
  colnames(mydata)=c("obstime","causesta","cif1","cif2","s")
  mydata$trueclass=ifelse((mydata$obstime<=predicttime & mydata$causesta==1), 1, ifelse((mydata$obstime<=predicttime & mydata$causesta==2),2,ifelse((mydata$obstime>predicttime),3,NA))) 
  # for those censored before predicttime, trueclass is NA
  cpmdata=mydata[order(mydata$obstime),] # order by obseved time
  cpmdataobs<-na.omit(cpmdata) # remove those censored before predicttime
  n.nmiss=length(cpmdataobs$trueclass) # number of subjects not missing
  
  obsclass1=which(cpmdataobs$trueclass==1)                      ## T<t0 and event=1 and T<C
  obsclass2=which(cpmdataobs$trueclass==2)                      ## T<t0 and event=2 and T<C
  obsclass3=which(cpmdataobs$trueclass==3)                      ## T>t0
  
  nclass1=length(obsclass1)
  nclass2=length(obsclass2)
  nclass3=length(obsclass3)
  
  if(nclass1==0 | nclass2==0 | nclass3==0){
    return(0)
  }
  
  Tevent1=cpmdataobs$obstime[obsclass1] # event time for event=1
  Tevent2=cpmdataobs$obstime[obsclass2] # event time for event=2
  Tevent3=rep(predicttime, nclass3) # censored at predicttime
  Teventobs=c(Tevent1, Tevent2, Tevent3)
  obsevent.npts=length(Teventobs)
  
  # estimate survival function of censoring
  dim.time = length(cpmdata$obstime)
  censoring.index=which(cpmdata$causesta==0)
  censoring.survival=NULL # KM estimator of censoring
  censoring.time=unique(sort(cpmdata$obstime[censoring.index])) # sorted unique time of censoring
  censoring.npts=length(censoring.time)
  
  timepool=matrix(rep(cpmdata$obstime, each=censoring.npts),censoring.npts, dim.time) # each row is a copy of obstime, censoring.npts copies/rows
  statuspool=matrix(rep(cpmdata$causesta, each=censoring.npts),censoring.npts, dim.time) # each row is a copy of causesta, censoring.npts copies/rows
  censor.Y=apply(ifelse(timepool>=censoring.time, 1,0), 1, sum) # at risk set
  censor.d=apply(ifelse(timepool==censoring.time & statuspool==0, 1,0),1,sum) # event
  censoring.survival=cumprod(1-censor.d/censor.Y)
  censoring.time=c(0, censoring.time)
  censoring.survival=c(1,censoring.survival)
  
  censortimepool=matrix(rep(censoring.time, each=obsevent.npts), obsevent.npts, length(censoring.time))
  # each row is a copy of censoring.time, obsevent.npts number of rows
  
  max_index = function(uniquetime.event_time) 
  {
    vec.len=length(uniquetime.event_time)
    uniquetime = uniquetime.event_time[1:vec.len-1]
    event_time = uniquetime.event_time[vec.len]
    if(event_time<=min(uniquetime))
      return(1)
    else
      return(max(which(uniquetime<event_time)))
  }
  # for n event times, from the first n-1 points, find the last one that is <= the nth time
  max_index_noevent = function(uniquetime.event_time) 
  {
    vec.len=length(uniquetime.event_time)
    uniquetime = uniquetime.event_time[1:vec.len-1]
    event_time = uniquetime.event_time[vec.len]
    if(sum(uniquetime>=event_time)==0)
      return(vec.len-1)
    else
      return(min(which(uniquetime>=event_time)))
    
  }
  # for n event times, from the first n-1 times, find the first one that is >= the nth time.
  
  eventatcensordist.index = apply(cbind(censortimepool, Teventobs), 1, max_index)
  probcen.event=censoring.survival[eventatcensordist.index]
  
  probcen.event1=probcen.event[1:nclass1]
  probcen.event2=probcen.event[(nclass1+1):(nclass1+nclass2)]
  noeventatcensordist.index = apply(cbind(censortimepool, Teventobs), 1, max_index_noevent)
  probcen.noevent=censoring.survival[noeventatcensordist.index]
  probcen.event3=probcen.noevent[(nclass1+nclass2+1):(nclass1+nclass2+nclass3)]
  probcensorweig=kronecker(probcen.event1%*%t(probcen.event2),probcen.event3)
  
  class1comp=cpmdataobs[obsclass1,]
  class2comp=cpmdataobs[obsclass2,]
  class3comp=cpmdataobs[obsclass3,]
  
  class1comp$dist1=(class1comp$cif1-1)^2+(class1comp$cif2)^2+(class1comp$s)^2
  class1comp$dist2=(class1comp$cif1)^2+(class1comp$cif2-1)^2+(class1comp$s)^2
  class1comp$dist3=(class1comp$cif1)^2+(class1comp$cif2)^2+(class1comp$s-1)^2
  
  class2comp$dist1=(class2comp$cif1-1)^2+(class2comp$cif2)^2+(class2comp$s)^2
  class2comp$dist2=(class2comp$cif1)^2+(class2comp$cif2-1)^2+(class2comp$s)^2
  class2comp$dist3=(class2comp$cif1)^2+(class2comp$cif2)^2+(class2comp$s-1)^2
  
  class3comp$dist1=(class3comp$cif1-1)^2+(class3comp$cif2)^2+(class3comp$s)^2
  class3comp$dist2=(class3comp$cif1)^2+(class3comp$cif2-1)^2+(class3comp$s)^2
  class3comp$dist3=(class3comp$cif1)^2+(class3comp$cif2)^2+(class3comp$s-1)^2
  
  Identcla3=rep(1, nclass3)
  Identcla12=matrix(1, nclass1, nclass2)
  class1m1=kronecker(matrix(rep(class1comp$dist1,nclass2), nclass1, nclass2), Identcla3)
  class1m2=kronecker(matrix(rep(class1comp$dist2,nclass2), nclass1, nclass2), Identcla3)
  class1m3=kronecker(matrix(rep(class1comp$dist3,nclass2), nclass1, nclass2), Identcla3)
  # inner: each column is a copy of y[obsclass1], nclass2 copies.
  # stack Identcla3 copies of the inner matrix
  class2m1=kronecker(matrix(rep(class2comp$dist1,each=nclass1), nclass1, nclass2), Identcla3)
  class2m2=kronecker(matrix(rep(class2comp$dist2,each=nclass1), nclass1, nclass2), Identcla3)
  class2m3=kronecker(matrix(rep(class2comp$dist3,each=nclass1), nclass1, nclass2), Identcla3)
  # inner: each row is a copy of y[obsclass2], nclass2 copies
  # stack Identcla3 copies of inner matrix
  class3m1=kronecker(Identcla12,class3comp$dist1)
  class3m2=kronecker(Identcla12,class3comp$dist2)
  class3m3=kronecker(Identcla12,class3comp$dist3)
  # each matrix is nclass1*nclass2 with all elements equal to y[obsclass3][i]
  # stack all matrices
  m1=class1m1+class2m2+class3m3
  m2=class1m1+class2m3+class3m2
  m3=class1m2+class2m1+class3m3
  m4=class1m2+class2m3+class3m1
  m5=class1m3+class2m1+class3m2
  m6=class1m3+class2m2+class3m1
  
  
  Num.IPCW=sum(as.numeric(m1<=m2 & m1<=m3 & m1<=m4 & m1<=m5 & m1<=m6)/probcensorweig)
  
  event1n=rep(1, nclass1)
  event2n=rep(1, nclass2)
  event3n=rep(1, nclass3)
  Demon.IPCW=sum(kronecker(event1n%*%t(event2n),event3n)/probcensorweig)
  
  weightedHUM=Num.IPCW/Demon.IPCW
  return(weightedHUM)
}

vus.true=function(TT, epsilon, cif1, cif2, surv, predicttime){
  class=epsilon
  class[TT>predicttime]=3
  cpmdata=data.frame(cbind(TT, epsilon, class, cif1, cif2, surv))
  cpmdata=cpmdata[order(TT),]
  colnames(cpmdata)=c("TT", "epsilon", "class", "cif1", "cif2", "s")
  
  class1comp=cpmdata[cpmdata$class==1,]
  class2comp=cpmdata[cpmdata$class==2,]
  class3comp=cpmdata[cpmdata$class==3,]
  nclass1=nrow(class1comp)
  nclass2=nrow(class2comp)
  nclass3=nrow(class3comp)
  
  class1comp$dist1=(class1comp$cif1-1)^2+(class1comp$cif2)^2+(class1comp$s)^2
  class1comp$dist2=(class1comp$cif1)^2+(class1comp$cif2-1)^2+(class1comp$s)^2
  class1comp$dist3=(class1comp$cif1)^2+(class1comp$cif2)^2+(class1comp$s-1)^2
  
  class2comp$dist1=(class2comp$cif1-1)^2+(class2comp$cif2)^2+(class2comp$s)^2
  class2comp$dist2=(class2comp$cif1)^2+(class2comp$cif2-1)^2+(class2comp$s)^2
  class2comp$dist3=(class2comp$cif1)^2+(class2comp$cif2)^2+(class2comp$s-1)^2
  
  class3comp$dist1=(class3comp$cif1-1)^2+(class3comp$cif2)^2+(class3comp$s)^2
  class3comp$dist2=(class3comp$cif1)^2+(class3comp$cif2-1)^2+(class3comp$s)^2
  class3comp$dist3=(class3comp$cif1)^2+(class3comp$cif2)^2+(class3comp$s-1)^2
  
  Identcla3=rep(1, nclass3)
  Identcla12=matrix(1, nclass1, nclass2)
  class1m1=kronecker(matrix(rep(class1comp$dist1,nclass2), nclass1, nclass2), Identcla3)
  class1m2=kronecker(matrix(rep(class1comp$dist2,nclass2), nclass1, nclass2), Identcla3)
  class1m3=kronecker(matrix(rep(class1comp$dist3,nclass2), nclass1, nclass2), Identcla3)
  # inner: each column is a copy of y[obsclass1], nclass2 copies.
  # stack Identcla3 copies of the inner matrix
  class2m1=kronecker(matrix(rep(class2comp$dist1,each=nclass1), nclass1, nclass2), Identcla3)
  class2m2=kronecker(matrix(rep(class2comp$dist2,each=nclass1), nclass1, nclass2), Identcla3)
  class2m3=kronecker(matrix(rep(class2comp$dist3,each=nclass1), nclass1, nclass2), Identcla3)
  # inner: each row is a copy of y[obsclass2], nclass2 copies
  # stack Identcla3 copies of inner matrix
  class3m1=kronecker(Identcla12,class3comp$dist1)
  class3m2=kronecker(Identcla12,class3comp$dist2)
  class3m3=kronecker(Identcla12,class3comp$dist3)
  # each matrix is nclass1*nclass2 with all elements equal to y[obsclass3][i]
  # stack all matrices
  m1=class1m1+class2m2+class3m3
  m2=class1m1+class2m3+class3m2
  m3=class1m2+class2m1+class3m3
  m4=class1m2+class2m3+class3m1
  m5=class1m3+class2m1+class3m2
  m6=class1m3+class2m2+class3m1
  
  
  Num.IPCW=sum(as.numeric(m1<=m2 & m1<=m3 & m1<=m4 & m1<=m5 & m1<=m6))
  
  event1n=rep(1, nclass1)
  event2n=rep(1, nclass2)
  event3n=rep(1, nclass3)
  Demon.IPCW=sum(kronecker(event1n%*%t(event2n),event3n))
  
  weightedHUM=Num.IPCW/Demon.IPCW
  
  return(weightedHUM)
}
