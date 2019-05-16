rm(list = ls())
library("survival")
library(MASS)
library(ggplot2)
library(Hmisc)
library(powerSurvEpi)
library(blockrand)
setwd("U:/Biostatistical Services/Private/Internship/StageJuliaLanger/Memoire_texte")
ovarian<- read.csv("Ovariand.csv", header = TRUE)
ovarian$TRIAL<-as.factor(ovarian$TRIAL)
summary(ovarian)

#Discretisation
b<-c(-1,0.9,1.9,8)
label<-c("ECOG 0","ECOG 1","ECOG 2 or worse")
ovarian$PS<-cut(ovarian$PS,breaks=b,labels=label)
summary(ovarian$PS)
PS2<-as.numeric(ovarian$PS)

b<-c(-1,0.9,1.9,8)
label<-c("no residual disease","<= 2 cm","> 2 cm")
ovarian$RESDISEA<-cut(ovarian$RESDISEA,breaks=b,labels=label)
summary(ovarian$RESDISEA)
RESDISEA2<-as.numeric(ovarian$RESDISEA)

b<-c(0,1,2,3)
label<-c("Good","Intermediate","Poor")
ovarian$GRADE<-cut(ovarian$GRADE,breaks=b,labels=label)
summary(ovarian$GRADE)
GRADE2<-as.numeric(ovarian$GRADE)

b<-c(0,2,3,4)
label<-c("Stage I, II","Stage III","Stage IV")
ovarian$STAGE<-cut(ovarian$STAGE,breaks=b,labels=label)
summary(ovarian$STAGE)
STAGE2<-as.numeric(ovarian$STAGE)

#Creation of the variable SCORE
SCORE<-0.48 * RESDISEA2 + 0.37* PS2 + 0.016 * ovarian$AGE
+ 0.14 *GRADE2 + 0.24 *STAGE2
ovarian<-cbind(ovarian,SCORE)
ovarian <- ovarian[order(ovarian$SCORE) , ]

b<-c(0,1.95 , 2.355 , 2.673 ,2.85,10)
label<-c("Best","Good","Intermediate","Bad","Worst")
ovarian$SCORE<-cut(ovarian$SCORE,breaks=b,labels=label)
table(ovarian$SCORE)

#Discretisation of the age
b<-c(17,54,64,100)
label<-c("<55", "55-65",">=65")
ovarian$AGE<-cut(ovarian$AGE,breaks=b,labels=label)
summary(ovarian$AGE)

#Creation of the variable CENTER

set.seed(40)
CENTER<-round(rweibull(1198,1.1,33)/10,digits=0)
ovarian<-cbind(ovarian,CENTER)
table(ovarian$CENTER)

#When sample size = 50
#set.seed(40)
#CENTER<-round(runif(1198,0,4),digits=0)
#ovarian<-cbind(ovarian,CENTER)
#table(ovarian$CENTER)

#Creation randomvariable

set.seed(20)
RANDOM<-round(runif(1198,1,4),digits=0)
ovarian<-cbind(ovarian,RANDOM)
head(ovarian)

#Distribution of the survival times
fitdistr(ovarian$TIME[ovarian$TRT==2 & ovarian$TIME>0], densfun="weibull")

plot(density(ovarian$TIME[ovarian$TRT==1]), main="Density of the survival time of the control group vs Weibull (1,3.48)")
lines(dweibull(seq(0:0.0005:100),1,3.48),col="red")
legend(x=8,y=0.2,legend=c("Control group", "Weibull"),col=c("black", "red"), lty=1:2)

plot(density(ovarian$TIME[ovarian$TRT==2]), main="Density of the survival time of the experimental group vs Weibull (1,4)")
lines(dweibull(seq(0:0.0005:100),1,4),col="red")
legend(x=8,y=0.2,legend=c("Experimental group", "Weibull"),col=c("black", "red"), lty=1:2)

#Kaplan-Meier
library("survival", lib.loc="~/R/win-library/3.3")
surv.ov1 = Surv( time = ovarian$TIME[ovarian$TRT==1] , event = ovarian$STATUS[ovarian$TRT==1])
fit1 = survfit ( surv.ov1 ~ 1)
surv.ov2 = Surv( time = ovarian$TIME[ovarian$TRT==2] , event = ovarian$STATUS[ovarian$TRT==2])
fit2 = survfit ( surv.ov2 ~ 1)

weib<-function(t,lambda){

exp(-t/lambda)

}


plot ( fit1 ,conf.int=FALSE, xlab = " Time " , ylab = "Survival function" , main = "Kaplan - Meier curve ",mark.time=T)
lines(weib(seq(0:0.0005:100),3.48),col="red",lwd=2)
legend(x=8,y=0.9,legend=c("Control group", "Weibull survival function"),col=c("black", "red"), lty=1:2)


plot ( fit2 ,conf.int=FALSE, xlab = " Time " , ylab = "Survival function" , main = "Kaplan - Meier curve ",mark.time=T)
lines(weib(seq(0:0.0005:100),4),col="red",lwd=2)
legend(x=8,y=0.9,legend=c("Experimental group", "Weibull survival function"),col=c("black", "red"), lty=1:2)


#Creation of a time with a not significative HR
HRold<-0.8502
for (j in 1:1198){
  if (ovarian$TRT[j]<1.9) {ovarian$TIME[j]=ovarian$TIME[j]/HRold}
  else {ovarian$TIME[j]=ovarian$TIME[j]}
}



#Creation of the variable TIME2
TIME2<-ovarian$TIME
ovarian<-cbind(ovarian,TIME2)

#Dataset without missing values and with only the prognostic factors
library(tidyr)
ovarian <- ovarian[ , c(1,4,5,6,8,10,14,21,23,24,25,26)]
ovarian<-ovarian%>%drop_na(c(AGE,RESDISEA,GRADE,STAGE,PS))



##################################################################################################
#################################Complete randomisation###########################################
##################################################################################################

N<-500
M<-1000
nsimu<-500
HR<-0.7
ovariansimple<-matrix(nrow=N,ncol=12)
pvalue<-numeric(M)
pvalue1<-numeric(M)
pvalue2<-numeric(M)
pvalue3<-numeric(M)
pvalue4<-numeric(M)
pvalue5<-numeric(M)
pvaluepow<-numeric(M)
pvaluepow1<-numeric(M)
pvaluepow2<-numeric(M)
pvaluepow3<-numeric(M)
pvaluepow4<-numeric(M)
pvaluepow5<-numeric(M)
coeff<-numeric(M)
coeff1<-numeric(M)
coeff2<-numeric(M)
coeff3<-numeric(M)
coeff4<-numeric(M)
coeff5<-numeric(M)


##############################Power of the test################################################
set.seed(5239)

for (i in 1:M){
  TRT<-round(runif(N, min=0, max=1), digits=0)
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  ovariansimple<-cbind(ovariansimple,TRT)
  for (j in 1:N){
    if (ovariansimple$TRT[j]<0.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
  }
  surv.ov2 = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
  fit = survfit ( surv.ov2 ~ 1)
  #Logrank
  surv2<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT,data=ovariansimple)
  pvaluepow[i]<-1 - pchisq(surv2$chisq, 1)
  #For the MSE
  cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT , data = ovariansimple )
  coeff[i]<- exp(cox2$coefficients)
  #Stratified Logrank (1 factor)
  strat21<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
  pvaluepow1[i]<- 1 - pchisq(strat21$chisq, 1)
  #For the MSE
  cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA) , data = ovariansimple )
  coeff1[i]<- exp(cox2$coefficients)
  #Stratified Logrank (3 factors)
  strat22<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
  pvaluepow2[i]<- 1 - pchisq(strat22$chisq, 1)
  #For the MSE
  cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS) , data = ovariansimple )
  coeff2[i]<- exp(cox2$coefficients)
  #Stratified Logrank (5 factors)
  strat23<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
  pvaluepow3[i]<- 1 - pchisq(strat23$chisq, 1)
  #For the MSE
  cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE) , data = ovariansimple )
  coeff3[i]<- exp(cox2$coefficients)
  #Stratified Logrank (Score)
  strat24<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
  pvaluepow4[i]<- 1 - pchisq(strat24$chisq, 1)
  #For the MSE
  cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE) , data = ovariansimple )
  coeff4[i]<- exp(cox2$coefficients)
  #Stratified Logrank (5 factors + Score + Center)
  strat25<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
  pvaluepow5[i]<- 1 - pchisq(strat25$chisq, 1)
  #For the MSE
  cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER) , data = ovariansimple )
  coeff5[i]<- exp(cox2$coefficients)
  
}

#Logrank
for (j in 1:M){
  
  if (pvaluepow[j]<0.05) {pvaluepow[j]=1}
  else {
    pvaluepow[j]=0
  }
  
}

powsr<-mean(pvaluepow)

#Stratified Logrank (1 factor)
for (j in 1:M){
  
  if (pvaluepow1[j]<0.05) {pvaluepow1[j]=1}
  else {
    pvaluepow1[j]=0
  }
  
}

powsr1<-mean(pvaluepow1)

#Stratified Logrank (3 factors) 
for (j in 1:M){
  
  if (pvaluepow2[j]<0.05) {pvaluepow2[j]=1}
  else {
    pvaluepow2[j]=0
  }
  
}

powsr2<-mean(pvaluepow2)

#Stratified Logrank (5 factors)
for (j in 1:M){
  
  if (pvaluepow3[j]<0.05) {pvaluepow3[j]=1}
  else {
    pvaluepow3[j]=0
  }
  
}

powsr3<-mean(pvaluepow3)

#Stratified Logrank (Score)
for (j in 1:M){
  
  if (pvaluepow4[j]<0.05) {pvaluepow4[j]=1}
  else {
    pvaluepow4[j]=0
  }
  
}

powsr4<-mean(pvaluepow4)

#Stratified Logrank (5 factors+ Score + Center)
for (j in 1:M){
  
  if (pvaluepow5[j]<0.05) {pvaluepow5[j]=1}
  else {
    pvaluepow5[j]=0
  }
  
}

powsr5<-mean(pvaluepow5)

##############################Size of the test#################################################
set.seed(5239)
for (i in 1:M){
TRT<-round(runif(N, min=0, max=1), digits=0)
ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
ovariansimple<-cbind(ovariansimple,TRT)
surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
fit = survfit ( surv.ov ~ 1)
#Stratified Logrank (1 factor)
strat1<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
pvalue1[i]<- 1 - pchisq(strat1$chisq, 1)
#Stratified Logrank (3 factors)
strat2<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
pvalue2[i]<- 1 - pchisq(strat2$chisq, 1)
#Stratified Logrank (5 factors)
strat3<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
pvalue3[i]<- 1 - pchisq(strat3$chisq, 1)
#Stratified Logrank (Score)
strat4<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
pvalue4[i]<- 1 - pchisq(strat4$chisq, 1)
#Stratified Logrank (5 factors + Score + Center)
strat5<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
pvalue5[i]<- 1 - pchisq(strat5$chisq, 1)
#Logrank
surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
pvalue[i]<-1 - pchisq(surv$chisq, 1)
}

#Logrank
for (j in 1:M){
  
  if (pvalue[j]<0.05) {pvalue[j]=1}
  else {
    pvalue[j]=0
  }
  
}

stsr<-mean(pvalue)
#Stratified Logrank (1 factor)
for (j in 1:M){
  
  if (pvalue1[j]<0.05) {pvalue1[j]=1}
  else {
    pvalue1[j]=0
  }
  
}

stsr1<-mean(pvalue1)

#Stratified Logrank (3 factors) 
for (j in 1:M){
  
  if (pvalue2[j]<0.05) {pvalue2[j]=1}
  else {
    pvalue2[j]=0
  }
  
}

stsr2<-mean(pvalue2)

#Stratified Logrank (5 factors)
for (j in 1:M){
  
  if (pvalue3[j]<0.05) {pvalue3[j]=1}
  else {
    pvalue3[j]=0
  }
  
}

stsr3<-mean(pvalue3)

#Stratified Logrank (Score)
for (j in 1:M){
  
  if (pvalue4[j]<0.05) {pvalue4[j]=1}
  else {
    pvalue4[j]=0
  }
  
}

stsr4<-mean(pvalue4)

#Stratified Logrank (5 factors+ Score + Center)
for (j in 1:M){
  
  if (pvalue5[j]<0.05) {pvalue5[j]=1}
  else {
    pvalue5[j]=0
  }
  
}

stsr5<-mean(pvalue5)


##############################MSE################################################################

MSE<-numeric(M)
MSE1<-numeric(M)
MSE2<-numeric(M)
MSE3<-numeric(M)
MSE4<-numeric(M)
MSE5<-numeric(M)

#No stratification
for (j in 1:M){
  MSE[j]<-(coeff[j]-HR)^2
}

MSEfinal<-mean(MSE)

#Stratification with 1 factor

for (j in 1:M){
  MSE1[j]<-(coeff1[j]-HR)^2
}

MSEfinal1<-mean(MSE1)

#Stratification with 3 factors

for (j in 1:M){
  MSE2[j]<-(coeff2[j]-HR)^2
}

MSEfinal2<-mean(MSE2)

#Stratification with 5 factors

for (j in 1:M){
  MSE3[j]<-(coeff3[j]-HR)^2
}

MSEfinal3<-mean(MSE3)

#Stratification with the score

for (j in 1:M){
  MSE4[j]<-(coeff4[j]-HR)^2
}

MSEfinal4<-mean(MSE4)

#Stratification with the score + 5 factors + the center

for (j in 1:M){
  MSE5[j]<-(coeff5[j]-HR)^2
}

MSEfinal5<-mean(MSE5)


###############################################################################################
##############################Rerandomization##################################################
###############################################################################################

##############################Size of the test#################################################
#Logrank  

set.seed(4444)

pval.sim=stat.sim=NULL

pvalrerandotest<-rep( 0, len=M)
pvalrerando<-numeric(M)

for (i in 1:M){
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  TRT<-round(runif(N, min=0, max=1), digits=0)
  ovariansimple<-cbind(ovariansimple,TRT)
  surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (k in 1:nsimu){
    randogroup[,1+k]<-round(runif(N, min=0, max=1), digits=0)
  }
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
 
  for (l in 1:nsimu){

      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerando[i]<-a/nsimu 
    }



for (l in 1:M){
  if (pvalrerando[l]<0.05) {pvalrerandotest[l]=1}
  else {pvalrerandotest[l]=0}
}


pvalrerand<-mean(pvalrerandotest)


#Stratified Logrank (1 factor)

set.seed(4444)
pval.sim=stat.sim=NULL

pvalrerando1<-numeric(M)
pvalrerandotest1<-rep( 0, len=M)

for (i in 1:M){
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  TRT<-round(runif(N, min=0, max=1), digits=0)
  ovariansimple<-cbind(ovariansimple,TRT)
  surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (k in 1:nsimu){
    randogroup[,1+k]<-round(runif(N, min=0, max=1), digits=0)
  }
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
  
  for (l in 1:nsimu){
    
    surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
    p<-1 - pchisq(surv2$chisq, 1)
    s<-surv2$chisq
    pval.sim=c(pval.sim,p)
    stat.sim=c(stat.sim,s)
    
  }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerando1[i]<-a/nsimu 
}



for (l in 1:M){
  if (pvalrerando1[l]<0.05) {pvalrerandotest1[l]=1}
  else {pvalrerandotest1[l]=0}
}


pvalrerand1<-mean(pvalrerandotest1)

#Stratified Logrank (3 factors)

set.seed(4444)

pvalrerando2<-numeric(M)
pvalrerandotest2<-rep( 0, len=M)

pval.sim=stat.sim=NULL

for (i in 1:M){
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  TRT<-round(runif(N, min=0, max=1), digits=0)
  ovariansimple<-cbind(ovariansimple,TRT)
  surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (k in 1:nsimu){
    randogroup[,1+k]<-round(runif(N, min=0, max=1), digits=0)
  }
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
  
  for (l in 1:nsimu){
    
    surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
    p<-1 - pchisq(surv2$chisq, 1)
    s<-surv2$chisq
    pval.sim=c(pval.sim,p)
    stat.sim=c(stat.sim,s)
    
  }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerando2[i]<-a/nsimu 
}


for (l in 1:M){
  if (pvalrerando2[l]<0.05) {pvalrerandotest2[l]=1}
  else {pvalrerandotest2[l]=0}
}


pvalrerand2<-mean(pvalrerandotest2)

#Stratified Logrank (5 factors)


set.seed(4444)

pvalrerando3<-numeric(M)
pvalrerandotest3<-rep( 0, len=M)

pval.sim=stat.sim=NULL

for (i in 1:M){
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  TRT<-round(runif(N, min=0, max=1), digits=0)
  ovariansimple<-cbind(ovariansimple,TRT)
  surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (k in 1:nsimu){
    randogroup[,1+k]<-round(runif(N, min=0, max=1), digits=0)
  }
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
  
  for (l in 1:nsimu){
    
    surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE),data=randogroup)
    p<-1 - pchisq(surv2$chisq, 1)
    s<-surv2$chisq
    pval.sim=c(pval.sim,p)
    stat.sim=c(stat.sim,s)
    
  }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerando3[i]<-a/nsimu 
}


for (l in 1:M){
  if (pvalrerando3[l]<0.05) {pvalrerandotest3[l]=1}
  else {pvalrerandotest3[l]=0}
}


pvalrerand3<-mean(pvalrerandotest3)

#Stratified Logrank (Score)



pvalrerando4<-numeric(M)
pvalrerandotest4<-rep( 0, len=M)

pval.sim=stat.sim=NULL

set.seed(4444)
for (i in 1:M){
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  TRT<-round(runif(N, min=0, max=1), digits=0)
  ovariansimple<-cbind(ovariansimple,TRT)
  surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (k in 1:nsimu){
    randogroup[,1+k]<-round(runif(N, min=0, max=1), digits=0)
  }
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
  
  for (l in 1:nsimu){
    
    surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
    p<-1 - pchisq(surv2$chisq, 1)
    s<-surv2$chisq
    pval.sim=c(pval.sim,p)
    stat.sim=c(stat.sim,s)
    
  }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerando4[i]<-a/nsimu 
}


for (l in 1:M){
  if (pvalrerando4[l]<0.05) {pvalrerandotest4[l]=1}
  else {pvalrerandotest4[l]=0}
}


pvalrerand4<-mean(pvalrerandotest4)

#Stratified Logrank (5 factors + Score + Center)

pvalrerando5<-numeric(M)
pvalrerandotest5<-rep( 0, len=M)

pval.sim=stat.sim=NULL
set.seed(4444)
for (i in 1:M){
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  TRT<-round(runif(N, min=0, max=1), digits=0)
  ovariansimple<-cbind(ovariansimple,TRT)
  surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (k in 1:nsimu){
    randogroup[,1+k]<-round(runif(N, min=0, max=1), digits=0)
  }
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
  
  for (l in 1:nsimu){
    
    surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE)+strata(randogroup$SCORE)+strata(randogroup$CENTER),data=randogroup)
    p<-1 - pchisq(surv2$chisq, 1)
    s<-surv2$chisq
    pval.sim=c(pval.sim,p)
    stat.sim=c(stat.sim,s)
    
  }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerando5[i]<-a/nsimu 
}



for (l in 1:M){
  if (pvalrerando5[l]<0.05) {pvalrerandotest5[l]=1}
  else {pvalrerandotest5[l]=0}
}


pvalrerand5<-mean(pvalrerandotest5)

##############################Power of the test################################################

#Logrank  




pvalrerandopow<-numeric(M)
pvalrerando<-numeric(M)
pvalrerandotestpow<-numeric(M)

pval.sim=stat.sim=NULL
set.seed(4452)

for (i in 1:M){
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  TRT<-round(runif(N, min=0, max=1), digits=0)
  ovariansimple<-cbind(ovariansimple,TRT)
  for (j in 1:N){
    if (ovariansimple$TRT[j]<0.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
  }
  surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (m in 1:nsimu){
    randogroup[,1+m]<-round(runif(N, min=0, max=1), digits=0)
  }
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  for (j in 1:length(randogroup[,1])){
    for (k in 1:nsimu){
      if (randogroup[j,k+1]<0.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
    else {randogroup$TIME2[j]=randogroup$TIME[j]}
    }
  }
  surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
  
  for (l in 1:nsimu){

    survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
    p<-1 - pchisq(survpow2$chisq, 1)
    s<-survpow2$chisq
    pval.sim=c(pval.sim,p)
    stat.sim=c(stat.sim,s)
    
  }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerandopow[i]<-a/nsimu 
}


for (l in 1:M){
  if (pvalrerandopow[l]<0.05) {pvalrerandotestpow[l]=1}
  else {pvalrerandotestpow[l]=0}
}




pvalrerandlogpow<-mean(pvalrerandotestpow)


#Stratified Logrank (1 factor)

pvalrerandopow1<-numeric(M)
pvalrerandotestpow1<-numeric(M)

pval.sim=stat.sim=NULL
set.seed(4452)
for (i in 1:M){
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  TRT<-round(runif(N, min=0, max=1), digits=0)
  ovariansimple<-cbind(ovariansimple,TRT)
  for (j in 1:N){
    if (ovariansimple$TRT[j]<0.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
  }
  surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (m in 1:nsimu){
    randogroup[,1+m]<-round(runif(N, min=0, max=1), digits=0)
  }
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  for (j in 1:length(randogroup[,1])){
    for (k in 1:nsimu){
      if (randogroup[j,k+1]<0.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
      else {randogroup$TIME2[j]=randogroup$TIME[j]}
    }
  }
  surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
  
  for (l in 1:nsimu){
    
    survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
    p<-1 - pchisq(survpow2$chisq, 1)
    s<-survpow2$chisq
    pval.sim=c(pval.sim,p)
    stat.sim=c(stat.sim,s)
    
  }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerandopow1[i]<-a/nsimu 
}


for (l in 1:M){
  if (pvalrerandopow1[l]<0.05) {pvalrerandotestpow1[l]=1}
  else {pvalrerandotestpow1[l]=0}
}


pvalrerandlogpow1<-mean(pvalrerandotestpow1)



#Stratified Logrank (3 factors)


pvalrerandopow2<-numeric(M)
pvalrerandotestpow2<-numeric(M)


pval.sim=stat.sim=NULL
set.seed(4452)
for (i in 1:M){
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  TRT<-round(runif(N, min=0, max=1), digits=0)
  ovariansimple<-cbind(ovariansimple,TRT)
  for (j in 1:N){
    if (ovariansimple$TRT[j]<0.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
  }
  surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (m in 1:nsimu){
    randogroup[,1+m]<-round(runif(N, min=0, max=1), digits=0)
  }
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  for (j in 1:length(randogroup[,1])){
    for (k in 1:nsimu){
      if (randogroup[j,k+1]<0.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
      else {randogroup$TIME2[j]=randogroup$TIME[j]}
    }
  }
  surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
  
  for (l in 1:nsimu){
    
    survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
    p<-1 - pchisq(survpow2$chisq, 1)
    s<-survpow2$chisq
    pval.sim=c(pval.sim,p)
    stat.sim=c(stat.sim,s)
    
  }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerandopow2[i]<-a/nsimu 
}



for (l in 1:M){
  if (pvalrerandopow2[l]<0.05) {pvalrerandotestpow2[l]=1}
  else {pvalrerandotestpow2[l]=0}
}


pvalrerandlogpow2<-mean(pvalrerandotestpow2)


#Stratified Logrank (5 factors)


pvalrerandopow3<-numeric(M)
pvalrerandotestpow3<-numeric(M)


pval.sim=stat.sim=NULL
set.seed(4452)
for (i in 1:M){
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  TRT<-round(runif(N, min=0, max=1), digits=0)
  ovariansimple<-cbind(ovariansimple,TRT)
  for (j in 1:N){
    if (ovariansimple$TRT[j]<0.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
  }
  surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (m in 1:nsimu){
    randogroup[,1+m]<-round(runif(N, min=0, max=1), digits=0)
  }
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  for (j in 1:length(randogroup[,1])){
    for (k in 1:nsimu){
      if (randogroup[j,k+1]<0.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
      else {randogroup$TIME2[j]=randogroup$TIME[j]}
    }
  }
  surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
  
  for (l in 1:nsimu){
    
    survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE),data=randogroup)
    p<-1 - pchisq(survpow2$chisq, 1)
    s<-survpow2$chisq
    pval.sim=c(pval.sim,p)
    stat.sim=c(stat.sim,s)
    
  }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerandopow3[i]<-a/nsimu 
}



for (l in 1:M){
  if (pvalrerandopow3[l]<0.05) {pvalrerandotestpow3[l]=1}
  else {pvalrerandotestpow3[l]=0}
}


pvalrerandlogpow3<-mean(pvalrerandotestpow3)

#Stratified Logrank (Score)


pvalrerandopow4<-numeric(M)
pvalrerandotestpow4<-numeric(M)

pval.sim=stat.sim=NULL
set.seed(4452)
for (i in 1:M){
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  TRT<-round(runif(N, min=0, max=1), digits=0)
  ovariansimple<-cbind(ovariansimple,TRT)
  for (j in 1:N){
    if (ovariansimple$TRT[j]<0.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
  }
  surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (m in 1:nsimu){
    randogroup[,1+m]<-round(runif(N, min=0, max=1), digits=0)
  }
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  for (j in 1:length(randogroup[,1])){
    for (k in 1:nsimu){
      if (randogroup[j,k+1]<0.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
    else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
  }
  surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
  
  for (l in 1:nsimu){
    
    survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
    p<-1 - pchisq(survpow2$chisq, 1)
    s<-survpow2$chisq
    pval.sim=c(pval.sim,p)
    stat.sim=c(stat.sim,s)
    
  }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerandopow4[i]<-a/nsimu 
}



for (l in 1:M){
  if (pvalrerandopow4[l]<0.05) {pvalrerandotestpow4[l]=1}
  else {pvalrerandotestpow4[l]=0}
}


pvalrerandlogpow4<-mean(pvalrerandotestpow4)

#Stratified Logrank (Score + 5 factors + Center)


pvalrerandopow5<-numeric(M)
pvalrerandotestpow5<-numeric(M)

pval.sim=stat.sim=NULL

for (i in 1:M){
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  TRT<-round(runif(N, min=0, max=1), digits=0)
  ovariansimple<-cbind(ovariansimple,TRT)
  for (j in 1:N){
    if (ovariansimple$TRT[j]<0.9 ) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
  }
  surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE)+strata(ovariansimple$CENTER)+strata(ovariansimple$SCORE),data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (m in 1:nsimu){
    randogroup[,1+m]<-round(runif(N, min=0, max=1), digits=0)
  }
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  for (j in 1:length(randogroup[,1])){
    for (k in 1:nsimu){
      if (randogroup[j,k+1]<0.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
     else {randogroup$TIME2[j]=randogroup$TIME[j]}
    }
  }
  surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
  
  for (l in 1:nsimu){
    
    survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$CENTER)+strata(randogroup$SCORE),data=randogroup)
    p<-1 - pchisq(survpow2$chisq, 1)
    s<-survpow2$chisq
    pval.sim=c(pval.sim,p)
    stat.sim=c(stat.sim,s)
    
  }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerandopow5[i]<-a/nsimu 
}


for (l in 1:M){
  if (pvalrerandopow5[l]<0.05) {pvalrerandotestpow5[l]=1}
  else {pvalrerandotestpow5[l]=0}
}


pvalrerandlogpow5<-mean(pvalrerandotestpow5)


##################################################################################################
#################################Permuted blocks with the center##################################
##################################################################################################

#################################Power############################################################

set.seed(5239)
ovariansimple<-matrix(nrow=N,ncol=13)
pvalue<-numeric(M)
pvalue1<-numeric(M)
pvalue2<-numeric(M)
pvalue3<-numeric(M)
pvalue4<-numeric(M)
pvalue5<-numeric(M)
pvaluepow<-numeric(M)
pvaluepow1<-numeric(M)
pvaluepow2<-numeric(M)
pvaluepow3<-numeric(M)
pvaluepow4<-numeric(M)
pvaluepow5<-numeric(M)
coeff<-numeric(M)
coeff1<-numeric(M)
coeff2<-numeric(M)
coeff3<-numeric(M)
coeff4<-numeric(M)
coeff5<-numeric(M)

Factor1<-c(0:16,20)
Factor2<-c("no residual disease","<= 2 cm","> 2 cm")
Factor3<-c("<55","55-65",">=65")
Factor4<-c("ECOG 0", "ECOG 1", "ECOG 2 or worse")
Factor5<-c("Good","Intermediate","Poor")
Factor6<-c("Stage III", "Stage IV", "Stage I, II")
Grid<-Factor1
nbrrand<-numeric(length(Grid))
nbrrandsimple<-numeric(length(Grid))

for(m in 1:length(Grid)){
  nbrrand[m]<-nrow(ovarian[ovarian$CENTER==Grid[m] , ])
  
}

set.seed(5239)
for (i in 1:M){
  Rando1fact1<-NULL
  TRTCENTER<-NULL
  
  for (m in 1:length(Grid)){
    Rando1fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
    
    Rando1fact<-Rando1fact[1:nbrrand[m],]
    
    Rando1fact1<-rbind(Rando1fact1,Rando1fact)
  }
  
  
  ovarian <- ovarian[order(ovarian$CENTER) , ]
  TRTCENTER<-na.omit(Rando1fact1)
  TRT<-TRTCENTER$treatment
  TRT<-as.numeric(TRT)
  
  ovarian<-cbind(ovarian[,1:12],TRT)
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  for (j in 1:N){
    if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
  }
  surv.ov2 = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
  fit = survfit ( surv.ov2 ~ 1)
  #Logrank
  surv2<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT,data=ovariansimple)
  pvaluepow[i]<-1 - pchisq(surv2$chisq, 1)
  #For the MSE
  cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT , data = ovariansimple )
  coeff[i]<- exp(cox2$coefficients)
  #Stratified Logrank (1 factor)
  strat21<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
  pvaluepow1[i]<- 1 - pchisq(strat21$chisq, 1)
  #For the MSE
  cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA) , data = ovariansimple )
  coeff1[i]<- exp(cox2$coefficients)
  #Stratified Logrank (3 factors)
  strat22<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
  pvaluepow2[i]<- 1 - pchisq(strat22$chisq, 1)
  #For the MSE
  cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS) , data = ovariansimple )
  coeff2[i]<- exp(cox2$coefficients)
  #Stratified Logrank (5 factors)
  strat23<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
  pvaluepow3[i]<- 1 - pchisq(strat23$chisq, 1)
  #For the MSE
  cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE) , data = ovariansimple )
  coeff3[i]<- exp(cox2$coefficients)
  #Stratified Logrank (Score)
  strat24<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
  pvaluepow4[i]<- 1 - pchisq(strat24$chisq, 1)
  #For the MSE
  cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE) , data = ovariansimple )
  coeff4[i]<- exp(cox2$coefficients)
  #Stratified Logrank (5 factors + Score + Center)
  strat25<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
  pvaluepow5[i]<- 1 - pchisq(strat25$chisq, 1)
  #For the MSE
  cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER) , data = ovariansimple )
  coeff5[i]<- exp(cox2$coefficients)
  
}


#Logrank
for (j in 1:M){
  
  if (pvaluepow[j]<0.05) {pvaluepow[j]=1}
  else {
    pvaluepow[j]=0
  }
  
}

pbpow<-mean(pvaluepow)

#Stratified Logrank (1 factor)
for (j in 1:M){
  
  if (pvaluepow1[j]<0.05) {pvaluepow1[j]=1}
  else {
    pvaluepow1[j]=0
  }
  
}

pbpow1<-mean(pvaluepow1)

#Stratified Logrank (3 factors) 
for (j in 1:M){
  
  if (pvaluepow2[j]<0.05) {pvaluepow2[j]=1}
  else {
    pvaluepow2[j]=0
  }
  
}

pbpow2<-mean(pvaluepow2)

#Stratified Logrank (5 factors)
for (j in 1:M){
  
  if (pvaluepow3[j]<0.05) {pvaluepow3[j]=1}
  else {
    pvaluepow3[j]=0
  }
  
}

pbpow3<-mean(pvaluepow3)

#Stratified Logrank (Score)
for (j in 1:M){
  
  if (pvaluepow4[j]<0.05) {pvaluepow4[j]=1}
  else {
    pvaluepow4[j]=0
  }
  
}

pbpow4<-mean(pvaluepow4)

#Stratified Logrank (5 factors+ Score + Center)
for (j in 1:M){
  
  if (pvaluepow5[j]<0.05) {pvaluepow5[j]=1}
  else {
    pvaluepow5[j]=0
  }
  
}

pbpow5<-mean(pvaluepow5)

##############################Size of the test#################################################

set.seed(5239)

for (i in 1:M){
  
  Rando1fact1<-NULL
  TRTCENTER<-NULL
  
  for (m in 1:length(Grid)){
    Rando1fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
    
    Rando1fact<-Rando1fact[1:nbrrand[m],]
    
    Rando1fact1<-rbind(Rando1fact1,Rando1fact)
  }
  
  
  ovarian <- ovarian[order(ovarian$CENTER) , ]
  TRTCENTER<-na.omit(Rando1fact1)
  TRT<-TRTCENTER$treatment
  TRT<-as.numeric(TRT)
  
  ovarian<-cbind(ovarian[,1:12],TRT)
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
  fit = survfit ( surv.ov ~ 1)
  #Stratified Logrank (1 factor)
  strat1<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
  pvalue1[i]<- 1 - pchisq(strat1$chisq, 1)
  #Stratified Logrank (3 factors)
  strat2<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
  pvalue2[i]<- 1 - pchisq(strat2$chisq, 1)
  #Stratified Logrank (5 factors)
  strat3<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
  pvalue3[i]<- 1 - pchisq(strat3$chisq, 1)
  #Stratified Logrank (Score)
  strat4<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
  pvalue4[i]<- 1 - pchisq(strat4$chisq, 1)
  #Stratified Logrank (5 factors + Score + Center)
  strat5<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
  pvalue5[i]<- 1 - pchisq(strat5$chisq, 1)
  #Logrank
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
  pvalue[i]<-1 - pchisq(surv$chisq, 1)
}

#Logrank
for (j in 1:M){
  
  if (pvalue[j]<0.05) {pvalue[j]=1}
  else {
    pvalue[j]=0
  }
  
}

pbs<-mean(pvalue)
#Stratified Logrank (1 factor)
for (j in 1:M){
  
  if (pvalue1[j]<0.05) {pvalue1[j]=1}
  else {
    pvalue1[j]=0
  }
  
}

pbs1<-mean(pvalue1)

#Stratified Logrank (3 factors) 
for (j in 1:M){
  
  if (pvalue2[j]<0.05) {pvalue2[j]=1}
  else {
    pvalue2[j]=0
  }
  
}

pbs2<-mean(pvalue2)

#Stratified Logrank (5 factors)
for (j in 1:M){
  
  if (pvalue3[j]<0.05) {pvalue3[j]=1}
  else {
    pvalue3[j]=0
  }
  
}

pbs3<-mean(pvalue3)

#Stratified Logrank (Score)
for (j in 1:M){
  
  if (pvalue4[j]<0.05) {pvalue4[j]=1}
  else {
    pvalue4[j]=0
  }
  
}

pbs4<-mean(pvalue4)

#Stratified Logrank (5 factors+ Score + Center)
for (j in 1:M){
  
  if (pvalue5[j]<0.05) {pvalue5[j]=1}
  else {
    pvalue5[j]=0
  }
  
}

pbs5<-mean(pvalue5)


##############################MSE################################################################

MSE<-numeric(M)
MSE1<-numeric(M)
MSE2<-numeric(M)
MSE3<-numeric(M)
MSE4<-numeric(M)
MSE5<-numeric(M)

#No stratification
for (j in 1:M){
  MSE[j]<-(coeff[j]-HR)^2
}

MSEfinalpb<-mean(MSE)

#Stratification with 1 factor

for (j in 1:M){
  MSE1[j]<-(coeff1[j]-HR)^2
}

MSEfinalpb1<-mean(MSE1)

#Stratification with 3 factors

for (j in 1:M){
  MSE2[j]<-(coeff2[j]-HR)^2
}

MSEfinalpb2<-mean(MSE2)

#Stratification with 5 factors

for (j in 1:M){
  MSE3[j]<-(coeff3[j]-HR)^2
}

MSEfinalpb3<-mean(MSE3)

#Stratification with the score

for (j in 1:M){
  MSE4[j]<-(coeff4[j]-HR)^2
}

MSEfinalpb4<-mean(MSE4)

#Stratification with the score + 5 factors + the center

for (j in 1:M){
  MSE5[j]<-(coeff5[j]-HR)^2
}

MSEfinalpb5<-mean(MSE5)

###############################################################################################
##############################Rerandomization##################################################
###############################################################################################

##############################Size of the test#################################################
#Logrank  

set.seed(4444)

pvalrerandotest<-rep( 0, len=M)
pvalrerando<-numeric(M)

pval.sim=stat.sim=NULL

for (i in 1:M){
  Rando1fact1<-NULL
  TRTCENTER<-NULL
  
  for (m in 1:length(Grid)){
    Rando1fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
    
    Rando1fact<-Rando1fact[1:nbrrand[m],]
    
    Rando1fact1<-rbind(Rando1fact1,Rando1fact)
  }
  
  
  ovarian <- ovarian[order(ovarian$CENTER) , ]
  TRTCENTER<-na.omit(Rando1fact1)
  TRT<-TRTCENTER$treatment
  TRT<-as.numeric(TRT)
  
  ovarian<-cbind(ovarian[,1:12],TRT)
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  ovariansimple <- ovariansimple[order(ovariansimple$CENTER) , ]
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (k in 1:nsimu){
    Rando1fact1<-NULL
    TRTCENTER<-NULL
    for(m in 1:length(Grid)){
      nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid[m] , ])
      
    }
    for (m in 1:length(Grid)){
      Rando1fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando1fact<-Rando1fact[1:nbrrandsimple[m],]
      
      Rando1fact1<-rbind(Rando1fact1,Rando1fact)
    }
    
    TRTCENTER<-na.omit(Rando1fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    randogroup[,1+k]<-TRT
  }
  
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
  
  for (l in 1:nsimu){
    
    surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
    p<-1 - pchisq(surv2$chisq, 1)
    s<-surv2$chisq
    pval.sim=c(pval.sim,p)
    stat.sim=c(stat.sim,s)
    
  }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerando[i]<-a/nsimu 
}



for (l in 1:M){
  if (pvalrerando[l]<0.05) {pvalrerandotest[l]=1}
  else {pvalrerandotest[l]=0}
}


pvalrerandpb<-mean(pvalrerandotest)
pvalrerandpb


#Stratified Logrank (1 factor)

set.seed(4444)

pvalrerando1<-numeric(M)
pvalrerandotest1<-rep( 0, len=M)


pval.sim=stat.sim=NULL

for (i in 1:M){
  Rando1fact1<-NULL
  TRTCENTER<-NULL
  
  for (m in 1:length(Grid)){
    Rando1fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
    
    Rando1fact<-Rando1fact[1:nbrrand[m],]
    
    Rando1fact1<-rbind(Rando1fact1,Rando1fact)
  }
  
  
  ovarian <- ovarian[order(ovarian$CENTER) , ]
  TRTCENTER<-na.omit(Rando1fact1)
  TRT<-TRTCENTER$treatment
  TRT<-as.numeric(TRT)
  
  ovarian<-cbind(ovarian[,1:12],TRT)
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  ovariansimple <- ovariansimple[order(ovariansimple$CENTER) , ]
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (k in 1:nsimu){
    Rando1fact1<-NULL
    TRTCENTER<-NULL
    for(m in 1:length(Grid)){
      nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid[m] , ])
      
    }
    for (m in 1:length(Grid)){
      Rando1fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando1fact<-Rando1fact[1:nbrrandsimple[m],]
      
      Rando1fact1<-rbind(Rando1fact1,Rando1fact)
    }
    
    TRTCENTER<-na.omit(Rando1fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    randogroup[,1+k]<-TRT
  }
  
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
  
  for (l in 1:nsimu){
    
    surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
    p<-1 - pchisq(surv2$chisq, 1)
    s<-surv2$chisq
    pval.sim=c(pval.sim,p)
    stat.sim=c(stat.sim,s)
    
  }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerando1[i]<-a/nsimu 
}



for (l in 1:M){
  if (pvalrerando1[l]<0.05) {pvalrerandotest1[l]=1}
  else {pvalrerandotest1[l]=0}
}


pvalrerandpb1<-mean(pvalrerandotest1)
pvalrerandpb1

#Stratified Logrank (3 factors)

set.seed(4444)

pvalrerando2<-numeric(M)
pvalrerandotest2<-rep( 0, len=M)


pval.sim=stat.sim=NULL

for (i in 1:M){
  Rando1fact1<-NULL
  TRTCENTER<-NULL
  
  for (m in 1:length(Grid)){
    Rando1fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
    
    Rando1fact<-Rando1fact[1:nbrrand[m],]
    
    Rando1fact1<-rbind(Rando1fact1,Rando1fact)
  }
  
  
  ovarian <- ovarian[order(ovarian$CENTER) , ]
  TRTCENTER<-na.omit(Rando1fact1)
  TRT<-TRTCENTER$treatment
  TRT<-as.numeric(TRT)
  
  ovarian<-cbind(ovarian[,1:12],TRT)
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  ovariansimple <- ovariansimple[order(ovariansimple$CENTER) , ]
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
  for (k in 1:nsimu){
    Rando1fact1<-NULL
    TRTCENTER<-NULL
    for(m in 1:length(Grid)){
      nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid[m] , ])
      
    }
    for (m in 1:length(Grid)){
      Rando1fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando1fact<-Rando1fact[1:nbrrandsimple[m],]
      
      Rando1fact1<-rbind(Rando1fact1,Rando1fact)
    }
    
    TRTCENTER<-na.omit(Rando1fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    randogroup[,1+k]<-TRT
  }
  
  colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
  colnames(ovariansimple)[1]="subject.id"
  randogroup=merge(randogroup,ovariansimple,by="subject.id")
  surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
  
  for (l in 1:nsimu){
    
    surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
    p<-1 - pchisq(surv2$chisq, 1)
    s<-surv2$chisq
    pval.sim=c(pval.sim,p)
    stat.sim=c(stat.sim,s)
    
  }
  sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
  a<-0
  for (j in 1 : nsimu)
  {
    if (sim[j,4]>sim[j,2])
    {
      a<-a+1
    }
  }
  pvalrerando2[i]<-a/nsimu 
}



for (l in 1:M){
  if (pvalrerando2[l]<0.05) {pvalrerandotest2[l]=1}
  else {pvalrerandotest2[l]=0}
}


pvalrerandpb2<-mean(pvalrerandotest2)
pvalrerandpb2

#Stratified Logrank (5 factors)


set.seed(4444)

pvalrerando3<-numeric(M)
pvalrerandotest3<-rep( 0, len=M)


pval.sim=stat.sim=NULL

for (i in 1:M){
  Rando1fact1<-NULL
  TRTCENTER<-NULL
  
  for (m in 1:length(Grid)){
    Rando1fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
    
    Rando1fact<-Rando1fact[1:nbrrand[m],]
    
    Rando1fact1<-rbind(Rando1fact1,Rando1fact)
  }
  
  
  ovarian <- ovarian[order(ovarian$CENTER) , ]
  TRTCENTER<-na.omit(Rando1fact1)
  TRT<-TRTCENTER$treatment
  TRT<-as.numeric(TRT)
  
  ovarian<-cbind(ovarian[,1:12],TRT)
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
  surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE),data=ovariansimple)
  p.obs<-1 - pchisq(surv$chisq, 1)
  s.obs<-surv$chisq
  ovariansimple<-ovariansimple[,1:12]
  ovariansimple <- ovariansimple[order(ovariansimple$CENTER) , ]
  randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando1fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid)){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid[m] , ])
        
      }
      for (m in 1:length(Grid)){
        Rando1fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando1fact<-Rando1fact[1:nbrrandsimple[m],]
        
        Rando1fact1<-rbind(Rando1fact1,Rando1fact)
      }
      
      TRTCENTER<-na.omit(Rando1fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando3[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando3[l]<0.05) {pvalrerandotest3[l]=1}
    else {pvalrerandotest3[l]=0}
  }
  
  
  pvalrerandpb3<-mean(pvalrerandotest3)
  pvalrerandpb3
  
  #Stratified Logrank (Score)
  
  set.seed(4444)
  
  pvalrerando4<-numeric(M)
  pvalrerandotest4<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando1fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid)){
      Rando1fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando1fact<-Rando1fact[1:nbrrand[m],]
      
      Rando1fact1<-rbind(Rando1fact1,Rando1fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER) , ]
    TRTCENTER<-na.omit(Rando1fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando1fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid)){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid[m] , ])
        
      }
      for (m in 1:length(Grid)){
        Rando1fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando1fact<-Rando1fact[1:nbrrandsimple[m],]
        
        Rando1fact1<-rbind(Rando1fact1,Rando1fact)
      }
      
      TRTCENTER<-na.omit(Rando1fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando4[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando4[l]<0.05) {pvalrerandotest4[l]=1}
    else {pvalrerandotest4[l]=0}
  }
  
  
  pvalrerandpb4<-mean(pvalrerandotest4)
  pvalrerandpb4
  #Stratified Logrank (5 factors + Score + Center)
  
  set.seed(4444)
  
  pvalrerando5<-numeric(M)
  pvalrerandotest5<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando1fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid)){
      Rando1fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando1fact<-Rando1fact[1:nbrrand[m],]
      
      Rando1fact1<-rbind(Rando1fact1,Rando1fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER) , ]
    TRTCENTER<-na.omit(Rando1fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE)+strata(ovariansimple$CENTER)+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando1fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid)){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid[m] , ])
        
      }
      for (m in 1:length(Grid)){
        Rando1fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando1fact<-Rando1fact[1:nbrrandsimple[m],]
        
        Rando1fact1<-rbind(Rando1fact1,Rando1fact)
      }
      
      TRTCENTER<-na.omit(Rando1fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$CENTER)+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando5[l]<0.05) {pvalrerandotest5[l]=1}
    else {pvalrerandotest5[l]=0}
  }
  
  
  pvalrerandpb5<-mean(pvalrerandotest5)

  
  ##############################Power of the test################################################
  
  #Logrank  
  
  set.seed(4452)
  
  pvalrerandopow<-numeric(M)
  pvalrerandotestpow<-numeric(M)

  pval.sim=stat.sim=NULL
 set.seed(4452)

  
  for (i in 1:M){
    Rando1fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid)){
      Rando1fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando1fact<-Rando1fact[1:nbrrand[m],]
      
      Rando1fact1<-rbind(Rando1fact1,Rando1fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER) , ]
    TRTCENTER<-na.omit(Rando1fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando1fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid)){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid[m] , ])
        
      }
      for (m in 1:length(Grid)){
        Rando1fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando1fact<-Rando1fact[1:nbrrandsimple[m],]
        
        Rando1fact1<-rbind(Rando1fact1,Rando1fact)
      }
      
      TRTCENTER<-na.omit(Rando1fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow[l]<0.05) {pvalrerandotestpow[l]=1}
    else {pvalrerandotestpow[l]=0}
  }
  
  
  pvalrerandpbpow<-mean(pvalrerandotestpow)
  pvalrerandpbpow

  
  #Stratified Logrank (1 factor)
  
  pvalrerandopow1<-numeric(M)
  pvalrerandotestpow1<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
 set.seed(4452)
  for (i in 1:M){
    Rando1fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid)){
      Rando1fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando1fact<-Rando1fact[1:nbrrand[m],]
      
      Rando1fact1<-rbind(Rando1fact1,Rando1fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER) , ]
    TRTCENTER<-na.omit(Rando1fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando1fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid)){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid[m] , ])
        
      }
      for (m in 1:length(Grid)){
        Rando1fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando1fact<-Rando1fact[1:nbrrandsimple[m],]
        
        Rando1fact1<-rbind(Rando1fact1,Rando1fact)
      }
      
      TRTCENTER<-na.omit(Rando1fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow1[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow1[l]<0.05) {pvalrerandotestpow1[l]=1}
    else {pvalrerandotestpow1[l]=0}
  }
  
  
  pvalrerandpbpow1<-mean(pvalrerandotestpow1)
  pvalrerandpbpow1

  
  #Stratified Logrank (3 factors)
  
  pvalrerandopow2<-numeric(M)
  pvalrerandotestpow2<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
   set.seed(4452)
  for (i in 1:M){
    Rando1fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid)){
      Rando1fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando1fact<-Rando1fact[1:nbrrand[m],]
      
      Rando1fact1<-rbind(Rando1fact1,Rando1fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER) , ]
    TRTCENTER<-na.omit(Rando1fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando1fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid)){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid[m] , ])
        
      }
      for (m in 1:length(Grid)){
        Rando1fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando1fact<-Rando1fact[1:nbrrandsimple[m],]
        
        Rando1fact1<-rbind(Rando1fact1,Rando1fact)
      }
      
      TRTCENTER<-na.omit(Rando1fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
     p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow2[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow2[l]<0.05) {pvalrerandotestpow2[l]=1}
    else {pvalrerandotestpow2[l]=0}
  }
  
  
  pvalrerandpbpow2<-mean(pvalrerandotestpow2)
  pvalrerandpbpow2

  #Stratified Logrank (5 factors)
  
  pvalrerandopow3<-numeric(M)
  pvalrerandotestpow3<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
   set.seed(4452)
  for (i in 1:M){
    Rando1fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid)){
      Rando1fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando1fact<-Rando1fact[1:nbrrand[m],]
      
      Rando1fact1<-rbind(Rando1fact1,Rando1fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER) , ]
    TRTCENTER<-na.omit(Rando1fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando1fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid)){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid[m] , ])
        
      }
      for (m in 1:length(Grid)){
        Rando1fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando1fact<-Rando1fact[1:nbrrandsimple[m],]
        
        Rando1fact1<-rbind(Rando1fact1,Rando1fact)
      }
      
      TRTCENTER<-na.omit(Rando1fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow3[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow3[l]<0.05) {pvalrerandotestpow3[l]=1}
    else {pvalrerandotestpow3[l]=0}
  }
  
  
  pvalrerandpbpow3<-mean(pvalrerandotestpow3)
  pvalrerandpbpow3
  
#Stratified Logrank (Score)
  set.seed(4444)
  pvalrerandopow4<-numeric(M)
  pvalrerandotestpow4<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando1fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid)){
      Rando1fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando1fact<-Rando1fact[1:nbrrand[m],]
      
      Rando1fact1<-rbind(Rando1fact1,Rando1fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER) , ]
    TRTCENTER<-na.omit(Rando1fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9 ) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando1fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid)){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid[m] , ])
        
      }
      for (m in 1:length(Grid)){
        Rando1fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando1fact<-Rando1fact[1:nbrrandsimple[m],]
        
        Rando1fact1<-rbind(Rando1fact1,Rando1fact)
      }
      
      TRTCENTER<-na.omit(Rando1fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9 ) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
     p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow4[i]<-a/nsimu 
  }
  
  for (l in 1:M){
    if (pvalrerandopow4[l]<0.05) {pvalrerandotestpow4[l]=1}
    else {pvalrerandotestpow4[l]=0}
  }
  
  
  pvalrerandpbpow4<-mean(pvalrerandotestpow4)
  pvalrerandpbpow4
  
    
  #Stratified Logrank (Score + 5 factors + Center)
  
  pvalrerandopow5<-numeric(M)
  pvalrerandotestpow5<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
   set.seed(4444)
  for (i in 1:M){
    Rando1fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid)){
      Rando1fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando1fact<-Rando1fact[1:nbrrand[m],]
      
      Rando1fact1<-rbind(Rando1fact1,Rando1fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER) , ]
    TRTCENTER<-na.omit(Rando1fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando1fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid)){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid[m] , ])
        
      }
      for (m in 1:length(Grid)){
        Rando1fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando1fact<-Rando1fact[1:nbrrandsimple[m],]
        
        Rando1fact1<-rbind(Rando1fact1,Rando1fact)
      }
      
      TRTCENTER<-na.omit(Rando1fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
       else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$SCORE)+strata(randogroup$CENTER),data=randogroup)
     p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow5[l]<0.05) {pvalrerandotestpow5[l]=1}
    else {pvalrerandotestpow5[l]=0}
  }
  
  
  pvalrerandpbpow5<-mean(pvalrerandotestpow5)
  
  
  ##################################################################################################
  #################################Permuted blocks with the center + 1 factor#######################
  ##################################################################################################
  
  #################################Power############################################################
  
  set.seed(5239)
  ovariansimple<-matrix(nrow=N,ncol=13)
  pvalue<-numeric(M)
  pvalue1<-numeric(M)
  pvalue2<-numeric(M)
  pvalue3<-numeric(M)
  pvalue4<-numeric(M)
  pvalue5<-numeric(M)
  pvaluepow<-numeric(M)
  pvaluepow1<-numeric(M)
  pvaluepow2<-numeric(M)
  pvaluepow3<-numeric(M)
  pvaluepow4<-numeric(M)
  pvaluepow5<-numeric(M)
  coeff<-numeric(M)
  coeff1<-numeric(M)
  coeff2<-numeric(M)
  coeff3<-numeric(M)
  coeff4<-numeric(M)
  coeff5<-numeric(M)
  Factor1<-c(0:16,20)
  Factor2<-c("no residual disease","<= 2 cm","> 2 cm")
  Grid<-expand.grid(Factor1,Factor2)
  Grid <- Grid[order(Grid$Var1,Grid$Var2) , ]
  nbrrand<-numeric(length(Grid[,1]))
  nbrrandsimple<-numeric(length(Grid[,1]))
  
  for(m in 1:length(Grid[,1])){
    nbrrand[m]<-nrow(ovarian[ovarian$CENTER==Grid$Var1[m] & ovarian$RESDISEA==Grid$Var2[m] , ])
    
  }
  
  for (i in 1:M){
    Rando2fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando2fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando2fact<-Rando2fact[1:nbrrand[m],]
      
      Rando2fact1<-rbind(Rando2fact1,Rando2fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA) , ]
    TRTCENTER<-na.omit(Rando2fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov2 = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov2 ~ 1)
    #Logrank
    surv2<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT,data=ovariansimple)
    pvaluepow[i]<-1 - pchisq(surv2$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT , data = ovariansimple )
    coeff[i]<- exp(cox2$coefficients)
    #Stratified Logrank (1 factor)
    strat21<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvaluepow1[i]<- 1 - pchisq(strat21$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA) , data = ovariansimple )
    coeff1[i]<- exp(cox2$coefficients)
    #Stratified Logrank (3 factors)
    strat22<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvaluepow2[i]<- 1 - pchisq(strat22$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS) , data = ovariansimple )
    coeff2[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors)
    strat23<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvaluepow3[i]<- 1 - pchisq(strat23$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE) , data = ovariansimple )
    coeff3[i]<- exp(cox2$coefficients)
    #Stratified Logrank (Score)
    strat24<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvaluepow4[i]<- 1 - pchisq(strat24$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE) , data = ovariansimple )
    coeff4[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors + Score + Center)
    strat25<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvaluepow5[i]<- 1 - pchisq(strat25$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER) , data = ovariansimple )
    coeff5[i]<- exp(cox2$coefficients)
    
  }
  
  
  #Logrank
  for (j in 1:M){
    
    if (pvaluepow[j]<0.05) {pvaluepow[j]=1}
    else {
      pvaluepow[j]=0
    }
    
  }
  
  pb1pow<-mean(pvaluepow)
  
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvaluepow1[j]<0.05) {pvaluepow1[j]=1}
    else {
      pvaluepow1[j]=0
    }
    
  }
  
  pb1pow1<-mean(pvaluepow1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvaluepow2[j]<0.05) {pvaluepow2[j]=1}
    else {
      pvaluepow2[j]=0
    }
    
  }
  
  pb1pow2<-mean(pvaluepow2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvaluepow3[j]<0.05) {pvaluepow3[j]=1}
    else {
      pvaluepow3[j]=0
    }
    
  }
  
  pb1pow3<-mean(pvaluepow3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvaluepow4[j]<0.05) {pvaluepow4[j]=1}
    else {
      pvaluepow4[j]=0
    }
    
  }
  
  pb1pow4<-mean(pvaluepow4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvaluepow5[j]<0.05) {pvaluepow5[j]=1}
    else {
      pvaluepow5[j]=0
    }
    
  }
  
  pb1pow5<-mean(pvaluepow5)

  
  ##############################Size of the test#################################################
  
  set.seed(5239)
  
  for (i in 1:M){
    Rando2fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando2fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando2fact<-Rando2fact[1:nbrrand[m],]
      
      Rando2fact1<-rbind(Rando2fact1,Rando2fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA) , ]
    TRTCENTER<-na.omit(Rando2fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov ~ 1)
    #Stratified Logrank (1 factor)
    strat1<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvalue1[i]<- 1 - pchisq(strat1$chisq, 1)
    #Stratified Logrank (3 factors)
    strat2<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvalue2[i]<- 1 - pchisq(strat2$chisq, 1)
    #Stratified Logrank (5 factors)
    strat3<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvalue3[i]<- 1 - pchisq(strat3$chisq, 1)
    #Stratified Logrank (Score)
    strat4<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvalue4[i]<- 1 - pchisq(strat4$chisq, 1)
    #Stratified Logrank (5 factors + Score + Center)
    strat5<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvalue5[i]<- 1 - pchisq(strat5$chisq, 1)
    #Logrank
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    pvalue[i]<-1 - pchisq(surv$chisq, 1)
  }
  
  #Logrank
  for (j in 1:M){
    
    if (pvalue[j]<0.05) {pvalue[j]=1}
    else {
      pvalue[j]=0
    }
    
  }
  
  pb1s<-mean(pvalue)
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvalue1[j]<0.05) {pvalue1[j]=1}
    else {
      pvalue1[j]=0
    }
    
  }
  
  pb1s1<-mean(pvalue1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvalue2[j]<0.05) {pvalue2[j]=1}
    else {
      pvalue2[j]=0
    }
    
  }
  
  pb1s2<-mean(pvalue2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvalue3[j]<0.05) {pvalue3[j]=1}
    else {
      pvalue3[j]=0
    }
    
  }
  
  pb1s3<-mean(pvalue3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvalue4[j]<0.05) {pvalue4[j]=1}
    else {
      pvalue4[j]=0
    }
    
  }
  
  pb1s4<-mean(pvalue4)
  
  #Stratified Logrank (Score + 5 factors + Center)
  for (j in 1:M){
    
    if (pvalue5[j]<0.05) {pvalue5[j]=1}
    else {
      pvalue5[j]=0
    }
    
  }
  
  pb1s5<-mean(pvalue5)
  
  test<-c(pb1s,pb1s1,pb1s2,pb1s3,pb1s4,pb1s5)
  
  barplot(test)
  
  ##############################MSE################################################################
  
  MSE<-numeric(M)
  MSE1<-numeric(M)
  MSE2<-numeric(M)
  MSE3<-numeric(M)
  MSE4<-numeric(M)
  MSE5<-numeric(M)
  
  #No stratification
  for (j in 1:M){
    MSE[j]<-(coeff[j]-HR)^2
  }
  
  MSEfinalpb10<-mean(MSE)
  
  #Stratification with 1 factor
  
  for (j in 1:M){
    MSE1[j]<-(coeff1[j]-HR)^2
  }
  
  MSEfinalpb11<-mean(MSE1)
  
  #Stratification with 3 factors
  
  for (j in 1:M){
    MSE2[j]<-(coeff2[j]-HR)^2
  }
  
  MSEfinalpb12<-mean(MSE2)
  
  #Stratification with 5 factors
  
  for (j in 1:M){
    MSE3[j]<-(coeff3[j]-HR)^2
  }
  
  MSEfinalpb13<-mean(MSE3)
  
  #Stratification with the score
  
  for (j in 1:M){
    MSE4[j]<-(coeff4[j]-HR)^2
  }
  
  MSEfinalpb14<-mean(MSE4)
  
  #Stratification with the score + 5 factors + the center
  
  for (j in 1:M){
    MSE5[j]<-(coeff5[j]-HR)^2
  }
  
  MSEfinalpb15<-mean(MSE5)

  ###############################################################################################
  ##############################Rerandomization##################################################
  ###############################################################################################
  
  ##############################Size of the test#################################################
  #Logrank  
  
  set.seed(4444)
  
  pvalrerandotest<-rep( 0, len=M)
  pvalrerando<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando2fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando2fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando2fact<-Rando2fact[1:nbrrand[m],]
      
      Rando2fact1<-rbind(Rando2fact1,Rando2fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA) , ]
    TRTCENTER<-na.omit(Rando2fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando2fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando2fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando2fact<-Rando2fact[1:nbrrandsimple[m],]
        
        Rando2fact1<-rbind(Rando2fact1,Rando2fact)
      }
      
      TRTCENTER<-na.omit(Rando2fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando[l]<0.05) {pvalrerandotest[l]=1}
    else {pvalrerandotest[l]=0}
  }
  
  
  pvalrerandpb10<-mean(pvalrerandotest)
  
  
  #Stratified Logrank (1 factor)
  
  set.seed(4444)
  
  pvalrerando1<-numeric(M)
  pvalrerandotest1<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando2fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando2fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando2fact<-Rando2fact[1:nbrrand[m],]
      
      Rando2fact1<-rbind(Rando2fact1,Rando2fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA) , ]
    TRTCENTER<-na.omit(Rando2fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando2fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando2fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando2fact<-Rando2fact[1:nbrrandsimple[m],]
        
        Rando2fact1<-rbind(Rando2fact1,Rando2fact)
      }
      
      TRTCENTER<-na.omit(Rando2fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando1[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando1[l]<0.05) {pvalrerandotest1[l]=1}
    else {pvalrerandotest1[l]=0}
  }
  
  
  pvalrerandpb11<-mean(pvalrerandotest1)
  
  #Stratified Logrank (3 factors)
  
  set.seed(4444)
  
  pvalrerando2<-numeric(M)
  pvalrerandotest2<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando2fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando2fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando2fact<-Rando2fact[1:nbrrand[m],]
      
      Rando2fact1<-rbind(Rando2fact1,Rando2fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA) , ]
    TRTCENTER<-na.omit(Rando2fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando2fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando2fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando2fact<-Rando2fact[1:nbrrandsimple[m],]
        
        Rando2fact1<-rbind(Rando2fact1,Rando2fact)
      }
      
      TRTCENTER<-na.omit(Rando2fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando2[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando2[l]<0.05) {pvalrerandotest2[l]=1}
    else {pvalrerandotest2[l]=0}
  }
  
  
  pvalrerandpb12<-mean(pvalrerandotest2)
  
  #Stratified Logrank (5 factors)
  
  
  set.seed(4444)
  
  pvalrerando3<-numeric(M)
  pvalrerandotest3<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando2fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando2fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando2fact<-Rando2fact[1:nbrrand[m],]
      
      Rando2fact1<-rbind(Rando2fact1,Rando2fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA) , ]
    TRTCENTER<-na.omit(Rando2fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando2fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando2fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando2fact<-Rando2fact[1:nbrrandsimple[m],]
        
        Rando2fact1<-rbind(Rando2fact1,Rando2fact)
      }
      
      TRTCENTER<-na.omit(Rando2fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando3[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando3[l]<0.05) {pvalrerandotest3[l]=1}
    else {pvalrerandotest3[l]=0}
  }
  
  
  pvalrerandpb13<-mean(pvalrerandotest3)
  
  #Stratified Logrank (Score)
  
  set.seed(4444)
  
  pvalrerando4<-numeric(M)
  pvalrerandotest4<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando2fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando2fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando2fact<-Rando2fact[1:nbrrand[m],]
      
      Rando2fact1<-rbind(Rando2fact1,Rando2fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA) , ]
    TRTCENTER<-na.omit(Rando2fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando2fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando2fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando2fact<-Rando2fact[1:nbrrandsimple[m],]
        
        Rando2fact1<-rbind(Rando2fact1,Rando2fact)
      }
      
      TRTCENTER<-na.omit(Rando2fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando4[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando4[l]<0.05) {pvalrerandotest4[l]=1}
    else {pvalrerandotest4[l]=0}
  }
  
  
  pvalrerandpb14<-mean(pvalrerandotest4)
  
  #Stratified Logrank (5 factors + Score + Center)
  
  set.seed(4444)
  
  pvalrerando5<-numeric(M)
  pvalrerandotest5<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando2fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando2fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando2fact<-Rando2fact[1:nbrrand[m],]
      
      Rando2fact1<-rbind(Rando2fact1,Rando2fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA) , ]
    TRTCENTER<-na.omit(Rando2fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE)+strata(ovariansimple$CENTER)+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando2fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando2fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando2fact<-Rando2fact[1:nbrrandsimple[m],]
        
        Rando2fact1<-rbind(Rando2fact1,Rando2fact)
      }
      
      TRTCENTER<-na.omit(Rando2fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$CENTER)+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando5[l]<0.05) {pvalrerandotest5[l]=1}
    else {pvalrerandotest5[l]=0}
  }
  
  
  pvalrerandpb15<-mean(pvalrerandotest5)
  
  ##############################Power of the test################################################
  
  #Logrank  
  
  set.seed(4452)
  
  pvalrerandopow<-numeric(M)
  pvalrerandotestpow<-numeric(M)
 
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando2fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando2fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando2fact<-Rando2fact[1:nbrrand[m],]
      
      Rando2fact1<-rbind(Rando2fact1,Rando2fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA) , ]
    TRTCENTER<-na.omit(Rando2fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando2fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando2fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando2fact<-Rando2fact[1:nbrrandsimple[m],]
        
        Rando2fact1<-rbind(Rando2fact1,Rando2fact)
      }
      
      TRTCENTER<-na.omit(Rando2fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow[l]<0.05) {pvalrerandotestpow[l]=1}
    else {pvalrerandotestpow[l]=0}
  }
  
  
  pvalrerandpb1pow<-mean(pvalrerandotestpow)
  pvalrerandpb1pow
  
 
  
  #Stratified Logrank (1 factor)
  
  pvalrerandopow1<-numeric(M)
  pvalrerandotestpow1<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
    set.seed(4452)
  for (i in 1:M){
    Rando2fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando2fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando2fact<-Rando2fact[1:nbrrand[m],]
      
      Rando2fact1<-rbind(Rando2fact1,Rando2fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA) , ]
    TRTCENTER<-na.omit(Rando2fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando2fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando2fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando2fact<-Rando2fact[1:nbrrandsimple[m],]
        
        Rando2fact1<-rbind(Rando2fact1,Rando2fact)
      }
      
      TRTCENTER<-na.omit(Rando2fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow1[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow1[l]<0.05) {pvalrerandotestpow1[l]=1}
    else {pvalrerandotestpow1[l]=0}
  }
  
  
  pvalrerandpb1pow1<-mean(pvalrerandotestpow1)
  pvalrerandpb1pow1
  

  #Stratified Logrank (3 factors)
  
  pvalrerandopow2<-numeric(M)
  pvalrerandotestpow2<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
    set.seed(4452)
  for (i in 1:M){
    Rando2fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando2fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando2fact<-Rando2fact[1:nbrrand[m],]
      
      Rando2fact1<-rbind(Rando2fact1,Rando2fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA) , ]
    TRTCENTER<-na.omit(Rando2fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando2fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando2fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando2fact<-Rando2fact[1:nbrrandsimple[m],]
        
        Rando2fact1<-rbind(Rando2fact1,Rando2fact)
      }
      
      TRTCENTER<-na.omit(Rando2fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow2[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow2[l]<0.05) {pvalrerandotestpow2[l]=1}
    else {pvalrerandotestpow2[l]=0}
  }
  
  
  pvalrerandpb1pow2<-mean(pvalrerandotestpow2)
  pvalrerandpb1pow2

  
  #Stratified Logrank (5 factors)
  
  pvalrerandopow3<-numeric(M)
  pvalrerandotestpow3<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
    set.seed(4452)
  for (i in 1:M){
    Rando2fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando2fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando2fact<-Rando2fact[1:nbrrand[m],]
      
      Rando2fact1<-rbind(Rando2fact1,Rando2fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA) , ]
    TRTCENTER<-na.omit(Rando2fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando2fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando2fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando2fact<-Rando2fact[1:nbrrandsimple[m],]
        
        Rando2fact1<-rbind(Rando2fact1,Rando2fact)
      }
      
      TRTCENTER<-na.omit(Rando2fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow3[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow3[l]<0.05) {pvalrerandotestpow3[l]=1}
    else {pvalrerandotestpow3[l]=0}
  }
  
  
  pvalrerandpb1pow3<-mean(pvalrerandotestpow3)
  pvalrerandpb1pow3
  
  #Stratified Logrank (Score)
  
  pvalrerandopow4<-numeric(M)
  pvalrerandotestpow4<-numeric(M)
  
  pval.sim=stat.sim=NULL
    set.seed(4452)
  for (i in 1:M){
    Rando2fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando2fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando2fact<-Rando2fact[1:nbrrand[m],]
      
      Rando2fact1<-rbind(Rando2fact1,Rando2fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA) , ]
    TRTCENTER<-na.omit(Rando2fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9 ) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando2fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando2fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando2fact<-Rando2fact[1:nbrrandsimple[m],]
        
        Rando2fact1<-rbind(Rando2fact1,Rando2fact)
      }
      
      TRTCENTER<-na.omit(Rando2fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9 ) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
       p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow4[i]<-a/nsimu 
  }
  
  for (l in 1:M){
    if (pvalrerandopow4[l]<0.05) {pvalrerandotestpow4[l]=1}
    else {pvalrerandotestpow4[l]=0}
  }
  
  
  pvalrerandpb1pow4<-mean(pvalrerandotestpow4)
  pvalrerandpb1pow4
  
  #Stratified Logrank (Score + 5 factors + Center)
  
  pvalrerandopow5<-numeric(M)
  pvalrerandotestpow5<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
    set.seed(4452)
  for (i in 1:M){
    
    Rando2fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando2fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando2fact<-Rando2fact[1:nbrrand[m],]
      
      Rando2fact1<-rbind(Rando2fact1,Rando2fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA) , ]
    TRTCENTER<-na.omit(Rando2fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando2fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando2fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando2fact<-Rando2fact[1:nbrrandsimple[m],]
        
        Rando2fact1<-rbind(Rando2fact1,Rando2fact)
      }
      
      TRTCENTER<-na.omit(Rando2fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$SCORE)+strata(randogroup$CENTER),data=randogroup)
     p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow5[l]<0.05) {pvalrerandotestpow5[l]=1}
    else {pvalrerandotestpow5[l]=0}
  }
  
  
  pvalrerandpb1pow5<-mean(pvalrerandotestpow5)
  
  ##################################################################################################
  #################################Permuted blocks with the center + 3 factors######################
  ##################################################################################################
  set.seed(5239)
  ovariansimple<-matrix(nrow=N,ncol=13)
  pvalue<-numeric(M)
  pvalue1<-numeric(M)
  pvalue2<-numeric(M)
  pvalue3<-numeric(M)
  pvalue4<-numeric(M)
  pvalue5<-numeric(M)
  pvaluepow<-numeric(M)
  pvaluepow1<-numeric(M)
  pvaluepow2<-numeric(M)
  pvaluepow3<-numeric(M)
  pvaluepow4<-numeric(M)
  pvaluepow5<-numeric(M)
  coeff<-numeric(M)
  coeff1<-numeric(M)
  coeff2<-numeric(M)
  coeff3<-numeric(M)
  coeff4<-numeric(M)
  coeff5<-numeric(M)
  Factor1<-c(0:16,20)
  Factor2<-c("no residual disease","<= 2 cm","> 2 cm")
  Factor3<-c("<55","55-65",">=65")
  Factor4<-c("ECOG 0", "ECOG 1", "ECOG 2 or worse")
  Grid<-expand.grid(Factor1,Factor2,Factor3,Factor4)
  Grid <- Grid[order(Grid$Var1,Grid$Var2,Grid$Var3,Grid$Var4) , ]
  nbrrand<-numeric(length(Grid[,1]))
  nbrrandsimple<-numeric(length(Grid[,1]))

##########################################Distribution in the strata##########################################################

Griddist<-expand.grid(Factor2,Factor3,Factor4)
Griddist <- Griddist[order(Griddist$Var1,Griddist$Var2,Griddist$Var3) , ]
nbrrandsimp<-numeric(length(Griddist[,1]))
distr<-matrix(nrow=length(Griddist[,1]), ncol=M)

##Distribution in the strata
for (i in 1:M){
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
for(m in 1:length(Griddist[,1])){
  nbrrandsimp[m]<-nrow(ovariansimple[ovariansimple$RESDISEA==Griddist$Var1[m] & ovariansimple$AGE==Griddist$Var2[m] 
                           & ovariansimple$PS==Griddist$Var3[m] , ])
  
}
  distr[,i]<-as.vector(nbrrandsimp)
}

mean_distr<-rowMeans(distr)
mean_distr<-sort(mean_distr,decreasing = TRUE)
plot(mean_distr,type="h",main="Distribution of 500 patients strata from 3 factors",xlab="Number of strata",ylab="Number of patients",lwd=3)
#################################Power###############################################
  
  for(m in 1:length(Grid[,1])){
    nbrrand[m]<-nrow(ovarian[ovarian$CENTER==Grid$Var1[m] & ovarian$RESDISEA==Grid$Var2[m] & ovarian$AGE==Grid$Var3[m] 
                             & ovarian$PS==Grid$Var4[m] , ])
    
  }
  
    set.seed(5239)
  for (i in 1:M){
    Rando3fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando3fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando3fact<-Rando3fact[1:nbrrand[m],]
      
      Rando3fact1<-rbind(Rando3fact1,Rando3fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS) , ]
    TRTCENTER<-na.omit(Rando3fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov2 = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov2 ~ 1)
    #Logrank
    surv2<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT,data=ovariansimple)
    pvaluepow[i]<-1 - pchisq(surv2$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT , data = ovariansimple )
    coeff[i]<- exp(cox2$coefficients)
    #Stratified Logrank (1 factor)
    strat21<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvaluepow1[i]<- 1 - pchisq(strat21$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA) , data = ovariansimple )
    coeff1[i]<- exp(cox2$coefficients)
    #Stratified Logrank (3 factors)
    strat22<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvaluepow2[i]<- 1 - pchisq(strat22$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS) , data = ovariansimple )
    coeff2[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors)
    strat23<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvaluepow3[i]<- 1 - pchisq(strat23$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE) , data = ovariansimple )
    coeff3[i]<- exp(cox2$coefficients)
    #Stratified Logrank (Score)
    strat24<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvaluepow4[i]<- 1 - pchisq(strat24$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE) , data = ovariansimple )
    coeff4[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors + Score + Center)
    strat25<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvaluepow5[i]<- 1 - pchisq(strat25$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER) , data = ovariansimple )
    coeff5[i]<- exp(cox2$coefficients)
    
  }
  
  
  #Logrank
  for (j in 1:M){
    
    if (pvaluepow[j]<0.05) {pvaluepow[j]=1}
    else {
      pvaluepow[j]=0
    }
    
  }
  
  pb3pow<-mean(pvaluepow)
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvaluepow1[j]<0.05) {pvaluepow1[j]=1}
    else {
      pvaluepow1[j]=0
    }
    
  }
  
  pb3pow1<-mean(pvaluepow1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvaluepow2[j]<0.05) {pvaluepow2[j]=1}
    else {
      pvaluepow2[j]=0
    }
    
  }
  
  pb3pow2<-mean(pvaluepow2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvaluepow3[j]<0.05) {pvaluepow3[j]=1}
    else {
      pvaluepow3[j]=0
    }
    
  }
  
  pb3pow3<-mean(pvaluepow3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvaluepow4[j]<0.05) {pvaluepow4[j]=1}
    else {
      pvaluepow4[j]=0
    }
    
  }
  
  pb3pow4<-mean(pvaluepow4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvaluepow5[j]<0.05) {pvaluepow5[j]=1}
    else {
      pvaluepow5[j]=0
    }
    
  }
  
  pb3pow5<-mean(pvaluepow5)
  
  ##############################Size of the test#################################################
  
  set.seed(5239)
  
  for (i in 1:M){
    Rando3fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando3fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando3fact<-Rando3fact[1:nbrrand[m],]
      
      Rando3fact1<-rbind(Rando3fact1,Rando3fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS) , ]
    TRTCENTER<-na.omit(Rando3fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov ~ 1)
    #Stratified Logrank (1 factor)
    strat1<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvalue1[i]<- 1 - pchisq(strat1$chisq, 1)
    #Stratified Logrank (3 factors)
    strat2<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvalue2[i]<- 1 - pchisq(strat2$chisq, 1)
    #Stratified Logrank (5 factors)
    strat3<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvalue3[i]<- 1 - pchisq(strat3$chisq, 1)
    #Stratified Logrank (Score)
    strat4<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvalue4[i]<- 1 - pchisq(strat4$chisq, 1)
    #Stratified Logrank (5 factors + Score + Center)
    strat5<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvalue5[i]<- 1 - pchisq(strat5$chisq, 1)
    #Logrank
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    pvalue[i]<-1 - pchisq(surv$chisq, 1)
  }
  
  #Logrank
  for (j in 1:M){
    
    if (pvalue[j]<0.05) {pvalue[j]=1}
    else {
      pvalue[j]=0
    }
    
  }
  
  pb3s<-mean(pvalue)
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvalue1[j]<0.05) {pvalue1[j]=1}
    else {
      pvalue1[j]=0
    }
    
  }
  
  pb3s1<-mean(pvalue1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvalue2[j]<0.05) {pvalue2[j]=1}
    else {
      pvalue2[j]=0
    }
    
  }
  
  pb3s2<-mean(pvalue2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvalue3[j]<0.05) {pvalue3[j]=1}
    else {
      pvalue3[j]=0
    }
    
  }
  
  pb3s3<-mean(pvalue3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvalue4[j]<0.05) {pvalue4[j]=1}
    else {
      pvalue4[j]=0
    }
    
  }
  
  pb3s4<-mean(pvalue4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvalue5[j]<0.05) {pvalue5[j]=1}
    else {
      pvalue5[j]=0
    }
    
  }
  
  pb3s5<-mean(pvalue5)
  ##############################MSE################################################################
  
  MSE<-numeric(M)
  MSE1<-numeric(M)
  MSE2<-numeric(M)
  MSE3<-numeric(M)
  MSE4<-numeric(M)
  MSE5<-numeric(M)
  
  #No stratification
  for (j in 1:M){
    MSE[j]<-(coeff[j]-HR)^2
  }
  
  MSEfinalpb30<-mean(MSE)
  
  #Stratification with 1 factor
  
  for (j in 1:M){
    MSE1[j]<-(coeff1[j]-HR)^2
  }
  
  MSEfinalpb31<-mean(MSE1)
  
  #Stratification with 3 factors
  
  for (j in 1:M){
    MSE2[j]<-(coeff2[j]-HR)^2
  }
  
  MSEfinalpb32<-mean(MSE2)
  
  #Stratification with 5 factors
  
  for (j in 1:M){
    MSE3[j]<-(coeff3[j]-HR)^2
  }
  
  MSEfinalpb33<-mean(MSE3)
  
  #Stratification with the score
  
  for (j in 1:M){
    MSE4[j]<-(coeff4[j]-HR)^2
  }
  
  MSEfinalpb34<-mean(MSE4)
  
  #Stratification with the score + 5 factors + the center
  
  for (j in 1:M){
    MSE5[j]<-(coeff5[j]-HR)^2
  }
  
  MSEfinalpb35<-mean(MSE5)
  
  MSEvector<-c(MSEfinalpb30,MSEfinalpb31,MSEfinalpb32,MSEfinalpb33,MSEfinalpb34,abs(MSEfinalpb35))
  
  barplot(MSEvector)
  
  ###############################################################################################
  ##############################Rerandomization##################################################
  ###############################################################################################
  
  ##############################Size of the test#################################################
  #Logrank  
  
  set.seed(4444)
  
  pvalrerandotest<-rep( 0, len=M)
  pvalrerando<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando3fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando3fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando3fact<-Rando3fact[1:nbrrand[m],]
      
      Rando3fact1<-rbind(Rando3fact1,Rando3fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS) , ]
    TRTCENTER<-na.omit(Rando3fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando3fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando3fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando3fact<-Rando3fact[1:nbrrandsimple[m],]
        
        Rando3fact1<-rbind(Rando3fact1,Rando3fact)
      }
      
      TRTCENTER<-na.omit(Rando3fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando[l]<0.05) {pvalrerandotest[l]=1}
    else {pvalrerandotest[l]=0}
  }
  
  
  pvalrerandpb30<-mean(pvalrerandotest)
  pvalrerandpb30
  
  
  #Stratified Logrank (1 factor)
  
  set.seed(4444)
  
  pvalrerando1<-numeric(M)
  pvalrerandotest1<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando3fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando3fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando3fact<-Rando3fact[1:nbrrand[m],]
      
      Rando3fact1<-rbind(Rando3fact1,Rando3fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS) , ]
    TRTCENTER<-na.omit(Rando3fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando3fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando3fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando3fact<-Rando3fact[1:nbrrandsimple[m],]
        
        Rando3fact1<-rbind(Rando3fact1,Rando3fact)
      }
      
      TRTCENTER<-na.omit(Rando3fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando1[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando1[l]<0.05) {pvalrerandotest1[l]=1}
    else {pvalrerandotest1[l]=0}
  }
  
  
  pvalrerandpb31<-mean(pvalrerandotest1)
  pvalrerandpb31
  
  #Stratified Logrank (3 factors)
  
  set.seed(4444)
  
  pvalrerando2<-numeric(M)
  pvalrerandotest2<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando3fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando3fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando3fact<-Rando3fact[1:nbrrand[m],]
      
      Rando3fact1<-rbind(Rando3fact1,Rando3fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS) , ]
    TRTCENTER<-na.omit(Rando3fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando3fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando3fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando3fact<-Rando3fact[1:nbrrandsimple[m],]
        
        Rando3fact1<-rbind(Rando3fact1,Rando3fact)
      }
      
      TRTCENTER<-na.omit(Rando3fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando2[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando2[l]<0.05) {pvalrerandotest2[l]=1}
    else {pvalrerandotest2[l]=0}
  }
  
  
  pvalrerandpb32<-mean(pvalrerandotest2)
  pvalrerandpb32
  
  #Stratified Logrank (5 factors)
  
  
  set.seed(4444)
  
  pvalrerando3<-numeric(M)
  pvalrerandotest3<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando3fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando3fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando3fact<-Rando3fact[1:nbrrand[m],]
      
      Rando3fact1<-rbind(Rando3fact1,Rando3fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS) , ]
    TRTCENTER<-na.omit(Rando3fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando3fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando3fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando3fact<-Rando3fact[1:nbrrandsimple[m],]
        
        Rando3fact1<-rbind(Rando3fact1,Rando3fact)
      }
      
      TRTCENTER<-na.omit(Rando3fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando3[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando3[l]<0.05) {pvalrerandotest3[l]=1}
    else {pvalrerandotest3[l]=0}
  }
  
  
  pvalrerandopb33<-mean(pvalrerandotest3)
  
  #Stratified Logrank (Score)
  
  set.seed(4444)
  
  pvalrerando4<-numeric(M)
  pvalrerandotest4<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando3fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando3fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando3fact<-Rando3fact[1:nbrrand[m],]
      
      Rando3fact1<-rbind(Rando3fact1,Rando3fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS) , ]
    TRTCENTER<-na.omit(Rando3fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando3fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando3fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando3fact<-Rando3fact[1:nbrrandsimple[m],]
        
        Rando3fact1<-rbind(Rando3fact1,Rando3fact)
      }
      
      TRTCENTER<-na.omit(Rando3fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando4[i]<-a/nsimu 
  }
  
  
 
  for (l in 1:M){
    if (pvalrerando4[l]<0.05) {pvalrerandotest4[l]=1}
    else {pvalrerandotest4[l]=0}
  }
  
  
  pvalrerandpb34<-mean(pvalrerandotest4)
  pvalrerandpb34
  
  #Stratified Logrank (5 factors + Score + Center)
  
  set.seed(4444)
  
  pvalrerando5<-numeric(M)
  pvalrerandotest5<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando3fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando3fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando3fact<-Rando3fact[1:nbrrand[m],]
      
      Rando3fact1<-rbind(Rando3fact1,Rando3fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS) , ]
    TRTCENTER<-na.omit(Rando3fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE)+strata(ovariansimple$CENTER)+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando3fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando3fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando3fact<-Rando3fact[1:nbrrandsimple[m],]
        
        Rando3fact1<-rbind(Rando3fact1,Rando3fact)
      }
      
      TRTCENTER<-na.omit(Rando3fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$CENTER)+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando5[l]<0.05) {pvalrerandotest5[l]=1}
    else {pvalrerandotest5[l]=0}
  }
  
  
  pvalrerandpb35<-mean(pvalrerandotest5)

  ##############################Power of the test################################################
  
  #Logrank  
  
  set.seed(4452)
  
  pvalrerandopow<-numeric(M)
  pvalrerandotestpow<-numeric(M)

  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando3fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando3fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando3fact<-Rando3fact[1:nbrrand[m],]
      
      Rando3fact1<-rbind(Rando3fact1,Rando3fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS) , ]
    TRTCENTER<-na.omit(Rando3fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando3fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando3fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando3fact<-Rando3fact[1:nbrrandsimple[m],]
        
        Rando3fact1<-rbind(Rando3fact1,Rando3fact)
      }
      
      TRTCENTER<-na.omit(Rando3fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow[l]<0.05) {pvalrerandotestpow[l]=1}
    else {pvalrerandotestpow[l]=0}
  }
  
  
  pvalrerandpb3pow<-mean(pvalrerandotestpow)
  pvalrerandpb3pow
  
  
  
  #Stratified Logrank (1 factor)
  
  pvalrerandopow1<-numeric(M)
  pvalrerandotestpow1<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
    set.seed(4452)
  for (i in 1:M){
    Rando3fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando3fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando3fact<-Rando3fact[1:nbrrand[m],]
      
      Rando3fact1<-rbind(Rando3fact1,Rando3fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS) , ]
    TRTCENTER<-na.omit(Rando3fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando3fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando3fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando3fact<-Rando3fact[1:nbrrandsimple[m],]
        
        Rando3fact1<-rbind(Rando3fact1,Rando3fact)
      }
      
      TRTCENTER<-na.omit(Rando3fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow1[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow1[l]<0.05) {pvalrerandotestpow1[l]=1}
    else {pvalrerandotestpow1[l]=0}
  }
  
  
  pvalrerandpb3pow1<-mean(pvalrerandotestpow1)
  pvalrerandpb3pow1
  
  
  #Stratified Logrank (3 factors)
  
  pvalrerandopow2<-numeric(M)
  pvalrerandotestpow2<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  set.seed(4444)
  
  for (i in 1:M){
    Rando3fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando3fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando3fact<-Rando3fact[1:nbrrand[m],]
      
      Rando3fact1<-rbind(Rando3fact1,Rando3fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS) , ]
    TRTCENTER<-na.omit(Rando3fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando3fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando3fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando3fact<-Rando3fact[1:nbrrandsimple[m],]
        
        Rando3fact1<-rbind(Rando3fact1,Rando3fact)
      }
      
      TRTCENTER<-na.omit(Rando3fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow2[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow2[l]<0.05) {pvalrerandotestpow2[l]=1}
    else {pvalrerandotestpow2[l]=0}
  }
  
  
  pvalrerandpb3pow2<-mean(pvalrerandotestpow2)
  pvalrerandpb3pow2
  
  #Stratified Logrank (5 factors)
  
  pvalrerandopow3<-numeric(M)
  pvalrerandotestpow3<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
    set.seed(4452)
  for (i in 1:M){
    Rando3fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando3fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando3fact<-Rando3fact[1:nbrrand[m],]
      
      Rando3fact1<-rbind(Rando3fact1,Rando3fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS) , ]
    TRTCENTER<-na.omit(Rando3fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando3fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando3fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando3fact<-Rando3fact[1:nbrrandsimple[m],]
        
        Rando3fact1<-rbind(Rando3fact1,Rando3fact)
      }
      
      TRTCENTER<-na.omit(Rando3fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow3[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow3[l]<0.05) {pvalrerandotestpow3[l]=1}
    else {pvalrerandotestpow3[l]=0}
  }
  
  
  pvalrerandpb3pow3<-mean(pvalrerandotestpow3)
  
  
  
  
  
  #Stratified Logrank (Score + 5 factors + Center)
  
  pvalrerandopow5<-numeric(M)
  pvalrerandotestpow5<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
    set.seed(4452)
  for (i in 1:M){
    
    Rando3fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando3fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando3fact<-Rando3fact[1:nbrrand[m],]
      
      Rando3fact1<-rbind(Rando3fact1,Rando3fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS) , ]
    TRTCENTER<-na.omit(Rando3fact1)
    
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
      }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando3fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando3fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando3fact<-Rando3fact[1:nbrrandsimple[m],]
        
        Rando3fact1<-rbind(Rando3fact1,Rando3fact)
      }
      
      TRTCENTER<-na.omit(Rando3fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
       else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$SCORE)+strata(randogroup$CENTER),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow5[l]<0.05) {pvalrerandotestpow5[l]=1}
    else {pvalrerandotestpow5[l]=0}
  }
  
  
  pvalrerandpb3pow5<-mean(pvalrerandotestpow5)

    ##################################################################################################
  #################################Permuted blocks with the center + 5 factors######################
  ##################################################################################################
  
  set.seed(5239)
  ovariansimple<-matrix(nrow=N,ncol=13)
  pvalue<-numeric(M)
  pvalue1<-numeric(M)
  pvalue2<-numeric(M)
  pvalue3<-numeric(M)
  pvalue4<-numeric(M)
  pvalue5<-numeric(M)
  pvaluepow<-numeric(M)
  pvaluepow1<-numeric(M)
  pvaluepow2<-numeric(M)
  pvaluepow3<-numeric(M)
  pvaluepow4<-numeric(M)
  pvaluepow5<-numeric(M)
  coeff<-numeric(M)
  coeff1<-numeric(M)
  coeff2<-numeric(M)
  coeff3<-numeric(M)
  coeff4<-numeric(M)
  coeff5<-numeric(M)
  Factor1<-c(0:16,20)
  Factor2<-c("no residual disease","<= 2 cm","> 2 cm")
  Factor3<-c("<55","55-65",">=65")
  Factor4<-c("ECOG 0", "ECOG 1", "ECOG 2 or worse")
  Factor5<-c("Good","Intermediate","Poor")
  Factor6<-c("Stage III", "Stage IV", "Stage I, II")
  Grid<-expand.grid(Factor1,Factor2,Factor3,Factor4,Factor5,Factor6)
  Grid <- Grid[order(Grid$Var1,Grid$Var2,Grid$Var3,Grid$Var4,Grid$Var5,Grid$Var6) , ]
  nbrrand<-numeric(length(Grid[,1]))
  nbrrandsimple<-numeric(length(Grid[,1]))
  ######################################Distribution in the strata####################################################

##5factors

Griddist2<-expand.grid(Factor2,Factor3,Factor4,Factor5,Factor6)
Griddist2 <- Griddist2[order(Griddist2$Var1,Griddist2$Var2,Griddist2$Var3,Griddist2$Var4,Griddist2$Var5) , ]
nbrrandsimp2<-numeric(length(Griddist2[,1]))
distr2<-matrix(nrow=length(Griddist2[,1]), ncol=M)

##Distribution in the strata
for (i in 1:M){
  ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
  for(m in 1:length(Griddist2[,1])){
    nbrrandsimp2[m]<-nrow(ovariansimple[ovariansimple$RESDISEA==Griddist2$Var1[m] & ovariansimple$AGE==Griddist2$Var2[m] 
                                 & ovariansimple$PS==Griddist2$Var3[m]& ovariansimple$GRADE==Griddist2$Var4[m] 
                                 & ovariansimple$STAGE==Griddist2$Var5[m] , ])
    
  }
  distr2[,i]<-as.vector(nbrrandsimp2)
}

mean_distr2<-rowMeans(distr2)
mean_distr2<-sort(mean_distr2,decreasing = TRUE)
plot(mean_distr2,type="h",main="Distribution of 500 patients strata from 5 factors",xlab="Number of strata",ylab="Number of patients",lwd=3)
#################################Power############################################################


  for(m in 1:length(Grid[,1])){
    nbrrand[m]<-nrow(ovarian[ovarian$CENTER==Grid$Var1[m] & ovarian$RESDISEA==Grid$Var2[m] & ovarian$AGE==Grid$Var3[m] 
                             & ovarian$PS==Grid$Var4[m] & ovarian$GRADE==Grid$Var5[m] & ovarian$STAGE==Grid$Var6[m] , ])
    
  }
  start.time <- Sys.time()

  set.seed(5239)
  for (i in 1:M){
    Rando5fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando5fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando5fact<-Rando5fact[1:nbrrand[m],]
      
      Rando5fact1<-rbind(Rando5fact1,Rando5fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS,ovarian$GRADE,ovarian$STAGE) , ]
    TRTCENTER<-na.omit(Rando5fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov2 = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov2 ~ 1)
    #Logrank
    surv2<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT,data=ovariansimple)
    pvaluepow[i]<-1 - pchisq(surv2$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT , data = ovariansimple )
    coeff[i]<- exp(cox2$coefficients)
    #Stratified Logrank (1 factor)
    strat21<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvaluepow1[i]<- 1 - pchisq(strat21$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA) , data = ovariansimple )
    coeff1[i]<- exp(cox2$coefficients)
    #Stratified Logrank (3 factors)
    strat22<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvaluepow2[i]<- 1 - pchisq(strat22$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS) , data = ovariansimple )
    coeff2[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors)
    strat23<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvaluepow3[i]<- 1 - pchisq(strat23$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE) , data = ovariansimple )
    coeff3[i]<- exp(cox2$coefficients)
    #Stratified Logrank (Score)
    strat24<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvaluepow4[i]<- 1 - pchisq(strat24$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE) , data = ovariansimple )
    coeff4[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors + Score + Center)
    strat25<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvaluepow5[i]<- 1 - pchisq(strat25$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER) , data = ovariansimple )
    coeff5[i]<- exp(cox2$coefficients)
    
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
  
  #Logrank
  for (j in 1:M){
    
    if (pvaluepow[j]<0.05) {pvaluepow[j]=1}
    else {
      pvaluepow[j]=0
    }
    
  }
  
  pb5pow<-mean(pvaluepow)
  pb5pow
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvaluepow1[j]<0.05) {pvaluepow1[j]=1}
    else {
      pvaluepow1[j]=0
    }
    
  }
  
  pb5pow1<-mean(pvaluepow1)
  pb5pow1
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvaluepow2[j]<0.05) {pvaluepow2[j]=1}
    else {
      pvaluepow2[j]=0
    }
    
  }
  
  pb5pow2<-mean(pvaluepow2)
  pb5pow2
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvaluepow3[j]<0.05) {pvaluepow3[j]=1}
    else {
      pvaluepow3[j]=0
    }
    
  }
  
  pb5pow3<-mean(pvaluepow3)
  pb5pow3
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvaluepow4[j]<0.05) {pvaluepow4[j]=1}
    else {
      pvaluepow4[j]=0
    }
    
  }
  
  pb5pow4<-mean(pvaluepow4)
  pb5pow4
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvaluepow5[j]<0.05) {pvaluepow5[j]=1}
    else {
      pvaluepow5[j]=0
    }
    
  }
  
  pb5pow5<-mean(pvaluepow5)
  
  ##############################Size of the test#################################################
  
  set.seed(5239)
  
  for (i in 1:M){
    Rando5fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando5fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando5fact<-Rando5fact[1:nbrrand[m],]
      
      Rando5fact1<-rbind(Rando5fact1,Rando5fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS,ovarian$GRADE,ovarian$STAGE) , ]
    TRTCENTER<-na.omit(Rando5fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov ~ 1)
    #Stratified Logrank (1 factor)
    strat1<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvalue1[i]<- 1 - pchisq(strat1$chisq, 1)
    #Stratified Logrank (3 factors)
    strat2<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvalue2[i]<- 1 - pchisq(strat2$chisq, 1)
    #Stratified Logrank (5 factors)
    strat3<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvalue3[i]<- 1 - pchisq(strat3$chisq, 1)
    #Stratified Logrank (Score)
    strat4<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvalue4[i]<- 1 - pchisq(strat4$chisq, 1)
    # #Stratified Logrank (5 factors + Score + Center)
    # strat5<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    # pvalue5[i]<- 1 - pchisq(strat5$chisq, 1)
    #Logrank
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    pvalue[i]<-1 - pchisq(surv$chisq, 1)
  }
  
  #Logrank
  for (j in 1:M){
    
    if (pvalue[j]<0.05) {pvalue[j]=1}
    else {
      pvalue[j]=0
    }
    
  }
  
  pb5s<-mean(pvalue)
  pb5s
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvalue1[j]<0.05) {pvalue1[j]=1}
    else {
      pvalue1[j]=0
    }
    
  }
  
  pb5s1<-mean(pvalue1)
  pb5s1
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvalue2[j]<0.05) {pvalue2[j]=1}
    else {
      pvalue2[j]=0
    }
    
  }
  
  pb5s2<-mean(pvalue2)
  pb5s2
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvalue3[j]<0.05) {pvalue3[j]=1}
    else {
      pvalue3[j]=0
    }
    
  }
  
  pb5s3<-mean(pvalue3)
  pb5s3
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvalue4[j]<0.05) {pvalue4[j]=1}
    else {
      pvalue4[j]=0
    }
    
  }
  
  pb5s4<-mean(pvalue4)
  pb5s4
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvalue5[j]<0.05) {pvalue5[j]=1}
    else {
      pvalue5[j]=0
    }
    
  }
  
  pb5s5<-mean(pvalue5)
  ##############################MSE################################################################
  
  MSE<-numeric(M)
  MSE1<-numeric(M)
  MSE2<-numeric(M)
  MSE3<-numeric(M)
  MSE4<-numeric(M)
  MSE5<-numeric(M)
  
  #No stratification
  for (j in 1:M){
    MSE[j]<-(coeff[j]-HR)^2
  }
  
  MSEfinalpb50<-mean(MSE)
  
  #Stratification with 1 factor
  
  for (j in 1:M){
    MSE1[j]<-(coeff1[j]-HR)^2
  }
  
  MSEfinalpb51<-mean(MSE1)
  
  #Stratification with 3 factors
  
  for (j in 1:M){
    MSE2[j]<-(coeff2[j]-HR)^2
  }
  
  MSEfinalpb52<-mean(MSE2)
  
  #Stratification with 5 factors
  
  for (j in 1:M){
    MSE3[j]<-(coeff3[j]-HR)^2
  }
  
  MSEfinalpb53<-mean(MSE3)
  
  #Stratification with the score
  
  for (j in 1:M){
    MSE4[j]<-(coeff4[j]-HR)^2
  }
  
  MSEfinalpb54<-mean(MSE4)
  
  #Stratification with the score + 5 factors + the center
  
  for (j in 1:M){
    MSE5[j]<-(coeff5[j]-HR)^2
  }
  
  MSEfinalpb55<-mean(MSE5)
  
  ###############################################################################################
  ##############################Rerandomization##################################################
  ###############################################################################################
  
  ##############################Size of the test#################################################
  #Logrank  
  
  set.seed(4444)
  nsimu=250
  M=500
  pvalrerandotest<-rep( 0, len=M)
  pvalrerando<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando5fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando5fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando5fact<-Rando5fact[1:nbrrand[m],]
      
      Rando5fact1<-rbind(Rando5fact1,Rando5fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS,ovarian$GRADE,ovarian$STAGE) , ]
    TRTCENTER<-na.omit(Rando5fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS,ovariansimple$GRADE,ovariansimple$STAGE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando5fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] & ovariansimple$GRADE==Grid$Var5[m] & ovariansimple$STAGE==Grid$Var6[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando5fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando5fact<-Rando5fact[1:nbrrandsimple[m],]
        
        Rando5fact1<-rbind(Rando5fact1,Rando5fact)
      }
      
      TRTCENTER<-na.omit(Rando5fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando[l]<0.05) {pvalrerandotest[l]=1}
    else {pvalrerandotest[l]=0}
  }
  
  
  pvalrerandpb50<-mean(pvalrerandotest)
  
  
  #Stratified Logrank (1 factor)
  
  set.seed(4444)
  
  pvalrerando1<-numeric(M)
  pvalrerandotest1<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando5fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando5fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando5fact<-Rando5fact[1:nbrrand[m],]
      
      Rando5fact1<-rbind(Rando5fact1,Rando5fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS,ovarian$GRADE,ovarian$STAGE) , ]
    TRTCENTER<-na.omit(Rando5fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS,ovariansimple$GRADE,ovariansimple$STAGE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando5fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] & ovariansimple$GRADE==Grid$Var5[m] & ovariansimple$STAGE==Grid$Var6[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando5fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando5fact<-Rando5fact[1:nbrrandsimple[m],]
        
        Rando5fact1<-rbind(Rando5fact1,Rando5fact)
      }
      
      TRTCENTER<-na.omit(Rando5fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando1[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando1[l]<0.05) {pvalrerandotest1[l]=1}
    else {pvalrerandotest1[l]=0}
  }
  
  
  pvalrerandpb51<-mean(pvalrerandotest1)
  
  #Stratified Logrank (3 factors)
  
  set.seed(4444)
  
  pvalrerando2<-numeric(M)
  pvalrerandotest2<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando5fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando5fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando5fact<-Rando5fact[1:nbrrand[m],]
      
      Rando5fact1<-rbind(Rando5fact1,Rando5fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS,ovarian$GRADE,ovarian$STAGE) , ]
    TRTCENTER<-na.omit(Rando5fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS,ovariansimple$GRADE,ovariansimple$STAGE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando5fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] & ovariansimple$GRADE==Grid$Var5[m] & ovariansimple$STAGE==Grid$Var6[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando5fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando5fact<-Rando5fact[1:nbrrandsimple[m],]
        
        Rando5fact1<-rbind(Rando5fact1,Rando5fact)
      }
      
      TRTCENTER<-na.omit(Rando5fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando2[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando2[l]<0.05) {pvalrerandotest2[l]=1}
    else {pvalrerandotest2[l]=0}
  }
  
  
  pvalrerandpb52<-mean(pvalrerandotest2)
  
  #Stratified Logrank (5 factors)
  
  
  set.seed(4444)
  
  pvalrerando3<-numeric(M)
  pvalrerandotest3<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando5fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando5fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando5fact<-Rando5fact[1:nbrrand[m],]
      
      Rando5fact1<-rbind(Rando5fact1,Rando5fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS,ovarian$GRADE,ovarian$STAGE) , ]
    TRTCENTER<-na.omit(Rando5fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS,ovariansimple$GRADE,ovariansimple$STAGE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando5fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] & ovariansimple$GRADE==Grid$Var5[m] & ovariansimple$STAGE==Grid$Var6[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando5fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando5fact<-Rando5fact[1:nbrrandsimple[m],]
        
        Rando5fact1<-rbind(Rando5fact1,Rando5fact)
      }
      
      TRTCENTER<-na.omit(Rando5fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando3[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando3[l]<0.05) {pvalrerandotest3[l]=1}
    else {pvalrerandotest3[l]=0}
  }
  
  
  pvalrerandpb53<-mean(pvalrerandotest3)
  
  #Stratified Logrank (Score)
  
  set.seed(4444)
  
  pvalrerando4<-numeric(M)
  pvalrerandotest4<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando5fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando5fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando5fact<-Rando5fact[1:nbrrand[m],]
      
      Rando5fact1<-rbind(Rando5fact1,Rando5fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS,ovarian$GRADE,ovarian$STAGE) , ]
    TRTCENTER<-na.omit(Rando5fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS,ovariansimple$GRADE,ovariansimple$STAGE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando5fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] & ovariansimple$GRADE==Grid$Var5[m] & ovariansimple$STAGE==Grid$Var6[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando5fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando5fact<-Rando5fact[1:nbrrandsimple[m],]
        
        Rando5fact1<-rbind(Rando5fact1,Rando5fact)
      }
      
      TRTCENTER<-na.omit(Rando5fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando4[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando4[l]<0.05) {pvalrerandotest4[l]=1}
    else {pvalrerandotest4[l]=0}
  }
  
  
  pvalrerandpb54<-mean(pvalrerandotest4)
  
  #Stratified Logrank (5 factors + Score + Center)
  
  set.seed(4444)
  
  pvalrerando5<-numeric(M)
  pvalrerandotest5<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando5fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando5fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando5fact<-Rando5fact[1:nbrrand[m],]
      
      Rando5fact1<-rbind(Rando5fact1,Rando5fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS,ovarian$GRADE,ovarian$STAGE) , ]
    TRTCENTER<-na.omit(Rando5fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE)+strata(ovariansimple$CENTER)+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS,ovariansimple$GRADE,ovariansimple$STAGE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando5fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] & ovariansimple$GRADE==Grid$Var5[m] & ovariansimple$STAGE==Grid$Var6[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando5fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando5fact<-Rando5fact[1:nbrrandsimple[m],]
        
        Rando5fact1<-rbind(Rando5fact1,Rando5fact)
      }
      
      TRTCENTER<-na.omit(Rando5fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$CENTER)+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando5[l]<0.05) {pvalrerandotest5[l]=1}
    else {pvalrerandotest5[l]=0}
  }
  
  
  pvalrerandpb55<-mean(pvalrerandotest5)
  
  ##############################Power of the test################################################
  
  #Logrank  
  
  set.seed(4452)
  
  pvalrerandopow<-numeric(M)
  pvalrerandotestpow<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Rando5fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando5fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando5fact<-Rando5fact[1:nbrrand[m],]
      
      Rando5fact1<-rbind(Rando5fact1,Rando5fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS,ovarian$GRADE,ovarian$STAGE) , ]
    TRTCENTER<-na.omit(Rando5fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS,ovariansimple$GRADE,ovariansimple$STAGE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando5fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] & ovariansimple$GRADE==Grid$Var5[m] & ovariansimple$STAGE==Grid$Var6[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando5fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando5fact<-Rando5fact[1:nbrrandsimple[m],]
        
        Rando5fact1<-rbind(Rando5fact1,Rando5fact)
      }
      
      TRTCENTER<-na.omit(Rando5fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow[l]<0.05) {pvalrerandotestpow[l]=1}
    else {pvalrerandotestpow[l]=0}
  }
  
  
  pvalrerandpb5pow<-mean(pvalrerandotestpow)
  
  
  
  #Stratified Logrank (1 factor)
  
  pvalrerandopow1<-numeric(M)
  pvalrerandotestpow1<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  set.seed(4452)
  for (i in 1:M){
    Rando5fact1<-NULL
    TRTCENTER<-NULL
    
    set.seed(4444)
    
    for (m in 1:length(Grid[,1])){
      Rando5fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando5fact<-Rando5fact[1:nbrrand[m],]
      
      Rando5fact1<-rbind(Rando5fact1,Rando5fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS,ovarian$GRADE,ovarian$STAGE) , ]
    TRTCENTER<-na.omit(Rando5fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS,ovariansimple$GRADE,ovariansimple$STAGE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando5fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] & ovariansimple$GRADE==Grid$Var5[m] & ovariansimple$STAGE==Grid$Var6[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando5fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando5fact<-Rando5fact[1:nbrrandsimple[m],]
        
        Rando5fact1<-rbind(Rando5fact1,Rando5fact)
      }
      
      TRTCENTER<-na.omit(Rando5fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow1[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow1[l]<0.05) {pvalrerandotestpow1[l]=1}
    else {pvalrerandotestpow1[l]=0}
  }
  
  
  pvalrerandpb5pow1<-mean(pvalrerandotestpow1)
  
  
  #Stratified Logrank (3 factors)
  
  pvalrerandopow2<-numeric(M)
  pvalrerandotestpow2<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    Rando5fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando5fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando5fact<-Rando5fact[1:nbrrand[m],]
      
      Rando5fact1<-rbind(Rando5fact1,Rando5fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS,ovarian$GRADE,ovarian$STAGE) , ]
    TRTCENTER<-na.omit(Rando5fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS,ovariansimple$GRADE,ovariansimple$STAGE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando5fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] & ovariansimple$GRADE==Grid$Var5[m] & ovariansimple$STAGE==Grid$Var6[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando5fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando5fact<-Rando5fact[1:nbrrandsimple[m],]
        
        Rando5fact1<-rbind(Rando5fact1,Rando5fact)
      }
      
      TRTCENTER<-na.omit(Rando5fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow2[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow2[l]<0.05) {pvalrerandotestpow2[l]=1}
    else {pvalrerandotestpow2[l]=0}
  }
  
  
  pvalrerandpb5pow2<-mean(pvalrerandotestpow2)
  
  #Stratified Logrank (5 factors)
  
  pvalrerandopow3<-numeric(M)
  pvalrerandotestpow3<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    Rando5fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando5fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando5fact<-Rando5fact[1:nbrrand[m],]
      
      Rando5fact1<-rbind(Rando5fact1,Rando5fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS,ovarian$GRADE,ovarian$STAGE) , ]
    TRTCENTER<-na.omit(Rando5fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS,ovariansimple$GRADE,ovariansimple$STAGE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando5fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] & ovariansimple$GRADE==Grid$Var5[m] & ovariansimple$STAGE==Grid$Var6[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando5fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando5fact<-Rando5fact[1:nbrrandsimple[m],]
        
        Rando5fact1<-rbind(Rando5fact1,Rando5fact)
      }
      
      TRTCENTER<-na.omit(Rando5fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow3[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow3[l]<0.05) {pvalrerandotestpow3[l]=1}
    else {pvalrerandotestpow3[l]=0}
  }
  
  
  pvalrerandpb5pow3<-mean(pvalrerandotestpow3)
  
  #Stratified Logrank (Score)
  
  pvalrerandopow4<-numeric(M)
  pvalrerandotestpow4<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    Rando5fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando5fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando5fact<-Rando5fact[1:nbrrand[m],]
      
      Rando5fact1<-rbind(Rando5fact1,Rando5fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS,ovarian$GRADE,ovarian$STAGE) , ]
    TRTCENTER<-na.omit(Rando5fact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS,ovariansimple$GRADE,ovariansimple$STAGE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando5fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] & ovariansimple$GRADE==Grid$Var5[m] & ovariansimple$STAGE==Grid$Var6[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando5fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando5fact<-Rando5fact[1:nbrrandsimple[m],]
        
        Rando5fact1<-rbind(Rando5fact1,Rando5fact)
      }
      
      TRTCENTER<-na.omit(Rando5fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
       else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow4[i]<-a/nsimu 
  }
  
  for (l in 1:M){
    if (pvalrerandopow4[l]<0.05) {pvalrerandotestpow4[l]=1}
    else {pvalrerandotestpow4[l]=0}
  }
  
  
  pvalrerandpb5pow4<-mean(pvalrerandotestpow4)
  
  
  
  #Stratified Logrank (Score + 5 factors + Center)
  
  pvalrerandopow5<-numeric(M)
  pvalrerandotestpow5<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    
    Rando5fact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Rando5fact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Rando5fact<-Rando5fact[1:nbrrand[m],]
      
      Rando5fact1<-rbind(Rando5fact1,Rando5fact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$RESDISEA,ovarian$AGE,ovarian$PS,ovarian$GRADE,ovarian$STAGE) , ]
    TRTCENTER<-na.omit(Rando5fact1)
    
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$RESDISEA,ovariansimple$AGE,ovariansimple$PS,ovariansimple$GRADE,ovariansimple$STAGE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Rando5fact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$RESDISEA==Grid$Var2[m] & ovariansimple$AGE==Grid$Var3[m] 
                                             & ovariansimple$PS==Grid$Var4[m] & ovariansimple$GRADE==Grid$Var5[m] & ovariansimple$STAGE==Grid$Var6[m] , ])
        
      }
      for (m in 1:length(Grid[,1])){
        Rando5fact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Rando5fact<-Rando5fact[1:nbrrandsimple[m],]
        
        Rando5fact1<-rbind(Rando5fact1,Rando5fact)
      }
      
      TRTCENTER<-na.omit(Rando5fact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
       else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$SCORE)+strata(randogroup$CENTER),data=randogroup)
     p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow5[l]<0.05) {pvalrerandotestpow5[l]=1}
    else {pvalrerandotestpow5[l]=0}
  }
  
  
  pvalrerandpb5pow5<-mean(pvalrerandotestpow5)
  
  ##################################################################################################
  #################################Permuted blocks with the center and the score####################
  ##################################################################################################
  
  #################################Power############################################################
  nsimu=500
  M=1000
  set.seed(5239)
  ovariansimple<-matrix(nrow=N,ncol=13)
  pvalue<-numeric(M)
  pvalue1<-numeric(M)
  pvalue2<-numeric(M)
  pvalue3<-numeric(M)
  pvalue4<-numeric(M)
  pvalue5<-numeric(M)
  pvaluepow<-numeric(M)
  pvaluepow1<-numeric(M)
  pvaluepow2<-numeric(M)
  pvaluepow3<-numeric(M)
  pvaluepow4<-numeric(M)
  pvaluepow5<-numeric(M)
  coeff<-numeric(M)
  coeff1<-numeric(M)
  coeff2<-numeric(M)
  coeff3<-numeric(M)
  coeff4<-numeric(M)
  coeff5<-numeric(M)
  Factor1<-c(0:16,20)
  Factor2<-c("Best","Good","Intermediate","Bad","Worst")
  Grid<-expand.grid(Factor1,Factor2)
  Grid <- Grid[order(Grid$Var1,Grid$Var2) , ]
  nbrrand<-numeric(length(Grid[,1]))
  nbrrandsimple<-numeric(length(Grid[,1]))
  
  for(m in 1:length(Grid[,1])){
    nbrrand[m]<-nrow(ovarian[ovarian$CENTER==Grid$Var1[m] & ovarian$SCORE==Grid$Var2[m] , ])
    
  }
  
    set.seed(5239)
  for (i in 1:M){
    
    Randoscofact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Randoscofact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Randoscofact<-Randoscofact[1:nbrrand[m],]
      
      Randoscofact1<-rbind(Randoscofact1,Randoscofact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$SCORE) , ]
    TRTCENTER<-na.omit(Randoscofact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov2 = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov2 ~ 1)
    #Logrank
    surv2<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT,data=ovariansimple)
    pvaluepow[i]<-1 - pchisq(surv2$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT , data = ovariansimple )
    coeff[i]<- exp(cox2$coefficients)
    #Stratified Logrank (1 factor)
    strat21<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvaluepow1[i]<- 1 - pchisq(strat21$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA) , data = ovariansimple )
    coeff1[i]<- exp(cox2$coefficients)
    #Stratified Logrank (3 factors)
    strat22<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvaluepow2[i]<- 1 - pchisq(strat22$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS) , data = ovariansimple )
    coeff2[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors)
    strat23<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvaluepow3[i]<- 1 - pchisq(strat23$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE) , data = ovariansimple )
    coeff3[i]<- exp(cox2$coefficients)
    #Stratified Logrank (Score)
    strat24<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvaluepow4[i]<- 1 - pchisq(strat24$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE) , data = ovariansimple )
    coeff4[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors + Score + Center)
    strat25<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvaluepow5[i]<- 1 - pchisq(strat25$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER) , data = ovariansimple )
    coeff5[i]<- exp(cox2$coefficients)
    
  }
  
  
  #Logrank
  for (j in 1:M){
    
    if (pvaluepow[j]<0.05) {pvaluepow[j]=1}
    else {
      pvaluepow[j]=0
    }
    
  }
  
  pbscopow<-mean(pvaluepow)
  
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvaluepow1[j]<0.05) {pvaluepow1[j]=1}
    else {
      pvaluepow1[j]=0
    }
    
  }
  
  pbscopow1<-mean(pvaluepow1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvaluepow2[j]<0.05) {pvaluepow2[j]=1}
    else {
      pvaluepow2[j]=0
    }
    
  }
  
  pbscopow2<-mean(pvaluepow2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvaluepow3[j]<0.05) {pvaluepow3[j]=1}
    else {
      pvaluepow3[j]=0
    }
    
  }
  
  pbscopow3<-mean(pvaluepow3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvaluepow4[j]<0.05) {pvaluepow4[j]=1}
    else {
      pvaluepow4[j]=0
    }
    
  }
  
  pbscopow4<-mean(pvaluepow4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvaluepow5[j]<0.05) {pvaluepow5[j]=1}
    else {
      pvaluepow5[j]=0
    }
    
  }
  
  pbscopow5<-mean(pvaluepow5)

  ##############################Size of the test#################################################
  
  set.seed(5239)
  
  for (i in 1:M){
    Randoscofact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Randoscofact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Randoscofact<-Randoscofact[1:nbrrand[m],]
      
      Randoscofact1<-rbind(Randoscofact1,Randoscofact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$SCORE) , ]
    TRTCENTER<-na.omit(Randoscofact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov ~ 1)
    #Stratified Logrank (1 factor)
    strat1<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvalue1[i]<- 1 - pchisq(strat1$chisq, 1)
    #Stratified Logrank (3 factors)
    strat2<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvalue2[i]<- 1 - pchisq(strat2$chisq, 1)
    #Stratified Logrank (5 factors)
    strat3<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvalue3[i]<- 1 - pchisq(strat3$chisq, 1)
    #Stratified Logrank (Score)
    strat4<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvalue4[i]<- 1 - pchisq(strat4$chisq, 1)
    #Stratified Logrank (5 factors + Score + Center)
    strat5<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvalue5[i]<- 1 - pchisq(strat5$chisq, 1)
    #Logrank
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    pvalue[i]<-1 - pchisq(surv$chisq, 1)
  }
  
  #Logrank
  for (j in 1:M){
    
    if (pvalue[j]<0.05) {pvalue[j]=1}
    else {
      pvalue[j]=0
    }
    
  }
  
  pbscos<-mean(pvalue)
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvalue1[j]<0.05) {pvalue1[j]=1}
    else {
      pvalue1[j]=0
    }
    
  }
  
  pbscos1<-mean(pvalue1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvalue2[j]<0.05) {pvalue2[j]=1}
    else {
      pvalue2[j]=0
    }
    
  }
  
  pbscos2<-mean(pvalue2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvalue3[j]<0.05) {pvalue3[j]=1}
    else {
      pvalue3[j]=0
    }
    
  }
  
  pbscos3<-mean(pvalue3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvalue4[j]<0.05) {pvalue4[j]=1}
    else {
      pvalue4[j]=0
    }
    
  }
  
  pbscos4<-mean(pvalue4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvalue5[j]<0.05) {pvalue5[j]=1}
    else {
      pvalue5[j]=0
    }
    
  }
  
  pbscos5<-mean(pvalue5)
  
  ##############################MSE################################################################
  
  MSE<-numeric(M)
  MSE1<-numeric(M)
  MSE2<-numeric(M)
  MSE3<-numeric(M)
  MSE4<-numeric(M)
  MSE5<-numeric(M)
  
  #No stratification
  for (j in 1:M){
    MSE[j]<-(coeff[j]-HR)^2
  }
  
  MSEfinalpbsco<-mean(MSE)
  
  #Stratification with 1 factor
  
  for (j in 1:M){
    MSE1[j]<-(coeff1[j]-HR)^2
  }
  
  MSEfinalpbsco1<-mean(MSE1)
  
  #Stratification with 3 factors
  
  for (j in 1:M){
    MSE2[j]<-(coeff2[j]-HR)^2
  }
  
  MSEfinalpbsco2<-mean(MSE2)
  
  #Stratification with 5 factors
  
  for (j in 1:M){
    MSE3[j]<-(coeff3[j]-HR)^2
  }
  
  MSEfinalpbsco3<-mean(MSE3)
  
  #Stratification with the score
  
  for (j in 1:M){
    MSE4[j]<-(coeff4[j]-HR)^2
  }
  
  MSEfinalpbsco4<-mean(MSE4)
  
  #Stratification with the score + 5 factors + the center
  
  for (j in 1:M){
    MSE5[j]<-(coeff5[j]-HR)^2
  }
  
  MSEfinalpbsco5<-mean(MSE5)
   
  
  ###############################################################################################
  ##############################Rerandomization##################################################
  ###############################################################################################
  
  ##############################Size of the test#################################################
  #Logrank  
  
  set.seed(4444)
  
  pvalrerandotest<-rep( 0, len=M)
  pvalrerando<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Randoscofact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Randoscofact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Randoscofact<-Randoscofact[1:nbrrand[m],]
      
      Randoscofact1<-rbind(Randoscofact1,Randoscofact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$SCORE) , ]
    TRTCENTER<-na.omit(Randoscofact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$SCORE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Randoscofact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$SCORE==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Randoscofact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Randoscofact<-Randoscofact[1:nbrrandsimple[m],]
        
        Randoscofact1<-rbind(Randoscofact1,Randoscofact)
      }
      
      TRTCENTER<-na.omit(Randoscofact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando[l]<0.05) {pvalrerandotest[l]=1}
    else {pvalrerandotest[l]=0}
  }
  
  
  pvalrerandopbsco<-mean(pvalrerandotest)
  
  
  #Stratified Logrank (1 factor)
  
  set.seed(4444)
  
  pvalrerando1<-numeric(M)
  pvalrerandotest1<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Randoscofact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Randoscofact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Randoscofact<-Randoscofact[1:nbrrand[m],]
      
      Randoscofact1<-rbind(Randoscofact1,Randoscofact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$SCORE) , ]
    TRTCENTER<-na.omit(Randoscofact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$SCORE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Randoscofact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$SCORE==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Randoscofact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Randoscofact<-Randoscofact[1:nbrrandsimple[m],]
        
        Randoscofact1<-rbind(Randoscofact1,Randoscofact)
      }
      
      TRTCENTER<-na.omit(Randoscofact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando1[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando1[l]<0.05) {pvalrerandotest1[l]=1}
    else {pvalrerandotest1[l]=0}
  }
  
  
  pvalrerandpbsco1<-mean(pvalrerandotest1)
  
  #Stratified Logrank (3 factors)
  
  set.seed(4444)
  
  pvalrerando2<-numeric(M)
  pvalrerandotest2<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Randoscofact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Randoscofact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Randoscofact<-Randoscofact[1:nbrrand[m],]
      
      Randoscofact1<-rbind(Randoscofact1,Randoscofact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$SCORE) , ]
    TRTCENTER<-na.omit(Randoscofact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$SCORE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Randoscofact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$SCORE==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Randoscofact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Randoscofact<-Randoscofact[1:nbrrandsimple[m],]
        
        Randoscofact1<-rbind(Randoscofact1,Randoscofact)
      }
      
      TRTCENTER<-na.omit(Randoscofact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando2[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando2[l]<0.05) {pvalrerandotest2[l]=1}
    else {pvalrerandotest2[l]=0}
  }
  
  
  pvalrerandpbsco2<-mean(pvalrerandotest2)
  
  #Stratified Logrank (5 factors)
  
  
  set.seed(4444)
  
  pvalrerando3<-numeric(M)
  pvalrerandotest3<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Randoscofact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Randoscofact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Randoscofact<-Randoscofact[1:nbrrand[m],]
      
      Randoscofact1<-rbind(Randoscofact1,Randoscofact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$SCORE) , ]
    TRTCENTER<-na.omit(Randoscofact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$SCORE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Randoscofact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$SCORE==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Randoscofact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Randoscofact<-Randoscofact[1:nbrrandsimple[m],]
        
        Randoscofact1<-rbind(Randoscofact1,Randoscofact)
      }
      
      TRTCENTER<-na.omit(Randoscofact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando3[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando3[l]<0.05) {pvalrerandotest3[l]=1}
    else {pvalrerandotest3[l]=0}
  }
  
  
  pvalrerandpbsco3<-mean(pvalrerandotest3)
  
  #Stratified Logrank (Score)
  
  set.seed(4444)
  
  pvalrerando4<-numeric(M)
  pvalrerandotest4<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Randoscofact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Randoscofact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Randoscofact<-Randoscofact[1:nbrrand[m],]
      
      Randoscofact1<-rbind(Randoscofact1,Randoscofact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$SCORE) , ]
    TRTCENTER<-na.omit(Randoscofact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$SCORE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Randoscofact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$SCORE==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Randoscofact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Randoscofact<-Randoscofact[1:nbrrandsimple[m],]
        
        Randoscofact1<-rbind(Randoscofact1,Randoscofact)
      }
      
      TRTCENTER<-na.omit(Randoscofact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando4[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando4[l]<0.05) {pvalrerandotest4[l]=1}
    else {pvalrerandotest4[l]=0}
  }
  
  
  pvalrerandpbsco4<-mean(pvalrerandotest4)
  
  #Stratified Logrank (5 factors + Score + Center)
  
  set.seed(4444)
  
  pvalrerando5<-numeric(M)
  pvalrerandotest5<-rep( 0, len=M)
  
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Randoscofact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Randoscofact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Randoscofact<-Randoscofact[1:nbrrand[m],]
      
      Randoscofact1<-rbind(Randoscofact1,Randoscofact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$SCORE) , ]
    TRTCENTER<-na.omit(Randoscofact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE)+strata(ovariansimple$CENTER)+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$SCORE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Randoscofact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$SCORE==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Randoscofact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Randoscofact<-Randoscofact[1:nbrrandsimple[m],]
        
        Randoscofact1<-rbind(Randoscofact1,Randoscofact)
      }
      
      TRTCENTER<-na.omit(Randoscofact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$CENTER)+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando5[l]<0.05) {pvalrerandotest5[l]=1}
    else {pvalrerandotest5[l]=0}
  }
  
  
  pvalrerandpbsco5<-mean(pvalrerandotest5)
  pvalrerandpbsco5

  ##############################Power of the test################################################
  
  #Logrank  
  
  set.seed(4452)
  
  pvalrerandopow<-numeric(M)
  pvalrerandotestpow<-numeric(M)
 
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    Randoscofact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Randoscofact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Randoscofact<-Randoscofact[1:nbrrand[m],]
      
      Randoscofact1<-rbind(Randoscofact1,Randoscofact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$SCORE) , ]
    TRTCENTER<-na.omit(Randoscofact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$SCORE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Randoscofact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$SCORE==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Randoscofact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Randoscofact<-Randoscofact[1:nbrrandsimple[m],]
        
        Randoscofact1<-rbind(Randoscofact1,Randoscofact)
      }
      
      TRTCENTER<-na.omit(Randoscofact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow[l]<0.05) {pvalrerandotestpow[l]=1}
    else {pvalrerandotestpow[l]=0}
  }
  
  
  pvalrerandpbscopow<-mean(pvalrerandotestpow)
  pvalrerandpbscopow
  
  
  
  #Stratified Logrank (1 factor)
  
  pvalrerandopow1<-numeric(M)
  pvalrerandotestpow1<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    Randoscofact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Randoscofact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Randoscofact<-Randoscofact[1:nbrrand[m],]
      
      Randoscofact1<-rbind(Randoscofact1,Randoscofact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$SCORE) , ]
    TRTCENTER<-na.omit(Randoscofact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$SCORE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Randoscofact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$SCORE==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Randoscofact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Randoscofact<-Randoscofact[1:nbrrandsimple[m],]
        
        Randoscofact1<-rbind(Randoscofact1,Randoscofact)
      }
      
      TRTCENTER<-na.omit(Randoscofact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow1[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow1[l]<0.05) {pvalrerandotestpow1[l]=1}
    else {pvalrerandotestpow1[l]=0}
  }
  
  
  pvalrerandpbscopow1<-mean(pvalrerandotestpow1)
  pvalrerandpbscopow1
  
  
  #Stratified Logrank (3 factors)
  
  pvalrerandopow2<-numeric(M)
  pvalrerandotestpow2<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    Randoscofact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Randoscofact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Randoscofact<-Randoscofact[1:nbrrand[m],]
      
      Randoscofact1<-rbind(Randoscofact1,Randoscofact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$SCORE) , ]
    TRTCENTER<-na.omit(Randoscofact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$SCORE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Randoscofact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$SCORE==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Randoscofact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Randoscofact<-Randoscofact[1:nbrrandsimple[m],]
        
        Randoscofact1<-rbind(Randoscofact1,Randoscofact)
      }
      
      TRTCENTER<-na.omit(Randoscofact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow2[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow2[l]<0.05) {pvalrerandotestpow2[l]=1}
    else {pvalrerandotestpow2[l]=0}
  }
  
  
  pvalrerandpbscopow2<-mean(pvalrerandotestpow2)
  pvalrerandpbscopow2
  
  #Stratified Logrank (5 factors)
  
  pvalrerandopow3<-numeric(M)
  pvalrerandotestpow3<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    Randoscofact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Randoscofact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Randoscofact<-Randoscofact[1:nbrrand[m],]
      
      Randoscofact1<-rbind(Randoscofact1,Randoscofact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$SCORE) , ]
    TRTCENTER<-na.omit(Randoscofact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$SCORE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Randoscofact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$SCORE==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Randoscofact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Randoscofact<-Randoscofact[1:nbrrandsimple[m],]
        
        Randoscofact1<-rbind(Randoscofact1,Randoscofact)
      }
      
      TRTCENTER<-na.omit(Randoscofact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow3[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow3[l]<0.05) {pvalrerandotestpow3[l]=1}
    else {pvalrerandotestpow3[l]=0}
  }
  
  
  pvalrerandpbscopow3<-mean(pvalrerandotestpow3)
  pvalrerandpbscopow3
  
  #Stratified Logrank (Score)
  
  pvalrerandopow4<-numeric(M)
  pvalrerandotestpow4<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    Randoscofact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Randoscofact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Randoscofact<-Randoscofact[1:nbrrand[m],]
      
      Randoscofact1<-rbind(Randoscofact1,Randoscofact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$SCORE) , ]
    TRTCENTER<-na.omit(Randoscofact1)
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$SCORE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Randoscofact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$SCORE==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Randoscofact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Randoscofact<-Randoscofact[1:nbrrandsimple[m],]
        
        Randoscofact1<-rbind(Randoscofact1,Randoscofact)
      }
      
      TRTCENTER<-na.omit(Randoscofact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
       else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow4[i]<-a/nsimu 
  }
  
  for (l in 1:M){
    if (pvalrerandopow4[l]<0.05) {pvalrerandotestpow4[l]=1}
    else {pvalrerandotestpow4[l]=0}
  }
  
  
  pvalrerandpbscopow4<-mean(pvalrerandotestpow4)
  pvalrerandpbscopow4
  
  
  
  #Stratified Logrank (Score + 5 factors + Center)
  
  pvalrerandopow5<-numeric(M)
  pvalrerandotestpow5<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    
    Randoscofact1<-NULL
    TRTCENTER<-NULL
    
    for (m in 1:length(Grid[,1])){
      Randoscofact<-blockrand(nbrrand[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
      
      Randoscofact<-Randoscofact[1:nbrrand[m],]
      
      Randoscofact1<-rbind(Randoscofact1,Randoscofact)
    }
    
    
    ovarian <- ovarian[order(ovarian$CENTER,ovarian$SCORE) , ]
    TRTCENTER<-na.omit(Randoscofact1)
    
    TRT<-TRTCENTER$treatment
    TRT<-as.numeric(TRT)
    
    ovarian<-cbind(ovarian[,1:12],TRT)
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    ovariansimple <- ovariansimple[order(ovariansimple$CENTER,ovariansimple$SCORE) , ]
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      Randoscofact1<-NULL
      TRTCENTER<-NULL
      for(m in 1:length(Grid[,1])){
        nbrrandsimple[m]<-nrow(ovariansimple[ovariansimple$CENTER==Grid$Var1[m] & ovariansimple$SCORE==Grid$Var2[m], ])
        
      }
      for (m in 1:length(Grid[,1])){
        Randoscofact<-blockrand(nbrrandsimple[m],num.levels = 2, levels = c("0","1"),block.sizes = c(2,2))
        
        Randoscofact<-Randoscofact[1:nbrrandsimple[m],]
        
        Randoscofact1<-rbind(Randoscofact1,Randoscofact)
      }
      
      TRTCENTER<-na.omit(Randoscofact1)
      TRT<-TRTCENTER$treatment
      TRT<-as.numeric(TRT)
      randogroup[,1+k]<-TRT
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
       else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$SCORE)+strata(randogroup$CENTER),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow5[l]<0.05) {pvalrerandotestpow5[l]=1}
    else {pvalrerandotestpow5[l]=0}
  }
  
  
  pvalrerandpbscopow5<-mean(pvalrerandotestpow5)
  pvalrerandpbscopow5
  
  
  
  ##################################################################################################
  #################################Minimization with the center#####################################
  ##################################################################################################
  
  ##################################Alogrithme minimization#########################################
  
  miniAlgo <- function(dataset,strata,nbgp,determinism.percent,rando.ratio,weight){
    
    nbpatients=nrow(dataset)
    treatment.arm=NULL
    
    for(pat.i in 1:nbpatients){
      
      B<-data.frame(NULL)
      for(k in 1:length(strata)){
        for(j in 1:nbgp){
          if(pat.i==1){
            B[k,j]=0
          } else {
            B[k,j]<-length( which( dataset[1:(pat.i-1),strata[k]] == dataset[pat.i,strata[k]] & treatment.arm[1:(pat.i-1)] == j ))
          }			
        }	
      }
      
      # compute Ai
      A=(1/rando.ratio)*apply(weight*B,2,sum)
      
      # Identify all treatments with the lowest value Ai
      ind.min=as.vector(which(A==min(A)))
      ind=as.vector(which(A!=min(A)))
      
      # determine whether the group is ind.min or ind
      R=runif(2)
      
      # Select randomly one of the treatments identified in the previous step
      if (length(ind)==0 | R[1]<determinism.percent) {
        if (length(ind.min)==1) {
          arm=ind.min
        } else {
          int<-levels(cut(0:1,dig.lab=2,breaks=length(ind.min)))
          inter<-cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", int) ),upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", int) ))
          inter[1]<-0
          for (i in 1:length(ind.min)) {
            if (inter[i,1]<R[2] & R[2]<=inter[i,2]){
              arm<-ind.min[i]
            }
          }
        }
      } else {
        if (length(ind)==1) {
          arm=ind
        } else {
          int<-levels(cut(0:1,dig.lab=2,breaks=length(ind)))
          inter<-cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", int) ),upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", int) ))
          inter[1]<-0
          for (i in 1:length(ind)) {
            if (inter[i,1]<R[2] & R[2]<=inter[i,2]){
              arm<-ind[i]
            }
          }
        }
      }
      
      treatment.arm=c(treatment.arm,arm)
    }
    
    return(treatment.arm)
  }
  
  ##############################Power of the test################################################
  ovariansimple<-matrix(nrow=N,ncol=12)
  pvalue<-numeric(M)
  pvalue1<-numeric(M)
  pvalue2<-numeric(M)
  pvalue3<-numeric(M)
  pvalue4<-numeric(M)
  pvalue5<-numeric(M)
  pvaluepow<-numeric(M)
  pvaluepow1<-numeric(M)
  pvaluepow2<-numeric(M)
  pvaluepow3<-numeric(M)
  pvaluepow4<-numeric(M)
  pvaluepow5<-numeric(M)
  coeff<-numeric(M)
  coeff1<-numeric(M)
  coeff2<-numeric(M)
  coeff3<-numeric(M)
  coeff4<-numeric(M)
  coeff5<-numeric(M)
  set.seed(5239)
  for (i in 1:M){
    
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    ovariansimple<-cbind(ovariansimple,TRT)
    
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov2 = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov2 ~ 1)
    #Logrank
    surv2<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT,data=ovariansimple)
    pvaluepow[i]<-1 - pchisq(surv2$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT , data = ovariansimple )
    coeff[i]<- exp(cox2$coefficients)
    #Stratified Logrank (1 factor)
    strat21<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvaluepow1[i]<- 1 - pchisq(strat21$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA) , data = ovariansimple )
    coeff1[i]<- exp(cox2$coefficients)
    #Stratified Logrank (3 factors)
    strat22<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvaluepow2[i]<- 1 - pchisq(strat22$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS) , data = ovariansimple )
    coeff2[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors)
    strat23<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvaluepow3[i]<- 1 - pchisq(strat23$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE) , data = ovariansimple )
    coeff3[i]<- exp(cox2$coefficients)
    #Stratified Logrank (Score)
    strat24<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvaluepow4[i]<- 1 - pchisq(strat24$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE) , data = ovariansimple )
    coeff4[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors + Score + Center)
    strat25<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvaluepow5[i]<- 1 - pchisq(strat25$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER) , data = ovariansimple )
    coeff5[i]<- exp(cox2$coefficients)
    
  }
  
  #Logrank
  for (j in 1:M){
    
    if (pvaluepow[j]<0.05) {pvaluepow[j]=1}
    else {
      pvaluepow[j]=0
    }
    
  }
  
  minipow<-mean(pvaluepow)
  minipow
  
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvaluepow1[j]<0.05) {pvaluepow1[j]=1}
    else {
      pvaluepow1[j]=0
    }
    
  }
  
  minipow1<-mean(pvaluepow1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvaluepow2[j]<0.05) {pvaluepow2[j]=1}
    else {
      pvaluepow2[j]=0
    }
    
  }
  
  minipow2<-mean(pvaluepow2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvaluepow3[j]<0.05) {pvaluepow3[j]=1}
    else {
      pvaluepow3[j]=0
    }
    
  }
  
  minipow3<-mean(pvaluepow3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvaluepow4[j]<0.05) {pvaluepow4[j]=1}
    else {
      pvaluepow4[j]=0
    }
    
  }
  
  minipow4<-mean(pvaluepow4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvaluepow5[j]<0.05) {pvaluepow5[j]=1}
    else {
      pvaluepow5[j]=0
    }
    
  }
  
  minipow5<-mean(pvaluepow5)
  
  ##############################Size of the test#################################################
  set.seed(5239)
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov ~ 1)
    #Stratified Logrank (1 factor)
    strat1<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvalue1[i]<- 1 - pchisq(strat1$chisq, 1)
    #Stratified Logrank (3 factors)
    strat2<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvalue2[i]<- 1 - pchisq(strat2$chisq, 1)
    #Stratified Logrank (5 factors)
    strat3<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvalue3[i]<- 1 - pchisq(strat3$chisq, 1)
    #Stratified Logrank (Score)
    strat4<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvalue4[i]<- 1 - pchisq(strat4$chisq, 1)
    #Stratified Logrank (5 factors + Score + Center)
    strat5<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvalue5[i]<- 1 - pchisq(strat5$chisq, 1)
    #Logrank
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    pvalue[i]<-1 - pchisq(surv$chisq, 1)
  }
  
  #Logrank
  for (j in 1:M){
    
    if (pvalue[j]<0.05) {pvalue[j]=1}
    else {
      pvalue[j]=0
    }
    
  }
  
  minis<-mean(pvalue)
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvalue1[j]<0.05) {pvalue1[j]=1}
    else {
      pvalue1[j]=0
    }
    
  }
  
  minis1<-mean(pvalue1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvalue2[j]<0.05) {pvalue2[j]=1}
    else {
      pvalue2[j]=0
    }
    
  }
  
  minis2<-mean(pvalue2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvalue3[j]<0.05) {pvalue3[j]=1}
    else {
      pvalue3[j]=0
    }
    
  }
  
  minis3<-mean(pvalue3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvalue4[j]<0.05) {pvalue4[j]=1}
    else {
      pvalue4[j]=0
    }
    
  }
  
  minis4<-mean(pvalue4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvalue5[j]<0.05) {pvalue5[j]=1}
    else {
      pvalue5[j]=0
    }
    
  }
  
  minis5<-mean(pvalue5)

  ##############################MSE################################################################
  
  MSE<-numeric(M)
  MSE1<-numeric(M)
  MSE2<-numeric(M)
  MSE3<-numeric(M)
  MSE4<-numeric(M)
  MSE5<-numeric(M)
  
  #No stratification
  for (j in 1:M){
    MSE[j]<-(coeff[j]-HR)^2
  }
  
  MSEfinalmini<-mean(MSE)
  
  #Stratification with 1 factor
  
  for (j in 1:M){
    MSE1[j]<-(coeff1[j]-HR)^2
  }
  
  MSEfinalmini1<-mean(MSE1)
  
  #Stratification with 3 factors
  
  for (j in 1:M){
    MSE2[j]<-(coeff2[j]-HR)^2
  }
  
  MSEfinalmini2<-mean(MSE2)
  
  #Stratification with 5 factors
  
  for (j in 1:M){
    MSE3[j]<-(coeff3[j]-HR)^2
  }
  
  MSEfinalmini3<-mean(MSE3)
  
  #Stratification with the score
  
  for (j in 1:M){
    MSE4[j]<-(coeff4[j]-HR)^2
  }
  
  MSEfinalmini4<-mean(MSE4)
  
  #Stratification with the score + 5 factors + the center
  
  for (j in 1:M){
    MSE5[j]<-(coeff5[j]-HR)^2
  }
  
  MSEfinalmini5<-mean(MSE5)
  
  ###############################################################################################
  ##############################Rerandomization##################################################
  ###############################################################################################
  
  ##############################Size of the test#################################################
  #Logrank  
  
  set.seed(4444)
  
  pval.sim=stat.sim=NULL
  
  pvalrerandotest<-rep( 0, len=M)
  pvalrerando<-numeric(M)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando[l]<0.05) {pvalrerandotest[l]=1}
    else {pvalrerandotest[l]=0}
  }
  
  
  pvalrerandmini<-mean(pvalrerandotest)
  pvalrerandmini
  
  
  #Stratified Logrank (1 factor)
  
  set.seed(4444)
  pval.sim=stat.sim=NULL
  
  pvalrerando1<-numeric(M)
  pvalrerandotest1<-rep( 0, len=M)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando1[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando1[l]<0.05) {pvalrerandotest1[l]=1}
    else {pvalrerandotest1[l]=0}
  }
  
  
  pvalrerandmini1<-mean(pvalrerandotest1)
  pvalrerandmini1
  
  #Stratified Logrank (3 factors)
  
  set.seed(4444)
  
  pvalrerando2<-numeric(M)
  pvalrerandotest2<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando2[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando2[l]<0.05) {pvalrerandotest2[l]=1}
    else {pvalrerandotest2[l]=0}
  }
  
  
  pvalrerandmini2<-mean(pvalrerandotest2)
  pvalrerandmini2
  
  #Stratified Logrank (5 factors)
  
  
  set.seed(4444)
  
  pvalrerando3<-numeric(M)
  pvalrerandotest3<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando3[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando3[l]<0.05) {pvalrerandotest3[l]=1}
    else {pvalrerandotest3[l]=0}
  }
  
  
  pvalrerandmini3<-mean(pvalrerandotest3)
  pvalrerandmini3
  
  #Stratified Logrank (Score)
  
  set.seed(4444)
  
  
  rerando4<-rep( 0, len=nsimu)
  pvalrerando4<-numeric(M)
  pvalrerandotest4<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando4[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando4[l]<0.05) {pvalrerandotest4[l]=1}
    else {pvalrerandotest4[l]=0}
  }
  
  
  pvalrerandmini4<-mean(pvalrerandotest4)
  pvalrerandmini4
  
  #Stratified Logrank (5 factors + Score + Center)
  
  set.seed(4444)
  
  pvalrerando5<-numeric(M)
  pvalrerandotest5<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE)+strata(randogroup$SCORE)+strata(randogroup$CENTER),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando5[l]<0.05) {pvalrerandotest5[l]=1}
    else {pvalrerandotest5[l]=0}
  }
  
  
  pvalrerandmini5<-mean(pvalrerandotest5)
  pvalrerandmini5
  
  sizererando<-c(pvalrerandmini,pvalrerandmini1,pvalrerandmini2,pvalrerandmini3,pvalrerandmini4,pvalrerandmini5)
  
  barplot(sizererando)
  
  ##############################Power of the test################################################
  
  #Logrank  
  
  set.seed(4452)
  
  
  pvalrerandopow<-numeric(M)
  pvalrerando<-numeric(M)

  pvalrerandotestpow<-numeric(M)

  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow[i]<-a/nsimu 
  }
  
  for (l in 1:M){
    if (pvalrerandopow[l]<0.05) {pvalrerandotestpow[l]=1}
    else {pvalrerandotestpow[l]=0}
  }
  
  
  pvalrerandminipow<-mean(pvalrerandotestpow)
  pvalrerandminipow

  
  #Stratified Logrank (1 factor)
  
  pvalrerandopow1<-numeric(M)
  pvalrerandotestpow1<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow1[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow1[l]<0.05) {pvalrerandotestpow1[l]=1}
    else {pvalrerandotestpow1[l]=0}
  }
  
  
  pvalrerandminipow1<-mean(pvalrerandotestpow1)
  pvalrerandminipow1

  
  
  #Stratified Logrank (3 factors)
  
  pvalrerandopow2<-numeric(M)
  pvalrerandotestpow2<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=0.3*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow2[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow2[l]<0.05) {pvalrerandotestpow2[l]=1}
    else {pvalrerandotestpow2[l]=0}
  }
  
  
  pvalrerandminipow2<-mean(pvalrerandotestpow2)
  pvalrerandminipow2

  
  #Stratified Logrank (5 factors)
  
  pvalrerandopow3<-numeric(M)
  pvalrerandotestpow3<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow3[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow3[l]<0.05) {pvalrerandotestpow3[l]=1}
    else {pvalrerandotestpow3[l]=0}
  }
  
  
  pvalrerandminipow3<-mean(pvalrerandotestpow3)
  pvalrerandminipow3
  
  #Stratified Logrank (Score)
  
  pvalrerandopow4<-numeric(M)
  pvalrerandotestpow4<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9 ) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
     }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow4[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow4[l]<0.05) {pvalrerandotestpow4[l]=1}
    else {pvalrerandotestpow4[l]=0}
  }
  
  
  pvalrerandminipow4<-mean(pvalrerandotestpow4)
  pvalrerandminipow4
  
  #Stratified Logrank (Score + 5 factors + Center)
  
  pvalrerandopow5<-numeric(M)
  pvalrerandotestpow5<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE)+strata(ovariansimple$CENTER)+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata="CENTER",nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=1)
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$CENTER)+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow5[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow5[l]<0.05) {pvalrerandotestpow5[l]=1}
    else {pvalrerandotestpow5[l]=0}
  }
  
  
  pvalrerandminipow5<-mean(pvalrerandotestpow5)
  pvalrerandminipow5
  
  
  ##################################################################################################
  #################################Minimization with the center and 1 factor########################
  ##################################################################################################
  
  ##############################Power of the test################################################
  ovariansimple<-matrix(nrow=N,ncol=12)
  pvalue<-numeric(M)
  pvalue1<-numeric(M)
  pvalue2<-numeric(M)
  pvalue3<-numeric(M)
  pvalue4<-numeric(M)
  pvalue5<-numeric(M)
  pvaluepow<-numeric(M)
  pvaluepow1<-numeric(M)
  pvaluepow2<-numeric(M)
  pvaluepow3<-numeric(M)
  pvaluepow4<-numeric(M)
  pvaluepow5<-numeric(M)
  coeff<-numeric(M)
  coeff1<-numeric(M)
  coeff2<-numeric(M)
  coeff3<-numeric(M)
  coeff4<-numeric(M)
  coeff5<-numeric(M)
  set.seed(5239)
  for (i in 1:M){
    
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov2 = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov2 ~ 1)
    #Logrank
    surv2<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT,data=ovariansimple)
    pvaluepow[i]<-1 - pchisq(surv2$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT , data = ovariansimple )
    coeff[i]<- exp(cox2$coefficients)
    #Stratified Logrank (1 factor)
    strat21<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvaluepow1[i]<- 1 - pchisq(strat21$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA) , data = ovariansimple )
    coeff1[i]<- exp(cox2$coefficients)
    #Stratified Logrank (3 factors)
    strat22<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvaluepow2[i]<- 1 - pchisq(strat22$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS) , data = ovariansimple )
    coeff2[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors)
    strat23<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvaluepow3[i]<- 1 - pchisq(strat23$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE) , data = ovariansimple )
    coeff3[i]<- exp(cox2$coefficients)
    #Stratified Logrank (Score)
    strat24<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvaluepow4[i]<- 1 - pchisq(strat24$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE) , data = ovariansimple )
    coeff4[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors + Score + Center)
    strat25<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvaluepow5[i]<- 1 - pchisq(strat25$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER) , data = ovariansimple )
    coeff5[i]<- exp(cox2$coefficients)
    
  }
  
  #Logrank
  for (j in 1:M){
    
    if (pvaluepow[j]<0.05) {pvaluepow[j]=1}
    else {
      pvaluepow[j]=0
    }
    
  }
  
  mini1pow<-mean(pvaluepow)
  
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvaluepow1[j]<0.05) {pvaluepow1[j]=1}
    else {
      pvaluepow1[j]=0
    }
    
  }
  
  mini1pow1<-mean(pvaluepow1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvaluepow2[j]<0.05) {pvaluepow2[j]=1}
    else {
      pvaluepow2[j]=0
    }
    
  }
  
  mini1pow2<-mean(pvaluepow2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvaluepow3[j]<0.05) {pvaluepow3[j]=1}
    else {
      pvaluepow3[j]=0
    }
    
  }
  
  mini1pow3<-mean(pvaluepow3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvaluepow4[j]<0.05) {pvaluepow4[j]=1}
    else {
      pvaluepow4[j]=0
    }
    
  }
  
  mini1pow4<-mean(pvaluepow4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvaluepow5[j]<0.05) {pvaluepow5[j]=1}
    else {
      pvaluepow5[j]=0
    }
    
  }
  
  mini1pow5<-mean(pvaluepow5)

  
  ##############################Size of the test#################################################
  set.seed(5239)
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov ~ 1)
    #Stratified Logrank (1 factor)
    strat1<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvalue1[i]<- 1 - pchisq(strat1$chisq, 1)
    #Stratified Logrank (3 factors)
    strat2<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvalue2[i]<- 1 - pchisq(strat2$chisq, 1)
    #Stratified Logrank (5 factors)
    strat3<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvalue3[i]<- 1 - pchisq(strat3$chisq, 1)
    #Stratified Logrank (Score)
    strat4<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvalue4[i]<- 1 - pchisq(strat4$chisq, 1)
    #Stratified Logrank (5 factors + Score + Center)
    strat5<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvalue5[i]<- 1 - pchisq(strat5$chisq, 1)
    #Logrank
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    pvalue[i]<-1 - pchisq(surv$chisq, 1)
  }
  
  #Logrank
  for (j in 1:M){
    
    if (pvalue[j]<0.05) {pvalue[j]=1}
    else {
      pvalue[j]=0
    }
    
  }
  
  mini1s<-mean(pvalue)
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvalue1[j]<0.05) {pvalue1[j]=1}
    else {
      pvalue1[j]=0
    }
    
  }
  
  mini1s1<-mean(pvalue1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvalue2[j]<0.05) {pvalue2[j]=1}
    else {
      pvalue2[j]=0
    }
    
  }
  
  mini1s2<-mean(pvalue2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvalue3[j]<0.05) {pvalue3[j]=1}
    else {
      pvalue3[j]=0
    }
    
  }
  
  mini1s3<-mean(pvalue3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvalue4[j]<0.05) {pvalue4[j]=1}
    else {
      pvalue4[j]=0
    }
    
  }
  
  mini1s4<-mean(pvalue4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvalue5[j]<0.05) {pvalue5[j]=1}
    else {
      pvalue5[j]=0
    }
    
  }
  
  mini1s5<-mean(pvalue5)

  
  ##############################MSE################################################################
  
  MSE<-numeric(M)
  MSE1<-numeric(M)
  MSE2<-numeric(M)
  MSE3<-numeric(M)
  MSE4<-numeric(M)
  MSE5<-numeric(M)
  
  #No stratification
  for (j in 1:M){
    MSE[j]<-(coeff[j]-HR)^2
  }
  
  MSEfinalmini10<-mean(MSE)
  
  #Stratification with 1 factor
  
  for (j in 1:M){
    MSE1[j]<-(coeff1[j]-HR)^2
  }
  
  MSEfinalmini11<-mean(MSE1)
  
  #Stratification with 3 factors
  
  for (j in 1:M){
    MSE2[j]<-(coeff2[j]-HR)^2
  }
  
  MSEfinalmini12<-mean(MSE2)
  
  #Stratification with 5 factors
  
  for (j in 1:M){
    MSE3[j]<-(coeff3[j]-HR)^2
  }
  
  MSEfinalmini13<-mean(MSE3)
  
  #Stratification with the score
  
  for (j in 1:M){
    MSE4[j]<-(coeff4[j]-HR)^2
  }
  
  MSEfinalmini14<-mean(MSE4)
  
  #Stratification with the score + 5 factors + the center
  
  for (j in 1:M){
    MSE5[j]<-(coeff5[j]-HR)^2
  }
  
  MSEfinalmini15<-mean(MSE5)
    
  ###############################################################################################
  ##############################Rerandomization##################################################
  ###############################################################################################
  
  ##############################Size of the test#################################################
  #Logrank  
  
  set.seed(4444)
  
  pval.sim=stat.sim=NULL
  
  pvalrerandotest<-rep( 0, len=M)
  pvalrerando<-numeric(M)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando[l]<0.05) {pvalrerandotest[l]=1}
    else {pvalrerandotest[l]=0}
  }
  
  
  pvalrerandmini10<-mean(pvalrerandotest)
  pvalrerandmini10
  
  
  #Stratified Logrank (1 factor)
  
  set.seed(4444)
  pval.sim=stat.sim=NULL
  
  pvalrerando1<-numeric(M)
  pvalrerandotest1<-rep( 0, len=M)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando1[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando1[l]<0.05) {pvalrerandotest1[l]=1}
    else {pvalrerandotest1[l]=0}
  }
  
  
  pvalrerandmini11<-mean(pvalrerandotest1)
  pvalrerandmini11
  
  #Stratified Logrank (3 factors)
  
  set.seed(4444)
  
  pvalrerando2<-numeric(M)
  pvalrerandotest2<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando2[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando2[l]<0.05) {pvalrerandotest2[l]=1}
    else {pvalrerandotest2[l]=0}
  }
  
  
  pvalrerandmini12<-mean(pvalrerandotest2)
  pvalrerandmini12
  
  #Stratified Logrank (5 factors)
  
  
  set.seed(4444)
  
  pvalrerando3<-numeric(M)
  pvalrerandotest3<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando3[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando3[l]<0.05) {pvalrerandotest3[l]=1}
    else {pvalrerandotest3[l]=0}
  }
  
  
  pvalrerandmini13<-mean(pvalrerandotest3)
  pvalrerandmini13
  
  #Stratified Logrank (Score)
  
  set.seed(4444)
  
  
  rerando4<-rep( 0, len=nsimu)
  pvalrerando4<-numeric(M)
  pvalrerandotest4<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando4[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando4[l]<0.05) {pvalrerandotest4[l]=1}
    else {pvalrerandotest4[l]=0}
  }
  
  
  pvalrerandmini14<-mean(pvalrerandotest4)
  pvalrerandmini14
  
  #Stratified Logrank (5 factors + Score + Center)
  
  set.seed(4444)
  
  pvalrerando5<-numeric(M)
  pvalrerandotest5<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE)+strata(randogroup$SCORE)+strata(randogroup$CENTER),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando5[l]<0.05) {pvalrerandotest5[l]=1}
    else {pvalrerandotest5[l]=0}
  }
  
  
  pvalrerandmini15<-mean(pvalrerandotest5)

  ##############################Power of the test################################################
  
  #Logrank  
  
  set.seed(4452)
  
  
  pvalrerandopow<-numeric(M)
  pvalrerando<-numeric(M)

  pvalrerandotestpow<-numeric(M)

  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow[i]<-a/nsimu 
  }
  
  for (l in 1:M){
    if (pvalrerandopow[l]<0.05) {pvalrerandotestpow[l]=1}
    else {pvalrerandotestpow[l]=0}
  }
  
  
  pvalrerandmini1pow<-mean(pvalrerandotestpow)
  pvalrerandmini1pow
  
  
  
  #Stratified Logrank (1 factor)
  
  pvalrerandopow1<-numeric(M)
  pvalrerandotestpow1<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow1[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow1[l]<0.05) {pvalrerandotestpow1[l]=1}
    else {pvalrerandotestpow1[l]=0}
  }
  
  
  pvalrerandmini1pow1<-mean(pvalrerandotestpow1)
  pvalrerandmini1pow1
  
  
  #Stratified Logrank (3 factors)
  
  pvalrerandopow2<-numeric(M)
  pvalrerandotestpow2<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow2[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow2[l]<0.05) {pvalrerandotestpow2[l]=1}
    else {pvalrerandotestpow2[l]=0}
  }
  
  
  pvalrerandmini1pow2<-mean(pvalrerandotestpow2)
  pvalrerandmini1pow2
  
  #Stratified Logrank (5 factors)
  
  pvalrerandopow3<-numeric(M)
  pvalrerandotestpow3<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow3[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow3[l]<0.05) {pvalrerandotestpow3[l]=1}
    else {pvalrerandotestpow3[l]=0}
  }
  
  
  pvalrerandmini1pow3<-mean(pvalrerandotestpow3)
  pvalrerandmini1pow3
  
  #Stratified Logrank (Score)
  
  pvalrerandopow4<-numeric(M)
  pvalrerandotestpow4<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
       else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow4[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow4[l]<0.05) {pvalrerandotestpow4[l]=1}
    else {pvalrerandotestpow4[l]=0}
  }
  
  
  pvalrerandmini1pow4<-mean(pvalrerandotestpow4)
  pvalrerandmini1pow4
  
  
  
  #Stratified Logrank (Score + 5 factors + Center)
  
  pvalrerandopow5<-numeric(M)
  pvalrerandotestpow5<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
     }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE)+strata(ovariansimple$CENTER)+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
       else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$CENTER)+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow5[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow5[l]<0.05) {pvalrerandotestpow5[l]=1}
    else {pvalrerandotestpow5[l]=0}
  }
  
  
  pvalrerandmini1pow5<-mean(pvalrerandotestpow5)
  pvalrerandmini1pow5
  
  ##################################################################################################
  #################################Minimization with the center and 3 factors########################
  ##################################################################################################
  
  ##############################Power of the test################################################
  ovariansimple<-matrix(nrow=N,ncol=12)
  pvalue<-numeric(M)
  pvalue1<-numeric(M)
  pvalue2<-numeric(M)
  pvalue3<-numeric(M)
  pvalue4<-numeric(M)
  pvalue5<-numeric(M)
  pvaluepow<-numeric(M)
  pvaluepow1<-numeric(M)
  pvaluepow2<-numeric(M)
  pvaluepow3<-numeric(M)
  pvaluepow4<-numeric(M)
  pvaluepow5<-numeric(M)
  coeff<-numeric(M)
  coeff1<-numeric(M)
  coeff2<-numeric(M)
  coeff3<-numeric(M)
  coeff4<-numeric(M)
  coeff5<-numeric(M)
  set.seed(5239)
  for (i in 1:M){
    
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov2 = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov2 ~ 1)
    #Logrank
    surv2<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT,data=ovariansimple)
    pvaluepow[i]<-1 - pchisq(surv2$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT , data = ovariansimple )
    coeff[i]<- exp(cox2$coefficients)
    #Stratified Logrank (1 factor)
    strat21<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvaluepow1[i]<- 1 - pchisq(strat21$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA) , data = ovariansimple )
    coeff1[i]<- exp(cox2$coefficients)
    #Stratified Logrank (3 factors)
    strat22<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvaluepow2[i]<- 1 - pchisq(strat22$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS) , data = ovariansimple )
    coeff2[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors)
    strat23<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvaluepow3[i]<- 1 - pchisq(strat23$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE) , data = ovariansimple )
    coeff3[i]<- exp(cox2$coefficients)
    #Stratified Logrank (Score)
    strat24<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvaluepow4[i]<- 1 - pchisq(strat24$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE) , data = ovariansimple )
    coeff4[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors + Score + Center)
    strat25<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvaluepow5[i]<- 1 - pchisq(strat25$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER) , data = ovariansimple )
    coeff5[i]<- exp(cox2$coefficients)
    
  }
  
  #Logrank
  for (j in 1:M){
    
    if (pvaluepow[j]<0.05) {pvaluepow[j]=1}
    else {
      pvaluepow[j]=0
    }
    
  }
  
  mini3pow<-mean(pvaluepow)
  
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvaluepow1[j]<0.05) {pvaluepow1[j]=1}
    else {
      pvaluepow1[j]=0
    }
    
  }
  
  mini3pow1<-mean(pvaluepow1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvaluepow2[j]<0.05) {pvaluepow2[j]=1}
    else {
      pvaluepow2[j]=0
    }
    
  }
  
  mini3pow2<-mean(pvaluepow2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvaluepow3[j]<0.05) {pvaluepow3[j]=1}
    else {
      pvaluepow3[j]=0
    }
    
  }
  
  mini3pow3<-mean(pvaluepow3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvaluepow4[j]<0.05) {pvaluepow4[j]=1}
    else {
      pvaluepow4[j]=0
    }
    
  }
  
  mini3pow4<-mean(pvaluepow4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvaluepow5[j]<0.05) {pvaluepow5[j]=1}
    else {
      pvaluepow5[j]=0
    }
    
  }
  
  mini3pow5<-mean(pvaluepow5)

  ##############################Size of the test#################################################
  set.seed(5239)
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov ~ 1)
    #Stratified Logrank (1 factor)
    strat1<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvalue1[i]<- 1 - pchisq(strat1$chisq, 1)
    #Stratified Logrank (3 factors)
    strat2<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvalue2[i]<- 1 - pchisq(strat2$chisq, 1)
    #Stratified Logrank (5 factors)
    strat3<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvalue3[i]<- 1 - pchisq(strat3$chisq, 1)
    #Stratified Logrank (Score)
    strat4<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvalue4[i]<- 1 - pchisq(strat4$chisq, 1)
    #Stratified Logrank (5 factors + Score + Center)
    strat5<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvalue5[i]<- 1 - pchisq(strat5$chisq, 1)
    #Logrank
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    pvalue[i]<-1 - pchisq(surv$chisq, 1)
  }
  
  #Logrank
  for (j in 1:M){
    
    if (pvalue[j]<0.05) {pvalue[j]=1}
    else {
      pvalue[j]=0
    }
    
  }
  
  mini3s<-mean(pvalue)
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvalue1[j]<0.05) {pvalue1[j]=1}
    else {
      pvalue1[j]=0
    }
    
  }
  
  mini3s1<-mean(pvalue1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvalue2[j]<0.05) {pvalue2[j]=1}
    else {
      pvalue2[j]=0
    }
    
  }
  
  mini3s2<-mean(pvalue2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvalue3[j]<0.05) {pvalue3[j]=1}
    else {
      pvalue3[j]=0
    }
    
  }
  
  mini3s3<-mean(pvalue3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvalue4[j]<0.05) {pvalue4[j]=1}
    else {
      pvalue4[j]=0
    }
    
  }
  
  mini3s4<-mean(pvalue4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvalue5[j]<0.05) {pvalue5[j]=1}
    else {
      pvalue5[j]=0
    }
    
  }
  
  mini3s5<-mean(pvalue5)

  ##############################MSE################################################################
  
  MSE<-numeric(M)
  MSE1<-numeric(M)
  MSE2<-numeric(M)
  MSE3<-numeric(M)
  MSE4<-numeric(M)
  MSE5<-numeric(M)
  
  #No stratification
  for (j in 1:M){
    MSE[j]<-(coeff[j]-HR)^2
  }
  
  MSEfinalmini30<-mean(MSE)
  
  #Stratification with 1 factor
  
  for (j in 1:M){
    MSE1[j]<-(coeff1[j]-HR)^2
  }
  
  MSEfinalmini31<-mean(MSE1)
  
  #Stratification with 3 factors
  
  for (j in 1:M){
    MSE2[j]<-(coeff2[j]-HR)^2
  }
  
  MSEfinalmini32<-mean(MSE2)
  
  #Stratification with 5 factors
  
  for (j in 1:M){
    MSE3[j]<-(coeff3[j]-HR)^2
  }
  
  MSEfinalmini33<-mean(MSE3)
  
  #Stratification with the score
  
  for (j in 1:M){
    MSE4[j]<-(coeff4[j]-HR)^2
  }
  
  MSEfinalmini34<-mean(MSE4)
  
  #Stratification with the score + 5 factors + the center
  
  for (j in 1:M){
    MSE5[j]<-(coeff5[j]-HR)^2
  }
  
  MSEfinalmini35<-mean(MSE5)

  ###############################################################################################
  ##############################Rerandomization##################################################
  ###############################################################################################
  
  ##############################Size of the test#################################################
  #Logrank  
  
  set.seed(4444)
  
  pval.sim=stat.sim=NULL
  
  pvalrerandotest<-rep( 0, len=M)
  pvalrerando<-numeric(M)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando[l]<0.05) {pvalrerandotest[l]=1}
    else {pvalrerandotest[l]=0}
  }
  
  
  pvalrerandmini30<-mean(pvalrerandotest)
  pvalrerandmini30
  
  
  #Stratified Logrank (1 factor)
  
  set.seed(4444)
  pval.sim=stat.sim=NULL
  
  pvalrerando1<-numeric(M)
  pvalrerandotest1<-rep( 0, len=M)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando1[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando1[l]<0.05) {pvalrerandotest1[l]=1}
    else {pvalrerandotest1[l]=0}
  }
  
  
  pvalrerandmini31<-mean(pvalrerandotest1)
  pvalrerandmini31
  
  #Stratified Logrank (3 factors)
  
  set.seed(4444)
  
  pvalrerando2<-numeric(M)
  pvalrerandotest2<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando2[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando2[l]<0.05) {pvalrerandotest2[l]=1}
    else {pvalrerandotest2[l]=0}
  }
  
  
  pvalrerandmini32<-mean(pvalrerandotest2)
  pvalrerandmini32
  
  #Stratified Logrank (5 factors)
  
  
  set.seed(4444)
  
  pvalrerando3<-numeric(M)
  pvalrerandotest3<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando3[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando3[l]<0.05) {pvalrerandotest3[l]=1}
    else {pvalrerandotest3[l]=0}
  }
  
  
  pvalrerandmini33<-mean(pvalrerandotest3)
  pvalrerandmini33
  
  #Stratified Logrank (Score)
  
  set.seed(4444)
  
  
  rerando4<-rep( 0, len=nsimu)
  pvalrerando4<-numeric(M)
  pvalrerandotest4<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando4[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando4[l]<0.05) {pvalrerandotest4[l]=1}
    else {pvalrerandotest4[l]=0}
  }
  
  
  pvalrerandmini34<-mean(pvalrerandotest4)
  pvalrerandmini34
  
  #Stratified Logrank (5 factors + Score + Center)
  
  set.seed(4444)
  
  pvalrerando5<-numeric(M)
  pvalrerandotest5<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE)+strata(randogroup$SCORE)+strata(randogroup$CENTER),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando5[l]<0.05) {pvalrerandotest5[l]=1}
    else {pvalrerandotest5[l]=0}
  }
  
  
  pvalrerandmini35<-mean(pvalrerandotest5)
  pvalrerandmini35

  ##############################Power of the test################################################
  
  #Logrank  
  
  set.seed(4452)
  
  
  pvalrerandopow<-numeric(M)
  pvalrerando<-numeric(M)
  pvalrerandotestpow<-numeric(M)
 
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow[i]<-a/nsimu 
  }
  
  for (l in 1:M){
    if (pvalrerandopow[l]<0.05) {pvalrerandotestpow[l]=1}
    else {pvalrerandotestpow[l]=0}
  }
  
  
  pvalrerandmini3pow<-mean(pvalrerandotestpow)
  pvalrerandmini3pow
  
  
  
  #Stratified Logrank (1 factor)
  
  pvalrerandopow1<-numeric(M)
  pvalrerandotestpow1<-numeric(M)
  
  pval.sim=stat.sim=NULL
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow1[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow1[l]<0.05) {pvalrerandotestpow1[l]=1}
    else {pvalrerandotestpow1[l]=0}
  }
  
  
  pvalrerandmini3pow1<-mean(pvalrerandotestpow1)
  pvalrerandmini3pow1
  
  
  #Stratified Logrank (3 factors)
  
  pvalrerandopow2<-numeric(M)
  pvalrerandotestpow2<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow2[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow2[l]<0.05) {pvalrerandotestpow2[l]=1}
    else {pvalrerandotestpow2[l]=0}
  }
  
  
  pvalrerandmini3pow2<-mean(pvalrerandotestpow2)
  pvalrerandmini3pow2
  
  #Stratified Logrank (5 factors)
  
  pvalrerandopow3<-numeric(M)
  pvalrerandotestpow3<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow3[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow3[l]<0.05) {pvalrerandotestpow3[l]=1}
    else {pvalrerandotestpow3[l]=0}
  }
  
  
  pvalrerandmini3pow3<-mean(pvalrerandotestpow3)
  pvalrerandmini3pow3
  
  #Stratified Logrank (Score)
  
  pvalrerandopow4<-numeric(M)
  pvalrerandotestpow4<-numeric(M)
  
  pval.sim=stat.sim=NULL
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
       else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow4[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow4[l]<0.05) {pvalrerandotestpow4[l]=1}
    else {pvalrerandotestpow4[l]=0}
  }
  
  
  pvalrerandmini3pow4<-mean(pvalrerandotestpow4)
  pvalrerandmini3pow4
  
  
  
  #Stratified Logrank (Score + 5 factors + Center)
  
  pvalrerandopow5<-numeric(M)
  pvalrerandotestpow5<-numeric(M)
  
  pval.sim=stat.sim=NULL
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
     }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE)+strata(ovariansimple$CENTER)+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
       else {randogroup$TIME2[j]=randogroup$TIME[j]}
}
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$CENTER)+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow5[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow5[l]<0.05) {pvalrerandotestpow5[l]=1}
    else {pvalrerandotestpow5[l]=0}
  }
  
  
  pvalrerandmini3pow5<-mean(pvalrerandotestpow5)
  pvalrerandmini3pow5

  
  ##################################################################################################
  #################################Minimization with the center and 5 factors########################
  ##################################################################################################
  
  ##############################Power of the test################################################
  ovariansimple<-matrix(nrow=N,ncol=12)
  pvalue<-numeric(M)
  pvalue1<-numeric(M)
  pvalue2<-numeric(M)
  pvalue3<-numeric(M)
  pvalue4<-numeric(M)
  pvalue5<-numeric(M)
  pvaluepow<-numeric(M)
  pvaluepow1<-numeric(M)
  pvaluepow2<-numeric(M)
  pvaluepow3<-numeric(M)
  pvaluepow4<-numeric(M)
  pvaluepow5<-numeric(M)
  coeff<-numeric(M)
  coeff1<-numeric(M)
  coeff2<-numeric(M)
  coeff3<-numeric(M)
  coeff4<-numeric(M)
  coeff5<-numeric(M)
  set.seed(5239)
  for (i in 1:M){
    
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov2 = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov2 ~ 1)
    #Logrank
    surv2<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT,data=ovariansimple)
    pvaluepow[i]<-1 - pchisq(surv2$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT , data = ovariansimple )
    coeff[i]<- exp(cox2$coefficients)
    #Stratified Logrank (1 factor)
    strat21<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvaluepow1[i]<- 1 - pchisq(strat21$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA) , data = ovariansimple )
    coeff1[i]<- exp(cox2$coefficients)
    #Stratified Logrank (3 factors)
    strat22<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvaluepow2[i]<- 1 - pchisq(strat22$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS) , data = ovariansimple )
    coeff2[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors)
    strat23<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvaluepow3[i]<- 1 - pchisq(strat23$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE) , data = ovariansimple )
    coeff3[i]<- exp(cox2$coefficients)
    #Stratified Logrank (Score)
    strat24<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvaluepow4[i]<- 1 - pchisq(strat24$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE) , data = ovariansimple )
    coeff4[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors + Score + Center)
    strat25<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvaluepow5[i]<- 1 - pchisq(strat25$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER) , data = ovariansimple )
    coeff5[i]<- exp(cox2$coefficients)
    
  }
  
  #Logrank
  for (j in 1:M){
    
    if (pvaluepow[j]<0.05) {pvaluepow[j]=1}
    else {
      pvaluepow[j]=0
    }
    
  }
  
  mini5pow<-mean(pvaluepow)
  
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvaluepow1[j]<0.05) {pvaluepow1[j]=1}
    else {
      pvaluepow1[j]=0
    }
    
  }
  
  mini5pow1-mean(pvaluepow1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvaluepow2[j]<0.05) {pvaluepow2[j]=1}
    else {
      pvaluepow2[j]=0
    }
    
  }
  
  mini5pow2<-mean(pvaluepow2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvaluepow3[j]<0.05) {pvaluepow3[j]=1}
    else {
      pvaluepow3[j]=0
    }
    
  }
  
  mini5pow3<-mean(pvaluepow3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvaluepow4[j]<0.05) {pvaluepow4[j]=1}
    else {
      pvaluepow4[j]=0
    }
    
  }
  
  mini5pow4<-mean(pvaluepow4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvaluepow5[j]<0.05) {pvaluepow5[j]=1}
    else {
      pvaluepow5[j]=0
    }
    
  }
  
  mini5pow5<-mean(pvaluepow5)
  
  ##############################Size of the test#################################################
  set.seed(5239)
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov ~ 1)
    #Stratified Logrank (1 factor)
    strat1<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvalue1[i]<- 1 - pchisq(strat1$chisq, 1)
    #Stratified Logrank (3 factors)
    strat2<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvalue2[i]<- 1 - pchisq(strat2$chisq, 1)
    #Stratified Logrank (5 factors)
    strat3<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvalue3[i]<- 1 - pchisq(strat3$chisq, 1)
    #Stratified Logrank (Score)
    strat4<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvalue4[i]<- 1 - pchisq(strat4$chisq, 1)
    #Stratified Logrank (5 factors + Score + Center)
    strat5<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvalue5[i]<- 1 - pchisq(strat5$chisq, 1)
    #Logrank
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    pvalue[i]<-1 - pchisq(surv$chisq, 1)
  }
  
  #Logrank
  for (j in 1:M){
    
    if (pvalue[j]<0.05) {pvalue[j]=1}
    else {
      pvalue[j]=0
    }
    
  }
  
  mini5s<-mean(pvalue)
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvalue1[j]<0.05) {pvalue1[j]=1}
    else {
      pvalue1[j]=0
    }
    
  }
  
  mini5s1<-mean(pvalue1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvalue2[j]<0.05) {pvalue2[j]=1}
    else {
      pvalue2[j]=0
    }
    
  }
  
  mini5s2<-mean(pvalue2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvalue3[j]<0.05) {pvalue3[j]=1}
    else {
      pvalue3[j]=0
    }
    
  }
  
  mini5s3<-mean(pvalue3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvalue4[j]<0.05) {pvalue4[j]=1}
    else {
      pvalue4[j]=0
    }
    
  }
  
  mini5s4<-mean(pvalue4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvalue5[j]<0.05) {pvalue5[j]=1}
    else {
      pvalue5[j]=0
    }
    
  }
  
  mini5s5<-mean(pvalue5)
  ##############################MSE###############################################################
  
  MSE<-numeric(M)
  MSE1<-numeric(M)
  MSE2<-numeric(M)
  MSE3<-numeric(M)
  MSE4<-numeric(M)
  MSE5<-numeric(M)
  
  #No stratification
  for (j in 1:M){
    MSE[j]<-(coeff[j]-HR)^2
  }
  
  MSEfinalmini50<-mean(MSE)
  
  #Stratification with 1 factor
  
  for (j in 1:M){
    MSE1[j]<-(coeff1[j]-HR)^2
  }
  
  MSEfinalmini51<-mean(MSE1)
  
  #Stratification with 3 factors
  
  for (j in 1:M){
    MSE2[j]<-(coeff2[j]-HR)^2
  }
  
  MSEfinalmini52<-mean(MSE2)
  
  #Stratification with 5 factors
  
  for (j in 1:M){
    MSE3[j]<-(coeff3[j]-HR)^2
  }
  
  MSEfinalmini53<-mean(MSE3)
  
  #Stratification with the score
  
  for (j in 1:M){
    MSE4[j]<-(coeff4[j]-HR)^2
  }
  
  MSEfinalmini54<-mean(MSE4)
  
  #Stratification with the score + 5 factors + the center
  
  for (j in 1:M){
    MSE5[j]<-(coeff5[j]-HR)^2
  }
  
  MSEfinalmini55<-mean(MSE5)

  ###############################################################################################
  ##############################Rerandomization##################################################
  ###############################################################################################
  
  ##############################Size of the test#################################################
  #Logrank  
  
  set.seed(4444)
  
  pval.sim=stat.sim=NULL
  
  pvalrerandotest<-rep( 0, len=M)
  pvalrerando<-numeric(M)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando[l]<0.05) {pvalrerandotest[l]=1}
    else {pvalrerandotest[l]=0}
  }
  
  
  pvalrerandmini50<-mean(pvalrerandotest)
  pvalrerandmini50
  
  
  #Stratified Logrank (1 factor)
  
  set.seed(4444)
  pval.sim=stat.sim=NULL
  
  pvalrerando1<-numeric(M)
  pvalrerandotest1<-rep( 0, len=M)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando1[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando1[l]<0.05) {pvalrerandotest1[l]=1}
    else {pvalrerandotest1[l]=0}
  }
  
  
  pvalrerandmini51<-mean(pvalrerandotest1)
  pvalrerandmini51
  
  #Stratified Logrank (3 factors)
  
  set.seed(4444)
  
  pvalrerando2<-numeric(M)
  pvalrerandotest2<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando2[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando2[l]<0.05) {pvalrerandotest2[l]=1}
    else {pvalrerandotest2[l]=0}
  }
  
  
  pvalrerandmini52<-mean(pvalrerandotest2)
  pvalrerandmini52
  
  #Stratified Logrank (5 factors)
  
  
  set.seed(4444)
  
  pvalrerando3<-numeric(M)
  pvalrerandotest3<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando3[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando3[l]<0.05) {pvalrerandotest3[l]=1}
    else {pvalrerandotest3[l]=0}
  }
  
  
  pvalrerandmini53<-mean(pvalrerandotest3)
  pvalrerandmini53
  
  #Stratified Logrank (Score)
  
  set.seed(4444)
  
  
  rerando4<-rep( 0, len=nsimu)
  pvalrerando4<-numeric(M)
  pvalrerandotest4<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando4[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando4[l]<0.05) {pvalrerandotest4[l]=1}
    else {pvalrerandotest4[l]=0}
  }
  
  
  pvalrerandmini54<-mean(pvalrerandotest4)
  pvalrerandmini54
  
  #Stratified Logrank (5 factors + Score + Center)
  
  set.seed(4444)
  
  pvalrerando5<-numeric(M)
  pvalrerandotest5<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE)+strata(randogroup$SCORE)+strata(randogroup$CENTER),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando5[l]<0.05) {pvalrerandotest5[l]=1}
    else {pvalrerandotest5[l]=0}
  }
  
  
  pvalrerandmini55<-mean(pvalrerandotest5)
  pvalrerandmini55
  
  ##############################Power of the test################################################
  
  #Logrank  
  
  set.seed(4452)
  
  
  pvalrerandopow<-numeric(M)
  pvalrerando<-numeric(M)

  pvalrerandotestpow<-numeric(M)

  
  pval.sim=stat.sim=NULL
   set.seed(4452)
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow[i]<-a/nsimu 
  }
  
  for (l in 1:M){
    if (pvalrerandopow[l]<0.05) {pvalrerandotestpow[l]=1}
    else {pvalrerandotestpow[l]=0}
  }
  
  
  pvalrerandmini5pow<-mean(pvalrerandotestpow)
  pvalrerandmini5pow
  
  
  
  #Stratified Logrank (1 factor)
  
  pvalrerandopow1<-numeric(M)
  pvalrerandotestpow1<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow1[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow1[l]<0.05) {pvalrerandotestpow1[l]=1}
    else {pvalrerandotestpow1[l]=0}
  }
  
  
  pvalrerandmini5pow1<-mean(pvalrerandotestpow1)
  pvalrerandmini5pow1
  
  
  #Stratified Logrank (3 factors)
  
  pvalrerandopow2<-numeric(M)
  pvalrerandotestpow2<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow2[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow2[l]<0.05) {pvalrerandotestpow2[l]=1}
    else {pvalrerandotestpow2[l]=0}
  }
  
  
  pvalrerandmini5pow2<-mean(pvalrerandotestpow2)
  pvalrerandmini5pow2
  
  #Stratified Logrank (5 factors)
  
  pvalrerandopow3<-numeric(M)
  pvalrerandotestpow3<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow3[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow3[l]<0.05) {pvalrerandotestpow3[l]=1}
    else {pvalrerandotestpow3[l]=0}
  }
  
  
  pvalrerandmini5pow3<-mean(pvalrerandotestpow3)
  pvalrerandmini5pow3
  
  #Stratified Logrank (Score)
  
  pvalrerandopow4<-numeric(M)
  pvalrerandotestpow4<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
     }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
       else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow4[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow4[l]<0.05) {pvalrerandotestpow4[l]=1}
    else {pvalrerandotestpow4[l]=0}
  }
  
  
  pvalrerandmini5pow4<-mean(pvalrerandotestpow4)
  pvalrerandmini5pow4
  
  
  
  #Stratified Logrank (Score + 5 factors + Center)
  
  pvalrerandopow5<-numeric(M)
  pvalrerandotestpow5<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
     }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE)+strata(ovariansimple$CENTER)+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","RESDISEA","AGE","PS","GRADE","STAGE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1,1,1,1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
       else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$CENTER)+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow5[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow5[l]<0.05) {pvalrerandotestpow5[l]=1}
    else {pvalrerandotestpow5[l]=0}
  }
  
  
  pvalrerandmini5pow5<-mean(pvalrerandotestpow5)
  pvalrerandmini5pow5
  
  
  ##################################################################################################
  #################################Minimization with the center and the score#######################
  ##################################################################################################
  
  ##############################Power of the test################################################
  ovariansimple<-matrix(nrow=N,ncol=12)
  pvalue<-numeric(M)
  pvalue1<-numeric(M)
  pvalue2<-numeric(M)
  pvalue3<-numeric(M)
  pvalue4<-numeric(M)
  pvalue5<-numeric(M)
  pvaluepow<-numeric(M)
  pvaluepow1<-numeric(M)
  pvaluepow2<-numeric(M)
  pvaluepow3<-numeric(M)
  pvaluepow4<-numeric(M)
  pvaluepow5<-numeric(M)
  coeff<-numeric(M)
  coeff1<-numeric(M)
  coeff2<-numeric(M)
  coeff3<-numeric(M)
  coeff4<-numeric(M)
  coeff5<-numeric(M)
  set.seed(5239)
  for (i in 1:M){
    
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov2 = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov2 ~ 1)
    #Logrank
    surv2<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT,data=ovariansimple)
    pvaluepow[i]<-1 - pchisq(surv2$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT , data = ovariansimple )
    coeff[i]<- exp(cox2$coefficients)
    #Stratified Logrank (1 factor)
    strat21<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvaluepow1[i]<- 1 - pchisq(strat21$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA) , data = ovariansimple )
    coeff1[i]<- exp(cox2$coefficients)
    #Stratified Logrank (3 factors)
    strat22<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvaluepow2[i]<- 1 - pchisq(strat22$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS) , data = ovariansimple )
    coeff2[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors)
    strat23<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvaluepow3[i]<- 1 - pchisq(strat23$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE) , data = ovariansimple )
    coeff3[i]<- exp(cox2$coefficients)
    #Stratified Logrank (Score)
    strat24<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvaluepow4[i]<- 1 - pchisq(strat24$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$SCORE) , data = ovariansimple )
    coeff4[i]<- exp(cox2$coefficients)
    #Stratified Logrank (5 factors + Score + Center)
    strat25<-survdiff(formula=surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvaluepow5[i]<- 1 - pchisq(strat25$chisq, 1)
    #For the MSE
    cox2 <- coxph ( surv.ov2 ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER) , data = ovariansimple )
    coeff5[i]<- exp(cox2$coefficients)
    
  }
  
  #Logrank
  for (j in 1:M){
    
    if (pvaluepow[j]<0.05) {pvaluepow[j]=1}
    else {
      pvaluepow[j]=0
    }
    
  }
  
  miniscopow<-mean(pvaluepow)
  
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvaluepow1[j]<0.05) {pvaluepow1[j]=1}
    else {
      pvaluepow1[j]=0
    }
    
  }
  
  miniscopow1<-mean(pvaluepow1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvaluepow2[j]<0.05) {pvaluepow2[j]=1}
    else {
      pvaluepow2[j]=0
    }
    
  }
  
  miniscopow2<-mean(pvaluepow2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvaluepow3[j]<0.05) {pvaluepow3[j]=1}
    else {
      pvaluepow3[j]=0
    }
    
  }
  
  miniscopow3<-mean(pvaluepow3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvaluepow4[j]<0.05) {pvaluepow4[j]=1}
    else {
      pvaluepow4[j]=0
    }
    
  }
  
  miniscopow4<-mean(pvaluepow4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvaluepow5[j]<0.05) {pvaluepow5[j]=1}
    else {
      pvaluepow5[j]=0
    }
    
  }
  
  miniscopow5<-mean(pvaluepow5)
  
  ##############################Size of the test#################################################
  set.seed(5239)
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    fit = survfit ( surv.ov ~ 1)
    #Stratified Logrank (1 factor)
    strat1<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    pvalue1[i]<- 1 - pchisq(strat1$chisq, 1)
    #Stratified Logrank (3 factors)
    strat2<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    pvalue2[i]<- 1 - pchisq(strat2$chisq, 1)
    #Stratified Logrank (5 factors)
    strat3<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    pvalue3[i]<- 1 - pchisq(strat3$chisq, 1)
    #Stratified Logrank (Score)
    strat4<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    pvalue4[i]<- 1 - pchisq(strat4$chisq, 1)
    #Stratified Logrank (5 factors + Score + Center)
    strat5<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    pvalue5[i]<- 1 - pchisq(strat5$chisq, 1)
    #Logrank
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    pvalue[i]<-1 - pchisq(surv$chisq, 1)
  }
  
  #Logrank
  for (j in 1:M){
    
    if (pvalue[j]<0.05) {pvalue[j]=1}
    else {
      pvalue[j]=0
    }
    
  }
  
  miniscos<-mean(pvalue)
  #Stratified Logrank (1 factor)
  for (j in 1:M){
    
    if (pvalue1[j]<0.05) {pvalue1[j]=1}
    else {
      pvalue1[j]=0
    }
    
  }
  
  miniscos1<-mean(pvalue1)
  
  #Stratified Logrank (3 factors) 
  for (j in 1:M){
    
    if (pvalue2[j]<0.05) {pvalue2[j]=1}
    else {
      pvalue2[j]=0
    }
    
  }
  
  miniscos2<-mean(pvalue2)
  
  #Stratified Logrank (5 factors)
  for (j in 1:M){
    
    if (pvalue3[j]<0.05) {pvalue3[j]=1}
    else {
      pvalue3[j]=0
    }
    
  }
  
  miniscos3<-mean(pvalue3)
  
  #Stratified Logrank (Score)
  for (j in 1:M){
    
    if (pvalue4[j]<0.05) {pvalue4[j]=1}
    else {
      pvalue4[j]=0
    }
    
  }
  
  miniscos4<-mean(pvalue4)
  
  #Stratified Logrank (5 factors+ Score + Center)
  for (j in 1:M){
    
    if (pvalue5[j]<0.05) {pvalue5[j]=1}
    else {
      pvalue5[j]=0
    }
    
  }
  
  miniscos5<-mean(pvalue5)
  
  ##############################MSE################################################################
  
  MSE<-numeric(M)
  MSE1<-numeric(M)
  MSE2<-numeric(M)
  MSE3<-numeric(M)
  MSE4<-numeric(M)
  MSE5<-numeric(M)
  
  #No stratification
  for (j in 1:M){
    MSE[j]<-(coeff[j]-HR)^2
  }
  
  MSEfinalminisco<-mean(MSE)
  
  #Stratification with 1 factor
  
  for (j in 1:M){
    MSE1[j]<-(coeff1[j]-HR)^2
  }
  
  MSEfinalminisco1<-mean(MSE1)
  
  #Stratification with 3 factors
  
  for (j in 1:M){
    MSE2[j]<-(coeff2[j]-HR)^2
  }
  
  MSEfinalminisco2<-mean(MSE2)
  
  #Stratification with 5 factors
  
  for (j in 1:M){
    MSE3[j]<-(coeff3[j]-HR)^2
  }
  
  MSEfinalminisco3<-mean(MSE3)
  
  #Stratification with the score
  
  for (j in 1:M){
    MSE4[j]<-(coeff4[j]-HR)^2
  }
  
  MSEfinalminisco4<-mean(MSE4)
  
  #Stratification with the score + 5 factors + the center
  
  for (j in 1:M){
    MSE5[j]<-(coeff5[j]-HR)^2
  }
  
  MSEfinalminisco5<-mean(MSE5)
  
  ###############################################################################################
  ##############################Rerandomization##################################################
  ###############################################################################################
  
  ##############################Size of the test#################################################
  #Logrank  
  
  set.seed(4444)
  
  pval.sim=stat.sim=NULL
  
  pvalrerandotest<-rep( 0, len=M)
  pvalrerando<-numeric(M)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando[l]<0.05) {pvalrerandotest[l]=1}
    else {pvalrerandotest[l]=0}
  }
  
  
  pvalrerandminisco<-mean(pvalrerandotest)
  
  
  #Stratified Logrank (1 factor)
  
  set.seed(4444)
  pval.sim=stat.sim=NULL
  
  pvalrerando1<-numeric(M)
  pvalrerandotest1<-rep( 0, len=M)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando1[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando1[l]<0.05) {pvalrerandotest1[l]=1}
    else {pvalrerandotest1[l]=0}
  }
  
  
  pvalrerandminisco1<-mean(pvalrerandotest1)
  
  #Stratified Logrank (3 factors)
  
  set.seed(4444)
  
  pvalrerando2<-numeric(M)
  pvalrerandotest2<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando2[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando2[l]<0.05) {pvalrerandotest2[l]=1}
    else {pvalrerandotest2[l]=0}
  }
  
  
  pvalrerandminisco2<-mean(pvalrerandotest2)
  
  #Stratified Logrank (5 factors)
  
  
  set.seed(4444)
  
  pvalrerando3<-numeric(M)
  pvalrerandotest3<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando3[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando3[l]<0.05) {pvalrerandotest3[l]=1}
    else {pvalrerandotest3[l]=0}
  }
  
  
  pvalrerandminisco3<-mean(pvalrerandotest3)
  
  #Stratified Logrank (Score)
  
  set.seed(4444)
  
  
  rerando4<-rep( 0, len=nsimu)
  pvalrerando4<-numeric(M)
  pvalrerandotest4<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando4[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerando4[l]<0.05) {pvalrerandotest4[l]=1}
    else {pvalrerandotest4[l]=0}
  }
  
  
  pvalrerandminisco4<-mean(pvalrerandotest4)
  
  #Stratified Logrank (5 factors + Score + Center)
  
  set.seed(4444)
  
  pvalrerando5<-numeric(M)
  pvalrerandotest5<-rep( 0, len=M)
  
  pval.sim=stat.sim=NULL
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    surv.ov = Surv( time = ovariansimple$TIME , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE)+strata(ovariansimple$SCORE)+strata(ovariansimple$CENTER),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (k in 1:nsimu){
      randogroup[,1+k]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    surv.ovrerand2 = Surv( time = randogroup$TIME , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      surv2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$GRADE)+strata(randogroup$STAGE)+strata(randogroup$SCORE)+strata(randogroup$CENTER),data=randogroup)
      p<-1 - pchisq(surv2$chisq, 1)
      s<-surv2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerando5[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerando5[l]<0.05) {pvalrerandotest5[l]=1}
    else {pvalrerandotest5[l]=0}
  }
  
  
  pvalrerandminisco5<-mean(pvalrerandotest5)

  ##############################Power of the test################################################
  
  #Logrank  
  
  set.seed(4452)
  
  
  pvalrerandopow<-numeric(M)
  pvalrerando<-numeric(M)

  pvalrerandotestpow<-numeric(M)
 
  
  pval.sim=stat.sim=NULL
  set.seed(4452)
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT,data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l],data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow[i]<-a/nsimu 
  }
  
  for (l in 1:M){
    if (pvalrerandopow[l]<0.05) {pvalrerandotestpow[l]=1}
    else {pvalrerandotestpow[l]=0}
  }
  
  
  pvalrerandminiscopow<-mean(pvalrerandotestpow)
  
  
  
  #Stratified Logrank (1 factor)
  
  pvalrerandopow1<-numeric(M)
  pvalrerandotestpow1<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow1[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow1[l]<0.05) {pvalrerandotestpow1[l]=1}
    else {pvalrerandotestpow1[l]=0}
  }
  
  
  pvalrerandminiscopow1<-mean(pvalrerandotestpow1)
  
  
  #Stratified Logrank (3 factors)
  
  pvalrerandopow2<-numeric(M)
  pvalrerandotestpow2<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow2[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow2[l]<0.05) {pvalrerandotestpow2[l]=1}
    else {pvalrerandotestpow2[l]=0}
  }
  
  
  pvalrerandminiscopow2<-mean(pvalrerandotestpow2)
  
  #Stratified Logrank (5 factors)
  
  pvalrerandopow3<-numeric(M)
  pvalrerandotestpow3<-numeric(M)
  
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
    }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$GRADE)+strata(ovariansimple$STAGE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow3[i]<-a/nsimu 

  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow3[l]<0.05) {pvalrerandotestpow3[l]=1}
    else {pvalrerandotestpow3[l]=0}
  }
  
  
  pvalrerandminiscopow3<-mean(pvalrerandotestpow3)
  
  #Stratified Logrank (Score)
  
  pvalrerandopow4<-numeric(M)
  pvalrerandotestpow4<-numeric(M)
  
  pval.sim=stat.sim=NULL
  
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
      }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
        else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow4[i]<-a/nsimu 
  }
  
  
  
  for (l in 1:M){
    if (pvalrerandopow4[l]<0.05) {pvalrerandotestpow4[l]=1}
    else {pvalrerandotestpow4[l]=0}
  }
  
  
  pvalrerandminiscopow4<-mean(pvalrerandotestpow4)
  
  
  
  #Stratified Logrank (Score + 5 factors + Center)
  
  pvalrerandopow5<-numeric(M)
  pvalrerandotestpow5<-numeric(M)
  
  pval.sim=stat.sim=NULL
  set.seed(4444)
  
  for (i in 1:M){
    ovariansimple<-ovarian[sample(1:nrow(ovarian),N),]
    ovariansimple<-ovariansimple[,1:12]
    TRT<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    ovariansimple<-cbind(ovariansimple,TRT)
    for (j in 1:N){
      if (ovariansimple$TRT[j]<1.9) {ovariansimple$TIME2[j]=HR*ovariansimple$TIME2[j]}
     }
    surv.ov = Surv( time = ovariansimple$TIME2 , event = ovariansimple$STATUS)
    surv<-survdiff(formula=surv.ov ~ ovariansimple$TRT+strata(ovariansimple$RESDISEA)+strata(ovariansimple$AGE)+strata(ovariansimple$PS)+strata(ovariansimple$STAGE)+strata(ovariansimple$GRADE)+strata(ovariansimple$CENTER)+strata(ovariansimple$SCORE),data=ovariansimple)
    p.obs<-1 - pchisq(surv$chisq, 1)
    s.obs<-surv$chisq
    ovariansimple<-ovariansimple[,1:12]
    
    randogroup<-data.frame(subject.id=ovariansimple[,"ï..RAND"])
    for (m in 1:nsimu){
      randogroup[,1+m]<-miniAlgo(dataset=ovariansimple,strata=c("CENTER","SCORE"),nbgp=2,determinism.percent=0.8,rando.ratio=c(1,1),weight=c(1,1))
    }
    colnames(randogroup)=c("subject.id",paste("r",1:nsimu,sep=""))
    colnames(ovariansimple)[1]="subject.id"
    randogroup=merge(randogroup,ovariansimple,by="subject.id")
    for (j in 1:length(randogroup[,1])){
      for (k in 1:nsimu){
        if (randogroup[j,k+1]<1.9) {randogroup$TIME2[j]=HR*randogroup$TIME[j]}
       else {randogroup$TIME2[j]=randogroup$TIME[j]}
      }
    }
    surv.ovrerand2 = Surv( time = randogroup$TIME2 , event = randogroup$STATUS)
    
    for (l in 1:nsimu){
      
      survpow2<-survdiff(formula=surv.ovrerand2 ~ randogroup[,1+l]+strata(randogroup$RESDISEA)+strata(randogroup$AGE)+strata(randogroup$PS)+strata(randogroup$STAGE)+strata(randogroup$GRADE)+strata(randogroup$CENTER)+strata(randogroup$SCORE),data=randogroup)
      p<-1 - pchisq(survpow2$chisq, 1)
      s<-survpow2$chisq
      pval.sim=c(pval.sim,p)
      stat.sim=c(stat.sim,s)
      
    }
    sim=data.frame(pval.obs=rep(p.obs,nsimu),stat.obs=rep(s.obs,nsimu),pval.sim,stat.sim)
    a<-0
    for (j in 1 : nsimu)
    {
      if (sim[j,4]>sim[j,2])
      {
        a<-a+1
      }
    }
    pvalrerandopow5[i]<-a/nsimu 
  }
  
  
  for (l in 1:M){
    if (pvalrerandopow5[l]<0.05) {pvalrerandotestpow5[l]=1}
    else {pvalrerandotestpow5[l]=0}
  }
  
  
  pvalrerandminiscopow5<-mean(pvalrerandotestpow5)

   
