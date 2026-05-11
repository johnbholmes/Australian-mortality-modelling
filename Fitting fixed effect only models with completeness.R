#Fixed effect only model with no ties to be used for comparison with models
#with latent variables, 22 October 2025



#Inputs required are: 
#y:        	The matrix of death rates. Please store as class matrix.
#iter:     	The number of iterations after burnin.
#burnin:   	The number of iterations to be discarded.
#X          Matrix of predictor variables. Must include intercept. 
#alpha, gamma	Parameter values for inverse-gamma priors on error variances.    
#Agegender:    age-gender categories: A vector indicating which death categories are age-gender groups, not cause of death groups.
#y_incomplete: a matrix of incomplete death counts from 12 months prior to current release. This should include doctor certified deaths (in first column) and all age/gender categories.
#n_incomplete: a matrix of 'complete' counts from current release matching the times of occurrence of incomplete counts. This should include doctor certified deaths and all age/gender categories.
#y_recent:     a matrix of incomplete death counts from the most recent 12 months in the current release.
#covid.index:  an indicator for which columns in X are covid age-standardised death rates. Needed for dealing with incompleteness in recent death counts.
#flu.index
#coroner.index:an indicator for column in y is coroner age-standardised death rates.


#Outputs are:
#beta:	a list with entry i corresponding to the estimated of observed fixed coefficient for death category j at iteration i.
#beta.total:   a list with entry i corresponding to the estimate of observed fixed coefficient summed over all death categories at iteration i.


#####################################################################
#####################################################################
fixedeffectsonly.withcompleteness<-function(y,iter,burnin,X,alpha,gamma,y_incomplete,n_incomplete,y_recent,covid.index,flu.index,coroner.index,age.gender){

  library(BayesLogit)
  y0<-y #Store the original y (the last year will get updates for incompleteness)
  n<-dim(y)[1] #observations 
  k<-dim(y)[2] #death categories.
  p<-dim(X)[2] #number of fixed coefficients.

  #Processing for incompleteness
  N_incomplete<-dim(n_incomplete)[1]
  mtimes<-N_incomplete:1
  X_incomplete<-cbind(1,mtimes,mtimes^2)
  n_incomplete[y_incomplete>n_incomplete]<-y_incomplete[y_incomplete>n_incomplete]
  kappa<-y_incomplete-0.5*n_incomplete #modified based on the results of Polson.
  p_incomplete<-dim(X_incomplete)[2]
  c_incomplete<-dim(y_incomplete)[2]
  beta_0logit<-matrix(rnorm(p_incomplete*c_incomplete),p_incomplete,c_incomplete)
  Xkappa<-t(X_incomplete)%*%kappa
  link_incomplete<-X_incomplete%*%beta_0logit
  
  #get covid rates for correcting for incompleteness in last year, needed for updating coroner rates. 
  covid.X<-X[n+1-mtimes,covid.index]
  #get Flu rates for correcting for incompleteness in last year, needed for updating coroner rates. 
  Flu.X<-X[n+1-mtimes,flu.index]
  #Define predictor matrix to be updated
  Xuse <- X
  
  
  tau<-rgamma(k,alpha,gamma)        	  #Initial error precision
  tau[tau<1e-8]<-1e-8               #lower inital value to avoid singularity.
  
  #Initial storage of estimates.
  beta.store           <-list()
  tau.store            <-list()   
  yhat.store           <-list()
  yimp.store           <-list() 
  yhat.allcause           <-matrix(0,n,iter)
  yhat.allagegender       <-matrix(0,n,iter)
  iter <- iter + burnin
   
  for(i in 1:iter){
    
    #Part zero: running logistic regression to estimate the probability of death registration. #Should be done parallel but this is not implemented.
    for(li in 1:c_incomplete){
      omega<-rpg(N_incomplete,h=as.numeric(n_incomplete[,li]),z=link_incomplete[,li]) #generate pg r.v
      err1<-rnorm(N_incomplete)*sqrt(omega) #generate random noise
      Xerr<-t(X_incomplete)%*%err1 
      XTXomega<-t(X_incomplete)%*%diag(omega)%*%X_incomplete #XWX needed for variance
      beta_0logit[,li]<-solve(XTXomega,Xkappa[,li]+Xerr) #estimate beta for category li.
    }
    
    #Saving results.
    link_incomplete<-X_incomplete%*%beta_0logit #estimate link
    prob_incomplete<-(1+exp(-link_incomplete))^(-1) #save probability estimates.  
    #Use probability to generate 'missing number of deaths'. #note this means the n is now the y
    missing.dead<-matrix(rnbinom(n=prod(dim(y_recent)),size=y_recent,prob=prob_incomplete),N_incomplete,c_incomplete)
    #Determine multiplicative factor for rates. (All doctor certified deaths (used for causes), and age-gender categories)
    mult.dead<-(missing.dead+y_recent)/y_recent
    #Determine multiplicative factor for rates for all rates.
    mult.dead.all<-rowSums(missing.dead[,-1]+y_recent[,-1])/rowSums(y_recent[,-1])
    #Construct all death rate with filling in for incompleteness
    all.death.rates<-rowSums(y0[n+1-mtimes,-c(age.gender)])+rowSums(covid.X)+rowSums(Flu.X)
    
    #Update covid rates. in last year
    Xuse[n+1-mtimes,covid.index]<-covid.X*mult.dead[,1]
    #Update covid rates. in last year
    Xuse[n+1-mtimes,flu.index]<-Flu.X*mult.dead[,1]
    
    #Update rates by all cause but coroner
    y[n+1-mtimes,-c(coroner.index,age.gender)]<-y0[n+1-mtimes,-c(coroner.index,age.gender)]*mult.dead[,1]  #update rates for doctor certified deaths by cause.
    #Update rate by all age-gender specific group
    y[n+1-mtimes,c(age.gender)]<-y0[n+1-mtimes,c(age.gender)]*mult.dead[,-1] #update age specific contribution to rate by age/gender group.
    #update coroner death rate. Note y has been updated already for doctor certified deaths at this point and so has covid.rates in X.
    y[n+1-mtimes,c(coroner.index)]<-all.death.rates*mult.dead.all-rowSums(y[n+1-mtimes,-c(coroner.index,age.gender)])-rowSums(Xuse[n+1-mtimes,c(flu.index,covid.index)])
    
    
    #############################################
    #Now run the actual Latent variable model for estimate processes affecting death rates.
    #Fit linear model
    
    #Estimate betas. 
    errmat<-matrix(rnorm(n*k),k,n)
    errmat<-errmat/sqrt(tau)
    errmat<-t(errmat)
    p1b<-crossprod(Xuse,y+errmat)
    XTXuse<-crossprod(Xuse)
    beta<-solve(XTXuse,p1b)
    
    #Update taus.
    yhat     <-Xuse%*%beta
    resid.est<-y-yhat
    SSEs     <-colSums(resid.est^2)
    tau      <-rgamma(k,0.5*n,0.5*SSEs)
    
    yhat.cause    <-rowSums(yhat[,-c(age.gender)])
    yhat.agegender<-rowSums(yhat[,c(age.gender)])
    
    #storing results.
    beta.store[[max(1,i-burnin)]]           <-beta 
    tau.store[[max(1,i-burnin)]]            <-tau   
    yhat.store[[max(1,i-burnin)]]           <-yhat
    yimp.store[[max(1,i-burnin)]]           <-y[n+1-mtimes,] 
    yhat.allcause[,max(1,i-burnin)]            <-yhat.cause
    yhat.allagegender[,max(1,i-burnin)]        <-yhat.agegender
  }
  
  x<-list(beta.store,tau.store,yhat.store,yimp.store,yhat.allcause,yhat.allagegender)
  names(x)<-c('beta','tau','yhat','yimp','yhat.allcause','yhat.allagegender')
  
  return(x)
}

#Flu included as predictor. split by pandemic to account for ascertainment changes
#This is used to confirm latent variables are needed.
system.time(test11<-fixedeffectsonly.withcompleteness(y=as.matrix(Deathcategory),iter=15000,burnin=15000,X=myX,alpha=0.1,gamma=0.1,y_incomplete=y_incomplete,n_incomplete=n_incomplete,y_recent=y_recent,covid.index=18:20,flu.index=16:17,coroner.index=11,age.gender=12:21))
system.time(test11.2<-fixedeffectsonly.withcompleteness(y=as.matrix(Deathcategory),iter=15000,burnin=15000,X=myX,alpha=0.1,gamma=0.1,y_incomplete=y_incomplete,n_incomplete=n_incomplete,y_recent=y_recent,covid.index=18:20,flu.index=16:17,coroner.index=11,age.gender=12:21))
system.time(test11.3<-fixedeffectsonly.withcompleteness(y=as.matrix(Deathcategory),iter=15000,burnin=15000,X=myX,alpha=0.1,gamma=0.1,y_incomplete=y_incomplete,n_incomplete=n_incomplete,y_recent=y_recent,covid.index=18:20,flu.index=16:17,coroner.index=11,age.gender=12:21))

save(test11,file='Storing results_fixedeffectonlyrun1May92026.RData')
rm(test11)
save(test11.2,file='Storing results_fixedeffectonlyrun2May92026.RData')
rm(test11.2)
save(test11.3,file='Storing results_fixedeffectonlyrun3May92026.RData')
rm(test11.3)

#Looking at r^2 for cause of death and age gender of death 
#separately before start of pandemic behavioural change and after behavioural change.

#Cor before 2020

#Focus only on post 2020 data.
ind.to.use2<-min(which(year > 2020.25)):(n-52)
rFluFEpost2020<-cor(rowSums(Deathcategory[ind.to.use2,1:11]),test11$yhat.allcause[ind.to.use,])
rFluFEpost2020a<-cor(rowSums(Deathcategory[ind.to.use2,12:21]),test11$yhat.allagegender[ind.to.use,])
#Focus on pre 2020 data.
ind.to.use<-1:(min(which(year > 2020.25))-1)
rFluFEpre2020<-cor(rowSums(Deathcategory[ind.to.use,1:11]),test11$yhat.allcause[ind.to.use,])
rFluFEpre2020a<-cor(rowSums(Deathcategory[ind.to.use,12:21]),test11$yhat.allagegender[ind.to.use,])

#Check eigenvalues of covariance/correlation matrix
trying<-Reduce('+',test11$yhat)
testcov<-cov(Deathcategory[ind.to.use,]-trying[ind.to.use,]/15000)
plot(eigen(cov2cor(testcov))$val)

load('Storing results_Oct15JanknotwithPneu_replace_FluPneu_2025run1.RData')

#Repeat this with models including latent variables.


