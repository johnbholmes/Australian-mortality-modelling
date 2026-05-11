#Modified to tie seasonal uplift and baseline curve (for non cancer, non coronial deaths) +tie on 17 April 2023.
#Completed 28 June 2023
#Modified on 22 July 2023 to include a correction for incomplete counts in the last year.
#Aim restrict the latent variable to the period where a Respiratory season is declared.
#plus a second year round latent variable to capture (ideally) temperature related processes.

#Running a Gibbs sampler for multivariable
#regression with one latent variable.

#To be used for age-standardised mortality data.


#Inputs required are: 
#y:        	The matrix of death rates. Please store as class matrix.
#fs:        Number of seasonal latent variable, fixed to 1.
#fy:        Number of year round latent variables
#iter:     	The number of iterations after burnin.
#burnin:   	The number of iterations to be discarded.
#X          Matrix of predictor variables. Must include intercept. 
#alphaL, gammaL  Parameter values for inverse-gamma priors on loading variances.
#alpha, gamma	Parameter values for inverse-gamma priors on error variances.    
#ind.epi   Indices of model matrix corresponding to epidemic variable (winter curves)
#ind.base  Indices of categories where the baseline curve should be tied together for estimation.
#ind.to.tie Variables where baseline curve is tied. 
#Agegender:    age-gender categories: A vector indicating which death categories are age-gender groups, not cause of death groups.
#y_incomplete: a matrix of incomplete death counts from 12 months prior to current release. This should include doctor certified deaths (in first column) and all age/gender categories.
#n_incomplete: a matrix of 'complete' counts from current release matching the times of occurrence of incomplete counts. This should include doctor certified deaths and all age/gender categories.
#y_recent:     a matrix of incomplete death counts from the most recent 12 months in the current release.
#covid.index:  an indicator for which columns in X are covid age-standardised death rates. Needed for dealing with incompleteness in recent death counts.
#coroner.index:an indicator for column in y is coroner age-standardised death rates.
#These will have ties as seasonal peaks differ by category and this is affecting stability of epidemic.
#There will be no seasonality for cancer or Coronor related deaths.


#Outputs are:
#Epidemic:      A list with entry i corresponding to the estimate of winter epidemic death rates in each category.
#var.unique:   	A list with entry i corresponding to the vector of estimated error variances at iteration i.
#Epidemic.total:	A list with entry i corresponding to the total epidemic death at week k at iteration i.
#Epi.frac:      A list with the estimated proportion of epidemic deaths occuring in each category.
#beta:	a list with entry i corresponding to the estimated of observed fixed coefficient for death category j at iteration i.
#beta.total:   a list with entry i corresponding to the estimate of observed fixed coefficient summed over all death categories at iteration i.


#####################################################################
#####################################################################

#Version Tie (Tie main effect of seasonal uplift by category (as timing of seasonality differs) Latent variable restricted to seasonal period:
MVR.fullTie2F.agegender.withcompleteness<-function(y,fs,fy,iter,burnin,X,alphaL,gammaL,alpha,gamma,ind.epi,ind.base,ind.to.tie,y_incomplete,n_incomplete,y_recent,covid.index,coroner.index,age.gender){
  Resp.ind<-X[,ind.epi[1]] #First epidemic indicator is always an intercept.
  Resp.ind<-as.numeric(Resp.ind)
  n.r     <-sum(Resp.ind) #number of observations that belong to the respiratory season.
  library(BayesLogit)
  y0<-y #Store the original y (the last year will get updates for incompleteness)
  n<-dim(y)[1] #observations 
  k<-dim(y)[2] #death categories.
  p<-dim(X)[2] #number of fixed coefficients.
  pse<-length(ind.epi) #number of seasonal coefficients. 
  xb.epi<-X[,ind.epi]%*%matrix(1,length(ind.epi),1) #Initial curve for epidemic effect.
  
  psb<-length(ind.to.tie) #number of baseline coefficient (one fixed at 1, two tied across categories)
  xb.base<-X[,ind.to.tie]%*%c(1,0.5,0.5) #Initial curve for 'true' seasonality.
  xb.base2<-X[,ind.to.tie]%*%c(0,0.5,0.5) #Initial curve for 'true' seasonality. in cancer, coroner and age-gender group.
  
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
  
  
  #xb.base will be dealt with separately.
  Xuse<-cbind(X[,-c(ind.epi,ind.to.tie)],xb.base,xb.epi) #X is modified for tying effects (All other causes)
  Xuse2<-cbind(X[,-c(ind.epi,ind.to.tie[-1])],xb.base2,xb.epi)   #X is modified for tying effects in baseline, but varying scale. (Cancer and Coroner, age gender)
  
  #extract indices which appear Xuse2 but not Xuse.
  ind.use1<-(1:p)[-c(ind.epi,ind.to.tie)] #index in 1
  ind.use2<-(1:p)[-c(ind.epi,ind.to.tie[-1])] #index in 2.
  ind.both<-which(ind.use2 %in% ind.use1)  #index of L2 that appears in L too.
  ind.diff<-(1:(p-pse))[-ind.both]         #index of L2 that does not appear in L.
  ind.covid1<-which(ind.use1 %in% covid.index) #indices of Xuse corresponding to covid rates.
  ind.covid2<-which(ind.use2 %in% covid.index) #indices of Xuse2 corresponding to covid rates.
  
  #seasonal uplift to vary by proportion across category. 
  
  coef.names<-colnames(X) #Names of coefficients.
  cat.names <-colnames(y) #Names of death categories.
  
  F<-matrix(rnorm(n*fs),n,fs)    		  #Initial Factors Seasonal
  F[Resp.ind %in% 0,]<-0
  Fy<-matrix(rnorm(n*fy),n,fy)    		#Initial Factors all year.
  pLv<-p-pse-psb+1+1            #Because pse and psb are replaced with one and one parameters in tied baseline categories.
  L<-matrix(0,fs+fy+pLv,k)    		  #Initial Lbeta	(Set to zero matrix.)	#This includes the proportions
  pLv2<-p-pse-(psb-1)+1+1       #Because pse is replaced with one parameter in non-tied baseline categories and sine, cosine fixed in position based on average.
  L2<-matrix(0,fs+fy+pLv2,k)    		#Initial Lbeta	(Set to zero matrix.)	#This includes the proportions
  
  
  
  tau<-rgamma(k,alpha,gamma)        	  #Initial error precision
  tau[tau<1e-8]<-1e-8               #lower inital value to avoid singularity.
  XF <- cbind(Xuse,F,Fy)            #The initial combining of X and F. for all other natural causes doctor certified.
  XF2 <- cbind(Xuse2,F,Fy)          #The initial combining of X and F. for Cancer and Coroner
  FL <-XF%*%L                    #Storing E(Y) in tied baseline categories.
  FL2 <-XF2%*%L2                 #Storing E(Y) in non-tied baseline categories.
  
  #Initial storage of estimates.
  Yearround.noise.cause     <-list()
  Yearround.noise.agegender <-list()
  Yearround.noise.total.cause    <-list()
  Yearround.noise.total.agegender<-list()
  Yearround.noise.frac     <-list()
  Epidemic.noise.cause     <-list()
  Epidemic.total.cause     <-list()
  Epidemic.noise.agegender <-list()
  Epidemic.total.agegender <-list()
  var.unique   <-list()
  Epidemic.noise.total.cause    <-list()
  Epidemic.all.total.cause      <-list()
  Epidemic.noise.total.agegender<-list()
  Epidemic.all.total.agegender  <-list()
  Epi.noise.frac      <-list()
  Epi.all.frac        <-list()
  beta                <-list()  
  beta.total.cause    <-list()
  beta.total.agegender<-list()
  base.curve          <-list()  
  yimp                <-list()
  
  ds<-rgamma(fs,alphaL,gammaL)	#Initial precision for factor loadings.
  dy<-rgamma(fy,alphaL,gammaL)
  ds<-min(ds,1e-8)
  dy<-min(dy,1e-8)
  alphad<-alphaL+0.5*k			#alpha parameters.
  alpha <-alpha+0.5*n				#alpha parameter for error variance.
  dc<-rgamma(1,alphaL,gammaL)	#Initial precision for factor curve (there is only one) loadings.
  db<-rgamma(1,alphaL,gammaL)	#Initial precision for baseline curve (there is only one) loadings.
  dc<-min(dc,1e-8)
  db<-min(db,1e-8)
  
  iter <- iter + burnin
  p1a<-matrix(0,(p-pse-psb+1+1+fs+fy),k)  #A zero matrix to add necessary extra noise for the latent factor (in last column)
  p1a2<-matrix(0,(p-pse-(psb-1)+1+1+fs+fy),k)  #A zero matrix to add necessary extra noise for the latent factor (in last column) when updating L2.
  beta.base<-rep(0,psb-1)        #starting put for coefficients for sin, cos in baseline curve.
  yuseF<-matrix(0,n.r,k)                #A zero matrix for initial store of residualised y needed for updating F (This will be over-written again and again).
  yuseFy<-matrix(0,n,k)                #A zero matrix for initial store of residualised y needed for updating F (This will be over-written again and again).
  yuse<-matrix(0,n,k)                 #Same, but for epidemic coefficients.
  yuseb<-matrix(0,n,k)                 #Same, but for baseline coefficients.
  
  
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
    all.death.rates<-rowSums(y0[n+1-mtimes,-c(age.gender)])+rowSums(covid.X)
    
    #Update covid rates. in last year
    Xuse[n+1-mtimes,ind.covid1]<-covid.X*mult.dead[,1]
    Xuse2[n+1-mtimes,ind.covid2]<-covid.X*mult.dead[,1]
    #Update rates by all cause but coroner
    y[n+1-mtimes,-c(coroner.index,age.gender)]<-y0[n+1-mtimes,-c(coroner.index,age.gender)]*mult.dead[,1]  #update rates for doctor certified deaths by cause.
    #Update rate by all age-gender specific group
    y[n+1-mtimes,c(age.gender)]<-y0[n+1-mtimes,c(age.gender)]*mult.dead[,-1] #update age specific contribution to rate by age/gender group.
    #update coroner death rate. Note y has been updated already for doctor certified deaths at this point and so has covid.rates in X.
    y[n+1-mtimes,c(coroner.index)]<-all.death.rates*mult.dead.all-rowSums(y[n+1-mtimes,-c(coroner.index,age.gender)])-rowSums(Xuse[n+1-mtimes,ind.covid1])
    
    
    #############################################
    #Now run the actual Latent variable model for estimate processes affecting death rates.
    
    
    #Part one: Update L, beta category by category. Strictly speaking should 
    FFT<-crossprod(XF)
    p1a[pLv+1:fs,]<-sqrt(ds)*matrix(rnorm(fs*k),fs,k)
    p1a[pLv+fs+1:fy,]<-sqrt(dy)*matrix(rnorm(fy*k),fy,k)
    
    p1a[pLv,]<-sqrt(dc)*matrix(rnorm(1*k),1,k) #Add error component for curve loading
    p1a[pLv-1,]<-sqrt(db)*matrix(rnorm(1*k),1,k) #Add error component for base loading
    p1<-t(XF)%*%t(matrix(rnorm(n*k),k,n)*sqrt(tau))+p1a
    p2<-t(XF)%*%y%*%diag(tau)
    pL<-p1+p2
    #Categories with ties for baseline seasonality. (Respiratory conditions, Circulatory conditions, Dementia, Diabetes, other non-cancer conditions)
    for(j in ind.base){
      varLinv<-FFT*tau[j]
      diag(varLinv)[pLv+1:fs]<-ds+diag(varLinv)[pLv+1:fs] #Add on extra tau for Loading.
      diag(varLinv)[pLv+fs+1:fy]<-dy+diag(varLinv)[pLv+fs+1:fy] #Add on extra tau for Loading.
      diag(varLinv)[pLv]<-dc+diag(varLinv)[pLv] #Add on extra tau for curve Loading.
      diag(varLinv)[pLv-1]<-db+diag(varLinv)[pLv-1] #Add on extra tau for baseline Loading.
      L[,j]<-solve(varLinv,pL[,j])         		#Update beta, L.
    }
    
    #Categories without ties for baseline seasonality. (So Cancer, Coroner, age group)
    FFT2<-crossprod(XF2)
    p1a2[pLv2+1:fs,]<-sqrt(ds)*matrix(rnorm(fs*k),fs,k)
    p1a2[pLv2+fs+1:fy,]<-sqrt(dy)*matrix(rnorm(fy*k),fy,k)
    
    p1a2[pLv2,]<-sqrt(dc)*matrix(rnorm(1*k),1,k) #Add error component for curve loading
    p12<-t(XF2)%*%t(matrix(rnorm(n*k),k,n)*sqrt(tau))+p1a2
    p22<-t(XF2)%*%y%*%diag(tau)
    pL2<-p12+p22
    for(j  in (1:k)[-ind.base]){
      varLinv2<-FFT2*tau[j]
      diag(varLinv2)[pLv2+1:fs]<-ds+diag(varLinv2)[pLv2+1:fs] #Add on extra tau for Loading.
      diag(varLinv2)[pLv2+fs+1:fy]<-dy+diag(varLinv2)[pLv2+fs+1:fy] #Add on extra tau for Loading.
      diag(varLinv2)[pLv2]<-dc+diag(varLinv2)[pLv2] #Add on extra tau for curve Loading.
      diag(varLinv2)[pLv2-1]<-db+diag(varLinv2)[pLv2-1] #Add on extra tau for baseline Loading.
      L2[,j]<-solve(varLinv2,pL2[,j])         		#Update beta, L.
    }    
    
    #Fill in L for categories with non-tied baseline (Not practice just dropping intercept to avoid issues. In the saved results, the intercept will appear)
    L[,-ind.base]<-L2[-ind.to.tie[1],-ind.base] #Not this appears before the tied base and tied epidemic curve effects.
    #This is needed for updating F later.
    
    #Fill in empty rows in L2 to get consistent result.
    #Fill in L for categories with non-tied baseline (This is filling in saving X\beta_{baseline curve} as single value)
    #This is adding the intercept in, and filling the rest.
    L2[-ind.to.tie[1],ind.base]<-L[,ind.base]
    L2[ind.to.tie[1],ind.base]<-L[(pLv-1),ind.base]    
    
    #Update d. for the factor 
    sumL2<-rowSums(L[-c(1:pLv),]^2)
    gammad<-gammaL+0.5*sumL2
    d<-rgamma(fs+fy,alphad,gammad)
    ds<-d[1:fs]
    dy<-d[fs+1:fy]
    
    #Update dc (precision for the combined epidemic curve loading (mean factor)).
    sumL2c<-sum(L[pLv,]^2)
    gammadc<-gammaL+0.5*sumL2c
    dc<-rgamma(1,alphad,gammadc)
    
    #Update db (precision for the combined baseline curve loading for tied baseline curves only.
    sumL2b<-sum(L[pLv-1,]^2)
    gammadb<-gammaL+0.5*sumL2b
    db<-rgamma(1,alphad,gammadb)
    
    #Update F.
    #Update F for seasonal variable.
    #residualised y for tied baseline categories.
    yuseF[,ind.base]<-y[Resp.ind %in% 1,ind.base]-Xuse[Resp.ind %in% 1,]%*%L[1:pLv,ind.base]-Fy[Resp.ind %in% 1,]%*%t(L[pLv+fs+1:fy,ind.base]) #As fy =1
    #residualised y for non-tied baseline categories.
    yuseF[,-ind.base]<-y[Resp.ind %in% 1,-ind.base]-Xuse2[Resp.ind %in% 1,]%*%L2[1:pLv2,-ind.base]-Fy[Resp.ind %in% 1,]%*%t(L2[pLv2+fs+1:fy,-ind.base])
    p1F<-diag(fs)+t(L[pLv+1:fs,])%*%diag(tau)%*%(L[pLv+1:fs,])
    errF<-matrix(rnorm(fs*n.r),fs,n.r)+L[pLv+1:fs,]%*%diag(sqrt(tau))%*%matrix(rnorm(k*n.r),k,n.r)
    p2F<-L[pLv+1:fs,]%*%diag(tau)%*%t(yuseF)+errF
    F[Resp.ind %in% 1,]<-t(solve(p1F,p2F))		
    
    #Update F for year-round variable.
    #residualised y for tied baseline categories.
    yuseFy[,ind.base]<-y[,ind.base]-Xuse%*%L[1:pLv,ind.base]-F%*%t(L[pLv+1:fs,ind.base])
    #residualised y for non-tied baseline categories.
    yuseFy[,-ind.base]<-y[,-ind.base]-Xuse2%*%L2[1:pLv2,-ind.base]-F%*%t(L2[pLv2+1:fs,-ind.base])
    p1F<-diag(fy)+t(L[pLv+fs+1:fy,])%*%diag(tau)%*%(L[pLv+fs+1:fy,])
    errF<-matrix(rnorm(fy*n),fy,n)+L[pLv+fs+1:fy,]%*%diag(sqrt(tau))%*%matrix(rnorm(k*n),k,n)
    p2F<-L[pLv+fs+1:fy,]%*%diag(tau)%*%t(yuseFy)+errF
    Fy<-t(solve(p1F,p2F))		
    
    
    XF <- cbind(Xuse,F,Fy)                             #combine X with new F.   tied baseline categories.
    XF2 <- cbind(Xuse2,F,Fy)                           #combine X with new F. non-tied baseline categories.
    err<-rep(0,k)
    #squared error for tied baseline categories.
    err[ind.base]<-colSums((y[,ind.base]-XF%*%L[,ind.base])^2)
    #squared error for non-tied baseline categories.
    err[-ind.base]<-colSums((y[,-ind.base]-XF2%*%L2[,-ind.base])^2) 
    tau<-rgamma(k,alpha,gamma+0.5*err) 						#Update tau 
    
    #Estimate the parameters of the combined variables (mean seasonal curve).
    XTXepi<-crossprod(X[,ind.epi])
    pi.k<-L[pLv,] #Category component to curve.
    #residualised y for tied baseline categories.
    yuse[,ind.base]<-y[,ind.base]-XF[,-pLv]%*%L[-pLv,ind.base]
    #residualised y for non-tied baseline categories.
    yuse[,-ind.base]<-y[,-ind.base]-XF2[,-pLv2]%*%L2[-pLv2,-ind.base]
    
    tyuse.k<-t(X[,ind.epi])%*%(yuse%*%(pi.k*tau))
    err.k <-t(X[,ind.epi])%*%matrix(rnorm(n),n,1)*sqrt(sum((pi.k^2)*tau))+matrix(rnorm(pse),pse,1)*sqrt(2)
    beta.epi<-solve(XTXepi*sum((pi.k^2)*tau)+2*diag(pse),tyuse.k+err.k)#prior for beta variance is 0.5, as diag(XX') = 2 always.
    xb.epi<-X[,ind.epi]%*%beta.epi
    
    
    #Estimate the parameters of the combined variables (this is xb.base) (mean baseline curve, fixed intercept at 1 for chronic conditions, 0 for cancer,coroner, age group (as intercepts can differ)).
    Lb<-L2[pLv2-1,]
    XTXb<-crossprod(X[,ind.to.tie[-1]])
    yuseb[,ind.base]<-y[,ind.base]-XF[,-c(pLv-1)]%*%L[-c(pLv-1),ind.base]-matrix(Lb[ind.base],n,length(ind.base),byrow = TRUE) #main effect. Lb is main multiplier.
    yuseb[,-ind.base]<-y[,-ind.base]-XF2[,-c(pLv2-1)]%*%L2[-c(pLv2-1),-ind.base]#Intercept is built in here.
    tyuse.kb<-t(X[,ind.to.tie[-1]])%*%(yuseb%*%(Lb*tau))
    err.kb <-t(X[,ind.to.tie[-1]])%*%matrix(rnorm(n),n,1)*sqrt(sum((Lb^2)*tau))#+matrix(rnorm(psb-1),psb-1,1)*sqrt(1)
    beta.base<-solve(XTXb*sum((Lb^2)*tau),tyuse.kb+err.kb)#+1*diag(psb-1)prior for beta variance is 0.5, as diag(XX') = 1 (sin^2 + cos^2) always.
    xb.base<-X[,ind.to.tie]%*%c(1,beta.base) #For chronic conditions
    xb.base2<-X[,ind.to.tie[-1]]%*%c(beta.base) #For cancer, coroner, and age-gender groups.
    
    #For models where baseline curves are tied, save new XF.
    Xuse<-cbind(X[,-c(ind.epi,ind.to.tie)],xb.base,xb.epi) #Update X matrix for new curve.
    Xuse2<-cbind(X[,-c(ind.epi,ind.to.tie[-1])],xb.base2,xb.epi) #Update X matrix for new curve.
    XF <- cbind(Xuse,F,Fy)                 #combine new X with  F (updated earlier).    
    XF2 <- cbind(Xuse2,F,Fy)                 #combine new X with F (updated  earlier).    
    
    #Fill in empty rows in L2 to get consistent result.
    
    #storing results.
    year.noise<-Fy%*%L[pLv+fs+1:fy,]
    Yearround.noise.cause[[max(1,i-burnin)]]     <- year.noise[,-age.gender]
    Yearround.noise.agegender[[max(1,i-burnin)]] <- year.noise[,age.gender]
    Yearround.noise.total.cause[[max(1,i-burnin)]]    <-rowSums(year.noise[,-age.gender])
    Yearround.noise.total.agegender[[max(1,i-burnin)]]<-rowSums(year.noise[,age.gender])
    frac.cause                       <-L[pLv+fs+1:fy,-age.gender]/sum(L[pLv+fs+1:fy,-age.gender])
    frac.age                         <-L[pLv+fs+1:fy,age.gender]/sum(L[pLv+fs+1:fy,age.gender])
    Yearround.noise.frac[[max(1,i-burnin)]]      <-c(frac.cause,frac.age)#colSums(epi.at.time)/sum(epi.at.time)
    
    
    epi.at.time                <-F%*%L[pLv+1:fs,]
    Epi.mean.at.time           <-Xuse[,pLv]%*%t(L[pLv,])
    Epidemic.noise.cause[[max(1,i-burnin)]] <-epi.at.time[,-age.gender]
    Epidemic.total.cause[[max(1,i-burnin)]] <-epi.at.time[,-age.gender]+Epi.mean.at.time[,-age.gender]
    Epidemic.noise.agegender[[max(1,i-burnin)]] <-epi.at.time[,age.gender]
    Epidemic.total.agegender[[max(1,i-burnin)]] <-epi.at.time[,age.gender]+Epi.mean.at.time[,age.gender]
    var.unique[[max(1,i-burnin)]]     <-1/tau
    Epidemic.noise.total.cause[[max(1,i-burnin)]] <-rowSums(epi.at.time[,-age.gender])
    Epidemic.all.total.cause[[max(1,i-burnin)]]   <-rowSums(epi.at.time[,-age.gender])+rowSums(Epi.mean.at.time[,-age.gender])
    Epidemic.noise.total.agegender[[max(1,i-burnin)]] <-rowSums(epi.at.time[,age.gender])
    Epidemic.all.total.agegender[[max(1,i-burnin)]]   <-rowSums(epi.at.time[,age.gender])+rowSums(Epi.mean.at.time[,age.gender])
    frac.cause                       <-L[pLv+1:fs,-age.gender]/sum(L[pLv+1:fs,-age.gender])
    frac.age                        <-L[pLv+1:fs,age.gender]/sum(L[pLv+1:fs,age.gender])
    Epi.noise.frac[[max(1,i-burnin)]]      <-c(frac.cause,frac.age)#colSums(epi.at.time)/sum(epi.at.time)
    frac.cause.all                  <-colSums(epi.at.time[,-age.gender]+Epi.mean.at.time[,-age.gender])/sum(epi.at.time[,-age.gender]+Epi.mean.at.time[,-age.gender])
    frac.age.all                    <-colSums(epi.at.time[,age.gender]+Epi.mean.at.time[,age.gender])/sum(epi.at.time[,age.gender]+Epi.mean.at.time[,age.gender])
    
    Epi.all.frac[[max(1,i-burnin)]]       <-c(frac.cause.all,frac.age.all)
    beta[[max(1,i-burnin)]]           <-L2[1:pLv2,] 
    beta.total.cause[[max(1,i-burnin)]]     <-rowSums(L2[1:pLv2,-age.gender])   
    beta.total.agegender[[max(1,i-burnin)]] <-rowSums(L2[1:pLv2,age.gender])  
    base.curve[[max(1,i-burnin)]]     <-xb.base  
    yimp[[max(1,i-burnin)]]     <-y[n+1-mtimes,] 
    
    
  }
  
  x<-list(Yearround.noise.cause,Yearround.noise.agegender,Yearround.noise.total.cause,Yearround.noise.total.agegender,Yearround.noise.frac,Epidemic.noise.cause,Epidemic.total.cause,Epidemic.noise.agegender,Epidemic.total.agegender,Epidemic.noise.total.cause,Epidemic.all.total.cause,Epidemic.noise.total.agegender,Epidemic.all.total.agegender,Epi.noise.frac,Epi.all.frac,var.unique,beta,beta.total.cause,beta.total.agegender,base.curve,yimp)
  names(x)<-c('Yearround.noise.cause','Yearround.noise.agegender','Yearround.noise.total.cause','Yearround.noise.total.agegender','Yearround.noise.frac','Epidemic.noise.cause','Epidemic.all.cause','Epidemic.noise.agegender','Epidemic.all.agegender','Epidemic.noise.total.cause','Epidemic.all.total.cause','Epidemic.noise.total.agegender','Epidemic.all.total.agegender','Epi.noise.frac','Epi.all.frac','var.unique','beta','beta.total.cause','beta.total.agegender','base.curve','yimp')
  
  return(x)
}



