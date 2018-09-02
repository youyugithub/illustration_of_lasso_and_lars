### new version

ndim<-10
A<-matrix(runif(ndim^2)*2-1,ncol=ndim) 
Sigma<-t(A)%*%A
# Sigma<-diag(ndim)
Sigma_inv<-solve(Sigma)
L<-solve(chol(Sigma))
max(abs(L%*%t(L)-solve(Sigma)))

x<-rnorm(ndim)
gammalist<-seq(0,5,0.1)
mumatrix<-matrix(NA,length(gammalist),ndim)
for(gammaidx in 1:length(gammalist)){
  muhatold<-rep(0,ndim)
  muhatnew<-x
  gamma<-gammalist[gammaidx]
  while(max(abs(muhatold-muhatnew))>1e-5){
    muhatold<-muhatnew
    l1<-2*(Sigma_inv%*%(muhatold)-x)+gamma*(muhatold/(abs(muhatold*x)+1e-5))
    l2<-2*Sigma_inv
    diag(l2)<-diag(l2)+gamma/(abs(x*muhatold)+1e-5)
    muhatnew<-muhatold-solve(l2,l1)
  }
  mumatrix[gammaidx,]<-muhatnew
}
matplot(mumatrix,type="l")

##########
## LARS ##
##########

set.seed(1000)

A<-matrix(runif(ndim^2)*2-1,ncol=ndim)
Sigma<-t(A)%*%A

ndim<-10
nobs<-30
beta_truth<-rep(0,ndim)
beta_truth[sample(1:ndim,6)]<-sample(c(-1,1),6,replace=TRUE)*(1:6/6)
X<-matrix(runif(nobs*ndim),nobs,ndim)%*%A
y<-X%*%beta_truth+rnorm(nobs,0,0.5)

## Pre-treatment: center and standardize data
X.mean<-apply(X,1,mean)
X.sd<-apply(X,1,sd)
y.mean<-mean(y)
XX<-scale(X)
yy<-scale(y,scale=FALSE)

## Let's start
beta_path<-matrix(NA,ndim,ndim)
# start with mu=0
muhat<-rep(0,nobs)
# start with beta=0
betahat<-rep(0,ndim)
# empty active set
active_set<-rep(FALSE,ndim)
# all-zero sign set
sign_set<-rep(0,ndim)


## First step: activate the first index
# correlation vector
cc<-t(XX)%*%(yy-muhat)
# update sign set
sign_set<-sign(cc)
# active set is initialized using the largest correlation
active_set[which.max(abs(cc))]<-TRUE

for(tempidx in 1:(ndim-1)){
  
  # correlation vector (help us check equality of correlations)
  cc<-t(XX)%*%(yy-muhat)
  # update sign set
  sign_set<-sign(cc)
  
  # Formula (2.4-2.7) in Efron et al (2004)
  one<-rep(1,sum(active_set))
  X_A<-sweep(as.matrix(XX[,active_set]),2,sign_set[active_set],"*")
  G_A<-t(X_A)%*%X_A
  a_A<-sqrt(t(one)%*%solve(G_A,one))
  w_A<-solve(G_A,one)%*%a_A
  u_A<-X_A%*%w_A
  # u_A is the direction of next step
  
  # Formula (2.11)
  aaa<-t(XX)%*%u_A
  
  # Formula (2.13): length of next step
  vneg<-((max(abs(cc))-as.vector(cc))/(as.vector(a_A)-aaa))
  vneg[vneg<=0]<-Inf
  vneg[active_set]<-Inf
  vpos<-((max(abs(cc))+as.vector(cc))/(as.vector(a_A)+aaa))
  vpos[vpos<=0]<-Inf
  vpos[active_set]<-Inf
  gamma_length<-min(pmin(vneg,vpos))
  # update muhat [formula (2.12)] and betahat
  # the formula for betahat is not found in Efron et al. (2004)
  muhat<-muhat+gamma_length*u_A
  betahat[active_set]<-betahat[active_set]+sign_set[active_set]*gamma_length*w_A #
  #beta_path[,tempidx]<-solve(t(XX)%*%XX,t(XX)%*%muhat)
  beta_path[,tempidx]<-betahat

  # activate the next index
  active_set[which.min(pmin(vneg,vpos))]<-TRUE
}

beta_path[,ndim]<-c(solve(t(XX)%*%XX,t(XX)%*%yy))
matplot(apply(abs(beta_path),2,sum),t(beta_path),type='b')
round(beta_path,4)

##################################
## LARS with LASSO modification ##
##################################

set.seed(4000)

A<-matrix(runif(ndim^2)*2-1,ncol=ndim)
Sigma<-t(A)%*%A

ndim<-10
nobs<-30
beta_truth<-rep(0,ndim)
beta_truth[sample(1:ndim,6)]<-sample(c(-1,1),6,replace=TRUE)*(1:6/6)
X<-matrix(runif(nobs*ndim),nobs,ndim)%*%A
y<-X%*%beta_truth+rnorm(nobs,0,0.5)

## Pre-treatment: center and standardize data
X.mean<-apply(X,1,mean)
X.sd<-apply(X,1,sd)
y.mean<-mean(y)
XX<-scale(X)
yy<-scale(y)

## Let's start
beta_path<-matrix(NA,ndim,0)
# start with mu=0
muhat<-rep(0,nobs)
# start with beta=0
betahat<-rep(0,ndim)
# empty active set
active_set<-rep(FALSE,ndim)
nonzero_set<-rep(FALSE,ndim)
# all-zero sign set
sign_set<-rep(0,ndim)


## First step: activate the first index
# correlation vector
cc<-t(XX)%*%(yy-muhat)
# update sign set
sign_set<-sign(cc)
# active set is initialized using the largest correlation
active_set[which.max(abs(cc))]<-TRUE
cnt<-0

while(!all(nonzero_set)){
  cnt<-cnt+1
  # correlation vector (help us check equality of correlations)
  cc<-t(XX)%*%(yy-muhat)
  # update sign set
  sign_set<-sign(cc)

  # Formula (2.4-2.7) in Efron et al (2004)
  one<-rep(1,sum(active_set))
  X_A<-sweep(as.matrix(XX[,active_set]),2,sign_set[active_set],"*")
  G_A<-t(X_A)%*%X_A
  a_A<-sqrt(t(one)%*%solve(G_A,one))
  w_A<-solve(G_A,one)%*%a_A
  u_A<-X_A%*%w_A
  # u_A is the direction of next step
  
  # Formula (2.11)
  aaa<-t(XX)%*%u_A
  
  # Formula (2.13): length of next step
  vneg<-((max(abs(cc))-as.vector(cc))/(as.vector(a_A)-aaa))
  vneg[vneg<=0]<-Inf
  vneg[active_set]<-Inf
  vpos<-((max(abs(cc))+as.vector(cc))/(as.vector(a_A)+aaa))
  vpos[vpos<=0]<-Inf
  vpos[active_set]<-Inf
  gamma_length<-min(pmin(vneg,vpos))

  # An implied step found in R package lars
  if(all(active_set))gamma_length<-max(abs(cc))/a_A
  
  # vector d defined right before formula (3.3)
  ddd<-rep(Inf,ndim)
  ddd[active_set]<-sign_set[active_set]*w_A
  # formula (3.4)
  gamma_temp<--betahat/ddd
  gamma_temp[gamma_temp<=0]<-Inf
  # formula (3.5)
  gamma_tilde<-min(gamma_temp)
  if(gamma_tilde<gamma_length){
    ## by Lasso modification
    # update muhat formula (3.6)
    muhat<-muhat+gamma_tilde*u_A
    # update betahat formula (3.3)
    betahat[active_set]<-betahat[active_set]+gamma_tilde*ddd[active_set] #sign_set[active_set]
    beta_path<-cbind(beta_path,betahat)
    nonzero_set<-active_set
    nonzero_set[which.min(gamma_temp)]<-FALSE
    active_set[which.min(gamma_temp)]<-FALSE
  }else{
    # if there is nothing to add
    if(all(active_set))break
    # otherwise proceed as usual
    # update muhat [formula (2.12)] and betahat
    # the formula for betahat is not found or possibly given in (3.3) in Efron et al. (2004)
    muhat<-muhat+gamma_length*u_A
    betahat[active_set]<-betahat[active_set]+sign_set[active_set]*gamma_length*w_A
    beta_path<-cbind(beta_path,betahat)
    nonzero_set<-active_set
    # activate the next index
    active_set[which.min(pmin(vneg,vpos))]<-TRUE
  }
}

dim(beta_path)

beta_path<-cbind(beta_path,c(solve(t(XX)%*%XX,t(XX)%*%yy)))
matplot(apply(abs(beta_path),2,sum),t(beta_path),type='b',pch="*")
matplot(t(beta_path),type='b',pch="*")
round(beta_path,4)

###########
## check ##
###########

library(lars)
amodel<-lars(XX,yy,type="lar")
plot(amodel)
max(abs(t(amodel$beta)[,-1]-beta_path))
dim(t(amodel$beta))
dim(beta_path)
matplot(apply(abs(beta_path),2,max),t(beta_path),type='b',pch="*")


amodel2<-lars(XX,yy,type="lasso",trace=TRUE)
plot(amodel2,breaks=F)
matplot(apply(abs(amodel2$beta),1,sum)[1:17],amodel2$beta[1:17,],type='b',pch="*")
matplot(apply(abs(beta_path),2,sum)[c(1:15,19)],t(beta_path)[c(1:15,19),],type='b',pch="*")
round(t(amodel2$beta)[,2:16]-beta_path[,1:15],3)
dim(amodel2$beta)
dim(beta_path)
image(round(t(amodel2$beta),4)==0)
image(round(beta_path,4)==0)
dim(amodel2$beta)
