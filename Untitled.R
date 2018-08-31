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

ndim<-10
nobs<-30
beta_truth<-rep(0,ndim)
beta_truth[sample(1:ndim,5)]<-1:5
X<-matrix(runif(nobs*ndim),nobs,ndim)
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
# empty active set
active_set<-rep(FALSE,ndim)

## First step
# correlation vector
cc<-t(XX)%*%(yy-muhat)
# active set is initialized using the largest correlation
active_set[which.max(abs(cc))]<-TRUE

for(tempidx in 1:(ndim-1)){

  # Formula (2.4-2.7) in Efron et al (2004)
  one<-rep(1,sum(active_set))
  X_A<-XX[,active_set]
  G_A<-t(X_A)%*%X_A
  a_A<-sqrt(t(one)%*%solve(G_A,one))
  w_A<-solve(G_A,one)%*%a_A
  u_A<-X_A%*%w_A
  # u_A is the direction of next step
  
  # Formula (2.11)
  aaa<-t(X)%*%u_A
  
  # Formula (2.13): length of next step
  vneg<-((max(abs(cc))-as.vector(cc))/(as.vector(a_A)-aaa))
  vneg[vneg<=0]<-Inf
  vneg[active_set]<-Inf
  vpos<-((max(abs(cc))+as.vector(cc))/(as.vector(a_A)+aaa))
  vpos[vpos<=0]<-Inf
  vpos[active_set]<-Inf
  gamma_length<-min(pmin(vneg,vpos))
  # activate the next index
  active_set[which.min(pmin(vneg,vpos))]<-TRUE
  # update muhat
  muhat<-muhat+gamma_length*u_A
  
  show(muhat)
  beta_path[,tempidx]<-round(solve(t(XX)%*%XX,t(XX)%*%muhat),15)
}

beta_path
matplot(apply(abs(beta_path),2,sum),t(beta_path),type='l')
