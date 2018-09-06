##################################
# Lasso using coordinate descent #
##################################

set.seed(10000)
ndim<-10
Sigma<-random_varcov_mat(ndim)
Sigma_inv<-solve(Sigma)
L<-chol(Sigma_inv)

xxx<-mvrnorm(1,rep(0,ndim),Sigma)
Lxx<-L%*%xxx
est<-rep(0,ndim)

S<-function(z,g){
  temp<-abs(z)-g
  temp[temp<0]<-0
  return(sign(z)*temp)
}

lstlam<-seq(1,0.01,-0.01)
betpat<-matrix(NA,ndim,length(lstlam))

est<-rep(0,ndim)
for(idxlam in 1:length(lstlam)){
  lam<-lstlam[idxlam]
  for(iter in 1:20){
    for(idx in 1:ndim){
      first<-sum(L[1:idx,idx]*(Lxx[1:idx]-L[1:idx,-idx]%*%est[-idx]))
      est[idx]<-S(first,lam/abs(Lxx[idx]))/sum(L[1:idx,idx]^2)
    }
  }
  betpat[,idxlam]<-est
}

matplot(t(betpat),type="l")

lstlam<-seq(1,0.01,-0.01)
betpat<-matrix(NA,ndim,length(lstlam))

est<-rep(0,ndim)
for(idxlam in 1:length(lstlam)){
  lam<-lstlam[idxlam]
  for(iter in 1:20){
    for(idx in 1:ndim){
      tmpest<-est
      tmpest[idx]<-0
      rrr<-Lxx-L%*%tmpest
      arg1<-t(L[,idx])%*%rrr
      
      est[idx]<-S(arg1,lam)/sum(L[,idx]^2)
    }
  }
  betpat[,idxlam]<-est
}

matplot(t(betpat),type="l")
temp<-glmnet(L,Lxx,intercept=FALSE,standardize=FALSE,lambda=lstlam/ndim)
matplot(t(temp$beta),type="l")

#    _________
#  /           \
# |  THE SAME!  |
#  \           /
#   ```````````
#