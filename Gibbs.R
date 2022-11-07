#iteration: number of iterations of the algorithm. 
##define initial values for parameters## 
#Packages "rgenoud" and "multinomRob" should be installed. 
library(rgenoud) 
library(multinomRob) 
k=2 
iteration=1000 
mix.new=mu.new= var.new=matrix(0,iteration,k) 
z.new=matrix(0,length(x),k) 
mix=mu=var=NULL 
##define hyperparameter for prior distibution### 
delta=mu.0= omega.0= betta.0=rep(1,k) 
tau2=rep(9,k) 
##generate random sample from the prior distributions## 
for(i in 1:k){ 
  mix[i]<-rgamma(n=1,shape=delta[i],rate=1) 
  mu[i]=rnorm(1,mu.0[i],sqrt(tau2[i])) 
  var[i]=1/rgamma(1,omega.0[i],betta.0[i]) 
} 
mix=mix/sum(mix) 
numer=matrix(0,nrow=length(x),ncol=k) 
## Iteration step##### 
for(it in 1:iteration) 
{ 
  ## find the latent variable z##### 
  for(i in 1:k) { 
    numer[,i]=(mix[i]*dnorm(x,mean=mu[i],sd=sqrt(var[i]))) 
  } 
  prob=numer/matrix(rep(rowSums(numer),k),ncol=k,byrow=F) 
  z=matrix(0,length(x),k) 
  for(j in 1:length(x)){ 
    z[j,]=t(rmultinomial(1,prob[j,])) 
  } 
  n.mix=apply(z,2,sum) 
  ## generate parameters from Posterior distribution #### 
  for(i in 1:k) { 
    mix[i]=rgamma(1,shape=delta[i]+n.mix[i],rate=1) 
    mu[i]=rnorm(1,(tau2[i]*sum(x[z[,i]==1])+mu.0[i]*var[i])/ 
                  (n.mix[i]*tau2[i]+var[i]),sqrt((var[i]*tau2[i])/ 
                                                   (n.mix[i]*tau2[i]+var[i]))) 
    var[i]=1/rgamma(1,shape=omega.0[i]+.5*n.mix[i],rate= 
                      betta.0[i] + .5*sum(z[,i]*(x-mu[i])^2)) 
  } 
  mix=mix/sum(mix) 
  ##Save### 
  z.new=z.new+z 
  mix.new[it,]=mix 
  mu.new[it,]=mu 
  var.new[it,]=var 
} 
##END OF ITERATION STAGE### 
##Compute mean between estimated parameter ##### 
apply(mix.new[(iteration/2):(iteration-10),],2,mean) 
apply(mu.new[(iteration/2):(iteration-10),],2,mean) 
apply(var.new[(iteration/2):(iteration-10),],2,mean)