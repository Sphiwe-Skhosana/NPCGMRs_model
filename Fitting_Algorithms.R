library(extraDistr)
library(np)

##A function to compute (yes, compute!) the complete-data log-likelihood for the NPCGMRs model
Q_function_non <- function(y, mu, prop, alpha, sigma2, eta,z_mat, v_mat) {
  n <- length(y)
  k <- length(eta)
  
  # ---- input checks ----
  stopifnot(
    is.matrix(mu),    nrow(mu)    == n, ncol(mu)    == k,
    is.matrix(prop),    nrow(prop)    == n, ncol(prop)    == k,
    is.matrix(sigma2),    nrow(sigma2)    == n, ncol(sigma2)    == k,
    is.matrix(z_mat), nrow(z_mat) == n, ncol(z_mat) == k,
    is.matrix(v_mat), nrow(v_mat) == n, ncol(v_mat) == k,
    length(alpha)  == k, all(alpha  > 0 & alpha  <= 1), all(sigma2 > 0),
    length(eta)    == k, all(eta    > 0)
  )
  
  Q <- 0
  
  for (j in 1:k) {
    
    z_j   <- z_mat[, j]            # E[Z_ij]
    v_j   <- v_mat[, j]            # E[V_ij | Z_ij = 1]
    mu_j  <- mu[, j]
    sd_j  <- sqrt(sigma2[,j])
    sd_j2 <- sqrt(eta[j]) * sd_j
    
    # ---- term 1: log mixing proportion ----
    Q <- Q + sum(z_j * log(prop[,j]))
    
    # ---- term 2: log alpha (weighted by z and v) ----
    Q <- Q + sum(z_j * v_j) * log(alpha[j])
    
    # ---- term 3: log(1 - alpha) (weighted by z and (1-v)) ----
    Q <- Q + sum(z_j * (1 - v_j)) * log(1 - alpha[j])
    
    # ---- term 4: log N(y; mu_j, sigma2_j) for good obs ----
    Q <- Q + sum(z_j * v_j * dnorm(y, mu_j, sd_j, log = TRUE))
    
    # ---- term 5: log N(y; mu_j, eta_j*sigma2_j) for bad obs ----
    Q <- Q + sum(z_j * (1 - v_j) * dnorm(y, mu_j, sd_j2, log = TRUE))
  }
  
  return(Q)
}

##A function to compute (again, yes compute!) the complete-data log-likelihood for the SPCGMRs model
Q_function_non_G <- function(y,mu,prop,sigma2,z_mat) {
  n <- length(y)
  k<- ncol(mu)
  # ---- input checks ----
  stopifnot(
    is.matrix(mu),    nrow(mu)    == n, ncol(mu)    == k,
    is.matrix(prop),    nrow(prop)    == n, ncol(prop)    == k,
    is.matrix(sigma2),    nrow(sigma2)    == n, ncol(sigma2)    == k,
    is.matrix(z_mat), nrow(z_mat) == n, ncol(z_mat) == k)
  Q <- 0
  for (j in 1:k) {
    z_j   <- z_mat[, j]            # E[Z_ij]
    mu_j  <- mu[, j]
    sd_j  <- sqrt(sigma2[,j])
    Q <- Q + sum(z_j * log(prop[,j]))
    Q <- Q + sum(z_j * dnorm(y, mu_j, sd_j, log = TRUE))
  }
  return(Q)
}

Q_function_semi_G <- function(y, mu, prop, sigma2,z_mat) {
  n <- length(y)
  k <- length(prop)
  # ---- input checks ----
  stopifnot(
    is.matrix(mu),    nrow(mu)    == n, ncol(mu)    == k,
    is.matrix(z_mat), nrow(z_mat) == n, ncol(z_mat) == k,
    length(sigma2) == k, all(sigma2 > 0)
  )
  Q <- 0
  
  for (j in 1:k) {
    z_j   <- z_mat[, j]            # E[Z_ij]
    mu_j  <- mu[, j]
    sd_j  <- sqrt(sigma2[j])
    Q <- Q + sum(z_j) * log(prop[j])
    Q <- Q + sum(z_j * dnorm(y, mu_j, sd_j, log = TRUE))
  }
  return(Q)
}


Q_function_semi <- function(y, mu, prop, alpha, sigma2, eta,z_mat, v_mat) {
  n <- length(y)
  k <- length(eta)
  # ---- input checks ----
  stopifnot(
    is.matrix(mu),    nrow(mu)    == n, ncol(mu)    == k,
    is.matrix(z_mat), nrow(z_mat) == n, ncol(z_mat) == k,
    is.matrix(v_mat), nrow(v_mat) == n, ncol(v_mat) == k,
    length(alpha)  == k, all(alpha  > 0 & alpha  <= 1),
    length(sigma2) == k, all(sigma2 > 0),
    length(eta)    == k, all(eta    > 0)
  )
  Q <- 0
  
  for (j in 1:k) {
    z_j   <- z_mat[, j]            # E[Z_ij]
    v_j   <- v_mat[, j]            # E[V_ij | Z_ij = 1]
    mu_j  <- mu[, j]
    sd_j  <- sqrt(sigma2[j])
    sd_j2 <- sqrt(eta[j]) * sd_j
    # ---- term 1: log mixing proportion ----
    Q <- Q + sum(z_j) * log(prop[j])
    # ---- term 2: log alpha (weighted by z and v) ----
    Q <- Q + sum(z_j * v_j) * log(alpha[j])
    # ---- term 3: log(1 - alpha) (weighted by z and (1-v)) ----
    Q <- Q + sum(z_j * (1 - v_j)) * log(1 - alpha[j])
    # ---- term 4: log N(y; mu_j, sigma2_j) for good obs ----
    Q <- Q + sum(z_j * v_j * dnorm(y, mu_j, sd_j, log = TRUE))
    # ---- term 5: log N(y; mu_j, eta_j*sigma2_j) for bad obs ----
    Q <- Q + sum(z_j * (1 - v_j) * dnorm(y, mu_j, sd_j2, log = TRUE))
  }
  return(Q)
}


##This function is definitely used to compute the hat matrix of the kernel linear smoother
hat.matrix=function(x,h){
  n=nrow(x)
  # Compute full hat matrix in one call
  ksum_full <- npksum(
    txdat  = data.frame(x = x),
    exdat  = data.frame(x = x),   # evaluate at ALL training points
    bws    = h,
    tydat  = diag(n),              # identity matrix → returns weight matrix
    leave.one.out = FALSE
  )
  # ksum_full$ksum is n x n (unnormalized)
  # Normalize each row by its sum
  H <- ksum_full$ksum / rowSums(ksum_full$ksum)
  return(H)
}

##Here we are calculating the information criteria for the non-parametric models
IC_non_parametric=function(x,bw,k,obs_llk,comp_llk){
  n=nrow(x)
  dfm=sum(diag(hat.matrix(x,bw)))
  dfn=(3*k-1)*dfm+2*k
  ICL=-2*comp_llk+log(n)*dfn
  BIC=-2*obs_llk+log(n)*dfn
  AIC=-2*obs_llk+2*dfn
  return(list(ICL=ICL,AIC=AIC,BIC=BIC,df=dfn))
}

##Here we are calculating the information criteria for the semi-parametric models

IC_semi_parametric=function(x,bw,k,obs_llk,comp_llk){
  n=nrow(x)
  dfm=sum(diag(hat.matrix(x,bw)))
  dfn=k*dfm+4*k-1
  ICL=-2*comp_llk+log(n)*dfn
  BIC=-2*obs_llk+log(n)*dfn
  AIC=-2*obs_llk+2*dfn
  return(list(ICL=ICL,AIC=AIC,BIC=BIC,df=dfn))
}

##This is a performance measure for the estimated non-parametric functions
RASE=function(est,true){
  library(combinat)
  n=nrow(est)
  k=ncol(est);idx=permn(k)
  out1=out2=NULL
  for(j in 1:factorial(k)){
    out1=c(out1,sqrt(sum(rowSums((est[,idx[[j]]]-true)^2))/n))
    out2=c(out2,sum(rowSums(est[,idx[[j]]]-true))/n)
  }
  return(round(c(min(out1),min(out2)),4))
}


###Kernel function
Kern<-function(x,x0,h){
  z=(x-x0)/h
  f=ifelse(abs(z)<=1,0.75*(1-z^2),0)/h
  #f=dnorm(z)/h
  out=f
  if(sum(f)>0){out=f/sum(f)};
  return(out)
}

###SimulatIng from a contaminated normal
rCNorm=function(n,mu,sigma,alpha,eta){
  z=sample(1:2,n,replace=T,prob=c(alpha,1-alpha))
  y=rep(0,n)
  for(i in 1:n){
    if(z[i]==1){
      y[i]=rnorm(1,mu,sigma)}else{
        y[i]=rnorm(1,mu,sqrt(eta)*sigma)
      }  }
  return(y)
}

##Contaminated normal denisty function
dCNorm=function(u,mu,sigma,alpha,eta){
  alpha*dnorm(u,mu,sigma)+(1-alpha)*dnorm(u,mu,sqrt(eta)*sigma)
}

eta_fun=function(eta,y,mu,sigma2,weig1,weig2){
  w1=weig1*(1-weig2)
  0.5*sum(w1*log(eta*sigma2)+(w1/(eta*sigma2))*(y-mu)^2)
}

alpha_fun=function(alpha,weig1,weig2){
  -sum(weig1*(weig2*log(alpha)+(1-weig2)*log(1-alpha)))
}

mixf_nonpar=function(y,prop,mu,sigma2,alpha,eta){
  n=length(y)
  k=ncol(prop)
  g=NULL
  for(j in 1:k){
    g=cbind(g,prop[,j]*dCNorm(y,mu[,j],sigma2[,j],alpha[j],eta[j]))
  }
  return(g)
}

mixf_nonpar1=function(y,prop,mu,sigma2,alpha,eta){
  n=length(y)
  k=ncol(prop)
  g=NULL
  for(j in 1:k){
    g=cbind(g,prop[,j]*dCNorm(y,mu[,j],sigma2[j],alpha[j],eta[j]))
  }
  return(g)
}

mixf_semipar=function(y,prop,mu,sigma2,alpha,eta){
  n=length(y)
  k=length(prop)
  g=NULL
  for(j in 1:k){
    g=cbind(g,prop[j]*dCNorm(y,mu[,j],sigma2[j],alpha[j],eta[j]))
  }
  return(g)
}

##A function to fit a contaminated Gaussian mixture of linear regressions
CN_mix_reg=function(y,x,k,init.model=NULL){
  n=length(y)
  x=as.matrix(x);p=ncol(x)
  z=cbind(1,x)
  ###Initialize
  initmu=init.model$initmu
  if(is.null(initmu)){Beta0=matrix(rnorm((p+1)*k,0,1),p+1,k);mu0=z%*%Beta0}else{mu0=initmu}
  initeta=init.model$initeta
  if(is.null(initeta)){eta0=runif(k,1,1e3)}else{eta0=initeta}
  initalpha=init.model$initalpha
  if(is.null(initalpha)){alpha0=runif(k,0.85,1)}else{alpha0=initalpha}
  initprop=init.model$initprop
  if(is.null(initprop)){prop0=rep(1/k,k)}else{prop0=initprop}
  initsigma2=init.model$initsigma2
  if(is.null(initsigma2)){sigma20=rgamma(k,2,1);}else{sigma20=initsigma2}
  ###
  LogLik0=sum(log(rowSums(mixf_semipar(y,prop0,mu0,sigma20,alpha0,eta0))))
  diff=1e6
  llk=NULL
  count=0
  while(diff>1e-10){
    ###E-step
    resp1=mixf_semipar(y,prop0,mu0,sigma20,alpha0,eta0);gn1=resp1/rowSums(resp1);
    resp2=sapply(1:k,function(j){(alpha0[j]*dnorm(y-mu0[,j],0,sqrt(sigma20[j])))/dCNorm(y-mu0[,j],0,sigma20[j],alpha0[j],eta0[j])});gn2=resp2
    ###CM-step1
    prop1=colSums(gn1)/n
    alpha1=Beta1=sigma21=NULL
    for(j in 1:k){
      prob=gn1[,j]*(gn2[,j]+((1-gn2[,j])/eta0[j]))
      W=diag(prob)
      alpha1=c(alpha1,optimize(alpha_fun,c(0,1),weig1=gn1[,j],weig2=gn2[,j])$minimum)
      #alpha1=c(alpha1,sum(gn1[,j]*gn2[,j])/sum(gn1[,j]))
      Beta=solve(t(z)%*%W%*%z)%*%t(z)%*%W%*%y
      Beta1=cbind(Beta1,Beta)
      sigma2=sum(prob*(y-z%*%Beta)^2)/sum(gn1[,j])
      sigma21=c(sigma21,sigma2)
    }
    mu1=sapply(1:k,function(j) z%*%Beta1[,j])
    ###CM-step2
    eta1=NULL
    for(j in 1:k){
      #a=sum(gn1[,j]*(1-gn2[,j]));b=sum((gn1[,j]*(1-gn2[,j]))*(y-z%*%Beta1[,j])^2)/sigma21[j]
      res=optimise(eta_fun,interval=c(1.02,1e2),y=y,mu=mu1[,j],sigma2=sigma21[j],weig1=gn1[,j],weig2=gn2[,j])
      if(res$minimum==1.01 | res$minimum==1e2){
        eta1=c(eta1,eta0[j])}else{
          eta1=c(eta1,res$minimum)
          #eta1=c(eta1,max(1.2,b/a))
        }
    }
    ##Checking for convergence
    LogLik1=sum(log(rowSums(mixf_semipar(y,prop1,mu1,sigma21,alpha1,eta1))))
    diff=abs(LogLik1-LogLik0)/abs(LogLik0)
    LogLik0=LogLik1
    eta0=eta1
    sigma20=sigma21
    prop0=prop1
    Beta0=Beta1
    mu0=mu1
    alpha0=alpha1
    llk=c(llk,LogLik1)
    count=count+1
    if(count==5e2) diff=1e-100 ##Set the max no. of iterations at 100
  }
  df=(4+p+1)*k-1
  BIC=-2*LogLik1+log(n)*df
  idx=order(prop1)
  prop1=prop1[idx];eta1=eta1[idx];alpha1=alpha1[idx];sigma21=sigma21[idx];Beta1=matrix(Beta1[,idx],p+1,k);
  mu1=sapply(1:k,function(j) z%*%Beta1[,j])
  resp1=mixf_semipar(y,prop1,Beta1,sigma21,alpha1,eta1);gn1=resp1/rowSums(resp1)
  resp2=sapply(1:k,function(j) (alpha1[j]*dnorm(y-z%*%Beta1[,j],0,sqrt(sigma21[j]))));gn2=resp2/colSums(resp2)
  z1=apply(gn1,1,which.max)
  z2=sapply(1:n,function(i) as.numeric(gn2[i,z1[i]]>0.5))
  out=cbind(z1,z2,ifelse(z2==1,"good","bad"))
  #plot.ts(llk)
  return(list(resp=gn1,z_clust=z1,z_out=z2,alpha=alpha1,out=out,eta=eta1,prop=prop1,sigma2=sigma21,Beta=Beta1,mu=mu1,LL=LogLik0,BIC=BIC,LLiter=llk))
}

##A function to fit a non-parametric Gaussian mixture of regressions model
Kernel_Mix_EM=function(x,y,xgrid=NULL,k,bw,init.model=NULL){
  n=length(y)
  x=as.matrix(x);p=ncol(x)
  z=cbind(1,x)
  if(is.null(xgrid)) xgrid=sapply(1:p,function(j) seq(min(x[,j]),max(x[,j]),length.out=100))
  ##Initial state
  initmu=init.model$initmu
  if(is.null(initmu)){Beta0=matrix(rnorm((p+1)*k,0,1),p+1,k);mu0=z%*%Beta0}else{mu0=initmu}
  initprop=init.model$initprop
  if(is.null(initprop)){pi0=matrix(rep(1/k,k),n,k,byrow=T)}else{pi0=matrix(initprop,n,k,byrow=T)}
  initsigma2=init.model$initsigma2
  if(is.null(initsigma2)){sigma20=matrix(rgamma(k,2,1),n,k,byrow=T)}else{sigma20=matrix(initsigma2,n,k,byrow=T)}
  ##
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[,j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[,j])))))
  diff=1e6
  count=0
  llk=NULL
  while(diff>1e-10){
    ##E-step
    g=sapply(1:k,function(j) pi0[,j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[,j])+1e-100);gn=g/rowSums(g)
    ##M-step
    mu1=pi1=sigma21=NULL
    for(j in 1:k){
      w=cbind(gn[,j])
      prop=npksum(txdat=x,exdat=xgrid,tydat=gn[,j],bws=bw)$ksum/npksum(txdat=x,exdat=xgrid,bws=bw)$ksum
      pi1=cbind(pi1,prop)
      mu=npksum(txdat=x,exdat=xgrid,tydat=y,weights=w,bws=bw)$ksum/npksum(txdat=x,exdat=xgrid,weights=w,bws=bw)$ksum
      mu1=cbind(mu1,mu)
      sigma21=cbind(sigma21,npksum(txdat=x,exdat=xgrid,weights=w,tydat=(y-mu1[,j])^2,bws=bw)$ksum/npksum(txdat=x,exdat=xgrid,weights=w,bws=bw)$ksum)
    }
    pi1=pi1#sapply(1:k,function(j) approx(xgrid,pi1[,j],x,rule=2)$y)
    sigma21=sigma21#sapply(1:k,function(j) approx(xgrid,sigma21[,j],x,rule=2)$y)
    pi1=pi1/rowSums(pi1)
    ##Evaluating convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[,j]*dnorm(y-mu1[,j],0,sqrt(sigma21)[,j])))))
    diff=abs(LogLik1-LogLik0)/abs(LogLik0)
    LogLik0=LogLik1
    mu0=mu1
    pi0=pi1
    sigma20=sigma21
    count=count+1
    llk=c(llk,LogLik0)
    if(count==1e3) diff=1e-100
  }
  comp_llk=Q_function_non_G(y,mu1,pi1,sigma21,gn)
  IC=IC_non_parametric(x,bw,k,LogLik0,comp_llk)
  res=order(colSums(mu1));mu1=mu1[,res];pi1=pi1[,res];sigma21=sigma21[,res]
  g=sapply(1:k,function(j) pi1[,j]*dnorm(y-mu1[,j],0,sqrt(sigma21[,j])));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  z=apply(gn,1,which.max)
  out=list(resp=gn,z_clust=z,mix.prop=pi1,mix.mu=mu1,mix.sigma2=sigma21,LL=LL1,LLiter=llk,IC=IC)
  return(out)
}

##Step1: A function to fit a non-parametric contaminated Gaussian mixture of regressions model using Algorithm 1
Algorithm1_step1=function(y,x,k,bw=h,init.model=NULL,xgrid=NULL,alpha.est=NULL,eta.est=NULL){
  n=length(y)
  x=as.matrix(x);p=ncol(x)
  z=cbind(1,x)
  if(is.null(xgrid)) xgrid=sapply(1:p,function(j) seq(min(x[,j]),max(x[,j]),length.out=100))
  ###Initialize
  initmu=init.model$mu
  if(is.null(initmu)){Beta0=matrix(rnorm((p+1)*k,0,1),p+1,k);mu0=z%*%Beta0}else{mu0=initmu}
  initprop=init.model$prop
  if(is.null(initprop)){prop0=matrix(rep(1/k,k),n,k,byrow=T)}else{prop0=matrix(initprop,n,k,byrow=T)}
  initsigma2=init.model$sigma2
  if(is.null(initsigma2)){sigma20=matrix(rgamma(k,2,1),n,k,byrow=T);}else{sigma20=matrix(initsigma2,n,k,byrow=T)}
  ###
  alpha=init.model$initalpha
  if(is.null(alpha)) alpha=runif(k,0.85,1)
  eta=init.model$initeta
  if(is.null(eta)) eta=runif(k,1,1e2)
  LogLik0 <- sum(log(rowSums(sapply(1:k, function(j) {
    prop0[,j]*(alpha[j]*dnorm(y, mu0[,j], sqrt(sigma20[j])) +
                 (1 - alpha[j]) * dnorm(y, mu0[,j], sqrt(eta[j]) * sqrt(sigma20[j])))}))))
  z <- matrix(0, n, k)
  v <- matrix(0, n, k)
  diff=1e6
  llk=NULL
  count=0
  while(diff>1e-10){
    ###E-step
    for (j in 1:k) {
      z[, j] <- prop0[,j] * (
        alpha[j] * dnorm(y, mu0[,j], sqrt(sigma20[,j])) +
          (1 - alpha[j]) * dnorm(y, mu0[,j], sqrt(eta[j]) * sqrt(sigma20[,j]))
      )+1e-100
    }
    z <- z / rowSums(z)
    
    # --- v-step ---
    for (j in 1:k) {
      denom <- alpha[j] * dnorm(y, mu0[,j], sqrt(sigma20[,j])) +
        (1 - alpha[j]) * dnorm(y, mu0[,j], sqrt(eta[j]) * sqrt(sigma20[,j]))
      v[, j] <- ifelse(denom > 0, (alpha[j] * dnorm(y, mu0[,j], sqrt(sigma20[,j]))) / denom,0.5)
    }
    ###CM-step1
    prop1=mu1=sigma21=NULL
    for(j in 1:k){
      prop=npksum(txdat=x,exdat=xgrid,tydat=cbind(z[,j]),bws=bw,ckertype = "epanechnikov")$ksum/npksum(txdat=x,exdat=xgrid,bws=bw,ckertype = "epanechnikov")$ksum
      prop1=cbind(prop1,prop)#multivariate_interpolate(x,xgrid,prop))
      #prob=gn1[,j]*(gn2[,j]+((1-gn2[,j])/eta0[j]))*Kh
      prob=cbind(z[,j]*(v[,j]+((1-v[,j])/eta[j])))
      mu=npksum(txdat=x,exdat=xgrid,tydat=y,weights=prob,bws=bw,ckertype = "epanechnikov")$ksum/npksum(txdat=x,exdat=xgrid,weights=prob,bws=bw,ckertype = "epanechnikov")$ksum
      mu1=cbind(mu1,mu)#multivariate_interpolate(x,xgrid,mu))
      sigma2=npksum(txdat=x,exdat=xgrid,tydat=(y-mu1[,j])^2,weights=prob,bws=bw,ckertype = "epanechnikov")$ksum/npksum(txdat=x,exdat=xgrid,weights=cbind(z[,j]),bws=bw,ckertype = "epanechnikov")$ksum
      sigma21=cbind(sigma21,sigma2)#multivariate_interpolate(x,xgrid,sigma2))
    }
    ##Checking for convergence
    LogLik1 <- sum(log(rowSums(sapply(1:k, function(j) {
      prop1[,j]*(alpha[j]*dnorm(y, mu1[,j], sqrt(sigma21[,j])) +
                   (1 - alpha[j]) * dnorm(y, mu1[,j], sqrt(eta[j]) * sqrt(sigma21[,j])))}))))
    diff=abs(LogLik1-LogLik0)/abs(LogLik0)
    LogLik0=LogLik1
    sigma20=sigma21
    prop0=prop1
    mu0=mu1
    llk=c(llk,LogLik1)
    count=count+1
    if(count==1e3) diff=1e-100 ##Set the max no. of iterations at 100
  }
  df=(4+p+1)*k-1
  BIC=-2*LogLik1+log(n)*df
  idx=order(colSums(mu1))
  prop1=prop1[,idx];sigma21=sigma21[,idx];mu1=mu1[,idx];
  resp1=sapply(1:k,function(j) 
    prop1[,j]*(alpha[j]*dnorm(y,mu1[,j],sqrt(sigma21[,j]))+(1-alpha[j])*dnorm(y,mu1[,j],sqrt(eta[j])*sqrt(sigma21[,j]))))
  gn1=resp1/rowSums(resp1)
  resp2=sapply(1:k,function(j){
    denom=alpha[j]*dnorm(y,mu1[,j],sqrt(sigma21[,j]))+(1-alpha[j])*dnorm(y,mu1[,j],sqrt(eta[j])*sqrt(sigma21[,j]))
    ifelse(denom > 0, (alpha[j] * dnorm(y, mu1[,j], sqrt(sigma21[,j]))) / denom,0.5)})
  gn2=resp2/rowSums(resp2)
  z1=apply(gn1,1,which.max)
  z2=sapply(1:n,function(i) as.numeric(gn2[i,z1[i]]>0.5))
  out=cbind(z1,z2,ifelse(z2==1,"good","bad"))
  return(list(resp=gn1,z_clust=z1,z_out=z2,out=out,prop=prop1,sigma2=sigma21,mu=mu1,LL=LogLik0,BIC=BIC,LLiter=llk))
}

##Step2: A function to fit a non-parametric contaminated Gaussian mixture of regressions model using Algorithm 1
Algorithm1_step2=function(y,x,k,bw,init.model=NULL,xgrid=NULL,mu.est=NULL,pi.est=NULL,sigma2.est=NULL){
  n=length(y)
  x=as.matrix(x);p=ncol(x)
  z=cbind(1,x)
  if(is.null(xgrid)) xgrid=sapply(1:p,function(j) seq(min(x[,j]),max(x[,j]),length.out=100))
  ###Initialize
  initeta=init.model$eta
  if(is.null(initeta)){eta0=runif(k,1,1e2)}else{eta0=initeta}
  initalpha=init.model$alpha
  if(is.null(initalpha)){alpha0=runif(k,0.85,1)}else{alpha0=initalpha}
  ###
  mu=mu.est
  prop=pi.est
  sigma2=sigma2.est
  LogLik0 <- sum(log(rowSums(sapply(1:k, function(j) {
    prop[,j]*(alpha0[j]*dnorm(y,mu[,j],sqrt(sigma2[j]))+(1-alpha0[j])*dnorm(y,mu[,j], sqrt(eta0[j]) * sqrt(sigma2[j])))}))))
  z <- matrix(0, n, k)
  v <- matrix(0, n, k)
  diff=1e6
  llk=NULL
  count=0
  while(diff>1e-10){
    ###E-step
    for (j in 1:k) {
      z[, j] <- prop[,j] * (
        alpha0[j] * dnorm(y, mu[,j], sqrt(sigma2[,j])) +
          (1 - alpha0[j]) * dnorm(y, mu[,j], sqrt(eta0[j]) * sqrt(sigma2[,j]))
      )
    }
    z <- z / rowSums(z)
    
    # --- v-step ---
    for (j in 1:k) {
      denom <- alpha0[j] * dnorm(y, mu[,j], sqrt(sigma2[,j])) +
        (1 - alpha0[j]) * dnorm(y, mu[,j], sqrt(eta0[j]) * sqrt(sigma2[,j]))
      v[, j] <- ifelse(denom > 0, (alpha0[j] * dnorm(y, mu[,j], sqrt(sigma2[,j]))) / denom,0.5)
    }
    ###CM-step1
    alpha1=NULL
    for(j in 1:k){
      #alpha1=c(alpha1,optimise(alpha_fun,c(0,1),weig1=z[,j],weig2=v[,j])$minimum)
      alpha1=c(alpha1,sum(z[,j]*v[,j])/sum(z[,j]))
    }
    ###CM-step2
    eta1=NULL
    for(j in 1:k){
      a=sum(z[,j]*(1-v[,j]));b=sum((z[,j]*(1-v[,j]))*(y-mu[,j])^2)/sigma2[,j]
      #res=optimise(eta_fun,interval=c(1,1e2),y=y,mu=mu[,j],sigma2=sigma2[,j],weig1=z[,j],weig2=v[,j])
      #eta1=c(eta1,res$minimum)
      eta1=c(eta1,max(1.01,b/a))
    }
    ##Checking for convergence
    LogLik1 <- sum(log(rowSums(sapply(1:k, function(j) {
      prop[,j]*(alpha1[j]*dnorm(y, mu[,j], sqrt(sigma2[,j])) +
                  (1 - alpha1[j]) * dnorm(y, mu[,j], sqrt(eta1[j]) * sqrt(sigma2[,j])))}))))
    diff=abs(LogLik1-LogLik0)/abs(LogLik0)
    LogLik0=LogLik1
    eta0=eta1
    alpha0=alpha1
    llk=c(llk,LogLik1)
    count=count+1
    if(count==1e3) diff=1e-100 ##Set the max no. of iterations at 100
  }
  comp_llk=Q_function_non(y,mu,prop,alpha1,sigma2,eta1,z,v);
  IC=IC_non_parametric(x,bw,k,LogLik0,comp_llk)
  idx=order(colSums(mu))
  eta1=eta1[idx];alpha1=alpha1[idx]
  resp1=sapply(1:k,function(j) 
    prop[,j]*(alpha1[j] * dnorm(y, mu[,j], sqrt(sigma2[,j]))+(1 - alpha1[j]) * dnorm(y, mu[,j], sqrt(eta1[j]) * sqrt(sigma2[,j]))
    ))
  gn1=resp1/rowSums(resp1)
  resp2=sapply(1:k,function(j){
    denom=alpha1[j] * dnorm(y, mu[,j], sqrt(sigma2[,j]))+(1-alpha1[j])*dnorm(y, mu[,j], sqrt(eta1[j])*sqrt(sigma2[,j]))
    ifelse(denom > 0, (alpha1[j] * dnorm(y, mu[,j], sqrt(sigma2[,j]))) / denom,0.5)})
  gn2=resp2/rowSums(resp2)
  z1=apply(gn1,1,which.max)
  z2=sapply(1:n,function(i) as.numeric(gn2[i,z1[i]]>0.5))
  out=cbind(z1,z2,ifelse(z2==1,"good","bad"))
  return(list(resp=gn1,z_clust=z1,z_out=z2,alpha=alpha1,out=out,eta=eta1,prop=prop,sigma2=sigma2,mu=mu,LL=LogLik0,IC=IC,LLiter=llk))
}

##A function to fit a non-parametric contaminated Gaussian mixture of regressions model using Algorithm 1
CN_mix_non_reg_alg1=function(y,x,k,bw=h,init.model=NULL,xgrid=NULL,alpha=NULL,eta=NULL){
  fit1=Algorithm1_step1(y,x,k,bw=bw,init.model = init.model,alpha.est=alpha,eta.est = eta,xgrid=x)
  fit2=Algorithm1_step2(y,x,k,bw=bw,init.model=init.model,mu.est=fit1$mu,pi.est=fit1$prop,sigma2.est=fit1$sigma2,xgrid=x)
  return(fit2)
}

##A function to fit a non-parametric contaminated Gaussian mixture of regressions model using Algorithm 2: ECM algorithm
CN_mix_non_reg=function(y,x,k,bw=h,init.model=NULL,xgrid=NULL){
  n=length(y)
  x=as.matrix(x);p=ncol(x)
  z=cbind(1,x)
  if(is.null(xgrid)) xgrid=sapply(1:p,function(j) seq(min(x[,j]),max(x[,j]),length.out=100))
  ###Initialize
  initmu=init.model$mu
  if(is.null(initmu)){Beta0=matrix(rnorm((p+1)*k,0,1),p+1,k);mu0=z%*%Beta0}else{mu0=initmu}
  initeta=init.model$eta
  if(is.null(initeta)){eta0=runif(k,1,1e2)}else{eta0=initeta}
  initalpha=init.model$alpha
  if(is.null(initalpha)){alpha0=runif(k,0.85,1)}else{alpha0=initalpha}
  initprop=init.model$prop
  if(is.null(initprop)){prop0=matrix(rep(1/k,k),n,k,byrow=T)}else{prop0=matrix(initprop,n,k,byrow=T)}
  initsigma2=init.model$sigma2
  if(is.null(initsigma2)){sigma20=matrix(rgamma(k,2,1),n,k,byrow=T);}else{sigma20=matrix(initsigma2,n,k,byrow=T)}
  ###
  LogLik0 <- sum(log(rowSums(sapply(1:k, function(j) {
    prop0[,j]*(alpha0[j]*dnorm(y, mu0[,j], sqrt(sigma20[j])) +
                 (1 - alpha0[j]) * dnorm(y, mu0[,j], sqrt(eta0[j]) * sqrt(sigma20[j])))}))))
  z <- matrix(0, n, k)
  v <- matrix(0, n, k)
  diff=1e6
  llk=NULL
  count=0
  while(diff>1e-10){
    ###E-step
    for (j in 1:k) {
      z[, j] <- prop0[,j] * (
        alpha0[j] * dnorm(y, mu0[,j], sqrt(sigma20[,j])) +
          (1 - alpha0[j]) * dnorm(y, mu0[,j], sqrt(eta0[j]) * sqrt(sigma20[,j]))
      )
    }
    z <- z / rowSums(z)
    
    # --- v-step ---
    for (j in 1:k) {
      denom <- alpha0[j] * dnorm(y, mu0[,j], sqrt(sigma20[,j])) +
        (1 - alpha0[j]) * dnorm(y, mu0[,j], sqrt(eta0[j]) * sqrt(sigma20[,j]))
      v[, j] <- ifelse(denom > 0, (alpha0[j] * dnorm(y, mu0[,j], sqrt(sigma20[,j]))) / denom,0.5)
    }
    ###CM-step1
    prop1=alpha1=mu1=sigma21=NULL
    for(j in 1:k){
      prop=npksum(txdat=x,exdat=xgrid,tydat=cbind(z[,j]),bws=bw,ckertype = "epanechnikov")$ksum/npksum(txdat=x,exdat=xgrid,bws=bw,ckertype = "epanechnikov")$ksum
      prop1=cbind(prop1,prop)#multivariate_interpolate(x,xgrid,prop))
      #prob=gn1[,j]*(gn2[,j]+((1-gn2[,j])/eta0[j]))*Kh
      prob=cbind(z[,j]*(v[,j]+((1-v[,j])/eta0[j])))
      #alpha1=c(alpha1,optimise(alpha_fun,c(0.8,1),weig1=z[,j],weig2=v[,j])$minimum)
      alpha1=c(alpha1,sum(z[,j]*v[,j])/sum(z[,j]))
      mu=npksum(txdat=x,exdat=xgrid,tydat=y,weights=prob,bws=bw,ckertype = "epanechnikov")$ksum/npksum(txdat=x,exdat=xgrid,weights=prob,bws=bw,ckertype = "epanechnikov")$ksum
      mu1=cbind(mu1,mu)#multivariate_interpolate(x,xgrid,mu))
      sigma2=npksum(txdat=x,exdat=xgrid,tydat=(y-mu1[,j])^2,weights=prob,bws=bw,ckertype = "epanechnikov")$ksum/npksum(txdat=x,exdat=xgrid,weights=cbind(z[,j]),bws=bw,ckertype = "epanechnikov")$ksum
      sigma21=cbind(sigma21,sigma2)#multivariate_interpolate(x,xgrid,sigma2))
    }
    #sigma21[which(sigma21==0)]=1e-6
    ###CM-step2
    eta1=NULL
    for(j in 1:k){
      a=sum(z[,j]*(1-v[,j]));b=sum((z[,j]*(1-v[,j]))*(y-mu1[,j])^2/sigma21[,j])
      #res=optimise(eta_fun,interval=c(1,1e2),y=y,mu=mu1[,j],sigma2=sigma21[,j],weig1=z[,j],weig2=v[,j])
      #eta1=c(eta1,res$minimum)
      eta1=c(eta1,max(1.01,b/a))
    }
    ##Checking for convergence
    LogLik1 <- sum(log(rowSums(sapply(1:k, function(j) {
      prop1[,j]*(alpha1[j]*dnorm(y, mu1[,j], sqrt(sigma21[,j])) +
                   (1 - alpha1[j]) * dnorm(y, mu1[,j], sqrt(eta1[j]) * sqrt(sigma21[,j])))})+1e-100)))
    diff=abs(LogLik1-LogLik0)/abs(LogLik0)
    if(is.na(diff)){diff=1e-20}else{
      LogLik0=LogLik1
      eta0=eta1
      sigma20=sigma21
      prop0=prop1
      mu0=mu1
      alpha0=alpha1
    }
    llk=c(llk,LogLik1)
    count=count+1
    if(count==1e3) diff=1e-100 ##Set the max no. of iterations at 100
  }
  comp_llk=Q_function_non(y,mu1,prop1,alpha1,sigma21,eta1,z,v)
  IC=IC_non_parametric(x,bw,k,LogLik0,comp_llk)
  idx=order(colSums(mu1))
  prop1=prop1[,idx];eta1=eta1[idx];alpha1=alpha1[idx];sigma21=sigma21[,idx];mu1=mu1[,idx];
  resp1=sapply(1:k,function(j) 
    prop1[,j]*(alpha1[j]*dnorm(y,mu1[,j],sqrt(sigma21[,j]))+(1-alpha1[j])*dnorm(y,mu1[,j],sqrt(eta1[j])*sqrt(sigma21[,j]))))
  gn1=resp1/rowSums(resp1)
  resp2=sapply(1:k,function(j){
    denom=alpha1[j]*dnorm(y, mu1[,j], sqrt(sigma21[,j]))+(1-alpha1[j])*dnorm(y,mu1[,j],sqrt(eta1[j])*sqrt(sigma21[,j]))
    ifelse(denom>0,(alpha1[j]*dnorm(y,mu1[,j],sqrt(sigma21[,j])))/denom,0.5)})
  gn2=resp2/rowSums(resp2)
  z1=apply(gn1,1,which.max)
  z2=sapply(1:n,function(i) as.numeric(gn2[i,z1[i]]>0.5))
  out=cbind(z1,z2,ifelse(z2==1,"good","bad"))
  return(list(resp=gn1,z_clust=z1,z_out=z2,alpha=alpha1,out=out,eta=eta1,prop=prop1,sigma2=sigma21,mu=mu1,LL=LogLik0,IC=IC,LLiter=llk))
}

##A function to fit a non-parametric contaminated Gaussian mixture of regressions model where the variance is a constant
#using Algorithm 2: ECM algorithm
CN_mix_non_reg1=function(y,x,k,bw=h,init.model=NULL,xgrid=NULL){
  n=length(y)
  x=as.matrix(x);p=ncol(x)
  z=cbind(1,x)
  if(is.null(xgrid)) xgrid=sapply(1:p,function(j) seq(min(x[,j]),max(x[,j]),length.out=100))
  ###Initialize
  initmu=init.model$mu
  if(is.null(initmu)){Beta0=matrix(rnorm((p+1)*k,0,1),p+1,k);mu0=z%*%Beta0}else{mu0=initmu}
  initeta=init.model$eta
  if(is.null(initeta)){eta0=runif(k,1,1e2)}else{eta0=initeta}
  initalpha=init.model$alpha
  if(is.null(initalpha)){alpha0=runif(k,0.85,1)}else{alpha0=initalpha}
  initprop=init.model$prop
  if(is.null(initprop)){prop0=matrix(rep(1/k,k),n,k,byrow=T)}else{prop0=matrix(initprop,n,k,byrow=T)}
  initsigma2=init.model$initsigma2
  if(is.null(initsigma2)){sigma20=rgamma(k,2,1)}else{sigma20=initsigma2}
  ###
  LogLik0 <- sum(log(rowSums(sapply(1:k, function(j) {
    prop0[,j]*(alpha0[j]*dnorm(y, mu0[,j], sqrt(sigma20[j])) +
                 (1 - alpha0[j]) * dnorm(y, mu0[,j], sqrt(eta0[j]) * sqrt(sigma20[j])))}))))
  
  z <- matrix(0, n, k)
  v <- matrix(0, n, k)
  diff=1e6
  llk=NULL
  count=0
  while(diff>1e-10){
    ###E-step
    for (j in 1:k) {
      z[, j] <- prop0[,j] * (
        alpha0[j] * dnorm(y, mu0[,j], sqrt(sigma20[j])) +
          (1 - alpha0[j]) * dnorm(y, mu0[,j], sqrt(eta0[j]) * sqrt(sigma20[j]))
      )
    }
    z <- z / rowSums(z)
    
    # --- v-step ---
    for (j in 1:k) {
      denom <- alpha0[j] * dnorm(y, mu0[,j], sqrt(sigma20[j])) +
        (1 - alpha0[j]) * dnorm(y, mu0[,j], sqrt(eta0[j]) * sqrt(sigma20[j]))
      v[, j] <- ifelse(denom > 0, (alpha0[j] * dnorm(y, mu0[,j], sqrt(sigma20[j]))) / denom,0.5)
    }
    ###CM-step1
    prop1=alpha1=mu1=sigma21=NULL
    for(j in 1:k){
      prop=npksum(txdat=x,exdat=xgrid,tydat=cbind(z[,j]),bws=bw,ckertype = "epanechnikov")$ksum/npksum(txdat=x,exdat=xgrid,bws=bw,ckertype = "epanechnikov")$ksum
      prop1=cbind(prop1,prop)#multivariate_interpolate(x,xgrid,prop))
      #prob=gn1[,j]*(gn2[,j]+((1-gn2[,j])/eta0[j]))*Kh
      prob=cbind(z[,j]*(v[,j]+((1-v[,j])/eta0[j])))
      #alpha1=c(alpha1,optimise(alpha_fun,c(0.8,1),weig1=z[,j],weig2=v[,j])$minimum)
      alpha1=c(alpha1,sum(gn1[,j]*gn2[,j])/sum(gn1[,j]))
      mu=npksum(txdat=x,exdat=xgrid,tydat=y,weights=prob,bws=bw,ckertype = "epanechnikov")$ksum/npksum(txdat=x,exdat=xgrid,weights=prob,bws=bw,ckertype = "epanechnikov")$ksum
      mu1=cbind(mu1,mu)#multivariate_interpolate(x,xgrid,mu))
      sigma2=sum(prob*(y-mu1[,j])^2)/sum(z[,j])
      sigma21=c(sigma21,sigma2)
    }
    ###CM-step2
    eta1=NULL
    for(j in 1:k){
      a=sum(gn1[,j]*(1-gn2[,j]));b=sum((gn1[,j]*(1-gn2[,j]))*(y-mu1[,j])^2)/sigma21[,j]
      #res=optimise(eta_fun,interval=c(1,1e2),y=y,mu=mu1[,j],sigma2=sigma21[j],weig1=z[,j],weig2=v[,j])
      #eta1=c(eta1,res$minimum)
      eta1=c(eta1,max(1.01,b/a))
    }
    ##Checking for convergence
    LogLik1 <- sum(log(rowSums(sapply(1:k, function(j) {
      prop1[,j]*(alpha1[j]*dnorm(y, mu1[,j], sqrt(sigma21[j])) +
                   (1 - alpha1[j]) * dnorm(y, mu1[,j], sqrt(eta1[j]) * sqrt(sigma21[j])))}))))
    diff=abs(LogLik1-LogLik0)/abs(LogLik0)
    LogLik0=LogLik1
    eta0=eta1
    sigma20=sigma21
    prop0=prop1
    mu0=mu1
    alpha0=alpha1
    llk=c(llk,LogLik1)
    count=count+1
    if(count==1e3) diff=1e-100 ##Set the max no. of iterations at 100
  }
  df=(4+p+1)*k-1
  BIC=-2*LogLik1+log(n)*df
  idx=order(colSums(mu1))
  prop1=prop1[,idx];eta1=eta1[idx];alpha1=alpha1[idx];sigma21=sigma21[idx];mu1=mu1[,idx];
  resp1=sapply(1:k,function(j) 
    prop1[,j]*(alpha1[j] * dnorm(y, mu1[,j], sqrt(sigma21[j]))+(1 - alpha1[j]) * dnorm(y, mu1[,j], sqrt(eta1[j])*sqrt(sigma21[j]))
    ))
  gn1=resp1/rowSums(resp1)
  resp2=sapply(1:k,function(j){
    denom=alpha1[j]*dnorm(y,mu1[,j],sqrt(sigma21[j]))+(1-alpha1[j])*dnorm(y, mu1[,j],sqrt(eta1[j])*sqrt(sigma21[j]))
    ifelse(denom>0,(alpha1[j]*dnorm(y, mu1[,j], sqrt(sigma21[j])))/denom,0.5)})
  gn2=resp2/rowSums(resp2)  
  z1=apply(gn1,1,which.max)
  z2=sapply(1:n,function(i) as.numeric(gn2[i,z1[i]]>0.5))
  out=cbind(z1,z2,ifelse(z2==1,"good","bad"))
  return(list(resp=gn1,z_clust=z1,z_out=z2,alpha=alpha1,out=out,eta=eta1,prop=prop1,sigma2=sigma21,mu=mu1,LL=LogLik0,BIC=BIC,LLiter=llk))
}

##A function to fit a semi-parametric Gaussian mixture of regressions model of Xiang and Yao (2018)
Kernel_Mix_EM_semi=function(x,y,k,bw,init.model=NULL,xgrid=NULL){
  n=length(y)
  x=as.matrix(x);p=ncol(x)
  z=cbind(1,x)
  if(is.null(xgrid)) xgrid=sapply(1:p,function(j) seq(min(x[,j]),max(x[,j]),length.out=100))
  ##Initial state
  initmu=init.model$initmu
  if(is.null(initmu)){Beta0=matrix(rnorm((p+1)*k,0,1),p+1,k);mu0=z%*%Beta0}else{mu0=initmu}
  initprop=init.model$initprop
  if(is.null(initprop)){pi0=rep(1/k,k)}else{pi0=initprop}
  initsigma2=init.model$initsigma2
  if(is.null(initsigma2)){sigma20=rgamma(k,2,1)}else{sigma20=initsigma2}
  ##
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[j])))))
  #Kh=multivariate_kernel(x,xgrid,h)
  diff=1e6
  count=0
  llk=NULL
  while(diff>1e-10){
    ##E-step
    g=sapply(1:k,function(j) pi0[j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[j])+1e-100);gn=g/rowSums(g)
    ##M-step
    pi1=colSums(gn)/n
    mu1=sigma21=NULL
    for(j in 1:k){
      #W=gn[,j]*Kh
      w=cbind(gn[,j])
      mu1=cbind(mu1,npksum(txdat=x,exdat=xgrid,tydat=y,weights=w,bws=bw,ckertype = "epanechnikov")$ksum/npksum(txdat=x,exdat=xgrid,weights=w,bws=bw,ckertype = "epanechnikov")$ksum)
      sigma21=cbind(sigma21,sum(gn[,j]*(y-mu1[,j])^2)/sum(gn[,j]))
    }
    ##Evaluating convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mu1[,j],0,sqrt(sigma21)[j])))))
    diff=abs(LogLik1-LogLik0)/abs(LogLik0)
    LogLik0=LogLik1
    mu0=mu1
    pi0=pi1
    sigma20=sigma21
    count=count+1
    llk=c(llk,LogLik0)
    if(count==1e3) diff=1e-100
  }
  comp_llk=Q_function_semi_G(y,mu1,pi1,sigma21,gn)
  IC=IC_semi_parametric(x,bw,k,LogLik0,comp_llk)
  idx=order(colSums(mu1))
  pi1=pi1[idx];sigma21=sigma21[idx];mu1=mu1[,idx];
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu1[,j],0,sqrt(sigma21[j])));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  z=apply(gn,1,which.max)
  out=list(resp=gn,z_clust=z,mix.prop=pi1,mix.mu=mu1,mix.sigma2=sigma21,LL=LL1,LLiter=llk,IC=IC)
  return(out)
}

##A function to fit a semi-parametric contaminated Gaussian mixture of regressions model proposed in Skhosana and Yao (2026)
CN_mix_semi_reg=function(y,x,k,bw,init.model=NULL,xgrid=NULL){
  n=length(y)
  x=as.matrix(x);p=ncol(x)
  z=cbind(1,x)
  if(is.null(xgrid)) xgrid=sapply(1:p,function(j) seq(min(x[,j]),max(x[,j]),length.out=100))
  ###Initialize
  initmu=init.model$initmu
  if(is.null(initmu)){Beta0=matrix(rnorm((p+1)*k,0,1),p+1,k);mu0=z%*%Beta0}else{mu0=initmu}
  initeta=init.model$initeta
  if(is.null(initeta)){eta0=runif(k,1,1e3)}else{eta0=initeta}
  initalpha=init.model$initalpha
  if(is.null(initalpha)){alpha0=runif(k,0.85,1)}else{alpha0=initalpha}
  initprop=init.model$initprop
  if(is.null(initprop)){prop0=rep(1/k,k)}else{prop0=initprop}
  initsigma2=init.model$initsigma2
  if(is.null(initsigma2)){sigma20=rgamma(k,2,1)}else{sigma20=initsigma2}
  ###
  LogLik0 <- sum(log(rowSums(sapply(1:k, function(j) {
    prop0[j]*(alpha0[j]*dnorm(y, mu0[,j], sqrt(sigma20[j])) +
                (1 - alpha0[j]) * dnorm(y, mu0[,j], sqrt(eta0[j]) * sqrt(sigma20[j])))}))))
  
  z <- matrix(0, n, k)
  v <- matrix(0, n, k)
  diff=1e6
  llk=NULL
  count=0
  while(diff>1e-10){
    ###E-step
    for (j in 1:k) {
      z[, j] <- prop0[j] * (
        alpha0[j] * dnorm(y, mu0[,j], sqrt(sigma20[j])) +
          (1 - alpha0[j]) * dnorm(y, mu0[,j], sqrt(eta0[j]) * sqrt(sigma20[j]))
      )
    }
    z <- z / rowSums(z)
    
    # --- v-step ---
    for (j in 1:k) {
      denom <- alpha0[j] * dnorm(y, mu0[,j], sqrt(sigma20[j])) +
        (1 - alpha0[j]) * dnorm(y, mu0[,j], sqrt(eta0[j]) * sqrt(sigma20[j]))
      v[, j] <- ifelse(denom > 0, (alpha0[j] * dnorm(y, mu0[,j], sqrt(sigma20[j]))) / denom,0.5)
    }
    ###CM-step1
    prop1=colSums(z)/n
    alpha1=mu1=sigma21=NULL
    for(j in 1:k){
      prob=z[,j]*(v[,j]+((1-v[,j])/eta0[j]))
      prob1=cbind(z[,j]*(v[,j]+((1-v[,j])/eta0[j])))
      alpha1=c(alpha1,optimise(alpha_fun,c(0.8,1),weig1=z[,j],weig2=v[,j])$minimum)
      #alpha1=c(alpha1,sum(z[,j]*v[,j])/sum(z[,j]))
      mu=npksum(txdat=x,exdat=xgrid,tydat=y,weights=prob1,bws=bw,ckertype = "epanechnikov")$ksum/npksum(txdat=x,exdat=xgrid,weights=prob1,bws=bw,ckertype = "epanechnikov")$ksum
      mu1=cbind(mu1,mu)#multivariate_interpolate(x,xgrid,mu))
      sigma2=sum(prob*(y-mu1[,j])^2)/sum(z[,j])
      sigma21=c(sigma21,sigma2)
    }
    ###CM-step2
    eta1=NULL
    for(j in 1:k){
      #a=sum(z[,j]*(1-v[,j]));b=sum((z[,j]*(1-v[,j]))*(y-mu1[,j])^2)/sigma21[j]
      res=optimise(eta_fun,interval=c(1,1e3),y=y,mu=mu1[,j],sigma2=sigma21[j],weig1=z[,j],weig2=v[,j])
      eta1=c(eta1,res$minimum)
      #eta1=c(eta1,max(1.01,b/a))
    }
    ##Checking for convergence
    LogLik1 <- sum(log(rowSums(sapply(1:k, function(j) {
      prop1[j]*(alpha1[j]*dnorm(y, mu1[,j], sqrt(sigma21[j]))+(1-alpha1[j])*dnorm(y,mu1[,j],sqrt(eta1[j])*sqrt(sigma21[j])))}))))
    diff=abs(LogLik1-LogLik0)/abs(LogLik0)
    LogLik0=LogLik1
    eta0=eta1
    sigma20=sigma21
    prop0=prop1
    mu0=mu1
    alpha0=alpha1
    llk=c(llk,LogLik1)
    count=count+1
    if(count==1e3) diff=1e-100 ##Set the max no. of iterations at 100
  }
  comp_llk=Q_function_semi(y,mu1,prop1,alpha1,sigma21,eta1,z,v)
  IC=IC_semi_parametric(x,bw,k,LogLik0,comp_llk)
  idx=order(colSums(mu1))
  prop1=prop1[idx];eta1=eta1[idx];alpha1=alpha1[idx];sigma21=sigma21[idx];mu1=mu1[,idx];
  resp1=sapply(1:k,function(j) 
    prop1[j]*(alpha1[j] * dnorm(y, mu1[,j], sqrt(sigma21[j]))+(1-alpha1[j])*dnorm(y,mu1[,j],sqrt(eta1[j])*sqrt(sigma21[j]))
    ))
  gn1=resp1/rowSums(resp1)
  resp2=sapply(1:k,function(j){
    denom=alpha1[j] * dnorm(y, mu1[,j], sqrt(sigma21[j]))+(1-alpha1[j])*dnorm(y, mu1[,j], sqrt(eta1[j])*sqrt(sigma21[j]))
    ifelse(denom>0,(alpha1[j]*dnorm(y, mu1[,j],sqrt(sigma21[j])))/denom,0.5)})
  gn2=resp2/rowSums(resp2)  
  z1=apply(gn1,1,which.max)
  z2=sapply(1:n,function(i) as.numeric(gn2[i,z1[i]]>0.5))
  out=cbind(z1,z2,ifelse(z2==1,"good","bad"))
  return(list(resp=gn1,z_clust=z1,z_out=z2,alpha=alpha1,out=out,eta=eta1,prop=prop1,sigma2=sigma21,mu=mu1,LL=LogLik0,IC=IC,LLiter=llk))
}

##Step 1: A function to perform the first step of Algorithm 3 for fitting a semi-parametric contaminated Gaussian mixture of regressions model
Algorithm1_semi_step1=function(y,x,k,bw=h,init.model=NULL,xgrid=NULL){
  n=length(y)
  x=as.matrix(x);p=ncol(x)
  z=cbind(1,x)
  if(is.null(xgrid)) xgrid=sapply(1:p,function(j) seq(min(x[,j]),max(x[,j]),length.out=100))
  ###Initialize
  initmu=init.model$initmu
  if(is.null(initmu)){Beta0=matrix(rnorm((p+1)*k,0,1),p+1,k);mu0=z%*%Beta0}else{mu0=initmu}
  ###
  alpha=init.model$initalpha
  if(is.null(alpha)) alpha=runif(k,0.85,1)
  eta=init.model$initeta
  if(is.null(eta)) eta=runif(k,1,1e2)
  sigma2=init.model$initsigma2
  if(is.null(sigma2)) sigma2=rgamma(k,2,1)
  prop=init.model$initprop
  if(is.null(prop)) prop=rep(1/k,k)
  
  LogLik0 <- sum(log(rowSums(sapply(1:k, function(j) {
    prop[j]*(alpha[j]*dnorm(y, mu0[,j], sqrt(sigma2[j])) +
               (1 - alpha[j]) * dnorm(y, mu0[,j], sqrt(eta[j]) * sqrt(sigma2[j])))}))))
  z <- matrix(0, n, k)
  v <- matrix(0, n, k)
  diff=1e6
  llk=NULL
  count=0
  time=NULL
  while(diff>1e-10){
    t0=Sys.time()
    ###E-step
    for (j in 1:k) {
      z[, j] <- prop[j] * (
        alpha[j] * dnorm(y, mu0[,j], sqrt(sigma2[j])) +
          (1 - alpha[j]) * dnorm(y, mu0[,j], sqrt(eta[j]) * sqrt(sigma2[j]))
      )+1e-100
    }
    z <- z / rowSums(z)
    
    # --- v-step ---
    for (j in 1:k) {
      denom <- alpha[j] * dnorm(y, mu0[,j], sqrt(sigma2[j])) +
        (1 - alpha[j]) * dnorm(y, mu0[,j], sqrt(eta[j]) * sqrt(sigma2[j]))
      v[, j] <- ifelse(denom > 0, (alpha[j] * dnorm(y, mu0[,j], sqrt(sigma2[j]))) / denom,0.5)
    }
    ###CM-step1
    prop1=mu1=sigma21=NULL
    for(j in 1:k){
      prob=cbind(z[,j]*(v[,j]+((1-v[,j])/eta[j])))
      mu=npksum(txdat=x,exdat=xgrid,tydat=y,weights=prob,bws=bw,ckertype = "epanechnikov")$ksum/npksum(txdat=x,exdat=xgrid,weights=prob,bws=bw,ckertype = "epanechnikov")$ksum
      mu1=cbind(mu1,mu)#multivariate_interpolate(x,xgrid,mu))
    }
    ##Checking for convergence
    LogLik1 <- sum(log(rowSums(sapply(1:k, function(j) {
      prop[j]*(alpha[j]*dnorm(y, mu1[,j], sqrt(sigma2[j])) +
                 (1 - alpha[j]) * dnorm(y, mu1[,j], sqrt(eta[j]) * sqrt(sigma2[j])))}))))
    diff=abs(LogLik1-LogLik0)/abs(LogLik0)
    LogLik0=LogLik1
    mu0=mu1
    llk=c(llk,LogLik1)
    count=count+1
    if(count==1e3) diff=1e-100 ##Set the max no. of iterations at 100
    time=c(time,difftime(Sys.time(),t0,units="mins"))
  }
  df=(4+p+1)*k-1
  BIC=-2*LogLik1+log(n)*df
  idx=order(colSums(mu1))
  prop=prop[idx];eta=eta[idx];alpha=alpha[idx];sigma2=sigma2[idx];mu1=mu1[,idx];
  resp1=sapply(1:k,function(j) 
    prop[j]*(alpha[j] * dnorm(y, mu1[,j], sqrt(sigma2[j]))+(1 - alpha[j]) * dnorm(y, mu1[,j], sqrt(eta[j]) * sqrt(sigma2[j]))
    ))
  gn1=resp1/rowSums(resp1)
  resp2=sapply(1:k,function(j){
    denom=alpha[j]*dnorm(y, mu1[,j], sqrt(sigma2[j]))+(1 - alpha[j]) * dnorm(y, mu1[,j], sqrt(eta[j]) * sqrt(sigma2[j]))
    ifelse(denom > 0, (alpha[j] * dnorm(y, mu1[,j], sqrt(sigma2[j]))) / denom,0.5)})
  gn2=resp2/rowSums(resp2)
  z1=apply(gn1,1,which.max)
  z2=sapply(1:n,function(i) as.numeric(gn2[i,z1[i]]>0.5))
  out=cbind(z1,z2,ifelse(z2==1,"good","bad"))
  return(list(resp=gn1,z_clust=z1,z_out=z2,out=out,prop=prop,sigma2=sigma2,mu=mu1,LL=LogLik0,BIC=BIC,LLiter=llk,time=time))
}

##Step 2: A function to perform the first step of Algorithm 3 for fitting a semi-parametric contaminated Gaussian mixture of regressions model
Algorithm1_semi_step2=function(y,x,k,bw=h,init.model=NULL,xgrid=NULL,mu.est=NULL){
  n=length(y)
  x=as.matrix(x);p=ncol(x)
  z=cbind(1,x)
  if(is.null(xgrid)) xgrid=sapply(1:p,function(j) seq(min(x[,j]),max(x[,j]),length.out=100))
  ###Initialize
  initeta=init.model$initeta
  if(is.null(initeta)){eta0=runif(k,1,1e3)}else{eta0=initeta}
  initalpha=init.model$initalpha
  if(is.null(initalpha)){alpha0=runif(k,0.85,1)}else{alpha0=initalpha}
  initprop=init.model$initprop
  if(is.null(initprop)){prop0=rep(1/k,k)}else{prop0=initprop}
  initsigma2=init.model$initsigma2
  if(is.null(initsigma2)){sigma20=rgamma(k,2,1)}else{sigma20=initsigma2}
  ###
  mu=mu.est
  LogLik0 <- sum(log(rowSums(sapply(1:k, function(j) {
    prop0[j]*(alpha0[j]*dnorm(y,mu[,j],sqrt(sigma20[j]))+(1-alpha0[j])*dnorm(y,mu[,j], sqrt(eta0[j]) * sqrt(sigma20[j])))}))))
  z <- matrix(0, n, k)
  v <- matrix(0, n, k)
  diff=1e6
  llk=NULL
  count=0
  time=NULL
  while(diff>1e-10){
    t0=Sys.time()
    ###E-step
    for (j in 1:k) {
      z[, j] <- prop0[j] * (
        alpha0[j] * dnorm(y, mu[,j], sqrt(sigma20[j])) +
          (1 - alpha0[j]) * dnorm(y, mu[,j], sqrt(eta0[j]) * sqrt(sigma20[j]))
      )
    }
    z <- z / rowSums(z)
    
    # --- v-step ---
    for (j in 1:k) {
      denom <- alpha0[j] * dnorm(y, mu[,j], sqrt(sigma20[j])) +
        (1 - alpha0[j]) * dnorm(y, mu[,j], sqrt(eta0[j]) * sqrt(sigma20[j]))
      v[, j] <- ifelse(denom > 0, (alpha0[j] * dnorm(y, mu[,j], sqrt(sigma20[j]))) / denom,0.5)
    }
    ###CM-step1
    prop1=colSums(z)/n
    alpha1=sigma21=NULL
    for(j in 1:k){
      prob=z[,j]*(v[,j]+((1-v[,j])/eta0[j]))
      #alpha1=c(alpha1,optimise(alpha_fun,c(0,1),weig1=z[,j],weig2=v[,j])$minimum)
      alpha1=c(alpha1,sum(z[,j]*v[,j])/sum(z[,j]))
      sigma2=sum(prob*(y-mu[,j])^2)/sum(z[,j])
      sigma21=c(sigma21,sigma2)
    }
    ###CM-step2
    eta1=NULL
    for(j in 1:k){
      a=sum(z[,j]*(1-v[,j]));b=sum((z[,j]*(1-v[,j]))*(y-mu[,j])^2)/sigma2[j]
      #res=optimise(eta_fun,interval=c(1,1e3),y=y,mu=mu[,j],sigma2=sigma20[j],weig1=z[,j],weig2=v[,j])
      #eta1=c(eta1,res$minimum)
      eta1=c(eta1,max(1.01,b/a))
    }
    ##Checking for convergence
    LogLik1 <- sum(log(rowSums(sapply(1:k, function(j) {
      prop1[j]*(alpha1[j]*dnorm(y, mu[,j], sqrt(sigma21[j])) +
                  (1 - alpha1[j]) * dnorm(y, mu[,j], sqrt(eta1[j]) * sqrt(sigma21[j])))}))))
    diff=abs(LogLik1-LogLik0)/abs(LogLik0)
    LogLik0=LogLik1
    eta0=eta1
    sigma20=sigma21
    prop0=prop1
    alpha0=alpha1
    llk=c(llk,LogLik1)
    count=count+1
    if(count==1e3) diff=1e-100 ##Set the max no. of iterations at 100
    time=c(time,difftime(Sys.time(),t0,units="mins"))
  }
  comp_llk=Q_function_semi(y,mu,prop1,alpha1,sigma21,eta1,z,v)
  IC=IC_semi_parametric(x,bw,k,LogLik0,comp_llk)
  idx=order(colSums(mu))
  prop1=prop1[idx];eta1=eta1[idx];alpha1=alpha1[idx];sigma21=sigma21[idx];mu=mu[,idx];
  resp1=sapply(1:k,function(j) 
    prop1[j]*(alpha1[j] * dnorm(y, mu[,j], sqrt(sigma21[j]))+(1 - alpha1[j]) * dnorm(y, mu[,j], sqrt(eta1[j]) * sqrt(sigma21[j]))
    ))
  gn1=resp1/rowSums(resp1)
  resp2=sapply(1:k,function(j){
    denom=alpha1[j] * dnorm(y, mu[,j], sqrt(sigma21[j]))+(1-alpha1[j])*dnorm(y, mu[,j], sqrt(eta1[j])*sqrt(sigma21[j]))
    ifelse(denom > 0, (alpha1[j] * dnorm(y, mu[,j], sqrt(sigma21[j]))) / denom,0.5)})
  gn2=resp2/rowSums(resp2)
  z1=apply(gn1,1,which.max)
  z2=sapply(1:n,function(i) as.numeric(gn2[i,z1[i]]>0.5))
  out=cbind(z1,z2,ifelse(z2==1,"good","bad"))
  return(list(resp=gn1,z_clust=z1,z_out=z2,alpha=alpha1,out=out,eta=eta1,prop=prop1,sigma2=sigma21,mu=mu,LL=LogLik0,IC=IC,LLiter=llk,time=time))
}

##A function to perform the first step of Algorithm 3 for fitting a semi-parametric contaminated Gaussian mixture of regressions model
CN_mix_semi_reg_alg1=function(y,x,k,bw=h,init.model=NULL,xgrid=NULL){
  fit1=Algorithm1_semi_step1(y,x,k,bw=bw,init.model = init.model,xgrid=x)
  fit2=Algorithm1_semi_step2(y,x,k,init.model=init.model,mu.est=fit1$mu,xgrid=x)
  n.iter=sum(length(fit1$LLiter),length(fit2$LLiter))
  time=sum(mean(fit1$time),mean(fit2$time))
  return(list(fit=fit2,n.iter=n.iter,time=time))
}

##A function to generate data from a variety of semi-parametric mixtures of regressions
Generate_data_semi=function(x,z,alpha,eta,mu,sigma2){
  if(!is.matrix(x)) x=as.matrix(x)
  n=nrow(x); p=ncol(x)
  ymat=matrix(0,nrow=n,5)
  xmat=array(0,dim=c(n,p,5)) 
  par(mfrow=c(2,2))
  ##Gaussian distributed errors##
  for(i in 1:n){
    ymat[i,1]=mu[i,z[i]]+rnorm(1,0,sqrt(sigma2[z[i]]))
  }
  ##t-distributed errors##
  for(i in 1:n){
    ymat[i,2]=mu[i,z[i]]+rt(1,5)
  }
  ##Contaminated normal errors##
  for(i in 1:n){
    ymat[i,3]=rCNorm(1,mu[i,z[i]],sigma2[z[i]],alpha[z[i]],eta[z[i]])
  }
  for(j in 1:3) xmat[,,j]=x
  ##Gaussian distributed errors with 5% data points substituted with (0,y), where y is generated from Uniform(10,15)##
  m=n*0.05
  idx=sample(1:n,m,replace=F)
  ymat[idx,4]=runif(m,10,15)
  xmat[idx,,4]=0.5
  for(i in setdiff(1:n,idx)){
    ymat[i,4]=mu[i,z[i]]+rnorm(1,0,sqrt(sigma2[z[i]]))
    xmat[i,,4]=x[i,]
  }
  ##Gaussian distributed errors with 5% y-values generated from the Uniform(-10,10)##
  m=n*0.05
  idx=sample(1:n,m,replace=F)
  ymat[idx,5]=runif(m,-10,10)
  xmat[idx,,5]=x[idx,]
  for(i in setdiff(1:n,idx)){
    ymat[i,5]=mu[i,z[i]]+rnorm(1,0,sqrt(sigma2[z[i]]))
    xmat[i,,5]=x[i,]
  }
  
  #for(j in 1:5) {plot(xmat[,j],ymat[,j])}
  return(list(ymat=ymat,xmat=xmat))
}

##A function to generate data from a variety of non-parametric mixtures of regressions
Generate_data_non=function(x,z,alpha,eta,mu,sigma2){
  if(!is.matrix(x)) x=as.matrix(x)
  n=nrow(x); p=ncol(x)
  ymat=matrix(0,nrow=n,5)
  xmat=array(0,dim=c(n,p,5)) 
  par(mfrow=c(2,2))
  ##Gaussian distributed errors##
  for(i in 1:n){
    ymat[i,1]=mu[i,z[i]]+rnorm(1,0,sqrt(sigma2[i,z[i]]))
  }
  ##t-distributed errors##
  for(i in 1:n){
    ymat[i,2]=mu[i,z[i]]+rt(1,5)
  }
  ##Contaminated normal errors##
  for(i in 1:n){
    ymat[i,3]=rCNorm(1,mu[i,z[i]],sigma2[i,z[i]],alpha[z[i]],eta[z[i]])
  }
  for(j in 1:3) xmat[,,j]=x
  ##Gaussian distributed errors with 5% data points substituted with (0,y), where y is generated from Uniform(10,15)##
  m=n*0.05
  idx=sample(1:n,m,replace=F)
  ymat[idx,4]=runif(m,10,15)
  xmat[idx,,4]=0.5
  for(i in setdiff(1:n,idx)){
    ymat[i,4]=mu[i,z[i]]+rnorm(1,0,sqrt(sigma2[i,z[i]]))
    xmat[i,,4]=x[i,]
  }
  ##Gaussian distributed errors with 5% y-values generated from the Uniform(-10,10)##
  m=n*0.05
  idx=sample(1:n,m,replace=F)
  ymat[idx,5]=runif(m,-10,10)
  xmat[idx,,5]=x[idx,]
  for(i in setdiff(1:n,idx)){
    ymat[i,5]=mu[i,z[i]]+rnorm(1,0,sqrt(sigma2[i,z[i]]))
    xmat[i,,5]=x[i,]
  }
  
  #for(j in 1:5) {plot(xmat[,j],ymat[,j])}
  return(list(ymat=ymat,xmat=xmat))
}

##A function to generate data to obtain initial values for the starting the fitting algorithms
initialize.model=function(x,y,k,method=NULL,true.init=NULL){
  n=length(y)
  x=as.matrix(x)
  BIC=1e6
  if(method=="1"){##Mixtures of Contaminated Gaussian regressions
    for(j in 1:10){
      m=list(BIC=1e6)
      try({m=CN_mix_reg(y,x,k=k,init.model=NULL)},silent=T)
      if(m$BIC<BIC){init.model=m;BIC=m$BIC}
    }
    init.model=list(initmu=m$mu,initsigma2=m$sigma2,initprop=m$prop,initeta=m$eta,initalpha=m$alpha)
  }
  if(method=="2"){##Mixtures of non-parametric Contaminated Gaussian regressions
    for(j in 1:1){
      m=list(BIC=1e6)
      try({m=Kernel_Mix_EM(x,y,k=k,bw=silverman_bandwidth(x),xgrid=x)},silent=T)
      if(m$IC$BIC<BIC){init.model=m;BIC=m$IC$BIC}
    }
    init.model=list(initmu=m$mu,initsigma2=m$sigma2,initprop=m$prop,initeta=m$eta,initalpha=m$alpha)
  }
  if(method=="3"){##True values
    m0=true.init
    init.model=list(initmu=m0$mu,initsigma2=m0$sigma2,initprop=m0$prop,initeta=m0$eta,initalpha=m0$alpha)
  }
  if(method=="4"){##Random starting point
    p=ncol(x)
    initeta=runif(k,1,1e3)
    initalpha=runif(k,0.85,1)
    initprop=rep(1/k,k)
    initsigma2=rgamma(k,2,1)
    Beta=matrix(rnorm((p+1)*k,0,1),p+1,k);initmu=cbind(1,x)%*%Beta
    init.model=list(initmu=initmu,initsigma2=initsigma2,initprop=initprop,initeta=initeta,initalpha=initalpha)
  }
  return(init.model)
}