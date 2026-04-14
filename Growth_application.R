library(np)
library(mixtools)
library(mclust)
library(dplyr)
data("oecdpanel")
attach(oecdpanel)

data=oecdpanel%>%filter(year==1965)%>%mutate(across(c(growth,initgdp,popgro,inv,humancap),standardizer))
attach(data)

X=cbind(initgdp,popgro,inv,humancap)
y=growth
n=nrow(X)

####
set.seed(135)
for(j in 0:1){
hist(data$growth[oecdpanel$oecd==j & oecdpanel$year==1965],main="",xlab="GDP growth rate per capita")
}

x=X

fit1=fit2=fit3=fit4=NULL
for(k in 2){
try({
init.model0=initialize.model(x,y,k,method=4,true.init = true.model)
h=rep(1.6,4)
fit0=regmixEM(y,x,k=k);df=k*(ncol(fit0$x))+2*k-1;BIC0=-2*fit0$loglik+log(length(x))*df
fit1=CN_mix_semi_reg(y,x,k=k,bw=h,init.model=init.model0,xgrid=x)
init.model00=list(initmu=fit1$mu,initeta=fit1$eta,initalpha=fit1$alpha,initsigma2=fit1$sigma2)
fit2=CN_mix_non_reg(y,x,k=k,bw=h,init.model=init.model0,xgrid=x)
fit3=Kernel_Mix_EM_semi(x,y,k=k,bw=h,init.model=init.model0,xgrid=x)
fit4=Kernel_Mix_EM(x,y,k=k,bw=h,init.model=init.model0,xgrid=x)
})

print(fit1$IC)
print(fit2$IC)
print(fit3$IC)
print(fit4$IC)
}

classError(data$oecd,ifelse(fit1$resp[,1]>=0.5,1,2))$errorRate
classError(data$oecd,ifelse(fit2$resp[,1]>=0.5,1,2))$errorRate
classError(data$oecd,ifelse(fit3$resp[,1]>=0.5,1,2))$errorRate
classError(data$oecd,ifelse(fit4$resp[,1]>=0.5,1,2))$errorRate
