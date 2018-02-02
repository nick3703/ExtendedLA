#REQUIRED LIBRARIES#############
library(maptools)
library(spdep)
library(zoo)
library(LaplacesDemon)
library(maptools)
library(Matrix)

#READ IN DATA################

N.mat=readMM("https://raw.githubusercontent.com/nick3703/ExpandedLA/master/chineigh.mtx")
dat.in.format.chi2<-read.csv("https://raw.githubusercontent.com/nick3703/ExpandedLA/master/Chi_Violent_Crimes.csv")
pop.in.format<-read.csv("https://raw.githubusercontent.com/nick3703/ExpandedLA/master/Chi_Population.csv")
dat<-dat.in.format.chi2$count


#Create Space and Time Neighborhoods
sp.size<-size<-length(unique(dat.in.format.chi2$community))
t.size<-num.obs<-length(unique(dat.in.format.chi2$wk))
spat.Neigh<- diag(1,t.size) %x% N.mat
time.Neigh<-Matrix(0,nrow=nrow(spat.Neigh),ncol=ncol(spat.Neigh),sparse=TRUE)
for(i in 1:(ncol(time.Neigh)-ncol(N.mat))){
  time.Neigh[i+nrow(N.mat),i]<-1
}


log.pop<-rep(log(pop.in.format$pop),(t.size))
wx<-read.csv("https://raw.githubusercontent.com/nick3703/ExpandedLA/master/Chi_Wx_2015.csv")
wx.as.covariate<-rep(wx$X32,each=sp.size)
wx.as.covariate<-(wx.as.covariate-mean(wx.as.covariate))/sqrt(var(wx.as.covariate))


iden2<-Matrix(0,nrow=(sp.size*t.size),ncol=(sp.size*t.size),sparse=T)
diag(iden2)<-1
iden<-Matrix(0,ncol=ncol(spat.Neigh),nrow=nrow(spat.Neigh))
diag(iden)<-1

M<-iden

y<-dat
########################FUNCTIONS###############################################################

#Standard Laplace Approximation
margs.chi<-function(params){
  t1=params[1]
  t2=params[2]
  eta=params[3]
  intercept=params[4]
  b2=params[5]
  b3=params[6]
  spat.Neigh.f<-t2*spat.Neigh
  b<- rep(intercept,(num.obs*size))+b2*wx.as.covariate+b3*log.pop
  diag(M)<-1/t1
  time.Neigh.f<-time.Neigh*eta
  Q<-Matrix((iden-spat.Neigh.f)%*%M,sparse=T)
  past.effect<-time.Neigh.f%*%y
  past.effect[1:size,]<-0
  
  #Functions to get first and second partial derivatives
  first.i<-function(z){
    (y*exp(z+b))/(exp(z+b)+past.effect)-exp(z+b)
  }
  second.i<-function(z){
    (y*exp(z+b))/(exp(z+b)+past.effect)-(exp(2*z+2*b)*y)/((exp(z+b)+past.effect)^2)-exp(z+b)
  }
  b.i<-function(z){
    first.i(z)-z*second.i(z)
  }
  c.i<-function(z){
    -second.i(z)
  }
  mu<-rep(0,(size*num.obs))
  mu.old<-mu
  mat<-as.numeric(c.i(mu))
  m<-Matrix(0,nrow=(size*num.obs),ncol=(size*num.obs))
  diag(m)<-mat
  Q.test<-Q+m
  mu<-rep(4,(size*num.obs))
  diff<-3
  count<-1
  #Optimization Step
  while(diff>.005 & count < 50 & diff<10000000){
    mu.old<-mu
    mat<-as.numeric(c.i(mu))
    m<-Matrix(0,nrow=(size*num.obs),ncol=(size*num.obs))
    diag(m)<-mat
    mu<-solve(Q+m,b.i(mu))
    diff<-sum((mu-mu.old)^2)
    count<-count+1
  }
  x.star<-mu
  d.star<-Matrix(Q+m,sparse=T)
  if(diff<.05){
    Q.det<-determinant(Q/(2*pi))
    det.d.star<-determinant(d.star/(2*pi))
    d.inv<-solve(d.star) #Only needed if effective number of obs calculation from Rue et. all is required or else can omit
    output<-list(log.lik=-1/2*t(x.star)%*%Q%*%(x.star)+sum(y*log(exp(x.star+b)+past.effect)-exp(x.star+b)-past.effect-lfactorial(y))-1/2*as.numeric(det.d.star$modulus)+1/2*as.numeric(Q.det$modulus),xi=x.star,sig=d.star,pd=(size*num.obs)-sum(diag(Q.test%*%d.inv)))
    return(as.numeric(output$log.lik))
  }
  else{
    #Something went wrong if log.likelihood is -300000
    return((log.lik=-300000))
  }
}


#Higher Order Laplace Approximation
margs.higher.order.full<-function(params,giicheck=0){
  t1=params[1]
  t2=params[2]
  eta=params[3]
  intercept=params[4]
  b2=params[5]
  b3=params[6]
  spat.Neigh.f<-t2*spat.Neigh
  b<- rep(intercept,(num.obs*size))+b2*wx.as.covariate+b3*log.pop
  diag(M)<-1/t1
  Q<-Matrix((iden-spat.Neigh.f)%*%M,sparse=T)
  dat.old<-time.Neigh%*%y
  dat.old[1:size,]<-0
  dat.ex<-y
  l.i<-function(z){
    -exp(z+b)-eta*dat.old+dat.ex*log(exp(z+b)+eta*dat.old)-lfactorial(dat.ex)
  }
  lam<-function(z){
    exp(z+b)+eta*dat.old
  }
  #Functions to get first through sixth derivatives
  first.i<-function(z){
    ((dat.ex*exp(z+b))/(exp(z+b)+eta*dat.old)-exp(z+b))
  }
  second.i<-function(z){
    ((dat.ex*exp(z+b))/(exp(z+b)+eta*dat.old)-(exp(2*z+2*b)*dat.ex)/((exp(z+b)+eta*dat.old)^2)-exp(z+b))
  }
  third.i<-function(z){
    -(exp(z+b)*(dat.ex/lam(z)-1)-3*exp(2*z+2*b)*(dat.ex/lam(z)^2)+2*exp(3*z+3*b)*(dat.ex/lam(z)^3))
  }
  fourth.i<-function(z){
    -(exp(z+b)*(dat.ex/lam(z)-1)-7*exp(2*z+2*b)*(dat.ex/lam(z)^2)+12*exp(3*z+3*b)*(dat.ex/lam(z)^3)-6*exp(4*z+4*b)*(dat.ex/lam(z)^4))
  }
  sixth.i<-function(z){
    -(exp(z+b)*(dat.ex/lam(z)-1)-31*exp(2*z+2*b)*(dat.ex/lam(z)^2)+180*exp(3*z+3*b)*(dat.ex/lam(z)^3)-438*exp(4*z+4*b)*(dat.ex/lam(z)^4)+
        408*exp(5*z+5*b)*(dat.ex/lam(z)^5)-120*exp(6*z+6*b)*(dat.ex/lam(z)^6))
  }
  b.i<-function(z){
    first.i(z)-z*second.i(z)
  }
  c.i<-function(z){
    -second.i(z)
  }
  mu<-rep(10,(size*num.obs))
  diff<-3
  count<-1
  #Same Optimization as in standard Laplace Approximation
  while(diff>.005 & count < 50 & diff<10000000){
    mu.old<-mu
    mat<-as.numeric(c.i(mu))
    m<-Matrix(0,nrow=(size*num.obs),ncol=(size*num.obs))
    diag(m)<-mat
    mu<-solve(Q+m,b.i(mu))
    diff<-sum((mu-mu.old)^2)
    count<-count+1
  }
  if(diff<.005){
  x.star<-mu
  V.j<-Matrix(0,nrow=size*num.obs,ncol=size*num.obs)
  place.hold<-Matrix(0,nrow=size*num.obs,ncol=size*num.obs)
  diag(place.hold)<- -second.i(mu)
  V.j<-Matrix(Q+place.hold,sparse=T)
  covs<-solve(V.j) #Invert sparse Matrix
  T4.v.term <- diag(covs)^2
  T4.h.term<- fourth.i(x.star)
  T6.v.term <- diag(covs)^4
  T6.h.term <- sixth.i(x.star)
  T3.h.term<- third.i(x.star)
  sigs<- diag(covs)
  third.term<-c()
  fourth.term<-c()
  priors=dnorm(intercept,0,1000,log=T)+dnorm(b2,0,1000,log=T)+dnorm(b3,0,1000,log=T)+dhalfcauchy(sqrt(t1),5,log=T)
  #Calculate additional terms needed in higher order Laplace Approximation
  for(o in 1:num.obs){
    h111.h222<-T3.h.term[((o-1)*size+1):(o*size)]%*%t(T3.h.term[((o-1)*size+1):(o*size)]) #Evaluation of all combinations of g_{iii}g_{jjj}
    h11.h22<-sigs[((o-1)*size+1):(o*size)]%*%t(sigs[((o-1)*size+1):(o*size)]) #Evaluation of g^{ii}g^{jj}
    third.term[o]<-sum(9/72*covs[((o-1)*size+1):(o*size),((o-1)*size+1):(o*size)]*h11.h22*h111.h222) #Evaluation of 9/72 g^{ii}g^{jj}g^{ij}g_{iii}g_{jjj} from (4.19)
    fourth.term[o]<-sum(6/72*h111.h222*covs[((o-1)*size+1):(o*size),((o-1)*size+1):(o*size)]^3) #Evaluation of 6/72 g^{ij}^3 g_{iii}g_{jjj}
  }
  output<- list(log.lik=1/2*as.numeric(determinant(Q/(2*pi))$modulus)-1/2*as.numeric(determinant(V.j/(2*pi))$modulus)+sum(l.i(mu))-1/2*t(mu)%*%Q%*%mu+sum(-1/8*T4.v.term*T4.h.term)+sum(third.term)+sum(fourth.term)+priors,
                r1s=sum(sum(-1/8*T4.v.term*T4.h.term)),r2s=sum(third.term),r3s=sum(fourth.term),r4s=sum(sum(-1/48*T6.v.term*T6.h.term)),gii=sigs)
  if(giicheck==0){  
  return(as.numeric(output$log.lik))
  }else{
    return(output$gii)
  }
  }else{
  return(log.lik=-300000)
}
}




#Newton-Raphson
ptm <- proc.time()

old<-matrix(c(.2, .183, .30,-6, .13, .50),nrow=6)


diff<-1
m<-1
q<-1
while(q <10 & (abs(diff)>.005)){
  if(diff>.005){
    c<-margs.higher.order.full(old)
    cen=c
    Hess.est=hessian(func=margs.higher.order.full,x=old,method.args=list(d=.001,r=2))
    g=grad(func=margs.higher.order.full,x=old,method.args=list(d=.001,r=2))
  }
  x.new<-old-1/m*solve(Hess.est)%*%g
  if(x.new[1,1]<0){x.new[1,1]=.1}
  if(x.new[2,1]<0){x.new[2,1]=.01}
  if(x.new[2,1]>.182){x.new[2,1]=.182}
  if(x.new[3,1]>1){x.new[3,1]=.9}
  if(x.new[3,1]<0){x.new[3,1]=.01}
  new.lik<-margs.higher.order.full(x.new)
  diff<-as.numeric(new.lik-cen)
  if(diff<0) {
    m<-2*m
    old<-old  }
  if(diff>=0) {
    old<-x.new
    m<-1
    q<-q+1
    print(old)
    print(cen)
    
    
  }
}

proc.time()-ptm

mode.gii=margs.higher.order.full(old,giicheck=1)

hist(mode.gii,xlab=expression(hat(g)^ii~' terms'),main=expression(hat(g)^ii~ 'Terms Calculated at Posterior Mode Using Extended LA'))

abline(v=1,col="red",lwd=3)



x.new-1.96*sqrt(diag(solve(-Hess.est)))
x.new+1.96*sqrt(diag(solve(-Hess.est)))






diff<-1
m<-1
q<-1
while(q <10 & (abs(diff)>.005)){
  if(diff>.005){
    cen<-margs.chi(old)
    Hess.est=hessian(func=margs.chi,x=old,method.args=list(d=.001,r=4))
    g=grad(func=margs.chi,x=old,method.args=list(d=.001,r=2))
    }
  x.new<-old-1/m*solve(Hess.est)%*%g
  if(x.new[1,1]<0){x.new[1,1]=.1}
  if(x.new[2,1]<0){x.new[2,1]=.01}
  if(x.new[2,1]>.182){x.new[2,1]=.182}
  if(x.new[3,1]>1){x.new[3,1]=.9}
  if(x.new[3,1]<0){x.new[3,1]=.01}
  new.lik<-margs.chi(x.new)
  diff<-as.numeric(new.lik-cen)
  if(diff<0) {
    m<-2*m
    old<-old  }
  if(diff>=0) {
    old<-x.new
    m<-1
    q<-q+1
    print(old)
    print(cen)
    
    
  }
}


x.new #Clearly the first parameter (corresponding to \sigma^2_{sp} is different)

x.new+1.96*sqrt(diag(-solve(Hess.est)))
x.new-1.96*sqrt(diag(-solve(Hess.est)))

#Minimal differences in other parameters
