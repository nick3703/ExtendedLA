library(spdep)
library(mvtnorm)
library(sparseMVN)
library(psych)
library(lattice)
library(spdep)
library(LaplacesDemon)

#Functions

##########Function to Fit########

#Size of Spatial Domain and Number of Temporal Observations
size<-100
num.obs<-100

#Data Structures Needed for Functions
mynb <- cell2nb(sqrt(size), sqrt(size), type="rook",torus=TRUE)
spat<-nb2mat(mynb,style="B")
spat<-Matrix(spat,sparse=T)
iden<-Matrix(0,nrow=num.obs,ncol=num.obs)
diag(iden)<-1
spat.Neigh<- iden %x% spat
diag(spat.Neigh)=0
spat<-Matrix(spat,sparse=T)
spat.Neigh<- diag(1,num.obs) %x% spat
time.Neigh<-Matrix(0,nrow=(size*num.obs),ncol=(size*num.obs))
for(i in 1:(ncol(time.Neigh)-size)){
  time.Neigh[i+size,i]<-1
}
iden<-Matrix(0,ncol=ncol(spat.Neigh),nrow=nrow(spat.Neigh))
diag(iden)<-1
M<-iden
#######################


margs.higher.order.full<-function(t1,t2,eta,intercept,include6=1){
  spat.Neigh.f<-t2*spat.Neigh
  b<- rep(intercept,(num.obs*size))
  diag(M)<-1/t1
  Q<-Matrix((iden-spat.Neigh.f)%*%M,sparse=T)
  dat.old<-time.Neigh%*%y
  dat.old[1:size,]<-attacks[,1]
  dat.ex<-y
  l.i<-function(z){
    -exp(z+b)-eta*dat.old+dat.ex*log(exp(z+b)+eta*dat.old)-lfactorial(dat.ex)
  }
  lam<-function(z){
    exp(z+b)+eta*dat.old
  }
  #Functions to get first and second partial derivatives
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
  while(diff>.05 & count < 50 & diff<10000000){
    mu.old<-mu
    mat<-as.numeric(c.i(mu))
    m<-Matrix(0,nrow=(size*num.obs),ncol=(size*num.obs))
    diag(m)<-mat
    mu<-solve(Q+m,b.i(mu))
    diff<-sum((mu-mu.old)^2)
    count<-count+1
  }
  x.star<-mu
  V.j<-Matrix(0,nrow=size*num.obs,ncol=size*num.obs)
  place.hold<-Matrix(0,nrow=size*num.obs,ncol=size*num.obs)
  diag(place.hold)<- -second.i(mu)
  V.j<-Matrix(Q+place.hold,sparse=T)
  covs<-solve(V.j)
  T4.v.term <- diag(covs)^2
  T4.h.term<- fourth.i(x.star)
  if(include6==1){
    T6.v.term=0
    T6.h.term=0
  }else{
    T6.v.term <- diag(covs)^4
    T6.h.term <- sixth.i(x.star) 
  }
  T3.h.term<- third.i(x.star)
  sigs<- diag(covs)
  third.term<-c()
  fourth.term<-c()
  priors=dnorm(intercept,0,100,log=TRUE)+dhalfcauchy(t1,scale=5,log=TRUE) #Uniform on eta and zeta
  for(o in 1:num.obs){
    h111.h222<-T3.h.term[((o-1)*size+1):(o*size)]%*%t(T3.h.term[((o-1)*size+1):(o*size)])
    h11.h22<-sigs[((o-1)*size+1):(o*size)]%*%t(sigs[((o-1)*size+1):(o*size)])
    third.term[o]<-sum(9/72*covs[((o-1)*size+1):(o*size),((o-1)*size+1):(o*size)]*h11.h22*h111.h222)
    fourth.term[o]<-sum(6/72*h111.h222*covs[((o-1)*size+1):(o*size),((o-1)*size+1):(o*size)]^3)
  }
  output<- list(log.lik=1/2*as.numeric(determinant(Q/(2*pi))$modulus)-1/2*as.numeric(determinant(V.j/(2*pi))$modulus)+sum(l.i(mu))-1/2*t(mu)%*%Q%*%mu+sum(-1/8*T4.v.term*T4.h.term)+sum(third.term)+sum(fourth.term)+sum(-1/48*T6.h.term*T6.v.term)+priors,
                r1s=sum(sum(-1/8*T4.v.term*T4.h.term)),r2s=sum(third.term),r3s=sum(fourth.term),r4s=sum(-1/48*T6.h.term*T6.v.term),k=diag(place.hold),l=V.j,m=covs,n=T4.h.term,q=T3.h.term,r=mu,s=(exp(mu+b)*(1/lam(mu)-1)),log.lik2=1/2*as.numeric(determinant(Q/(2*pi))$modulus)-1/2*as.numeric(determinant(V.j/(2*pi))$modulus)+sum(l.i(mu))-1/2*t(mu)%*%Q%*%mu)
  return(output)
}


#############Finite Difference Approximation to Hessian##############################


Hess.ex<-function(t1,t2,eta,b,t1step,t2step,etastep,bstep,cen){
  min<- -300001
  while(min < -299999){
    t1.for<-margs.higher.order.full(t1+t1step,t2,eta,b)
    t1.bac<-margs.higher.order.full(t1-t1step,t2,eta,b)
    t2.for<-margs.higher.order.full(t1,t2+t2step,eta,b)
    t2.bac<-margs.higher.order.full(t1,t2-t2step,eta,b)
    eta.for<-margs.higher.order.full(t1,t2,eta+etastep,b)
    eta.bac<-margs.higher.order.full(t1,t2,eta-etastep,b)
    b.for<-margs.higher.order.full(t1,t2,eta,b+bstep)
    b.bac<-margs.higher.order.full(t1,t2,eta,b-bstep)
    t1.for.t2.for<-margs.higher.order.full(t1+t1step,t2+t2step,eta,b)
    t1.for.t2.bac<-margs.higher.order.full(t1+t1step,t2-t2step,eta,b)
    t1.bac.t2.for<-margs.higher.order.full(t1-t1step,t2+t2step,eta,b)
    t1.bac.t2.bac<-margs.higher.order.full(t1-t1step,t2-t2step,eta,b)
    t1.for.eta.for<-margs.higher.order.full(t1+t1step,t2,eta+etastep,b)
    t1.for.eta.bac<-margs.higher.order.full(t1+t1step,t2,eta-etastep,b)
    t1.bac.eta.for<-margs.higher.order.full(t1-t1step,t2,eta+etastep,b)
    t1.bac.eta.bac<-margs.higher.order.full(t1-t1step,t2,eta-etastep,b)
    t1.for.b.for<-margs.higher.order.full(t1+t1step,t2,eta,b+bstep)
    t1.for.b.bac<-margs.higher.order.full(t1+t1step,t2,eta,b-bstep)
    t1.bac.b.for<-margs.higher.order.full(t1-t1step,t2,eta,b+bstep)
    t1.bac.b.bac<-margs.higher.order.full(t1-t1step,t2,eta,b-bstep)
    t2.for.eta.for<-margs.higher.order.full(t1,t2+t2step,eta+etastep,b)
    t2.for.eta.bac<-margs.higher.order.full(t1,t2+t2step,eta-etastep,b)
    t2.bac.eta.for<-margs.higher.order.full(t1,t2-t2step,eta+etastep,b)
    t2.bac.eta.bac<-margs.higher.order.full(t1,t2-t2step,eta-etastep,b)
    t2.for.b.for<-margs.higher.order.full(t1,t2+t2step,eta,b+bstep)
    t2.for.b.bac<-margs.higher.order.full(t1,t2+t2step,eta,b-bstep)
    t2.bac.b.for<-margs.higher.order.full(t1,t2-t2step,eta,b+bstep)
    t2.bac.b.bac<-margs.higher.order.full(t1,t2-t2step,eta,b-bstep)
    eta.for.b.for<-margs.higher.order.full(t1,t2,eta+etastep,b+bstep)
    eta.for.b.bac<-margs.higher.order.full(t1,t2,eta+etastep,b-bstep)
    eta.bac.b.for<-margs.higher.order.full(t1,t2,eta-etastep,b+bstep)
    eta.bac.b.bac<-margs.higher.order.full(t1,t2,eta-etastep,b-bstep)
    H<-Matrix(0,nrow=4,ncol=4,sparse=FALSE)
    H[1,1]<-(t1.for$log.lik-2*cen+t1.bac$log.lik)/(t1step^2)
    H[2,2]<-(t2.for$log.lik-2*cen+t2.bac$log.lik)/(t2step^2)
    H[3,3]<-(eta.for$log.lik-2*cen+eta.bac$log.lik)/(etastep^2)
    H[4,4]<-(b.for$log.lik-2*cen+b.bac$log.lik)/(bstep^2)
    H[1,2]<-H[2,1]<-(t1.for.t2.for$log.lik-t1.for.t2.bac$log.lik-t1.bac.t2.for$log.lik+t1.bac.t2.bac$log.lik)/(4*t2step*t1step)
    H[1,3]<-H[3,1]<-(t1.for.eta.for$log.lik-t1.for.eta.bac$log.lik-t1.bac.eta.for$log.lik+t1.bac.eta.bac$log.lik)/(4*etastep*t1step)
    H[1,4]<-H[4,1]<-(t1.for.b.for$log.lik-t1.for.b.bac$log.lik-t1.bac.b.for$log.lik+t1.bac.b.bac$log.lik)/(4*bstep*t1step)
    H[2,3]<-H[3,2]<-(t2.for.eta.for$log.lik-t2.for.eta.bac$log.lik-t2.bac.eta.for$log.lik+t2.bac.eta.bac$log.lik)/(4*t2step*etastep)
    H[2,4]<-H[4,2]<-(t2.for.b.for$log.lik-t2.for.b.bac$log.lik-t2.bac.b.for$log.lik+t2.bac.b.bac$log.lik)/(4*bstep*t2step)
    H[3,4]<-H[4,3]<-(eta.for.b.for$log.lik-eta.for.b.bac$log.lik-eta.bac.b.for$log.lik+eta.bac.b.bac$log.lik)/(4*bstep*etastep)
    min1<-min(as.numeric(t1.for$log.lik),as.numeric(t1.bac$log.lik),as.numeric(t2.for$log.lik),as.numeric(t2.bac$log.lik),as.numeric(eta.for$log.lik),as.numeric(eta.bac$log.lik))
    min2<-min(as.numeric(t1.for.t2.for$log.lik),as.numeric(t1.for.t2.bac$log.lik),as.numeric(t1.for.eta.for$log.lik),as.numeric(t1.for.eta.bac$log.lik),
              as.numeric(t1.bac.t2.for$log.lik),as.numeric(t1.bac.t2.bac$log.lik),
              as.numeric(t1.bac.eta.for$log.lik),as.numeric(t1.bac.eta.bac$log.lik),
              as.numeric(t2.for.eta.for$log.lik),as.numeric(t2.for.eta.bac$log.lik),
              as.numeric(t2.bac.eta.for$log.lik),as.numeric(t2.bac.eta.bac$log.lik))
    prevstep<-t1step
    t2step<-t2step/2
    t1step<-etastep<-prevstep/2
    #print(t1step)
    min<-min(min1,min2)
    output<-list(t1.p=(t1.for$log.lik-t1.bac$log.lik)/(2*t1step),t2.p=(t2.for$log.lik-t2.bac$log.lik)/(2*t2step),eta.p=(eta.for$log.lik-eta.bac$log.lik)/(2*etastep),b.p=(b.for$log.lik-b.bac$log.lik)/(2*bstep),H=H)
    
  }
  return(output)
}


#Parameters to Adjust
t1<-.7
t2<-.245

eta<-.1
etacross<-0
alpha<- rep(0,size)
#####
#Simulate Data
spat.Neigh.f<-t2*spat.Neigh
iden<-Matrix(0,ncol=ncol(spat.Neigh),nrow=nrow(spat.Neigh))
diag(iden)<-1
L<-iden
diag(L)<-1/t1
Q2<-Matrix(((iden-spat.Neigh.f)%*%L),sparse=T)
vecs<-c(0,1,rep(0,(size*num.obs-2)))
t<-solve(Q2,vecs)
c.Q<-Cholesky(Q2)
mu<-rep(alpha,(num.obs))
test<-rmvn.sparse(1,mu=mu,CH=c.Q,TRUE)
theta<-eta
counts<-num.obs
attacks<-matrix(NA,nrow=size,ncol=counts+1)
tot.resp<-matrix(NA,nrow=size,ncol=counts+1)
attacks[,1]<-0
tot.resp<-attacks2<-attacks
for(i in 1:counts){
  tot.resp[,i+1]<-exp(test[((i-1)*size+1):(i*size)])+theta*attacks[,i]+etacross*tot.resp[,i]
  attacks[,i+1]<-rpois(size,tot.resp[,i+1])
}

y<-rep(0,(size*num.obs))
count<-1
for(i in 2:(num.obs+1)){
  for(j in 1:size){
    y[count]<-attacks[j,i]
    count<-count+1
  }
}

dat<-y
#####################INITIAL VALUES FOR NEWTON-RAPHSON BASED ON SAMPLE MOMENTS


cors1<-c()

v2=c()
for(k in 1:size){
 cors1[k]=acf(attacks[k,],plot=F)[[1]][2]
}

fish.correct<-fisherz(cors1)
fish.mean<-mean(fish.correct)
eta.hat<-fisherz2r(fish.mean)
fish.correct<-fisherz(cors1)
fish.mean<-mean(fish.correct)
eta.hat<-fisherz2r(fish.mean)
var.expect<-1/(1-eta^2)*(exp(2*alpha)*(exp(2*t[2])-exp(t[2]))+exp(alpha)/(1-eta)*exp(t[2]/2))
v<-var(as.numeric(attacks))
a.hat<- -log(v*(1-eta.hat^2)-mean(attacks)+((1-eta.hat)*mean(attacks))^2)/2+2*log((1-eta.hat)*mean(attacks))
sig.11.hat<-log(mean(attacks)*(1-eta.hat))*2-2*a.hat
cor.mat<-cor(t(attacks)) #Ignore warning
spat.cor.mat<-lower.tri(spat)*spat*cor.mat
spat.cor.mat.vec<-as.numeric(spat.cor.mat)
spat.cor.mat.vec<-spat.cor.mat.vec[spat.cor.mat.vec!=0]
fish.correct.spat<-fisherz(spat.cor.mat.vec)
fish.spat.mean<-mean(fish.correct.spat)
spat.cor<-fisherz2r(fish.spat.mean)
diag(spat)=0
spat.f<-.248*spat
iden.mom<-Matrix(0,ncol=ncol(spat),nrow=nrow(spat))
diag(iden.mom)<-1
sig.12.hat<-log(spat.cor*(exp(2*sig.11.hat)-exp(sig.11.hat)+(exp(-a.hat)/(1-eta.hat))*exp(sig.11.hat/2))+exp(sig.11.hat))-sig.11.hat
sig.row.est=sig.12.hat*spat[1,]
sig.row.est[1]=sig.11.hat
Q.matrix<-Matrix(((iden.mom-spat.f)),sparse=T)
t<-solve(solve(Q.matrix),sig.row.est)
intercept.init=a.hat
eta.init=eta.hat
spat.init=.248
sig.init=t[1]
#####################################################


#Newton Raphson
old=matrix(c(sig.init,spat.init,eta.init-.1,intercept.init))

if(old[3,1]<0){old[3,1]=.01}

ptm=proc.time() #For Timing

diff<-1
m<-1
q<-1
while(q <10 & (abs(diff)>.005)){
  if(diff>.0005){
    c=margs.higher.order.full(old[1,1],old[2,1],old[3,1],old[4,1])
    cen=c$log.lik
    H<-Hess.ex(old[1,1],old[2,1],old[3,1],old[4,1],.001,.001,.001,.001,cen=cen)
    grad<-Matrix(c(as.numeric(H$t1.p),as.numeric(H$t2.p),as.numeric(H$eta.p),as.numeric(H$b.p)),nrow=4)
    Hess.est<-H$H
  }
  x.new<-old-1/m*solve(Hess.est)%*%grad
  if(x.new[1,1]<0){x.new[1,1]=.1}
  if(x.new[2,1]<0){x.new[2,1]=.01}
  if(x.new[2,1]>=.249){x.new[2,1]=.248}
  if(x.new[3,1]>1){x.new[3,1]=.9}
  if(x.new[3,1]<0){x.new[3,1]=.01}
  new.lik<-margs.higher.order.full(x.new[1],x.new[2],x.new[3],x.new[4])
  diff<-as.numeric(new.lik$log.lik-cen)
  if(diff<0) {
    m<-2*m
    old<-old  }
  if(diff>=0) {
    old<-x.new
    m<-1
    q<-q+1
    
  }
}

proc.time() - ptm

#############################################################################

old+1.96*sqrt(diag(solve(-Hess.est)))
old-1.96*sqrt(diag(solve(-Hess.est)))



