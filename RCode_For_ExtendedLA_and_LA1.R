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
margs.chi<-function(t1,t2,eta,intercept,b2,b3){
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
    return(output)
  }
  else{
    #Something went wrong if log.likelihood is -300000
    return(list(log.lik=-300000))
  }
}


#Higher Order Laplace Approximation
margs.higher.order.full<-function(t1,t2,eta,intercept,b2,b3){
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
  #Calculate additional terms needed in higher order Laplace Approximation
  for(o in 1:num.obs){
    h111.h222<-T3.h.term[((o-1)*size+1):(o*size)]%*%t(T3.h.term[((o-1)*size+1):(o*size)]) #Evaluation of all combinations of g_{iii}g_{jjj}
    h11.h22<-sigs[((o-1)*size+1):(o*size)]%*%t(sigs[((o-1)*size+1):(o*size)]) #Evaluation of g^{ii}g^{jj}
    third.term[o]<-sum(9/72*covs[((o-1)*size+1):(o*size),((o-1)*size+1):(o*size)]*h11.h22*h111.h222) #Evaluation of 9/72 g^{ii}g^{jj}g^{ij}g_{iii}g_{jjj} from (4.19)
    fourth.term[o]<-sum(6/72*h111.h222*covs[((o-1)*size+1):(o*size),((o-1)*size+1):(o*size)]^3) #Evaluation of 6/72 g^{ij}^3 g_{iii}g_{jjj}
  }
  output<- list(log.lik=1/2*as.numeric(determinant(Q/(2*pi))$modulus)-1/2*as.numeric(determinant(V.j/(2*pi))$modulus)+sum(l.i(mu))-1/2*t(mu)%*%Q%*%mu+sum(-1/8*T4.v.term*T4.h.term)+sum(third.term)+sum(fourth.term),
                r1s=sum(sum(-1/8*T4.v.term*T4.h.term)),r2s=sum(third.term),r3s=sum(fourth.term),r4s=sum(sum(-1/48*T6.v.term*T6.h.term)),gii=sigs)
  return(output)
  }else{
  return(list(log.lik=-300000))
}
}

#Optimization Step via Newton-Raphson
#Function to find Hessian (Could be sped up through BFGS)

Hess.full<-function(t1,t2,eta,b1,b2,b3,t1step,t2step,etastep,b1step,b2step,b3step,cen=cen){
  min<- -300000
  while(min < -200000){
    t1.for<-margs.higher.order.full(t1+t1step,t2,eta,b1,b2,b3)
    t1.bac<-margs.higher.order.full(t1-t1step,t2,eta,b1,b2,b3)
    t2.for<-margs.higher.order.full(t1,t2+t2step,eta,b1,b2,b3)
    t2.bac<-margs.higher.order.full(t1,t2-t2step,eta,b1,b2,b3)
    eta.for<-margs.higher.order.full(t1,t2,eta+etastep,b1,b2,b3)
    eta.bac<-margs.higher.order.full(t1,t2,eta-etastep,b1,b2,b3)
    b1.for<-margs.higher.order.full(t1,t2,eta,b1+b1step,b2,b3)
    b1.bac<-margs.higher.order.full(t1,t2,eta,b1-b1step,b2,b3)
    b2.for<-margs.higher.order.full(t1,t2,eta,b1,b2+b2step,b3)
    b2.bac<-margs.higher.order.full(t1,t2,eta,b1,b2-b2step,b3)
    b3.for<-margs.higher.order.full(t1,t2,eta,b1,b2,b3+b3step)
    b3.bac<-margs.higher.order.full(t1,t2,eta,b1,b2,b3-b3step)
    t1.for.t2.for<-margs.higher.order.full(t1+t1step,t2+t2step,eta,b1,b2,b3)
    t1.for.t2.bac<-margs.higher.order.full(t1+t1step,t2-t2step,eta,b1,b2,b3)
    t1.bac.t2.for<-margs.higher.order.full(t1-t1step,t2+t2step,eta,b1,b2,b3)
    t1.bac.t2.bac<-margs.higher.order.full(t1-t1step,t2-t2step,eta,b1,b2,b3)
    t1.for.eta.for<-margs.higher.order.full(t1+t1step,t2,eta+etastep,b1,b2,b3)
    t1.for.eta.bac<-margs.higher.order.full(t1+t1step,t2,eta-etastep,b1,b2,b3)
    t1.bac.eta.for<-margs.higher.order.full(t1-t1step,t2,eta+etastep,b1,b2,b3)
    t1.bac.eta.bac<-margs.higher.order.full(t1-t1step,t2,eta-etastep,b1,b2,b3)
    t1.for.b1.for<-margs.higher.order.full(t1+t1step,t2,eta,b1+b1step,b2,b3)
    t1.for.b1.bac<-margs.higher.order.full(t1+t1step,t2,eta,b1-b1step,b2,b3)
    t1.bac.b1.for<-margs.higher.order.full(t1-t1step,t2,eta,b1+b1step,b2,b3)
    t1.bac.b1.bac<-margs.higher.order.full(t1-t1step,t2,eta,b1-b1step,b2,b3)
    t1.for.b2.for<-margs.higher.order.full(t1+t1step,t2,eta,b1,b2+b2step,b3)
    t1.for.b2.bac<-margs.higher.order.full(t1+t1step,t2,eta,b1,b2-b2step,b3)
    t1.bac.b2.for<-margs.higher.order.full(t1-t1step,t2,eta,b1,b2+b2step,b3)
    t1.bac.b2.bac<-margs.higher.order.full(t1-t1step,t2,eta,b1,b2-b2step,b3)
    t1.for.b3.for<-margs.higher.order.full(t1+t1step,t2,eta,b1,b2,b3+b3step)
    t1.for.b3.bac<-margs.higher.order.full(t1+t1step,t2,eta,b1,b2,b3-b3step)
    t1.bac.b3.for<-margs.higher.order.full(t1-t1step,t2,eta,b1,b2,b3+b3step)
    t1.bac.b3.bac<-margs.higher.order.full(t1-t1step,t2,eta,b1,b2,b3-b3step)
    t2.for.eta.for<-margs.higher.order.full(t1,t2+t2step,eta+etastep,b1,b2,b3)
    t2.for.eta.bac<-margs.higher.order.full(t1,t2+t2step,eta-etastep,b1,b2,b3)
    t2.bac.eta.for<-margs.higher.order.full(t1,t2-t2step,eta+etastep,b1,b2,b3)
    t2.bac.eta.bac<-margs.higher.order.full(t1,t2-t2step,eta-etastep,b1,b2,b3)
    t2.for.b1.for<-margs.higher.order.full(t1,t2+t2step,eta,b1+b1step,b2,b3)
    t2.for.b1.bac<-margs.higher.order.full(t1,t2+t2step,eta,b1-b1step,b2,b3)
    t2.bac.b1.for<-margs.higher.order.full(t1,t2-t2step,eta,b1+b1step,b2,b3)
    t2.bac.b1.bac<-margs.higher.order.full(t1,t2-t2step,eta,b1-b1step,b2,b3)
    t2.for.b2.for<-margs.higher.order.full(t1,t2+t2step,eta,b1,b2+b2step,b3)
    t2.for.b2.bac<-margs.higher.order.full(t1,t2+t2step,eta,b1,b2-b2step,b3)
    t2.bac.b2.for<-margs.higher.order.full(t1,t2-t2step,eta,b1,b2+b2step,b3)
    t2.bac.b2.bac<-margs.higher.order.full(t1,t2-t2step,eta,b1,b2-b2step,b3)
    t2.for.b3.for<-margs.higher.order.full(t1,t2+t2step,eta,b1,b2,b3+b3step)
    t2.for.b3.bac<-margs.higher.order.full(t1,t2+t2step,eta,b1,b2,b3-b3step)
    t2.bac.b3.for<-margs.higher.order.full(t1,t2-t2step,eta,b1,b2,b3+b3step)
    t2.bac.b3.bac<-margs.higher.order.full(t1,t2-t2step,eta,b1,b2,b3-b3step)
    eta.for.b1.for<-margs.higher.order.full(t1,t2,eta+etastep,b1+b1step,b2,b3)
    eta.for.b1.bac<-margs.higher.order.full(t1,t2,eta+etastep,b1-b1step,b2,b3)
    eta.bac.b1.for<-margs.higher.order.full(t1,t2,eta-etastep,b1+b1step,b2,b3)
    eta.bac.b1.bac<-margs.higher.order.full(t1,t2,eta-etastep,b1-b1step,b2,b3)
    eta.for.b2.for<-margs.higher.order.full(t1,t2,eta+etastep,b1,b2+b2step,b3)
    eta.for.b2.bac<-margs.higher.order.full(t1,t2,eta+etastep,b1,b2-b2step,b3)
    eta.bac.b2.for<-margs.higher.order.full(t1,t2,eta-etastep,b1,b2+b2step,b3)
    eta.bac.b2.bac<-margs.higher.order.full(t1,t2,eta-etastep,b1,b2-b2step,b3)
    eta.for.b3.for<-margs.higher.order.full(t1,t2,eta+etastep,b1,b2,b3+b3step)
    eta.for.b3.bac<-margs.higher.order.full(t1,t2,eta+etastep,b1,b2,b3-b3step)
    eta.bac.b3.for<-margs.higher.order.full(t1,t2,eta-etastep,b1,b2,b3+b3step)
    eta.bac.b3.bac<-margs.higher.order.full(t1,t2,eta-etastep,b1,b2,b3-b3step)
    b1.for.b2.for<-margs.higher.order.full(t1,t2,eta,b1+b1step,b2+b2step,b3)
    b1.for.b2.bac<-margs.higher.order.full(t1,t2,eta,b1+b1step,b2-b2step,b3)
    b1.bac.b2.for<-margs.higher.order.full(t1,t2,eta,b1-b1step,b2+b2step,b3)
    b1.bac.b2.bac<-margs.higher.order.full(t1,t2,eta,b1-b1step,b2-b2step,b3)
    b1.for.b3.for<-margs.higher.order.full(t1,t2,eta,b1+b1step,b2,b3+b3step)
    b1.for.b3.bac<-margs.higher.order.full(t1,t2,eta,b1+b1step,b2,b3-b3step)
    b1.bac.b3.for<-margs.higher.order.full(t1,t2,eta,b1-b1step,b2,b3+b3step)
    b1.bac.b3.bac<-margs.higher.order.full(t1,t2,eta,b1-b1step,b2,b3-b3step)
    b2.for.b3.for<-margs.higher.order.full(t1,t2,eta,b1,b2+b2step,b3+b3step)
    b2.for.b3.bac<-margs.higher.order.full(t1,t2,eta,b1,b2+b2step,b3-b3step)
    b2.bac.b3.for<-margs.higher.order.full(t1,t2,eta,b1,b2-b2step,b3+b3step)
    b2.bac.b3.bac<-margs.higher.order.full(t1,t2,eta,b1,b2-b2step,b3-b3step)
    H<-Matrix(0,nrow=6,ncol=6,sparse=FALSE)
    H[1,1]<-(t1.for$log.lik-2*cen+t1.bac$log.lik)/(t1step^2)
    H[2,2]<-(t2.for$log.lik-2*cen+t2.bac$log.lik)/(t2step^2)
    H[3,3]<-(eta.for$log.lik-2*cen+eta.bac$log.lik)/(etastep^2)
    H[4,4]<-(b1.for$log.lik-2*cen+b1.bac$log.lik)/(b1step^2)
    H[5,5]<-(b2.for$log.lik-2*cen+b2.bac$log.lik)/(b2step^2)
    H[6,6]<-(b3.for$log.lik-2*cen+b3.bac$log.lik)/(b3step^2)
    H[1,2]<-H[2,1]<-(t1.for.t2.for$log.lik-t1.for.t2.bac$log.lik-t1.bac.t2.for$log.lik+t1.bac.t2.bac$log.lik)/(4*t2step*t1step)
    H[1,3]<-H[3,1]<-(t1.for.eta.for$log.lik-t1.for.eta.bac$log.lik-t1.bac.eta.for$log.lik+t1.bac.eta.bac$log.lik)/(4*etastep*t1step)
    H[1,4]<-H[4,1]<-(t1.for.b1.for$log.lik-t1.for.b1.bac$log.lik-t1.bac.b1.for$log.lik+t1.bac.b1.bac$log.lik)/(4*b1step*t1step)
    H[1,5]<-H[5,1]<-(t1.for.b2.for$log.lik-t1.for.b2.bac$log.lik-t1.bac.b2.for$log.lik+t1.bac.b2.bac$log.lik)/(4*b2step*t1step)
    H[1,6]<-H[6,1]<-(t1.for.b3.for$log.lik-t1.for.b3.bac$log.lik-t1.bac.b3.for$log.lik+t1.bac.b3.bac$log.lik)/(4*b3step*t1step)
    H[2,3]<-H[3,2]<-(t2.for.eta.for$log.lik-t2.for.eta.bac$log.lik-t2.bac.eta.for$log.lik+t2.bac.eta.bac$log.lik)/(4*t2step*etastep)
    H[2,4]<-H[4,2]<-(t2.for.b1.for$log.lik-t2.for.b1.bac$log.lik-t2.bac.b1.for$log.lik+t2.bac.b1.bac$log.lik)/(4*t2step*b1step)
    H[2,5]<-H[5,2]<-(t2.for.b2.for$log.lik-t2.for.b2.bac$log.lik-t2.bac.b2.for$log.lik+t2.bac.b2.bac$log.lik)/(4*t2step*b2step)
    H[2,6]<-H[6,2]<-(t2.for.b3.for$log.lik-t2.for.b3.bac$log.lik-t2.bac.b3.for$log.lik+t2.bac.b3.bac$log.lik)/(4*t2step*b3step)
    H[3,4]<-H[4,3]<-(eta.for.b1.for$log.lik-eta.for.b1.bac$log.lik-eta.bac.b1.for$log.lik+eta.bac.b1.bac$log.lik)/(4*etastep*b1step)
    H[3,5]<-H[5,3]<-(eta.for.b2.for$log.lik-eta.for.b2.bac$log.lik-eta.bac.b2.for$log.lik+eta.bac.b2.bac$log.lik)/(4*etastep*b2step)
    H[3,6]<-H[6,3]<-(eta.for.b3.for$log.lik-eta.for.b3.bac$log.lik-eta.bac.b3.for$log.lik+eta.bac.b3.bac$log.lik)/(4*etastep*b3step)
    H[4,5]<-H[5,4]<-(b1.for.b2.for$log.lik-b1.for.b2.bac$log.lik-b1.bac.b2.for$log.lik+b1.bac.b2.bac$log.lik)/(4*b2step*b1step)
    H[4,6]<-H[6,4]<-(b1.for.b3.for$log.lik-b1.for.b3.bac$log.lik-b1.bac.b3.for$log.lik+b1.bac.b3.bac$log.lik)/(4*b3step*b1step)
    H[5,6]<-H[6,5]<-(b2.for.b3.for$log.lik-b2.for.b3.bac$log.lik-b2.bac.b3.for$log.lik+b2.bac.b3.bac$log.lik)/(4*b2step*b3step)
    min1<-min(as.numeric(t1.for$log.lik),as.numeric(t1.bac$log.lik),as.numeric(t2.for$log.lik),as.numeric(t2.bac$log.lik),as.numeric(eta.for$log.lik),as.numeric(eta.bac$log.lik),
              as.numeric(b1.for$log.lik),as.numeric(b1.bac$log.lik),as.numeric(b2.for$log.lik),as.numeric(b2.bac$log.lik),
              as.numeric(b3.for$log.lik),as.numeric(b3.bac$log.lik))
    min2<-min(as.numeric(t1.for.t2.for$log.lik),as.numeric(t1.for.t2.bac$log.lik),as.numeric(t1.for.eta.for$log.lik),as.numeric(t1.for.eta.bac$log.lik),
              as.numeric(t1.bac.t2.for$log.lik),as.numeric(t1.bac.t2.bac$log.lik),
              as.numeric(t1.bac.eta.for$log.lik),as.numeric(t1.bac.eta.bac$log.lik),
              as.numeric(t2.for.eta.for$log.lik),as.numeric(t2.for.eta.bac$log.lik),
              as.numeric(t2.bac.eta.for$log.lik),as.numeric(eta.for.b1.for$log.lik),as.numeric(eta.for.b1.bac$log.lik),as.numeric(eta.bac.b1.for$log.lik),
              as.numeric(eta.bac.b1.bac$log.lik),as.numeric(b1.for.b2.for$log.lik),as.numeric(b1.for.b2.bac$log.lik),
              as.numeric(b1.bac.b2.for$log.lik),as.numeric(b1.bac.b2.bac$log.lik),as.numeric(b1.for.b3.for$log.lik),
              as.numeric(b1.for.b3.bac$log.lik),as.numeric(b1.bac.b3.for$log.lik),as.numeric(b1.bac.b3.bac$log.lik),
              as.numeric(b2.for.b3.for$log.lik),as.numeric(b2.for.b3.bac$log.lik),as.numeric(b2.bac.b3.for$log.lik),
              as.numeric(b2.bac.b3.bac$log.lik))
    prevstep<-t1step
    t2step<-t2step/2
    b1step<-b2step<-b3step<-t1step<-etastep<-prevstep/2
    min<-min(min1,min2) #If Step is takes outside of parameter space, this will decrease step and push us back in
  }
    output<-list(t1.p=(t1.for$log.lik-t1.bac$log.lik)/(2*t1step),t2.p=(t2.for$log.lik-t2.bac$log.lik)/(2*t2step),eta.p=(eta.for$log.lik-eta.bac$log.lik)/(2*etastep),b1.p=(b1.for$log.lik-b1.bac$log.lik)/(2*b1step),b2.p=(b2.for$log.lik-b2.bac$log.lik)/(2*b2step),b3.p=(b3.for$log.lik-b3.bac$log.lik)/(2*b3step),H=H)
    return(output)
}


#Newtom-Raphson

old<-matrix(c(.3, .181, .50,-5, .13, .50),nrow=6)

diff<-1
m<-1
q<-1
while(q <10 & (abs(diff)>.005)){
  if(diff>.005){
    cen<-margs.higher.order.full(old[1,1],old[2,1],old[3,1],old[4,1],old[5,1],old[6,1])$log.lik
    Hess.est<-Hess.full(old[1,1],old[2,1],old[3,1],old[4,1],old[5,1],old[6,1],.001,.001,.001,.001,.001,.001,cen=cen)
    H<-Hess.est$H
    grad<-Matrix(c(as.numeric(Hess.est$t1.p),as.numeric(Hess.est$t2.p),as.numeric(Hess.est$eta.p),as.numeric(Hess.est$b1.p),as.numeric(Hess.est$b2.p),as.numeric(Hess.est$b3.p)),nrow=6)
  }
  x.new<-old-1/m*solve(H)%*%grad
  if(x.new[1,1]<0){x.new[1,1]=.1}
  if(x.new[2,1]<0){x.new[2,1]=.01}
  if(x.new[2,1]>.182){x.new[2,1]=.182}
  if(x.new[3,1]>1){x.new[3,1]=.9}
  if(x.new[3,1]<0){x.new[3,1]=.01}
  new.lik<-margs.higher.order.full(x.new[1],x.new[2],x.new[3],x.new[4],x.new[5],x.new[6])
  diff<-as.numeric(new.lik$log.lik-cen)
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


x.new-1.96*sqrt(diag(solve(-H)))
x.new+1.96*sqrt(diag(solve(-H)))


#Estimated Covariance at posterior mode
solve(-H)


te=margs.higher.order.full(x.new[1,1],x.new[2,1],x.new[3,1],x.new[4,1],x.new[5,1],x.new[6,1])

hist(as.numeric(te$gii))  #Majority of g^{ii} terms are well behaved
abline(v=1,col="red",lwd=3)

te$log.lik #Log-likelihood

te$r4s #Addition of sixth order partials only minimally impacts log-likelihood

#Inference using Standard Laplace Approximation

#Function to find Hessian for Standard Laplace Approximation

Hess.LA1<-function(t1,t2,eta,b1,b2,b3,t1step,t2step,etastep,b1step,b2step,b3step,cen=cen){
  min<- -300000
  while(min < -200000){
    t1.for<-margs.chi(t1+t1step,t2,eta,b1,b2,b3)
    t1.bac<-margs.chi(t1-t1step,t2,eta,b1,b2,b3)
    t2.for<-margs.chi(t1,t2+t2step,eta,b1,b2,b3)
    t2.bac<-margs.chi(t1,t2-t2step,eta,b1,b2,b3)
    eta.for<-margs.chi(t1,t2,eta+etastep,b1,b2,b3)
    eta.bac<-margs.chi(t1,t2,eta-etastep,b1,b2,b3)
    b1.for<-margs.chi(t1,t2,eta,b1+b1step,b2,b3)
    b1.bac<-margs.chi(t1,t2,eta,b1-b1step,b2,b3)
    b2.for<-margs.chi(t1,t2,eta,b1,b2+b2step,b3)
    b2.bac<-margs.chi(t1,t2,eta,b1,b2-b2step,b3)
    b3.for<-margs.chi(t1,t2,eta,b1,b2,b3+b3step)
    b3.bac<-margs.chi(t1,t2,eta,b1,b2,b3-b3step)
    t1.for.t2.for<-margs.chi(t1+t1step,t2+t2step,eta,b1,b2,b3)
    t1.for.t2.bac<-margs.chi(t1+t1step,t2-t2step,eta,b1,b2,b3)
    t1.bac.t2.for<-margs.chi(t1-t1step,t2+t2step,eta,b1,b2,b3)
    t1.bac.t2.bac<-margs.chi(t1-t1step,t2-t2step,eta,b1,b2,b3)
    t1.for.eta.for<-margs.chi(t1+t1step,t2,eta+etastep,b1,b2,b3)
    t1.for.eta.bac<-margs.chi(t1+t1step,t2,eta-etastep,b1,b2,b3)
    t1.bac.eta.for<-margs.chi(t1-t1step,t2,eta+etastep,b1,b2,b3)
    t1.bac.eta.bac<-margs.chi(t1-t1step,t2,eta-etastep,b1,b2,b3)
    t1.for.b1.for<-margs.chi(t1+t1step,t2,eta,b1+b1step,b2,b3)
    t1.for.b1.bac<-margs.chi(t1+t1step,t2,eta,b1-b1step,b2,b3)
    t1.bac.b1.for<-margs.chi(t1-t1step,t2,eta,b1+b1step,b2,b3)
    t1.bac.b1.bac<-margs.chi(t1-t1step,t2,eta,b1-b1step,b2,b3)
    t1.for.b2.for<-margs.chi(t1+t1step,t2,eta,b1,b2+b2step,b3)
    t1.for.b2.bac<-margs.chi(t1+t1step,t2,eta,b1,b2-b2step,b3)
    t1.bac.b2.for<-margs.chi(t1-t1step,t2,eta,b1,b2+b2step,b3)
    t1.bac.b2.bac<-margs.chi(t1-t1step,t2,eta,b1,b2-b2step,b3)
    t1.for.b3.for<-margs.chi(t1+t1step,t2,eta,b1,b2,b3+b3step)
    t1.for.b3.bac<-margs.chi(t1+t1step,t2,eta,b1,b2,b3-b3step)
    t1.bac.b3.for<-margs.chi(t1-t1step,t2,eta,b1,b2,b3+b3step)
    t1.bac.b3.bac<-margs.chi(t1-t1step,t2,eta,b1,b2,b3-b3step)
    t2.for.eta.for<-margs.chi(t1,t2+t2step,eta+etastep,b1,b2,b3)
    t2.for.eta.bac<-margs.chi(t1,t2+t2step,eta-etastep,b1,b2,b3)
    t2.bac.eta.for<-margs.chi(t1,t2-t2step,eta+etastep,b1,b2,b3)
    t2.bac.eta.bac<-margs.chi(t1,t2-t2step,eta-etastep,b1,b2,b3)
    t2.for.b1.for<-margs.chi(t1,t2+t2step,eta,b1+b1step,b2,b3)
    t2.for.b1.bac<-margs.chi(t1,t2+t2step,eta,b1-b1step,b2,b3)
    t2.bac.b1.for<-margs.chi(t1,t2-t2step,eta,b1+b1step,b2,b3)
    t2.bac.b1.bac<-margs.chi(t1,t2-t2step,eta,b1-b1step,b2,b3)
    t2.for.b2.for<-margs.chi(t1,t2+t2step,eta,b1,b2+b2step,b3)
    t2.for.b2.bac<-margs.chi(t1,t2+t2step,eta,b1,b2-b2step,b3)
    t2.bac.b2.for<-margs.chi(t1,t2-t2step,eta,b1,b2+b2step,b3)
    t2.bac.b2.bac<-margs.chi(t1,t2-t2step,eta,b1,b2-b2step,b3)
    t2.for.b3.for<-margs.chi(t1,t2+t2step,eta,b1,b2,b3+b3step)
    t2.for.b3.bac<-margs.chi(t1,t2+t2step,eta,b1,b2,b3-b3step)
    t2.bac.b3.for<-margs.chi(t1,t2-t2step,eta,b1,b2,b3+b3step)
    t2.bac.b3.bac<-margs.chi(t1,t2-t2step,eta,b1,b2,b3-b3step)
    eta.for.b1.for<-margs.chi(t1,t2,eta+etastep,b1+b1step,b2,b3)
    eta.for.b1.bac<-margs.chi(t1,t2,eta+etastep,b1-b1step,b2,b3)
    eta.bac.b1.for<-margs.chi(t1,t2,eta-etastep,b1+b1step,b2,b3)
    eta.bac.b1.bac<-margs.chi(t1,t2,eta-etastep,b1-b1step,b2,b3)
    eta.for.b2.for<-margs.chi(t1,t2,eta+etastep,b1,b2+b2step,b3)
    eta.for.b2.bac<-margs.chi(t1,t2,eta+etastep,b1,b2-b2step,b3)
    eta.bac.b2.for<-margs.chi(t1,t2,eta-etastep,b1,b2+b2step,b3)
    eta.bac.b2.bac<-margs.chi(t1,t2,eta-etastep,b1,b2-b2step,b3)
    eta.for.b3.for<-margs.chi(t1,t2,eta+etastep,b1,b2,b3+b3step)
    eta.for.b3.bac<-margs.chi(t1,t2,eta+etastep,b1,b2,b3-b3step)
    eta.bac.b3.for<-margs.chi(t1,t2,eta-etastep,b1,b2,b3+b3step)
    eta.bac.b3.bac<-margs.chi(t1,t2,eta-etastep,b1,b2,b3-b3step)
    b1.for.b2.for<-margs.chi(t1,t2,eta,b1+b1step,b2+b2step,b3)
    b1.for.b2.bac<-margs.chi(t1,t2,eta,b1+b1step,b2-b2step,b3)
    b1.bac.b2.for<-margs.chi(t1,t2,eta,b1-b1step,b2+b2step,b3)
    b1.bac.b2.bac<-margs.chi(t1,t2,eta,b1-b1step,b2-b2step,b3)
    b1.for.b3.for<-margs.chi(t1,t2,eta,b1+b1step,b2,b3+b3step)
    b1.for.b3.bac<-margs.chi(t1,t2,eta,b1+b1step,b2,b3-b3step)
    b1.bac.b3.for<-margs.chi(t1,t2,eta,b1-b1step,b2,b3+b3step)
    b1.bac.b3.bac<-margs.chi(t1,t2,eta,b1-b1step,b2,b3-b3step)
    b2.for.b3.for<-margs.chi(t1,t2,eta,b1,b2+b2step,b3+b3step)
    b2.for.b3.bac<-margs.chi(t1,t2,eta,b1,b2+b2step,b3-b3step)
    b2.bac.b3.for<-margs.chi(t1,t2,eta,b1,b2-b2step,b3+b3step)
    b2.bac.b3.bac<-margs.chi(t1,t2,eta,b1,b2-b2step,b3-b3step)
    
    
    H<-Matrix(0,nrow=6,ncol=6,sparse=FALSE)
    H[1,1]<-(t1.for$log.lik-2*cen+t1.bac$log.lik)/(t1step^2)
    H[2,2]<-(t2.for$log.lik-2*cen+t2.bac$log.lik)/(t2step^2)
    H[3,3]<-(eta.for$log.lik-2*cen+eta.bac$log.lik)/(etastep^2)
    H[4,4]<-(b1.for$log.lik-2*cen+b1.bac$log.lik)/(b1step^2)
    H[5,5]<-(b2.for$log.lik-2*cen+b2.bac$log.lik)/(b2step^2)
    H[6,6]<-(b3.for$log.lik-2*cen+b3.bac$log.lik)/(b3step^2)
    
    
    H[1,2]<-H[2,1]<-(t1.for.t2.for$log.lik-t1.for.t2.bac$log.lik-t1.bac.t2.for$log.lik+t1.bac.t2.bac$log.lik)/(4*t2step*t1step)
    H[1,3]<-H[3,1]<-(t1.for.eta.for$log.lik-t1.for.eta.bac$log.lik-t1.bac.eta.for$log.lik+t1.bac.eta.bac$log.lik)/(4*etastep*t1step)
    H[1,4]<-H[4,1]<-(t1.for.b1.for$log.lik-t1.for.b1.bac$log.lik-t1.bac.b1.for$log.lik+t1.bac.b1.bac$log.lik)/(4*b1step*t1step)
    H[1,5]<-H[5,1]<-(t1.for.b2.for$log.lik-t1.for.b2.bac$log.lik-t1.bac.b2.for$log.lik+t1.bac.b2.bac$log.lik)/(4*b2step*t1step)
    H[1,6]<-H[6,1]<-(t1.for.b3.for$log.lik-t1.for.b3.bac$log.lik-t1.bac.b3.for$log.lik+t1.bac.b3.bac$log.lik)/(4*b3step*t1step)
    H[2,3]<-H[3,2]<-(t2.for.eta.for$log.lik-t2.for.eta.bac$log.lik-t2.bac.eta.for$log.lik+t2.bac.eta.bac$log.lik)/(4*t2step*etastep)
    H[2,4]<-H[4,2]<-(t2.for.b1.for$log.lik-t2.for.b1.bac$log.lik-t2.bac.b1.for$log.lik+t2.bac.b1.bac$log.lik)/(4*t2step*b1step)
    H[2,5]<-H[5,2]<-(t2.for.b2.for$log.lik-t2.for.b2.bac$log.lik-t2.bac.b2.for$log.lik+t2.bac.b2.bac$log.lik)/(4*t2step*b2step)
    H[2,6]<-H[6,2]<-(t2.for.b3.for$log.lik-t2.for.b3.bac$log.lik-t2.bac.b3.for$log.lik+t2.bac.b3.bac$log.lik)/(4*t2step*b3step)
    H[3,4]<-H[4,3]<-(eta.for.b1.for$log.lik-eta.for.b1.bac$log.lik-eta.bac.b1.for$log.lik+eta.bac.b1.bac$log.lik)/(4*etastep*b1step)
    H[3,5]<-H[5,3]<-(eta.for.b2.for$log.lik-eta.for.b2.bac$log.lik-eta.bac.b2.for$log.lik+eta.bac.b2.bac$log.lik)/(4*etastep*b2step)
    H[3,6]<-H[6,3]<-(eta.for.b3.for$log.lik-eta.for.b3.bac$log.lik-eta.bac.b3.for$log.lik+eta.bac.b3.bac$log.lik)/(4*etastep*b3step)
    H[4,5]<-H[5,4]<-(b1.for.b2.for$log.lik-b1.for.b2.bac$log.lik-b1.bac.b2.for$log.lik+b1.bac.b2.bac$log.lik)/(4*b2step*b1step)
    H[4,6]<-H[6,4]<-(b1.for.b3.for$log.lik-b1.for.b3.bac$log.lik-b1.bac.b3.for$log.lik+b1.bac.b3.bac$log.lik)/(4*b3step*b1step)
    H[5,6]<-H[6,5]<-(b2.for.b3.for$log.lik-b2.for.b3.bac$log.lik-b2.bac.b3.for$log.lik+b2.bac.b3.bac$log.lik)/(4*b2step*b3step)
    
    
    min1<-min(as.numeric(t1.for$log.lik),as.numeric(t1.bac$log.lik),as.numeric(t2.for$log.lik),as.numeric(t2.bac$log.lik),as.numeric(eta.for$log.lik),as.numeric(eta.bac$log.lik),
              as.numeric(b1.for$log.lik),as.numeric(b1.bac$log.lik),as.numeric(b2.for$log.lik),as.numeric(b2.bac$log.lik),
              as.numeric(b3.for$log.lik),as.numeric(b3.bac$log.lik))
    min2<-min(as.numeric(t1.for.t2.for$log.lik),as.numeric(t1.for.t2.bac$log.lik),as.numeric(t1.for.eta.for$log.lik),as.numeric(t1.for.eta.bac$log.lik),
              as.numeric(t1.bac.t2.for$log.lik),as.numeric(t1.bac.t2.bac$log.lik),
              as.numeric(t1.bac.eta.for$log.lik),as.numeric(t1.bac.eta.bac$log.lik),
              as.numeric(t2.for.eta.for$log.lik),as.numeric(t2.for.eta.bac$log.lik),
              as.numeric(t2.bac.eta.for$log.lik),as.numeric(eta.for.b1.for$log.lik),as.numeric(eta.for.b1.bac$log.lik),as.numeric(eta.bac.b1.for$log.lik),
              as.numeric(eta.bac.b1.bac$log.lik),as.numeric(b1.for.b2.for$log.lik),as.numeric(b1.for.b2.bac$log.lik),
              as.numeric(b1.bac.b2.for$log.lik),as.numeric(b1.bac.b2.bac$log.lik),as.numeric(b1.for.b3.for$log.lik),
              as.numeric(b1.for.b3.bac$log.lik),as.numeric(b1.bac.b3.for$log.lik),as.numeric(b1.bac.b3.bac$log.lik),
              as.numeric(b2.for.b3.for$log.lik),as.numeric(b2.for.b3.bac$log.lik),as.numeric(b2.bac.b3.for$log.lik),
              as.numeric(b2.bac.b3.bac$log.lik))
    prevstep<-t1step
    t2step<-t2step/2
    b1step<-b2step<-b3step<-t1step<-etastep<-prevstep/2
    min<-min(min1,min2)
  }
  output<-list(t1.p=(t1.for$log.lik-t1.bac$log.lik)/(2*t1step),t2.p=(t2.for$log.lik-t2.bac$log.lik)/(2*t2step),eta.p=(eta.for$log.lik-eta.bac$log.lik)/(2*etastep),b1.p=(b1.for$log.lik-b1.bac$log.lik)/(2*b1step),b2.p=(b2.for$log.lik-b2.bac$log.lik)/(2*b2step),b3.p=(b3.for$log.lik-b3.bac$log.lik)/(2*b3step),H=H)
  return(output)
  }






old<-matrix(c(.5,.18,.5,-6,.15,.5),nrow=6)

diff<-1
m<-1
q<-1
while(q <10 & (abs(diff)>.005)){
  if(diff>.005){
    cen<-margs.chi(old[1,1],old[2,1],old[3,1],old[4,1],old[5,1],old[6,1])$log.lik
    Hess.est<-Hess.LA1(old[1,1],old[2,1],old[3,1],old[4,1],old[5,1],old[6,1],.001,.001,.001,.001,.001,.001,cen=cen)
    H<-Hess.est$H
    grad<-Matrix(c(as.numeric(Hess.est$t1.p),as.numeric(Hess.est$t2.p),as.numeric(Hess.est$eta.p),as.numeric(Hess.est$b1.p),as.numeric(Hess.est$b2.p),as.numeric(Hess.est$b3.p)),nrow=6)
  }
  x.new<-old-1/m*solve(H)%*%grad
  if(x.new[1,1]<0){x.new[1,1]=.1}
  if(x.new[2,1]<0){x.new[2,1]=.01}
  if(x.new[2,1]>.182){x.new[2,1]=.182}
  if(x.new[3,1]>1){x.new[3,1]=.9}
  if(x.new[3,1]<0){x.new[3,1]=.01}
  new.lik<-margs.chi(x.new[1],x.new[2],x.new[3],x.new[4],x.new[5],x.new[6])
  diff<-as.numeric(new.lik$log.lik-cen)
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

x.new+1.96*sqrt(diag(-solve(H)))
x.new-1.96*sqrt(diag(-solve(H)))

#Minimal differences in other parameters
