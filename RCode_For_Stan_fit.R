
library(rstan)

library(spdep)
library(INLA)
library(zoo)
library(LaplacesDemon)
library(maptools)

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

past.effect<-time.Neigh%*%dat
past.effect[1:size,]<-0
Xs<-model.matrix(~1+wx.as.covariate+log.pop)



#Get Neighborhood Information
pairs<-triu(spat.Neigh)
pairs<-summary(pairs)[summary(pairs)$x==1,]


model2="functions {
real sparse_car_lpdf(vector phi, real tau, real zeta, 
int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
row_vector[n] phit_D; 
row_vector[n] phit_W; 
vector[n] ldet_terms;

phit_D = (phi .* D_sparse)';
phit_W = rep_row_vector(0, n);
for (i in 1:W_n) {
phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
}

for (i in 1:n) ldet_terms[i] =log1m(zeta * lambda[i]);
return 0.5 * (n * log(tau)
+ sum(ldet_terms)
- tau * (phit_D * phi - zeta * (phit_W * phi)));
}
}
data {
int<lower = 1> n;
int<lower = 1> p;
int<lower =1> W_n;
matrix[n, p] X;
int<lower = 0> y[n];
int W_sparse[W_n, 2];   // adjacency pairs
vector[n] D_sparse; // number of adjacent region pairs
vector[n] yold;
vector[n] lambda;

}
transformed data {

}
parameters {
vector[p] beta;
vector[n] phi;
real<lower = 0> stds;
real<lower = 0, upper = .185> zeta;
real<lower = 0, upper = 1> eta;
}
transformed parameters{
real<lower = 0> tau;
real<lower = 0> spvar;
tau = 1/stds^2;
spvar=stds^2;
}

model {
phi ~ sparse_car(tau, zeta, W_sparse, D_sparse, lambda, n, W_n);
beta ~ normal(0, 1000);
stds ~ cauchy(0, 5);

y ~ poisson(exp(X * beta + phi)+eta*yold);

}"


lams<-eigen(N.mat)$values


sp_d <- list(n = nrow(Xs),         # number of total observations
             p = ncol(Xs),         # number of coefficients
             X = Xs,               # design matrix
             y = y,                 # observed number of cases
             yold = as.numeric(past.effect),         #AR term
             W_n = sum(spat.Neigh)/2,    # number of neighbor pairs
             W_sparse = pairs[,1:2],    # sparse representation of adjacency matrix
             D_sparse=rep(1,t.size*sp.size),
             lambda=rep(lams,t.size)   #Eigenvalues
)               


m<-stan_model(model_code=model2)


niter=10000
nchains=3
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
sp_fit <- stan(model_code=model2, data = sp_d, 
               iter = niter, chains = nchains,verbose = FALSE,control = list(max_treedepth = 15))

print(sp_fit, pars = c('beta', 'spvar', 'zeta', 'eta','lp__'))

