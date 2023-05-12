# ExtendedLA
R Code for Chapter 4

SimulationandFittingExtendedLA.R is the R code for simulating from the equations 2-4 of the manuscript "Extended Laplace Approximation for Self-Exciting Spatio-Temporal Models of Count Data" and fitting it using the extended Laplace Approximation technique.  The file also contains moment based estimates for initializing a Newton Raphson algorithm.  The numerical approximation to the Hessian used here is based on finite differences.


RCode_for_ExtendedLA_and_LA1.R is the R code to fit the Chicago violent crime data using extended LA and LA(1).
RCode_for_Stan_fit.R is the R code to fit the Chicago violent crime data using Stan relying on precomputing eigenvalues.  The numerical approximation to the Hessian is based on the R function hessian() which uses the Richardson method.

RCode_for_Stan_fit is the R code to fit the MCMC based fit of the Chicago violent crime data. 

