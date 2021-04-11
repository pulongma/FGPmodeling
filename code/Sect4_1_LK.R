############################################################################
# This file contains an example code to implemente LK in Section 4.1 when
# the underlying true covariance model is the exponential covariance
############################################################################

rm(list=ls())

source("CRPS.R")

##### simulate GP with Deterministic function


library(geoR)
library(fields)
nx = 50
ny = 50
s1 = 0
s2 = 50
L0 = s2 - s1
xg = seq(s1, s2, length.out=nx)
yg = seq(s1, s2, length.out=ny)
locs = as.matrix(expand.grid(xg, yg))

n = dim(locs)[1]


distmat = rdist(locs, locs)
sig2=16
sig2eps = 1.6

### cov model 1
rho = Matern.cor.to.range(20, 0.5, cor.target = 0.05)
nu=0.5


mat_corr = Matern(distmat, range=rho, nu=nu)
mat_cov = sig2*mat_corr
L = t(chol(mat_cov))
#cov.EK = Matern(out$d, range=rho, nu=nu)



set.seed(123)
Y = c(L%*%rnorm(n))
#n_test = floor(0.1*n)
#ind_MAR = sample(seq(1, n), n_test)
ind_MAR = c(read.table("MAR_ind_sim2Dcode.txt")$V1)
n_test = length(ind_MAR)
locs_pred = locs[ind_MAR, ,drop=F]
Y_test = Y[ind_MAR]
Yobs = Y[-ind_MAR]
locs_obs = locs[-ind_MAR, ]



#################################################################
##### Lattice Kriging

library(LatticeKrig)


set.seed(357)

nsim = 30


LK.MSPE = rep(NA, nsim)
LK.MSPE = rep(NA, nsim)


tStart.LK <- proc.time()

for(i in 1:nsim){

Z = Y + sqrt(sig2eps)*rnorm(n)
Zobs = Z[-ind_MAR]
Z_test = Z[ind_MAR]

#LKinfo = LKrigSetup( x=as.matrix(locs), LKGeometry="LKInterval", nlevel=3,
#                    fixedFunction = NULL, 
#                    mean.neighbor=50, BasisFunction="WendlandFunction", 
#                    NC=50, verbose=FALSE)

# obj = LatticeKrig(as.matrix(locs_obs), Zobs, nlevel=3)
obj = LatticeKrig(as.matrix(locs_obs), Zobs, nlevel=3, findAwght=TRUE)


outputSim = LKrig.sim.conditional(obj, x.grid=as.matrix(locs_pred), M=30)

LK.se <- sqrt(apply(outputSim$g.draw, 1, "var"))
LK.fhat <- outputSim$ghat

resZ = LK.fhat - Z_test

LK.MSPE[i] = mean(resZ^2)


LK.CRPS[i] = mean(CRPS(Z_test, LK.fhat, LK.se))

}

tEnd.LK <- proc.time() - tStart.LK


# save(tEnd.LK, LK.MSPE, LK.MSPE, file="LK_nu05_EF20.RData")



