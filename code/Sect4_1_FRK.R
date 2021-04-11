############################################################################
# This file contains an example code to implemente FRK with 3 levels of 
# basis functions in Section 4.1 when the underlying true covariance model
# is the exponential covariance
############################################################################

rm(list=ls())

source("CRPS.R")


##### simulate GP with Deterministic function

library(geoR)
library(fields)
nx = 50
ny = 50
s1 = 0
s2 = 100
L0 = s2 - s1
xg = seq(s1, s2, length.out=nx)
yg = seq(s1, s2, length.out=ny)
locs = expand.grid(xg, yg)

n = dim(locs)[1]
distmat = as.matrix(rdist(locs))

sig2 = 16
sig2eps = 4
phi=3
mat_cov = sig2*exp(-distmat/phi)

L = t(chol(mat_cov))

Y = c(L%*%rnorm(n))

### hold out observations

## hold out a region [s1+0.25*L0,s1+0.75*L0]x[s1+0.25*L0,s1+0.75*L0]
flag_region = rep(FALSE, n)
locs_region = locs[flag_region==TRUE, ]
locs_remain = locs[flag_region==FALSE, ]
n_remain = dim(locs_remain)[1]
Y_remain = Y[flag_region==FALSE]
label = seq(1,n)
ind_MBD = label[flag_region==TRUE]
ind_remain = label[flag_region==FALSE]

nobs = floor(0.9*n_remain)  

set.seed(230)
indobs = sample(seq(1,n_remain), nobs)

Yobs = Y_remain[indobs]
locs_obs = locs_remain[indobs, ]

ind_MAR = ind_remain[-indobs]

locs_pred = locs[c(ind_MBD, ind_MAR), ]
Y_test = Y[c(ind_MBD, ind_MAR)]
test_MBD = seq(1, length(ind_MBD))
test_MAR = seq(length(ind_MBD)+1, length(Y_test))







#################################################################
## Load packages
library(FRK)
library(sp)
library(ggplot2)
library(gridExtra)
library(INLA)
library(splancs)






set.seed(123)
Ztemp = Y + sqrt(sig2eps)*rnorm(n)
df = data.frame(Lon=locs[ ,1], Lat=locs[ ,2], obs=Y)

## Make BAUs as SpatialPixels
BAUs <- df                            # assign BAUs                 # remove data from BAUs
BAUs$fs <- 1                          # set fs variation to unity
coordinates(BAUs)  <- ~Lon+Lat        # convert to SpatialPointsDataFrame
gridded(BAUs) <- TRUE                 # convert to SpatialPixelsDataFrame

## Make Data as SpatialPoints
dat <- subset(df,!is.na(obs))        # no missing data in data frame
coordinates(dat)  <- ~Lon+Lat         # Convert to SpatialPointsDataFrame

#### create basis functions at 2 resolutions 
basis <- auto_basis(plane(),          # we are on the plane
                    data = dat,       # data around which to make basis
                    #regular = 0,      # regular basis with 16 basis at coarsest resolution
                    nres = 3,         # 3 resolutions
                    scale_aperture = 1,   # aperture scaling of basis functions 
                    type="bisquare")

## save basis functions 
basis.df = as.data.frame(basis)
dim(basis.df[basis.df$res==1, ])
dim(basis.df[basis.df$res==2, ])
dim(basis.df[basis.df$res==3, ])

# 24+137+700

basis_info = as.matrix(basis.df)
# write.table(basis_info,
#            file='basis_info_3Res_2500grids.txt',
#            sep=" ", row.names=FALSE, col.names=FALSE)

# write.table(basis_info, 
#             file='basis_info_3Res_2500grids.txt', 
#             sep=" ", row.names=FALSE, col.names=FALSE)


########  begin FRK 
set.seed(123)

nsim = 30

tStart.FRK3r = vector(mode="list", nsim)
tEnd.FRK3r = vector(mode="list", nsim)
FRK.3rMSPE = rep(NA, nsim)
FRK.3rCRPS = rep(NA, nsim)



for(i in 1:nsim){

Z = Y + sqrt(sig2eps)*rnorm(n)
Z_remain = Z[flag_region==FALSE]
Zobs = Z_remain[indobs]
Z_test = Z[c(ind_MBD, ind_MAR)]

tStart.FRK3r[[i]] <- proc.time()

dat.obs = data.frame(Lon=locs_obs[ ,1], Lat=locs_obs[ ,2], Z=Zobs, std=sqrt(sig2eps))
coordinates(dat.obs)  <- ~Lon+Lat         # Convert to SpatialPointsDataFrame

## Estimate using ML
S <- FRK(f = Z~1,                       # formula for SRE model
         data = dat.obs,               # data
         basis = basis,               # Basis
         BAUs = BAUs,                 # BAUs
         est_error=0,
         K_type="unstructured",
         tol = 1e-3, print_lik=FALSE)                   # EM iterations


## Predict
BAUs_pred <- SRE.predict(S)           # predict over all BAUs
BAUs_pred_df <- data.frame(BAUs_pred) # convert to data frame

FRK.yhat = BAUs_pred_df$mu[c(ind_MBD, ind_MAR)]
FRK.3rMSPE[i] = mean((FRK.yhat-Y_test)^2)

BAUs_pred_df$sd_obs <- BAUs_pred_df$sd
FRK.yse = BAUs_pred_df$sd_obs[c(ind_MBD, ind_MAR)]
FRK.3rCRPS[i] = mean(CRPS(Y_test, FRK.yhat, FRK.yse))

tEnd.FRK3r[[i]] <- proc.time() - tStart.FRK3r[[i]]

}

# save(tEnd.FRK3r, FRK.3rMSPE, FRK.3rCRPS, file="FRK3r_results.RData")

