##### R2D2 #####
# This script runs the R2D2 code used in Van de Velde et al. (2021)

# Packages
install.packages("D:\\Users\\jpvdveld\\Documents\\PhD\\Code\\Packages\\R2D2_2.0.tar.gz", repos = NULL, type="source")
install.packages("R.matlab")

# Pad and libraries
rm(list=ls(all=TRUE))
setwd("E:\\Users\\jpvdveld\\Onderzoek\\Data\\")
library(R.matlab)
library(R2D2)

# Data 

data_ba.tmp <- readMat('1_biascorrection\\MPI-rcp45_threshold_qdm_results.mat')
data_ref.tmp <- readMat('0_original\\Uccle_xho.mat')

data_ba <- data_ba.tmp[1][[1]][[1]]
data_ba[[1]][6616, 4] <- 0
data_ref <- data_ref.tmp$xho

# Bias adjustment

seasons <- rbind(c(12,1,2),c(3,4,5), c(6,7,8), c(9,10,11))

r2d2_correction <- list()
r2d2_correction$E <- matrix(0, length(data_ba[[1]][,1]),3)
r2d2_correction$T <- matrix(0, length(data_ba[[1]][,1]),3)
r2d2_correction$P <- matrix(0, length(data_ba[[1]][,1]),3)

for(season in 1:4){
  
  r2d2_correction.E <- r2d2(data_ref[which(data_ref[,2] %in% seasons[season,]),4:6], 
                            data_ba[[1]][which(data_ba[[1]][,2] %in% seasons[season,]),4:6], 
                            icond = c(1), lag_search = 5, lag_keep = 3)
  r2d2_correction.T <- r2d2(data_ref[which(data_ref[,2] %in% seasons[season,]),4:6], 
                            data_ba[[1]][which(data_ba[[1]][,2] %in% seasons[season,]),4:6], 
                            icond = c(2), lag_search = 5, lag_keep = 3)
  r2d2_correction.P <- r2d2(data_ref[which(data_ref[,2] %in% seasons[season,]),4:6], 
                            data_ba[[1]][which(data_ba[[1]][,2] %in% seasons[season,]),4:6], 
                            icond = c(3), lag_search = 5, lag_keep = 3)
  
  r2d2_correction$E[which(data_ba[[1]][,2] %in% seasons[season,]),1:3] <- r2d2_correction.E$r2d2_bc
  r2d2_correction$T[which(data_ba[[1]][,2] %in% seasons[season,]),1:3] <- r2d2_correction.T$r2d2_bc
  r2d2_correction$P[which(data_ba[[1]][,2] %in% seasons[season,]),1:3] <- r2d2_correction.P$r2d2_bc
}


# Saving data in a Matlab-friendly way

writeMat('1_biascorrection\\MPI-rcp45_threshold_r2d2_results.mat', out = r2d2_correction)