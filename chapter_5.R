#==============================================================================#
# Chapter 5: Best Linear Unbiased Prediction of Breeding Value:
#               Multivariate Animal Models
#==============================================================================#

# Author: Austin Putz (putz (dot) austin (at) gmail (dot) com)
# Created: April 13, 2017
# Modified: April 13, 2017
# License: GPLv2

#==============================================================================#
# Description
#==============================================================================#

# Reference: Mrode 2014 (3rd Edition)

# Models
#   y1 = X1b1 + Z1a1 + e1
#   y2 = X2b2 + Z2a2 + e2

# Solve as
# y1 = | X1   0  | b1  + | Z1    0 | a1 + e1
# y2 = | 0   X2  | b2  + |  0   Z2 | a2 + e2

# This is ONLY for EQUAL design matrices, no missing traits for either

#==============================================================================#
# Libraries
#==============================================================================#

  rm(list=ls())

# load libraries
  library(pedigree)
  library(pedigreemm)

# source functions
# 1) createA function for A matrix creation
  source(paste0("/Users/austinputz/Documents/Programming/R/", 
          "Animal_Breeding/Gota_Morota/Pedigrees/createA.R"))

#==============================================================================#
# Libraries
#==============================================================================#
  
# data from Table 5.1
  data.5.1 <- data.frame(calf = 4:8,
  											 sex  = c("male","female","female","male","male"),
  											 WWG  = c(4.5, 2.9, 3.9, 3.5, 5.0),
  											 PWG  = c(6.8, 5.0, 6.8, 6.0, 7.5))
  print(data.5.1)
  
# set n1 and n2
  n  = nrow(data.5.1)
  n1 = length(data.5.1$WWG)
  n2 = length(data.5.1$PWG)
  
# pedigree from Table 5.1
  ped.5.1 <- data.frame(animal = 1:8,
  												sire   = c(0, 0, 0, 1, 3, 1, 4, 3),
  												dam    = c(0, 0, 0, 0, 2, 2, 5, 6))
  print(ped.5.1)
  
# set number in pedigree
  N <- nrow(ped.5.1)
  
#==============================================================================#
# Parameters
#==============================================================================#
  
# 1 = WWG (weaning weight gain)
# 2 = PWG (pre-weaning gain)
  
# Genetic and Residual covariance matrix
  G = matrix(c(20, 18, 18, 40), nrow=2)
  R = matrix(c(40, 11, 11, 30), nrow=2)
  
# get inverses of matrix
  Ginv = solve(G)
  Rinv = solve(R)
  
#==============================================================================#
# Set up matrices
#==============================================================================#
  
#----------------------------------------#
# Create relationship matrices
#----------------------------------------#
  
# get A
  A    <- createA(ped.5.1)
  Ainv <- solve(A)
  
#----------------------------------------#
# Subset Data
#----------------------------------------#
  
# subset to not missing WWG
  data.5.1.WWG <- data.5.1[!is.na(data.5.1$WWG), ]
    
# subset to those not missing PWG
  data.5.1.PWG <- data.5.1[!is.na(data.5.1$PWG), ]
  
#----------------------------------------#
# Response
#----------------------------------------#
  
# set Y's
  Y1 <- matrix(data.5.1.WWG$WWG, ncol=1)
  Y2 <- matrix(data.5.1.PWG$PWG, ncol=1)
  
#----------------------------------------#
# Fixed
#----------------------------------------#
  
# X1 and Z2
  X1 <- model.matrix(~ -1 + sex, data=data.5.1)
  X2 <- model.matrix(~ -1 + sex, data=data.5.1)
  
#----------------------------------------#
# Random
#----------------------------------------#
  
# Z1 and Z2
  Z1 <- matrix(0, nrow=n1, ncol=N)
  Z2 <- matrix(0, nrow=n2, ncol=N)

# fill in Z1
  for (i in 1:nrow(data.5.1.WWG)) {
    # get the calf number for 1:n1 in the data
    index = data.5.1.WWG[i, "calf"]
    Z1[i, index] = 1
  }
  print(Z1)
  
  rownames(Z1) <- data.5.1$calf
  colnames(Z1) <- ped.5.1$animal

# fill in Z2
  for (i in 1:nrow(data.5.1.PWG)) {
    # get the calf number for 1:n1 in the data
    index = data.5.1.PWG[i, "calf"]
    Z2[i, index] = 1
  }
  print(Z2)
  
  rownames(Z2) <- data.5.1$calf
  colnames(Z2) <- ped.5.1$animal
  
#----------------------------------------#
# Residual
#----------------------------------------#
  
# get inverses
  R11inv <- diag(5) * Rinv[1,1]
  R22inv <- diag(5) * Rinv[2,2]
  R12inv <- diag(5) * Rinv[1,2]
  R21inv <- t(R12inv)
  
#==============================================================================#
# Set up MME Solver
#==============================================================================#
  
#----------------------------------------#
# LHS
#----------------------------------------#
  
# first row of LHS MME
  LHS11 <- t(X1) %*% R11inv %*% X1
  LHS12 <- t(X1) %*% R12inv %*% X2
  LHS13 <- t(X1) %*% R11inv %*% Z1
  LHS14 <- t(X1) %*% R12inv %*% Z2

# second row of LHS MME
  LHS21 <- t(X2) %*% R21inv %*% X1
  LHS22 <- t(X2) %*% R22inv %*% X2
  LHS23 <- t(X2) %*% R21inv %*% Z1
  LHS24 <- t(X2) %*% R22inv %*% Z2

# third row of LHS MME
  LHS31 <- (t(Z1) %*% R11inv %*% X1)
  LHS32 <- (t(Z1) %*% R12inv %*% X2)
  LHS33 <- (t(Z1) %*% R11inv %*% Z1) + Ainv*Ginv[1,1]
  LHS34 <- (t(Z1) %*% R12inv %*% Z2) + Ainv*Ginv[1,2]

# third row of LHS MME
  LHS41 <- (t(Z2) %*% R21inv %*% X1)
  LHS42 <- (t(Z2) %*% R22inv %*% X2)
  LHS43 <- (t(Z2) %*% R21inv %*% Z1) + Ainv*Ginv[2,1]
  LHS44 <- (t(Z2) %*% R22inv %*% Z2) + Ainv*Ginv[2,2]

# build LHS MME
  LHS1 <- cbind(LHS11, LHS12, LHS13, LHS14)
  LHS2 <- cbind(LHS21, LHS22, LHS23, LHS24)
  LHS3 <- cbind(LHS31, LHS32, LHS33, LHS34)
  LHS4 <- cbind(LHS41, LHS42, LHS43, LHS44)

# build LHS MME
  LHS <- rbind(LHS1, LHS2, LHS3, LHS4)
  
# show LHS
  round(LHS, 3)
  
#----------------------------------------#
# RHS
#----------------------------------------#
  
# Set up RHS
  RHS11 <- (t(X1) %*% R11inv %*% Y1) + (t(X1) %*% R12inv %*% Y2)
  RHS21 <- (t(X2) %*% R21inv %*% Y1) + (t(X2) %*% R22inv %*% Y2)
  RHS31 <- (t(Z1) %*% R11inv %*% Y1) + (t(Z1) %*% R12inv %*% Y2)
  RHS41 <- (t(Z2) %*% R21inv %*% Y1) + (t(Z2) %*% R22inv %*% Y2)
  
# stack RHS
  RHS <- rbind(RHS11, RHS21, RHS31, RHS41)
  
# print RHS
  round(RHS, 3)
  
#----------------------------------------#
# Solve MME
#----------------------------------------#
  
# get solutions
  solutions <- solve(LHS) %*% RHS
  
# print solutions
  round(solutions, 3)
  
  
  
  
  
  
  
  













