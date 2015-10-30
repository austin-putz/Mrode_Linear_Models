#==============================================================================#
# Chapter 3: Best Linear Unbiased Prediction of Breeding Value:
#            Univariate Models with One Random Effect
#==============================================================================#

#------------------------------------------------------------------------------#
# Description
#------------------------------------------------------------------------------#

# Raphael A. Mrode
# Linear Models for the Prediction of Animal Breeding Values

# Author:   Austin Putz <putz[dot]austin[at]gmail[dot]com>
# Created:  Unknown
# Modified: 2015-07-31
# License:  GPLv2

#------------------------------------------------------------------------------#
# Libraries and functions
#------------------------------------------------------------------------------#

# install pedigreemm if not yet installed
  if (!("pedigreemm" %in% as.character(as.matrix(installed.packages())[, 1]))){
	  install.packages("pedigreemm")
  }

# pedigreemm
  library(pedigreemm)

#------------------------------------------------------------------------------#
# Basic model
#------------------------------------------------------------------------------#

# y = Xb + Za + e

# y = (nx1) vector of observations, n = # of records
# b = (px1) vector of fixed effects, p = # of levels of fixed effects
# a = (qx1) vector of random animal effects, q = # of levels of random effects (animals in pedigree)
# e = (nx1) vector of random residual effects
# X = (nxp) design matrix relating records to fixed effects
# Z = (nxq) design matrix relating records to random effects (animals)

# var(a) = G
# var(e) = R
# var(y) = V = ZGZ' + R

# MME (assuming homogeneous variance)
# Simplifies to:

# | X'X         X'Z       | b | = | X'y |
# | Z'X   Z'Z + A-1*alpha | a | = | Z'y |

#------------------------------------------------------------------------------#
# Example 3.1
#------------------------------------------------------------------------------#

# This example includes 5 records on 5 separate individual calves. 
# Pre-weaning gain was measured as the response. Fixed effects included
# sex with 3 males and 2 females. There are 3 ancestors in this data set.

#----------------------------------------#
# Data
#----------------------------------------#

# given parameters (page 37)
  sigmaE <- 40; sigmaE
  sigmaA <- 20; sigmaA
  alpha  <- sigmaE / sigmaA; alpha

# Data (page 37)
  calf   <- c(4:8)
  sex    <- c("male", "female", "female", "male", "male")
  WWG    <- c(4.5, 2.9, 3.9, 3.5, 5.0)

# Dataframe
  data.3.1 <- data.frame(calf, sex, WWG)
  print(data.3.1)
  rm(list=(c("calf", "sex", "WWG")))

#----------------------------------------#
# Pedigree
#----------------------------------------#

# set up pedigree
  calf    <- c(1:8)
  sire    <- c(NA, NA, NA, 1, 3, 1, 4, 3)
  dam     <- c(NA, NA, NA, NA, 2, 2, 5, 6)
  ped.3.1 <- data.frame(calf, sire, dam)
  print(ped.3.1)
  rm(list=c("calf","sire","dam"))
  
# editPed() to add parents to top of pedigree
  ped.edit <- editPed(sire=ped.3.1$sire, 
                      dam=ped.3.1$dam, 
                      label=ped.3.1$calf)
  print(ped.edit)

# pedigree() function to create pedigree S4 object
  ped.complete <- pedigree(sire= ped.edit$sire, 
                           dam= ped.edit$dam, 
                           label= ped.edit$label)
  print(ped.complete)

# create A matrix (3rd ed, page 23) (uses the matrix package, thus the "."s)
  A <- getA(ped.complete)
  print(A)
  
# A inverse
  round(solve(A), 3)

#----------------------------------------#
# Matrices
#----------------------------------------#

# set up X
  X <- with(data.3.1, model.matrix(WWG ~ sex - 1))
  print(X)

# NOTE: X'X  can be  X %*% t(X)  OR  crossprod(X)
# crossprod function is slightly faster

# view X'X
  crossprod(X)

# set up Z (n by n*) n* = number of animals in whole pedigree (8)
  Z <- matrix(0, nrow=nrow(data.3.1), ncol=nrow(ped.3.1))

# fill in Z
  for (i in 1:nrow(data.3.1)) {
    index = data.3.1[i, "calf" ]
    Z[i, index] = 1
  }
  print(Z)

# view Z'Z
  crossprod(Z)

# view X'Z
  crossprod(X, Z)

# view X'Y
  crossprod(X, data.3.1$WWG)

# view Z'Y
  crossprod(Z, data.3.1$WWG)

#----------------------------------------#
# Function to solve MME
#----------------------------------------#

# set up MME for PE solutions
  basicMME <- function(X, Z, A, y, alpha) {
    
	# Calculate blocks of LHS
    XpX = crossprod(X)
    ZpZ = crossprod(Z) + (solve(A) * alpha)
    XpZ = crossprod(X, Z)
    
	# Paste top and bottom together for LHS
    toprow    = cbind(XpX,    XpZ)
    bottomrow = cbind(t(XpZ), ZpZ)
    
	# Put top and bottom together for left hand side (LHS)
    LHS = rbind(toprow, bottomrow)
    
	# Elements for RHS
    XpY = crossprod(X, y)
    ZpY = crossprod(Z, y)
    
	# Calculate right hand side (RHS)
    RHS = rbind(XpY, ZpY) 
    
	# calculate solutions by direct inversion of LHS
    solutions = solve(LHS) %*% RHS
    
	# Return LHS, RHS, and solutions
    return(list(LHS=LHS, RHS=RHS, solutions=solutions))
    
  }

#----------------------------------------#
# LHS and RHS
#----------------------------------------#

# solve MME (LHS, RHS, and solutions)
  output.3.1 <- basicMME(X=X, Z=Z, A=A, y=data.3.1$WWG, alpha=alpha)

# get LHS
  output.3.1$LHS

# get inverse of LHS
  LHS.inv <- solve(output.3.1$LHS)
  round(LHS.inv, 3)

# get RHS
  output.3.1$RHS

#----------------------------------------#
# Solutions
#----------------------------------------#

# get solutions (page 39)
  sol.names          <- as.matrix(c("Female", "Male", "An_1_BV", "An_2_BV", "An_3_BV", 
                                    "An_4_BV", "An_5_BV", "An_6_BV", "An_7_BV", "An_8_BV" ), ncol=1)
  sols               <- as.matrix(round(output.3.1$solutions, 3))
  row.names(sols)    <- NULL
  solutions          <- as.data.frame(cbind(sol.names, sols))
  names(solutions)   <- c("Effect", "Solution")
  solutions$Effect   <- as.character(solutions$Effect)
  solutions$Solution <- as.numeric(as.character(solutions$Solution))
  print(solutions)

#----------------------------------------#
# Yield Deviations
#----------------------------------------#

# get Yield Deviations (YD) (page 40); YD = RESIDUAL = Y - Y_hat because they only have 1 observation
  YDs <- solve(crossprod(Z[1:5, 4:8])) %*% t(Z[1:5, 4:8]) %*% as.matrix(data.3.1$WWG - (X %*% matrix(solutions[1:2, 2])))
  print(YDs)

#----------------------------------------#
# Accuracy (3rd ed, page 44)
#----------------------------------------#

# PEV = C22.inv * sigmaE
# SEP = sqrt(PEV)
# rel = 1 - (C22.inv * alpha)
# acc = sqrt(rel)

# elements to make table on page 45
  animal    <- 1:8
  diagonals <- diag(LHS.inv[3:10, 3:10])  # don't take fixed diagonals
  r2        <- 1 - (diagonals * alpha)
  r         <- sqrt(r2)
  SEP       <- sqrt(diagonals * sigmaE)
  
# create table for r2, r, SEP
  accuracy.tab <- data.frame(animal, diagonals, r2, r, SEP)
  rm(list=c("animal", "diagonals", "r2", "r", "SEP"))
  round(accuracy.tab, 3)

#------------------------------------------------------------------------------#
# 3.2: Sire Model
#------------------------------------------------------------------------------#

# remove everything for next example
  rm(list=ls())

#####################################
# Caution: Work in progess
#####################################

# source functions
  source("/Users/austinputz/Documents/Animal_Breeding/Mrode_Examples/functions.R")
    # 1) createA function for A matrix creation

# read in data
  sex       <- c("male", "female", "female", "male", "male")
  sire      <- c(1, 3, 1, 4, 3)
  sire.sire <- c(0, 0, 0, 1, 0)
  dam.sire  <- rep(0, 5)
  WWG       <- c(4.5, 2.9, 3.9, 3.5, 5.0)

# create data frame
  data.3.2 <- data.frame(sex, sire, sire.sire, dam.sire, WWG)

# pedigree
  sire   <- c(1, 3, 4)
  sire.s <- c(0, 0, 1)
  sire.d <- c(0, 0, 0)
  ped    <- data.frame(sire, sire.s, sire.d)

# set up A matrix
  A     <- createA(s=ped$sire.s, d=ped$sire.d)
  A.inv <- solve(A)

# set up alpha
  sigma2.p <- 40 + 20              # 60
  sigma2.s <- 0.25 * 20            # 5
  sigma2.e <- sigma2.p - sigma2.s  # 55
  alpha    <- sigma2.e / sigma2.s  # 11

# set up X and Z matrices
  X             <- model.matrix(lm(WWG ~ sex -1, data=data.3.2))
  X             <- X[,c("sexmale","sexfemale")]
  data.3.2$sire <- as.factor(data.3.2$sire)
  Z             <- model.matrix(lm(WWG ~ sire -1, data=data.3.2))

# set up mixed model equations
  XpX <- crossprod(X)
  XpZ <- crossprod(X, Z)
  ZpX <- t(XpZ)
  ZpZ <- crossprod(Z)
  XpY <- crossprod(X, data.3.2$WWG)
  ZpY <- crossprod(Z, data.3.2$WWG)

# create large matrices
  LHS.top <- cbind(XpX, XpZ)
  LHS.bot <- cbind(ZpX, (ZpZ + (A.inv * alpha)) )  # remember to add A^-1 * alpha
  LHS     <- rbind(LHS.top, LHS.bot)
  RHS     <- rbind(XpY, ZpY)

# solutions
  solutions <- solve(LHS) %*% RHS
  solutions

#------------------------------------------------------------------------------#
# 3.6: Animal Model with Groups
#------------------------------------------------------------------------------#

#####################################
# Caution: Work in progess
#####################################

# Pedigree
  calf    <- 1:8
  sire    <- c(0,0,0,1,3,1,4,3)
  dam     <- c(0,0,0,0,2,2,5,6)
  ped.3.6 <- data.frame(calf, sire, dam)
  
# Remove individual vectors
  rm(list=c("calf","sire","dam"))

# Assume sire and dam are of different genetic merit!
  ped.3.6[1:3, 2]  <- 9        # Sires = G1 (group 1)
  ped.3.6[1:4, 3]  <- 10        # Dams  = G2 (group 2)

# # Recode as [G1 = 9] and [G2 = 10] (animals end at 8)
#   ped.3.6[ped.3.6 == "G1"] = 9
#   ped.3.6[ped.3.6 == "G2"] = 10
  
# Data
  calves  <- c(4:8)
  sex     <- c("male", "female", "female", "male", "male")
  WWG     <- c(4.5, 2.9, 3.9, 3.5, 5.0)
  data.3.6 <- data.frame(calves, sex, WWG)

# Remove individual vectors
  rm(list=c("calves", "sex", "WWG"))

# parameters
  n.obs    <- nrow(data.3.6)        # 5 observations
  n.an     <- nrow(ped.3.6)         # 8 animals
  n.groups <- max(ped.3.6) - nrow(ped.3.6)   # 2 groups

# set up X, Z, and Y matrices
  X <- model.matrix(lm(WWG ~ sex -1, data=data.3.6))
  Z <- matrix(0, nrow=(max(apply(ped.3.6, 2, max))), 
                 ncol=(max(apply(ped.3.6, 2, max))))
  Z.sub = diag(nrow(data.3.6))


























