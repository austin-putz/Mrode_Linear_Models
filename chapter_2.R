#==============================================================================#
# Chapter 2: Genetic Covariance Between Relatives
#==============================================================================#

#------------------------------------------------------------------------------#
# Description
#------------------------------------------------------------------------------#

# Raphael A. Mrode
# Linear Models for the Prediction of Animal Breeding Values

# Author:   Austin Putz <putz[dot]austin[at]gmail[dot]com>
# Created:  Unknown
# Modified: 2015-07-22
# License:  GPLv2

#------------------------------------------------------------------------------#
# install packages
#------------------------------------------------------------------------------#

# install pedigreemm if not yet installed
  if (!("pedigreemm" %in% as.character(as.matrix(installed.packages())[, 1]))) install.packages("pedigreemm")

#------------------------------------------------------------------------------#
# libraries
#------------------------------------------------------------------------------#

# Import library
	library(pedigreemm)

#------------------------------------------------------------------------------#
# Example pedigree (3rd ed, page 23)
#------------------------------------------------------------------------------#

# pedigree in table 2.1
	ped  <- data.frame(calf = as.character(c(3, 4, 5, 6)),
						sire = c(1, 1, 4, 5), 
						dam = c(2, NA, 3, 2))
	print(ped)

#------------------------------------------------------------------------------#
# Use pedigreemm package
#------------------------------------------------------------------------------#

# editPed() to add parents to top of pedigree
	ped.edit <- editPed(sire=ped$sire, dam=ped$dam, label=ped$calf)
	print(ped.edit)

# pedigree() function to create pedigree S4 object
	ped.complete <- pedigree(sire= ped.edit$sire, 
							dam= ped.edit$dam, 
							label= ped.edit$label)
	print(ped.complete)

# create A matrix (3rd ed, page 23)
# uses the matrix package, thus the "."s
	A     <- getA(ped.complete)
	print(A)
	
# create A inverse (3rd ed, page 27)
	A.inv <- getAInv(ped.complete)
	print(A.inv)

# get inbreeding coefficients (= 1 - diagonal)
	inbred.coefs <- inbreeding(ped.complete)
	print(inbred.coefs)

# D matrix for the A = TDT' equation
	diag(Dmat(ped.complete))

# L' matrix (3rd ed, page 29)
	relfactor(ped.complete)

























