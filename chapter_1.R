#==============================================================================#
# Chapter 1: Genetic Evaluation with Different Sources of Records
#==============================================================================#

#------------------------------------------------------------------------------#
# Description
#------------------------------------------------------------------------------#

# Raphael A. Mrode
# Linear Models for the Prediction of Animal Breeding Values

# Author:   Austin Putz <putz[dot]austin[at]gmail[dot]com>
# Created:  2015-08-02
# Modified: 2015-08-02
# License:  GPLv2

#------------------------------------------------------------------------------#
# 1.2
#------------------------------------------------------------------------------#

# Basic linear model
# 
# y_ij = mu_i + g_i + e_ij
#
# y_ij = record j of the ith animal
# mu_i = non-random fixed environmental effects
# g_i  = sum of additive, dominance, and epistatic genetic values
# e_ij = sum of random environmental effects
# 

#==============================================================================#
# 1.3 Breeding Value Prediction from the animals own performance
#==============================================================================#

#------------------------------------------------------------------------------#
# 1.3.1 Single record
#------------------------------------------------------------------------------#

calcR <- function(i, r2, var_y){
	
	# the r2 is equal to the heritability
	
	# calculate response to selection
	R = i * r2 * var_y
	
	# return variable to global env
	return(R)
	
}

#----------------------------------------#
# Example 1.1
#----------------------------------------#

# yearling weight = 320 kg
# herd mean = 250 kg
# h2 = heritability = 0.45
# 
# calculate her bredding value and accuracy

# Breeding value function for single records (SR)
  calcEBVSR <- function(h2, phen, avg){
	
	# calculate the breeding value
	a = h2 * (phen - avg)
	
	# return the variable to global env
	return(a)
	
  }

# Accuracy function for single records (SR)
  calcAccSR <- function(h2){
	
	# calculate the accuracy
	acc = sqrt(h2)
	
	# return the variable to the global env
	return(acc)
	
  }

# answer to breeding value
  calcEBVSR(h2=0.45, phen=320, avg=250)

# answer to accuracy
  calcAccSR(h2=0.45)




#------------------------------------------------------------------------------#
# 1.3.2 Repated records
#------------------------------------------------------------------------------#

# t     = (var(g) + var(pe)) / var(y)
# 1 - t = var(te) / var(y)

# Calculate t (repeatablity)
  calct <- function(var_g, var_pe, var_y){
	
	# calculate repeatability
	t = (var_g + var_pe) / var_y
	
	# return t
	return(t)
	
  }

# a_i           = b(y_avg_i - mu)
# b             = cov(a, y_avg) / var_y_avg
# cov(a, y_avg) = cov(a, g + pe + Sum te / n) = var_a
# var(y_avg)    = var(g) + var(pe) + var(te) / n
# var(t)        = (t + (1 - t) / n) * var_y
# b             = var_a / (t + (1 - t) / n) * var_y
#               = (n * h2) / (1 + (n - 1) * t)
# b depends on heritability, repeatability, and the number of records

# function for an EBV with repeated records (RR)
  calcEBVRR <- function(h2, t, n, avg_y, mu){
  	
  	# calculate b
  	b = (n * h2) / (1 + (n - 1) * t)
  	
  	# calculate the breeding value
  	a = b * (avg_y - mu)
  	
  	# return the EBV (a)
  	return(a)
  	
  }
  
# My own litter size example
  calcEBVRR(h2=0.1, t=0.15, n=4, avg_y=13, mu=10)

# accuracy of EBV is:
# r_ay    = cov(a, y_avg) / (sd_a * sd_y)
#         = var_a / (var_a * sqrt((t + (1 - t) / n) * var_y))
#         = h * sqrt((n / (1 + (n - 1) * t)))
#         = sqrt((n * h2) / (1 + (n - 1) * t))
#         = sqrt(b)

  calcAccRR <- function(h2, t, n){
  	
  	# calculate the accuracy for a repeated record
  	acc = sqrt((n * h2) / (1 + (n - 1) * t))
  	
  	# return the accuracy
  	return(acc)
  	
  }

# calculate an accuracy for a trait with:
  calcAccRR(h2 = 0.1, t = 0.15, n = 4)
  
# compare with single record accuracy
  calcAccSR(h2 = 0.1)
  
# That's a 66% increase in accuracy!

# Note: This will be more for traits with low repeatability
#       Repeatability generally follows heritability
#       i.e. if h2 is low, t is low
# t cannot be smaller than h2


#==============================================================================#
# 1.4 Breeding Value Prediction from Progeny Records
#==============================================================================#




#==============================================================================#
# 1.5 Breeding Value Prediction from Pedigree
#==============================================================================#




#==============================================================================#
# 1.6 Breeding Value Prediction for One Trait from Another
#==============================================================================#




#==============================================================================#
# 1.7 Selection Index
#==============================================================================#



#------------------------------------------------------------------------------#
# 1.7.1
#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
# 1.7.2
#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
# 1.7.3
#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
# 1.7.4
#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
# 1.7.5
#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
# 1.7.6
#------------------------------------------------------------------------------#










