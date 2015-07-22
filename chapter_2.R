#==============================================================================#
# Chapter 2 Genetic Covariance Between Relatives
#==============================================================================#

#----------------------------------------#
# libraries
#----------------------------------------#

# Import library
	library(pedigreemm)

#----------------------------------------#
# Example pedigree (3rd ed, page 23)
#----------------------------------------#

# pedigree in table 2.1
	ped  <- data.frame(calf = as.character(c(3, 4, 5, 6)),
						sire = c(1, 1, 4, 5), 
						dam = c(2, NA, 3, 2))
	print(ped)

#----------------------------------------#
# Use pedigreemm package
#----------------------------------------#

# editPed() to add parents to top of pedigree
	ped.edit <- editPed(sire=ped$sire, dam=ped$dam, label=ped$calf)
	print(ped.edit)

# pedigree() function to create pedigree S4 object
	ped.complete <- pedigree(sire= ped.edit $sire, 
							dam= ped.edit $dam, 
							label= ped.edit $label)
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
	Dmat(ped.complete)

# L' matrix (3rd ed, page 29)
	relfactor(ped.complete)

























