
## Return coefficient matrix (lambda) of n-th order Legendre polynomials 

## Arguments
## 	n      : order of polynomials
##	gengler: logical value. If TRUE, Genlger's scaling (1999) will be 
## 				applied. If not specified, TRUE is assumed.

## Literatures
## Mrode, R.A. 2005. Linear Models for the Prediction of Animal Breeding Values. 
## 		CAB International, Oxon, UK.
## Gengler, N. et. al. 1999. Estimation of (Co)variance Function Coefficients for 
## 		Test Day Yield with a Expectation-Maximization Restricted Maximum Likelihood 
##  		Algorithm.  Journal of Dairy Science.  82

## Author: Gota Morota <morota at wisc dot edu>
## Create: 31-Mar-2010
## Last-Modified: 2-Apr-2010
## License: GPLv3 or later

`legendre` <- function(n, gengler=FALSE){
	
	if (!(gengler %in% c(TRUE, FALSE))){
		warning("option for gengler not set correctly: will be set to FALSE")
		gengler=FALSE
	}
	
	# set order of matrix (order=n + 1)
	N <- n+1
	
	# set up matrix to fill with loop
	L <- matrix(0, nrow=N, ncol=N)
	
	# begin loop to fill the matrix
	for(i in (1:N)){
		
		if(i==1){
	 		L[i,i] <- 1
		}
		else if(i==2){
			L[i,i] <- 1
		}
		else  { # this will be i > 2
			tmp <- L[i-1, ] # this will print the line before it
		    tmp2 <- as.numeric()  # this will initiate a numeric vector
			tmp2 <- c(0, tmp[1:(N-1)])
			L[i,] <- (1/(i-2+1))*((2*(i-2) + 1)*tmp2 - (i-2)*L[i-2,])
		}
	}
	
	# Normalize
	for (j in (1:N)){	
		L[j, ] <- (sqrt( (2*(j-1) + 1)/2)  )*L[j, ]
	}
	
	# Gengler (1999)
	if (gengler==TRUE){
		L <- sqrt(2)*L
	}
	
	# return the matrix of coefficients
	return(t(L))
	
}







