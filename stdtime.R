

## Given time points covariate and order of fit for Legendre polynomials, 
## return matrix 'M' containing the polynomials of standardized time. 'M' is 
## order t (number of time points) by k (order of Legendre polynomials).

## Arguments 
##	t: a vector of time, age or days
## 	n: order of polynomials

## Literature: Mrode, R.A. 2005. Linear Models for the Prediction of Animal 
## Breeding Values. CAB International, Oxon, UK.

## Author: Gota Morota <morota at wisc dot edu>
## Create: 31-Mar-2010
## Last-Modified: 2-Apr-2010
## License: GPLv3 or later

`stdtime` <- function(t, n){
  
	N <- n+1
	M <- matrix(0, nrow=length(t), ncol=N)
	a <- -1 + 2*(t-t[1])/(t[length(t)] - t[1])
	M[, 1] <- 1
	
	for (i in 2:N){
		M[,i] <- a^(i-1)
	}
	
	return(M)
	
}



