###########################################################################
##
##   Copyright 2009, 2010, 2011 Charles Danko.
##
##   This program is part of the GRO-seq R package
##
##   GRO-seq is free software: you can redistribute it and/or modify it 
##   under the terms of the GNU General Public License as published by 
##   the Free Software Foundation, either version 3 of the License, or  
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful, but 
##   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
##   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
##   for more details.
##
##   You should have received a copy of the GNU General Public License along 
##   with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##########################################################################


########################################################################
##
##	Distribution Fitting Functions -- MLEfit given X.
##	Date: 2009-11-22
##
##	Presently, fits:
##		-- Gamma distribution  (RgammaMLE)
##		-- Normal distrubution (Rnorm)
##
##	Todo:
##	(1) Add fits for: 
##		-- Normal by MLE
##		-- Poisson
##		-- Negative binomial (??)
##
########################################################################
RgammaMLE <- function(X) {

	N <- as.real(NROW(X))
	sumxis <- as.real(sum(X))
	sumlogxis <- as.real(sum(log(X)))
	Fit <- .Call("RgammaMLE", N, sumxis, sumlogxis, PACKAGE = "groHMM")
	return(Fit)

}

Rnorm <- function(X) {

	returnList      <- list()
	returnList$mean <- mean(X)
	returnList$var  <- var(X)

	return(returnList)

}

################################
##
##  R interface to MLEFit for alpha*Normal+(1-alpha)*Exponential hybrid distribution.
##
################################
Rnorm.exp <- function(xi, wi=rep(1,NROW(xi)), guess=c(0.5, 0, 1, 1), tol=sqrt(.Machine$double.eps), maxit=10000) {
  	Fit <- .Call("RNormExpMLE", xi, wi, guess, tol, as.integer(maxit), PACKAGE = "groHMM")
}

################################
##
## weighted.var -- Computes the weighted varience, where the varience is weighted by some factor...
##
###############################

weighted.var <- function(x, w, na.rm = FALSE) {
    if (na.rm) {
        w <- w[i <- !is.na(x)]
        x <- x[i]
    }
    sum.w <- sum(w)
    sum.w2 <- sum(w^2)
    mean.w <- sum(x * w) / sum(w)
    (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm = na.rm)
}


