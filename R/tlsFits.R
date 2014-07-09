###########################################################################
##
##   Copyright 2009, 2010, 2011 Charles Danko and Minho Chae.
##
##   This program is part of the groHMM R package
##
##   groHMM is free software: you can redistribute it and/or modify it 
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


#####################################################
##
## tlsFits.R -- Total least squares fits for several functional forms, 
## including:
##  -- Linear; done using Demming regression (tlsDemming), and svd (tlsSvd).
##  -- LOESS; hacked approximation done by transforming data.  Rotats data so 
##      that the LOESS fit approximates a given angle (by default 45 degrees; 
##      assumes varience of two data types is equal). 
##
##
## WARNING -- FUNCTIONS CURRENTLY EXPERIMENTAL!! USE WITH EXTREME CATION!!
##
## WARNING -- CURRENTLY DEMMING AND SVD METHODS FOR FITTING TOTAL LEAST 
##          SQUARES DO NOT AGREE!!  SEE BELOW FOR MORE DISCUSSION.
##
#####################################################


#' A 'total least squares'-like hack for LOESS.  Works by rotating points 
#' 45 degrees, fitting LOESS, and rotating back.
#'
#' @param x X values.
#' @param y Y values.
#' @param theta Amount to rotate, sets the ratio of variences that are assumed 
#' by the hack.  Default: -pi/4 radians (45 degrees) for orthogonal regression.
#' @param span The LOESS span parameter.  Default: 1
#' @return List of input values and LOESS predictions.
#' @author Charles G. Danko
## Transform and fit LOESS!
## 
## By transforming X and Y by 45 degrees, we can fit using LOESS and it's 
## essentially "erros in variables" LOESS.
##
## Returns x and y values giving information on the fit...
tlsLoess <- function(x, y, theta= -pi/4, span= 1) {
 mx <- mean(x)
 my <- mean(y)
 x <- x-mx
 y <- y-my
 tx <- x*cos(theta)-y*sin(theta)
 ty <- x*sin(theta)+y*cos(theta)
 
 lo <- loess(ty~tx, family="gaussian", span = span)
 yp_loess  <- predict(lo,data.frame(x=tx))

 theta <- -theta
 x <- tx*cos(theta)-yp_loess*sin(theta)
 y <- tx*sin(theta)+yp_loess*cos(theta)
 x <- x+mx
 y <- y+my
 
 ord <- order(x)
 x <- x[ord]
 y <- y[ord]
 
 retVar <- list()
 retVar$x <- x
 retVar$y <- y
 
 return(retVar)
}

#' A 'total least squares' implementation using singular value demposition.
#'
#' @param x X values.
#' @param y Y values.
#' @return Parameters for the linear model Y~a*X+e.
#' @author Charles G. Danko
## linear total least squares by SVD.
## Test: 
## A <- as.matrix(rbind(c(1,1,1), c(2,-1,2), c(-1,4,3), c(4,2,1), c(3,-3,4)))
## b <- c(1,2,-1,4,8)
## 
tlsSvd <- function(x,y) {
 n <- NCOL(x)
 Cmat <- as.matrix(data.frame(x, y))
 s <- svd(Cmat)
 VAB <- s$v[c(1:n),(1+n):NCOL(s$v)]
 VBB <- s$v[(1+n):NROW(s$v),(1+n):NCOL(s$v)]
 X = -VAB/VBB
 int <- mean(y) - X * mean(x)
 return(c(int,X))
}

#' A 'total least squares' implementation using demming regression.
#'
#' @param x X values.
#' @param y Y values.
#' @param d Ratio of variences. Default: 1, for orthogonal regression.
#' @return Parameters for the linear model.
#' @author Charles G. Danko
## linear total least squares by Deming regression.
tlsDeming <- function(x,y,d=1) {
 sxx <- 1/(NROW(x)-1) * sum((x-mean(x))^2)
 sxy <- 1/(NROW(x)-1) * sum((x-mean(x))*(y-mean(y)))
 syy <- 1/(NROW(y)-1) * sum((y-mean(y))^2)
 X <- (syy-d*sxx+sqrt((syy-d*sxx)^2+4*d*(sxy^2)))/(2*sxy)
 int <- mean(y) - X * mean(x)
 return(c(int,X))
}
