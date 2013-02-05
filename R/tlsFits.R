#####################################################
##
## tlsFits.R -- Total least squares fits for several functional forms, including:
## 	-- Linear; done using Demming regression (tlsDemming), and svd (tlsSvd).
##  -- LOESS; hacked approximation done by transforming data.  Rotats data so that the LOESS fit approximates a 
##				given angle (by default 45 degrees; assumes varience of two data types is equal). 
##
##
## WARNING -- FUNCTIONS CURRENTLY EXPERIMENTAL!! USE WITH EXTREME CATION!!
##
## WARNING -- CURRENTLY DEMMING AND SVD METHODS FOR FITTING TOTAL LEAST SQUARES DO NOT AGREE!!  SEE BELOW FOR MORE DISCUSSION.
##
#####################################################


########################################################### 
## Transform and fit LOESS!
## 
## By transforming X and Y by 45 degrees, we can fit using LOESS and it's essentially "erros in variables" LOESS.
##
## Returns x and y values giving information on the fit...
##
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


## WARNING -- CURRENTLY DEMMING AND SVD METHODS FOR FITTING TOTAL LEAST SQUARES DO NOT AGREE!!
## This warning seems to be related to the normalization applied to the data!
##
## Differences in fits due to normalization of the input variables.
#> tlsSvd((expr[use]-mean(expr[use]))/var(expr[use]), (rate[use]-mean(rate[use]))/var(expr[use]))
#[1] 5.159783e-36 7.182321e-06
#> tlsDemming((expr[use]-mean(expr[use]))/var(expr[use]), (rate[use]-mean(rate[use]))/var(expr[use]))
#[1] 5.159783e-36 7.182321e-06
#> pcafit(data.frame((expr[use]-mean(expr[use]))/var(expr[use]), (rate[use]-mean(rate[use]))/var(expr[use])))
#$slope
#[1] 7.182321e-06
#  ...
##
## Conversely, the SVD results can be pretty closely reproduced by assuming that d = (var(rate[use])^2)/((var(expr[use]))^2)

 #####################################################################3
 ## linear total least squares by SVD.
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
## Test: 
## A <- as.matrix(rbind(c(1,1,1), c(2,-1,2), c(-1,4,3), c(4,2,1), c(3,-3,4)))
## b <- c(1,2,-1,4,8)
## 

 #####################################################################3
 ## linear total least squares by Demming regression.
tlsDemming <- function(x,y,d=1) {
 sxx <- 1/(NROW(x)-1) * sum((x-mean(x))^2)
 sxy <- 1/(NROW(x)-1) * sum((x-mean(x))*(y-mean(y)))
 syy <- 1/(NROW(y)-1) * sum((y-mean(y))^2)
 X <- (syy-d*sxx+sqrt((syy-d*sxx)^2+4*d*(sxy^2)))/(2*sxy)
 int <- mean(y) - X * mean(x)
 return(c(int,X))
}



