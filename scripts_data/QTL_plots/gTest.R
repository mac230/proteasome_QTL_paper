# Log-likelihood tests of independence & goodness of fit
# Does Williams' and Yates' correction
# does Monte Carlo simulation of p-values, via gtestsim.c
#
# G & q calculation from Sokal & Rohlf (1995) Biometry 3rd ed.
# TOI Yates' correction taken from Mike Camann's 2x2 G-test fn.
# GOF Yates' correction as described in Zar (2000)
# more stuff taken from ctest's chisq.test()
#
# ToDo:
# 1) Beautify
# 2) Add warnings for violations
# 3) Make appropriate corrections happen by default
#
# V3.3 Pete Hurd Sept 29 2001. phurd@ualberta.ca

g.test <- function(x, y = NULL, correct="none",
p = rep(1/length(x), length(x)), simulate.p.value = FALSE, B = 2000)
{
	DNAME <- deparse(substitute(x))
	if (is.data.frame(x)) x <- as.matrix(x)
	if (is.matrix(x)) {
		if (min(dim(x)) == 1) 
		x <- as.vector(x)
	}
	if (!is.matrix(x) && !is.null(y)) {
		if (length(x) != length(y)) 
		stop("x and y must have the same length")
		DNAME <- paste(DNAME, "and", deparse(substitute(y)))
		OK <- complete.cases(x, y)
		x <- as.factor(x[OK])
		y <- as.factor(y[OK])
		if ((nlevels(x) < 2) || (nlevels(y) < 2)) 
		stop("x and y must have at least 2 levels")
		x <- table(x, y)
	}
	if (any(x < 0) || any(is.na(x))) 
    stop("all entries of x must be nonnegative and finite")
	if ((n <- sum(x)) == 0) 
    stop("at least one entry of x must be positive")
#If x is matrix, do test of independence
	if (is.matrix(x)) {
#Test of Independence
		nrows<-nrow(x)
		ncols<-ncol(x)
		if (correct=="yates"){ # Do Yates' correction?
			if(dim(x)[1]!=2 || dim(x)[2]!=2) # check for 2x2 matrix
			stop("Yates' correction requires a 2 x 2 matrix")
			if((x[1,1]*x[2,2])-(x[1,2]*x[2,1]) > 0)
			{
			x[1,1] <- x[1,1] - 0.5
			x[2,2] <- x[2,2] - 0.5
			x[1,2] <- x[1,2] + 0.5
			x[2,1] <- x[2,1] + 0.5
			}
			else
			{
			x[1,1] <- x[1,1] + 0.5
			x[2,2] <- x[2,2] + 0.5
			x[1,2] <- x[1,2] - 0.5
			x[2,1] <- x[2,1] - 0.5
			}
			}
			
			sr <- apply(x,1,sum)
			sc <- apply(x,2,sum)
			E <- outer(sr,sc, "*")/n
			# are we doing a monte-carlo?
			# no monte carlo GOF?
			if (simulate.p.value){
			METHOD <- paste("Log likelihood ratio (G-test) test of independence\n\t with simulated p-value based on", B, "replicates")
			tmp <- .C("gtestsim", as.integer(nrows), as.integer(ncols),
			as.integer(sr), as.integer(sc), as.integer(n), as.integer(B),
			as.double(E), integer(nrows * ncols), double(n+1),
			integer(ncols), results=double(B), PACKAGE= "ctest")
			g <- 0
			for (i in 1:nrows){
			for (j in 1:ncols){
			if (x[i,j] != 0) g <- g + x[i,j] * log(x[i,j]/E[i,j])
			}
			}
			STATISTIC <- G <- 2 * g
			PARAMETER <- NA
			PVAL <- sum(tmp$results >= STATISTIC)/B
			}
			else {
			# no monte-carlo
			# calculate G
			g <- 0
			for (i in 1:nrows){
			for (j in 1:ncols){
			if (x[i,j] != 0) g <- g + x[i,j] * log(x[i,j]/E[i,j])
			}
			}
			q <- 1
			if (correct=="williams"){ # Do Williams' correction
			row.tot <- col.tot <- 0    
			for (i in 1:nrows){ row.tot <- row.tot + 1/(sum(x[i,])) }
			for (j in 1:ncols){ col.tot <- col.tot + 1/(sum(x[,j])) }
			q <- 1+ ((n*row.tot-1)*(n*col.tot-1))/(6*n*(ncols-1)*(nrows-1))
			}
			STATISTIC <- G <- 2 * g / q
			PARAMETER <- (nrow(x)-1)*(ncol(x)-1)
			PVAL <- 1-pchisq(STATISTIC,df=PARAMETER)
#			PVAL <- 1-exp(pchisq(STATISTIC,df=PARAMETER, log.p = TRUE))
			if(correct=="none")
			METHOD <- "Log likelihood ratio (G-test) test of independence without correction"
			if(correct=="williams")
			METHOD <- "Log likelihood ratio (G-test) test of independence with Williams' correction"
			if(correct=="yates")
			METHOD <- "Log likelihood ratio (G-test) test of independence with Yates' correction"
			}
			}
			else {
			# x is not a matrix, so we do Goodness of Fit
			METHOD <- "Log likelihood ratio (G-test) goodness of fit test"
			if (length(x) == 1) 
			stop("x must at least have 2 elements")
			if (length(x) != length(p)) 
			stop("x and p must have the same number of elements")
			E <- n * p
			
			if (correct=="yates"){ # Do Yates' correction
			if(length(x)!=2)
			stop("Yates' correction requires 2 data values")
			if ( (x[1]-E[1]) > 0.25) {
			x[1] <- x[1]-0.5
			x[2] <- x[2]+0.5
			}
			else if ( (E[1]-x[1]) > 0.25){
			x[1] <- x[1]+0.5
			x[2] <- x[2]-0.5
			}
			}
			names(E) <- names(x)
			g <- 0
			for (i in 1:length(x)){
			if (x[i] != 0) g <- g + x[i] * log(x[i]/E[i])
			}
			q <- 1
			if (correct=="williams"){ # Do Williams' correction
			q <- 1+(length(x)+1)/(6*n)
		}
		STATISTIC <- G <- 2*g/q
		PARAMETER <- length(x) - 1
#		PVAL <- pchisq(STATISTIC, PARAMETER, lower = FALSE, log.p = TRUE)
		PVAL <- pchisq(STATISTIC, PARAMETER, lower = FALSE)
	}
	names(STATISTIC) <- "Log likelihood ratio statistic (G)"
	names(PARAMETER) <- "X-squared df"
	names(PVAL) <- "p.value"
	structure(list(statistic=STATISTIC,parameter=PARAMETER,p.value=PVAL,
				   method=METHOD,data.name=DNAME, observed=x, expected=E),
			  class="htest")
}
