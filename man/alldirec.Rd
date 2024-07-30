\name{alldirec}
\alias{alldirec}
\alias{print.eqclist}
\title{
	Direct Equating Coefficients Between All Pairs of a List of Forms
}
\description{
	Calculates direct equating coefficients and standard errors using 
	IRT methods between all pairs of a list of forms.
}
\usage{
alldirec(mods, method = "mean-mean", all = FALSE, quadrature = TRUE, nq = 30, 
	direction = "both")
}
\arguments{
	\item{mods}{
	an object of the class \code{modIRT} containing 
	item parameter coefficients and their covariance matrix
	of the forms to be equated.
	}
	\item{method}{
	the equating method to be used. This should be one of
	"\code{mean-mean}", "\code{mean-gmean}", "\code{mean-sigma}",
	"\code{Haebara}" or "\code{Stocking-Lord}".
	}
	\item{all}{
	logical; if FALSE forms that do not have common items
	will not appear in the output.
	}
	\item{quadrature}{logical; if TRUE the Gauss-Hermite quadrature is used
	to approximate the integral in the function that is minimized in the 
	Haebara and Stocking-Lord methods.
	If FALSE the integral is replaced with a sum over 40
	equally spaced values ranging from -4 to 4 with an increment 
	of 0.05 and weights equal to one for all values.
	}
	\item{nq}{number of quadrature points used for the Gauss-Hermite quadrature
	if \code{quadrature} is TRUE
	}
	\item{direction}{a character string specifying the direction used to equate 
	pair of forms. If "\code{back}" newer forms are equated to older forms,
	if "\code{forward}" older forms are equated to newer forms, if "\code{both}"
	(the default) equating is performed in both directions. Form A is considered 
	to be newer than Form B if it is listed after Form B in \code{mods}.
	}
}
\value{
	An object of class \code{eqclist} consisting in a list with length 
	equal to the number of pairs of forms equated. Each component of
	the list is an object of class \code{eqc} returned by function
	\code{\link{direc}}.
}
\references{
	Battauz, M. (2015). equateIRT: An R Package for IRT Test Equating.
	\emph{Journal of Statistical Software}, \bold{68}, 1--22.
	
	Kolen, M.J. and Brennan, R.L. (2014). \emph{Test equating, scaling, 
	and linking: methods and practices}, 3rd ed., New York: Springer

	Ogasawara, H. (2000). Asymptotic standard errors of IRT equating coefficients 
	using moments. \emph{Economic Review (Otaru University of Commerce)}, 
	\bold{51}, 1--23.

	Ogasawara, H. (2001). Standard Errors of Item Response Theory Equating/Linking
	by Response Function Methods. \emph{Applied Psychological Measurement}, 
	\bold{25}, 53--67.
}
\author{
	Michela Battauz
}

\seealso{
	\code{\link{direc}}, \code{\link{eqc}}, \code{\link{itm}}, \code{\link{modIRT}}, \code{\link{summary.eqclist}}
}
\examples{

# the following code can be used to start from item response data

\dontrun{
library("mirt")
data("data2pl")

m1 <- mirt(data2pl[[1]], SE = TRUE)
m2 <- mirt(data2pl[[2]], SE = TRUE)
m3 <- mirt(data2pl[[3]], SE = TRUE)
m4 <- mirt(data2pl[[4]], SE = TRUE)
m5 <- mirt(data2pl[[5]], SE = TRUE)

mlist<- list(m1,m2,m3,m4,m5)
test <- paste("test", 1:5, sep = "")
mod2pl <- modIRT(est.mods = mlist, names = test, display = FALSE)

direclist2pl <- alldirec(mods = mod2pl, method = "Haebara")
summary(direclist2pl)
summary(direclist2pl$test2.test1)
}

# ===========================================================================
# the following code uses item parameter estimates previously obtained

# three-parameter logistic model
# direct equating coefficients using the "Stocking-Lord" method
data(est3pl)
test <- paste("test", 1:5, sep = "")
mod3pl <- modIRT(coef = est3pl$coef, var = est3pl$var, names = test, display = FALSE)
direclist3pl <- alldirec(mods = mod3pl, method = "Stocking-Lord")
summary(direclist3pl)
summary(direclist3pl$test1.test2)

# two-parameter logistic model
# direct equating coefficients using the "Haebara" method
data(est2pl)
test <- paste("test", 1:5, sep = "")
mod2pl <- modIRT(coef = est2pl$coef, var = est2pl$var, names = test, display = FALSE)
direclist2pl <- alldirec(mods = mod2pl, method = "Haebara")
summary(direclist2pl)
summary(direclist2pl$test1.test5)

# Rasch model
# direct equating coefficients using the "mean-mean" method
data(estrasch)
test <- paste("test", 1:5, sep = "")
modrasch <- modIRT(coef = estrasch$coef, var = estrasch$var, names = test, 
	display = FALSE)
direclistrasch <- alldirec(mods = modrasch, method = "mean-mean", all = TRUE)
summary(direclistrasch)
summary(direclistrasch$test5.test4)
}
