\name{alldirec}
\alias{alldirec}
\title{
	Direct equating coefficients between all pairs of a list of forms
}
\description{
	Calculates direct equating coefficients and standard errors using 
	IRT methods between all pairs of a list of forms.
}
\usage{
alldirec(mods, method = "mean-mean", all = FALSE, itmp, quadrature = TRUE, nq = 30)
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
	\item{itmp}{
	number of item parameters of the IRT model. This should be 1
	for the Rasch model, 2 for the two-parameter logistic model and 3 for 
	the three-parameter logistic model. It should be coherent with the item parameters
	contained in \code{mods}.
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
}
\value{
	An object of class \code{eqclist} consisting in a list with length 
	equal to the number of pairs of forms equated. Each component of
	the list is an object of class \code{eqc} returned by function
	\code{\link{direc}}.
}
\references{
	Kolen, M.J. and Brennan, R.L. (2004). \emph{Test equating, scaling, 
	and linking: methods and practices}, 2nd ed., New York: Springer

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
	\code{\link{direc}}, \code{\link{modIRT}}, \code{\link{summary.eqclist}}
}
\examples{
# three-parameter logistic model
# direct equating coefficients using the "Stocking-Lord" method
data(est3pl)
test <- paste("test", 1:5, sep = "")
mod3pl <- modIRT(coef = est3pl$coef, var = est3pl$var, names = test, ltparam = TRUE, lparam = TRUE)
direclist3pl <- alldirec(mods = mod3pl, method = "Stocking-Lord", itmp = 3)
summary(direclist3pl)
summary(direclist3pl$test1.test2)

# two-parameter logistic model
# direct equating coefficients using the "Haebara" method
data(est2pl)
test <- paste("test", 1:5, sep = "")
mod2pl <- modIRT(coef = est2pl$coef, var = est2pl$var, names = test, ltparam = TRUE)
direclist2pl <- alldirec(mods = mod2pl, method = "Haebara", itmp = 2)
summary(direclist2pl)
summary(direclist2pl$test1.test5)

# Rasch model
# direct equating coefficients using the "mean-mean" method
data(estrasch)
test <- paste("test", 1:5, sep = "")
modrasch <- modIRT(coef = estrasch$coef, var = estrasch$var, names = test, ltparam = TRUE)
direclistrasch <- alldirec(mods = modrasch, method = "mean-mean", all = TRUE, itmp = 1)
summary(direclistrasch)
summary(direclistrasch$test5.test4)
}
