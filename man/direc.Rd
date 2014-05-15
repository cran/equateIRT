\name{direc}
\alias{direc}

\title{
	Direct equating coefficients
}
\description{
	Calculates direct equating coefficients and standard errors using 
	IRT methods.
}
\usage{
direc(mod1, mod2, method = "mean-mean", suff1 = ".1", suff2 = ".2",
   D = 1, quadrature = TRUE, nq = 30)
}
\arguments{
	\item{mod1}{an object of the class \code{modIRT} containing 
	item parameter coefficients and their covariance matrix
	of the first form.
	}
	\item{mod2}{an object of the class \code{modIRT} containing 
	item parameter coefficients and their covariance matrix
	of the second form.
	}
	\item{method}{the equating method to be used. This should be one of
	"\code{mean-mean}", "\code{mean-gmean}", "\code{mean-sigma}",
	"\code{Haebara}" or "\code{Stocking-Lord}".
	}
	\item{suff1}{suffix to identify the first form to be equated.
	}
	\item{suff2}{suffix to identify the second form to be equated.
	}
	\item{D}{constant D of the IRT model used to estimate item 
	parameters. See below for more details.
	}
	\item{quadrature}{logical; if TRUE the Gauss-Hermite quadrature is used
	to approximate the integral in the function that is minimized in the 
	Haebara and Stocking-Lord methods.
	If FALSE the integral is replaced with a sum over 40
	equally spaced values ranging from -4 to 4 with an increment 
	of 0.05 and weights equal to one for all values.
	}
	\item{nq}{number of quadrature points used for the Gauss-Hermite quadrature
	if \code{quadrature} is TRUE.
	}
}
\details{
	Equating coefficients perform the conversion from the scale of the 
	first form to the scale of the second form.

	In the three-parameter logistic model  the probability of a 
	positive response on item \eqn{i} is
	\deqn{\pi_i = c_i + (1 - c_i) \frac{\exp[D a_i (\theta - b_i)]}{1 + 
	\exp[D a_i (\theta - b_i)]}.}{ \pi_i = c_i + (1 - c_i) *  
	exp\{D a_i (\theta - b_i)\} / [1 + exp\{D a_i (\theta - b_i)\}],}
	where \eqn{a_i} is the item discrimination parameter, \eqn{b_i} is 
	the item difficulty parameter, \eqn{c_i} is the item guessing parameter 
	and	\eqn{\theta} is the latent ability.
	The constant \eqn{D} can be specified using argument \code{D} of the 
	\code{direc} function.
	The two-parameter logistic model can be obtained by setting \eqn{c_i} 
	equal to zero, the one-parameter logistic model can be obtained 
	by setting \eqn{c_i} equal to zero and \eqn{a_i} costant across items,
	while the Rasch model can be obtained by setting \eqn{c_i} 
	equal to zero and \eqn{a_i} equal to 1.
}
\value{
	An object of class \code{eqc} with components
	\item{tab1}{item parameters of the first form.}
	\item{tab2}{item parameters of the second form.}
	\item{tab}{data frame containing item parameters of the first form (value1),
		the second form (value2), and the first form converted in the scale of the second
		form (value12).}
	\item{var12}{covariance matrix of item parameters of the first 
		and the second form.}
	\item{partial}{partial derivatives of equating coefficients A and B
		with respect to item parameters.}
	\item{A}{equating coefficient A.}
	\item{B}{equating coefficient B.}
	\item{varAB}{covariance matrix of the equating coefficients.}
	\item{commonitem}{list of length 1 containing the names of common item parameters.}
	\item{ni}{number of common items.}
	\item{forms}{names of equated forms.}
	\item{method}{the equating method used.}
	\item{itmp}{number of item parameters of the IRT model.}
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
	\code{\link{modIRT}}, \code{\link{summary.eqc}}
}
\examples{
# three-parameter logistic model
# direct equating coefficients between forms 1 and 2 using the Stocking-Lord method
data(est3pl)
test <- paste("test", 1:5, sep = "")
mod3pl <- modIRT(coef = est3pl$coef, var = est3pl$var, names = test, display = FALSE)
l12 <- direc(mod1 = mod3pl[1], mod2 = mod3pl[2], method = "Stocking-Lord")
summary(l12)

# two-parameter logistic model
# direct equating coefficients between forms 1 and 5 using the Haebara method
data(est2pl)
test <- paste("test", 1:5, sep = "")
mod2pl <- modIRT(coef = est2pl$coef, var = est2pl$var, names = test, display = FALSE)
l15 <- direc(mod1 = mod2pl[1], mod2 = mod2pl[5], method = "Haebara")
summary(l15)

# Rasch model
# direct equating coefficients between forms 5 and 4 using the mean-mean method
data(estrasch)
test <- paste("test", 1:5, sep = "")
modrasch <- modIRT(coef = estrasch$coef, var = estrasch$var, names = test, 
	display = FALSE)
l54 <- direc(mod1 = modrasch[5], mod2 = modrasch[4], method = "mean-mean")
summary(l54)
}

