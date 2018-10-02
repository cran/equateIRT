\name{modIRT}
\alias{modIRT}

\title{
	Estimated Coefficients and Covariance Matrix of IRT Models 
}
\description{
	Creates an object of the class \code{modIRT} containing estimated
	coefficients and covariance matrices of IRT models. 
	Rasch, one-parameter logistic, two-parameter logistic and three-parameter 
	logistic models	are included.
}
\usage{
modIRT(coef, var = NULL, names = NULL, ltparam = TRUE, lparam = TRUE, 
	display = TRUE, digits = 2)
}
\arguments{
	\item{coef}{list of matrices (one for each form) containing 
		item parameter estimates. Guessing, difficulty and discrimination 
		parameters should strictly be given in this order and they are 
		contained in different columns of the matrix.
		The names of the	rows of each matrix should be the names of the items.}
	\item{var}{list of matrices (one for each form) containing the 
		covariance matrix of item parameter estimates. They should be
		given in the same order of coefficients.}
	\item{names}{character vector containing the names of the forms. 
		This should have the same length of \code{coef} and \code{var}.
		If NULL, the names of the forms are assigned by 
		function \code{modIRT}.}
	\item{ltparam}{logical; if TRUE the latent trait parameterization is used
		for difficulty parameters and the \code{modIRT} function performs a 
		transformation of item parameters to return them in the usual IRT
		parameterization. Set to FALSE to avoid transformations. See below for more details.}
	\item{lparam}{logical; if TRUE the logistic parameterization is used
		for guessing parameters and the \code{modIRT} function performs a 
		transformation of item parameters to return them in the usual IRT
		parameterization. Set to FALSE to avoid transformations. See below for more details.}
	\item{display}{logical; if TRUE coefficients and standard errors are printed.}
	\item{digits}{integer indicating the number of decimal places to be used
		if \code{display} is \code{TRUE}.}
}
\details{
	\code{ltparam} and \code{lparam} refers the the parameterization used
	by the software used to estimate item parameters. The R package \code{ltm},
	and the programs IRTPRO and flexMIRT use these parameterizations.
	If \code{ltparam} is \code{TRUE} the latent trait parameterization is used. 
	Under this parameterization, the three-parameter logistic model is as follows
	\deqn{\pi_i = c_i + (1 - c_i) \frac{\exp(\beta_{1i} + \beta_{2i} z)}{1 + 
	\exp(\beta_{1i} + \beta_{2i} z)},}{ \pi_i = c_i + (1 - c_i) * \{exp(\beta_{1i} + \beta_{2i} z)\}/
	\{1 + exp(\beta_{1i} + \beta_{2i} z)\}, } where 
	\eqn{\pi_i} denotes the conditional probability of responding correctly to the \eqn{i}th item given \eqn{z}, 
	\eqn{c_i} denotes the guessing parameter, \eqn{\beta_{1i}} is the easiness parameter, 
	\eqn{\beta_{2i}} is the discrimination parameter, and \eqn{z} denotes the 
	latent ability.
	The two-parameter logistic model, the one-parameter logistic model 
	and the Rasch model present the same
	formulation. The two-parameter logistic model can be obtained 
	by setting \eqn{c_i} equal to zero,
	the one-parameter logistic model can be obtained 
	by setting \eqn{c_i} equal to zero and \eqn{\beta_{2i}} 
	costant across items,
	while the Rasch model can be obtained by setting \eqn{c_i} 
	equal to zero and \eqn{\beta_{2i}} equal to 1.
	
	If \code{lparam} is TRUE the guessing parameters are given under this
	parameterization
	\deqn{c_i = \frac{\exp(c_i^*)}{1+\exp(c_i^*)}. }{c_i = exp(c_i^*)/{1+exp(c_i^*)}. }

	The \code{modIRT} function returns parameter estimates
	under the usual IRT parameterization, that is,
	\deqn{\pi_i = c_i + (1 - c_i) \frac{\exp[D a_i (\theta - b_i)]}{1 + 
	\exp[D a_i (\theta - b_i)]},}{ \pi_i = c_i + (1 - c_i) *  
	exp\{D a_i (\theta - b_i)\} / [1 + exp\{D a_i (\theta - b_i)\}],}
	where \eqn{D a_i = \beta_{2i}}, \eqn{b_i = -\beta_{1i}/\beta_{2i}} and
	\eqn{\theta = z}.
	
	If \code{ltparam} or \code{lparam} are TRUE, the covariance matrix 
	is calculated using the delta method.
	
	If item parameters are already given under the usual IRT parameterization,
	arguments \code{ltparam} and \code{lparam} should be set to \code{FALSE}.
}
\value{
	An object of class \code{modIRT} consisting in a list with length equal
	to the number of forms containing lists with components
	\item{coefficients}{item parameter estimates.}
	\item{var}{covariance matrix of item parameter estimates.}
	\item{itmp}{number of item parameters of the IRT model. This is 1 for the 
	Rasch model, 2 for the one-parameter logistic model with constant discriminations, 
	2 for the	two-parameter logistic model and 3 for the three-parameter logistic
	model.}
}
\references{
	Battauz, M. (2015). equateIRT: An R Package for IRT Test Equating.
	\emph{Journal of Statistical Software}, \bold{68}, 1--22.
	
	Bartholomew, D., Knott, M. and Moustaki, I. (2011) 
	\emph{Latent Variable Models and Factor Analysis: 
	a Unified Approach}, 3rd ed. Wiley.
	
	Rizopoulos, D. (2006). ltm: an R package for latent variable modelling 
	and item response theory analyses. \emph{Journal of Statistical Software}, 
	\bold{17}, 1--25.
}
\author{
	Michela Battauz
}

\seealso{
	\code{\link{direc}}, \code{\link{import.ltm}}
}
\examples{
# three-parameter logistic model
data(est3pl)
test <- paste("test", 1:5, sep = "")
mod3pl <- modIRT(coef = est3pl$coef, var = est3pl$var, names = test, display = FALSE)

# two-parameter logistic model
data(est2pl)
test <- paste("test", 1:5, sep = "")
mod2pl <- modIRT(coef = est2pl$coef, var = est2pl$var, names = test, display = FALSE)

# Rasch model
data(estrasch)
test <- paste("test", 1:5, sep = "")
modrasch <- modIRT(coef = estrasch$coef, var = estrasch$var, names = test, 
	display = FALSE)

# one-parameter logistic model imported from the R package ltm
library(ltm)
mod1pl <- rasch(LSAT)
summary(mod1pl)
est.mod1pl <- import.ltm(mod1pl)
mod1pl.ltm <- modIRT(coef = list(est.mod1pl$coef), var = list(est.mod1pl$var), digits = 4)

}

