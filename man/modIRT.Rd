\name{modIRT}
\alias{modIRT}

\title{
	Estimated coefficients and covariance matrix of IRT models 
}
\description{
	Creates an object of the class \code{modIRT} containing estimated
	coefficients and covariance matrix of IRT models. 
	Rasch, two-parameter logistic and three-parameter logistic models
	are included.
}
\usage{
modIRT(coef, var, names, ltparam = FALSE, lparam = FALSE)
}
\arguments{
	\item{coef}{list of matrices (one for each form) containing 
		item parameter estimates. Guessing, difficulty and discrimination 
		parameters should strictly be given in this order and they are 
		contained in different columns of the matrix.
		The name of the	rows of each matrix should be the name of the item.}
	\item{var}{list of matrices (one for each form) containing the 
		covariance matrix of item parameter estimates. They should be
		given in the same order of coefficients.}
	\item{names}{character vector containing the names of the forms. 
		This should have the same length of \code{coef} and \code{var}.}
	\item{ltparam}{logical; if TRUE the latent trait parametrization is used
		for difficulty parameters. See below for more details.}
	\item{lparam}{logical; if TRUE the logistic parametrization is used
		for guessing parameters. See below for more details.}
}
\details{
	If \code{ltparam} is TRUE the latent trait parametrization is used 
	as in functions \code{rasch}, \code{ltm} and \code{tpm} of the \code{ltm} package.
	Under this parametrization, the three-parameter logistic model is as follows
	\deqn{\pi_i = c_i + (1 - c_i) \frac{\exp(\beta_{1i} + \beta_{2i} z)}{1 + 
	\exp(\beta_{1i} + \beta_{2i} z)},}{ \pi_i = c_i + (1 - c_i) * \{exp(\beta_{1i} + \beta_{2i} z)\}/
	\{1 + exp(\beta_{1i} + \beta_{2i} z)\}, } where 
	\eqn{\pi_i} denotes the conditional probability of responding correctly to the \eqn{i}th item given \eqn{z}, 
	\eqn{c_i} denotes the guessing parameter, \eqn{\beta_{1i}} is the easiness parameter, 
	\eqn{\beta_{2i}} is the discrimination parameter, and \eqn{z} denotes the 
	latent ability.
	The two-parameter logistic model and the Rasch model present the same
	formulation. The two-parameter logistic model can be obtained by setting \eqn{c_i} equal to zero,
	while the Rasch model can be obtained by setting \eqn{c_i} equal to zero and
	\eqn{\beta_{2i}} equal to 1.
	
	If \code{lparam} is TRUE the guessing parameters are given under this
	parametrization
	\deqn{c_i = \frac{\exp(c_i^*)}{1+\exp(c_i^*)}. }{c_i = exp(c_i^*)/{1+exp(c_i^*)}. }

	The \code{modIRT} function returns parameter estimates
	under the usual IRT parameterization, that is,
	\deqn{\pi_i = c_i + (1 - c_i) \frac{\exp[\beta_{2i} (z - \beta_{1i}^*)]}{1 + 
	\exp[\beta_{2i} (z - \beta_{1i}^*)]}.}{ \pi_i = c_i + (1 - c_i) *  exp\{\beta_{2i} (z - \beta_{1i}^*)\} / [1 + exp\{\beta_{2i} (z - \beta_{1i}^*)\}].}

	If \code{ltparam} or \code{lparam} are TRUE, the covariance matrix 
	is calculated using the delta method.
}
\value{
	An object of class \code{modIRT} consisting in a list with length equal
	to the number of forms containing lists with components
	\item{coef}{item parameter estimates.}
	\item{var}{covariance matrix of item parameter estimates.}
}
\references{
	Bartholomew, D., Knott, M. and Moustaki, I. (2011) \emph{Latent Variable Models and Factor Analysis: 
	a Unified Approach}, 3rd ed. Wiley.
	
	Rizopoulos, D. (2006). ltm: an R package for latent variable modelling 
	and item response theory analyses. \emph{Journal of Statistical Software}, 
	\bold{17}, 1--25.
}
\author{
	Michela Battauz
}

\seealso{
	\code{\link{direc}}
}
\examples{
# three-parameter logistic model
data(est3pl)
test <- paste("test", 1:5, sep = "")
mod3pl <- modIRT(coef = est3pl$coef, var = est3pl$var, names = test, ltparam = TRUE, lparam = TRUE)

# two-parameter logistic model
data(est2pl)
test <- paste("test", 1:5, sep = "")
mod2pl <- modIRT(coef = est2pl$coef, var = est2pl$var, names = test, ltparam = TRUE)

# Rasch model
data(estrasch)
test <- paste("test", 1:5, sep = "")
modrasch <- modIRT(coef = estrasch$coef, var = estrasch$var, names = test, ltparam = TRUE)
}

