\name{estrasch}
\alias{estrasch}
\docType{data}
\title{
	Item Parameter Estimates and Covariance Matrices of a Rasch Model
}
\description{
	This dataset includes item parameter estimates and covariance matrices of a
	Rasch model applied to 5 simulated datasets with common items.
	See details for more information on the linkage plan.
}
\usage{
data(estrasch)
}
\format{
	A list of length 2 with components:
	\describe{
	\item{\code{coef}}{a list of length 5 containing the matrices of item parameter 
		estimates. Each matrix presents 2 columns; the first column contains 
		difficulty parameters and the second column is equal to 1. 
		See details for information on the parameterization used. 
		Names of rows correspond to the	names of the items.}
	\item{\code{var}}{a list of length 5 containing the covariance matrices of 
		item parameter estimates.}
	}
}
\details{
	Every form is composed by 20 items and presents 10 items in common
	with adjacent forms. Furthermore, forms 1 and 5 present 10 common items.
	Use \code{\link{linkp}} to obtain a matrix with elements equal to the
	number of common items between different forms.

	Item parameters are given under the parameterization used in the \code{ltm} package. 
	Under this parameterization, the Rasch model is as follows
	\deqn{\pi_i = \frac{\exp(\beta_{1i} + z)}{1 + \exp(\beta_{1i} + z)},}{ \pi_i = \{exp(\beta_{1i} + z)\}/\{1 + exp(\beta_{1i} + z)\}, } 
	where \eqn{\pi_i} denotes the conditional probability of responding correctly to the \eqn{i}th item given \eqn{z}, 
    \eqn{\beta_{1i}} is the easiness parameter, and \eqn{z} denotes the 
    latent ability.
}
\author{
	Michela Battauz
}

\seealso{
	\code{\link{linkp}}, \code{\link{modIRT}}
}
\examples{
data(estrasch)
estrasch$coef
estrasch$var
linkp(coef = estrasch$coef)
}
\keyword{datasets}
