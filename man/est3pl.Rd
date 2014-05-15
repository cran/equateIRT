\name{est3pl}
\alias{est3pl}
\docType{data}
\title{
	Item parameter estimates and covariance matrices of a three-parameter logistic model
}
\description{
	This dataset includes item parameter estimates and covariance matrices of a
	three-parameter logistic model applied to 5 simulated datasets with common items.
	See details for more information on the linkage plan.
}
\usage{
data(est3pl)
}
\format{
	A list of length 2 with components:
	\describe{
	\item{\code{coef}}{a list of length 5 containing the matrices of item parameter 
		estimates. Each matrix presents 3 columns; the first column contains guessing
		parameters, the second column contains difficulty parameters and the third
		column contains discrimination parameters. See details for information on
		the parameterization used. Names of rows correspond to the
		names of the items.}
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
	Under this parameterization, the three-parameter logistic model is as follows
	\deqn{\pi_i = c_i + (1 - c_i) \frac{\exp(\beta_{1i} + \beta_{2i} z)}{1 + 
    \exp(\beta_{1i} + \beta_{2i} z)},}{ \pi_i = c_i + (1 - c_i) * \{exp(\beta_{1i} + \beta_{2i} z)\}/
	\{1 + exp(\beta_{1i} + \beta_{2i} z)\}, } where 
    \eqn{\pi_i} denotes the conditional probability of responding correctly to the \eqn{i}th item given \eqn{z}, 
    \eqn{c_i} denotes the guessing parameter, \eqn{\beta_{1i}} is the easiness parameter, 
    \eqn{\beta_{2i}} is the discrimination parameter, and \eqn{z} denotes the 
    latent ability.
	Furthermore, the guessing parameters are given under this parameterization
	\deqn{c_i = \frac{\exp(c_i^*)}{1+\exp(c_i^*)}. }{c_i = exp(c_i^*)/{1+exp(c_i^*)}. }
}
\author{
	Michela Battauz
}

\seealso{
	\code{\link{linkp}}, \code{\link{modIRT}}
}
\examples{
data(est3pl)
est3pl$coef
est3pl$var
linkp(coef = est3pl$coef)
}
\keyword{datasets}
