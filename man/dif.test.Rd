\name{dif.test}
\alias{dif.test}
\alias{print.dift}

\title{
	DIF Test
}
\description{
	Performs a Wald-type test for Differential Item Functioning detection.
}
\usage{
dif.test(est.mods = NULL, coef = NULL, var = NULL, names = NULL, 
  reference = NULL, method = "mean-mean",
  quadrature = TRUE, nq = 30, DIFtype = NULL, purification = FALSE, 
  signif.level = 0.05, trace = FALSE, maxiter = 30, anchor = NULL)
}

\arguments{
  \item{est.mods}{list of output objects from functions 
    \code{rasch}, \code{ltm}, or \code{tpm} 
		of the \code{ltm} package or from function \code{mirt} of the \code{mirt} package.}
	\item{coef}{list of matrices (one for each group) containing 
		the item parameter estimates. Guessing, difficulty and discrimination 
		parameters should strictly be given in this order and they are 
		contained in different columns of the matrix.
		The names of the	rows of each matrix should be the names of the items.}
	\item{var}{list of matrices (one for each group) containing the 
		covariance matrix of item parameter estimates. They should be
		given in the same order of coefficients.}
	\item{names}{character vector containing the names of the groups. 
		This should have the same length of \code{coef} and \code{var}.
		If \code{NULL}, the names of the groups will be "T1", "T2", ...}
	\item{reference}{
	  reference group. Can be specified by name or number.
	  The default is the first group.
	}
	\item{method}{the equating method to be used in function 
	  \code{\link{direc}} to convert the item parameters
	  to the scale of the reference group. This should be one of
  	"\code{mean-mean}", "\code{mean-gmean}", "\code{mean-sigma}",
  	"\code{Haebara}" or "\code{Stocking-Lord}".
	}
	\item{quadrature}{logical; if TRUE the Gauss-Hermite quadrature is used
  	in function \code{\link{direc}}
  	to approximate the integral in the function that is minimized in the 
  	Haebara and Stocking-Lord methods.
  	If FALSE the integral is replaced with a sum over 40
  	equally spaced values ranging from -4 to 4 with an increment 
  	of 0.05 and weights equal to one for all values.
	}
	\item{nq}{number of quadrature points used for the Gauss-Hermite quadrature
  	if \code{quadrature} is TRUE.
	}
	\item{DIFtype}{character indicating which parameters to test for DIF.
	  If \code{NULL} all parameters are tested for DIF.
	  Use "\code{beta1}" for \eqn{\beta_1}, "\code{beta2}" for \eqn{\beta_2},
	  "\code{beta3}" for \eqn{\beta_3}, "\code{beta12}" for \eqn{\beta_1}
	  and \eqn{\beta_2}, "\code{beta123}" for \eqn{\beta_1}, \eqn{\beta_2}
	  and \eqn{\beta_3}. See details.
	}
	\item{purification}{logical. if TRUE the procedure described in
	  Candell and Drasgow (1988) is applied.
	}
	\item{signif.level}{
	  significance level to use in the purification process.
	}
	\item{trace}{
	  logical. If TRUE tracing information is produced.
	}
	\item{maxiter}{
	  The maximum number of iterations in the purification process.
	}
	\item{anchor}{
	  Optional character vector containing the names of the items 
	  to use for equating. These should be items free of DIF.
	}
}
\details{
	The parameterization of the IRT model is that commonly used for estimation.
	Under this parameterization, the three-parameter logistic model is as follows
	\deqn{\pi_i = c_i + (1 - c_i) \frac{\exp(\beta_{1i} + \beta_{2i} z)}{1 + 
	\exp(\beta_{1i} + \beta_{2i} z)},}{ \pi_i = c_i + (1 - c_i) * 
	\{exp(\beta_{1i} + \beta_{2i} z)\}/
	\{1 + exp(\beta_{1i} + \beta_{2i} z)\}, } where 
	\eqn{\pi_i} denotes the conditional probability of responding correctly to 
	the \eqn{i}th item given \eqn{z}, 
	\eqn{c_i} denotes the guessing parameter, \eqn{\beta_{1i}} is the easiness parameter, 
	\eqn{\beta_{2i}} is the discrimination parameter, and \eqn{z} denotes the 
	latent ability.
  Furthermore, the guessing parameter is equal to
	\deqn{c_i = \frac{\exp(\beta_{3i})}{1+\exp(\beta_{3i})}}.
	
	The test verifies whether the item parameters 
	\eqn{\beta_{1i}, \beta_{2i}, \beta_{3i}}
	are invariant across two or more groups as explained in
	Battauz (2018).
}
\value{
	An object of class \code{dift} with components
	\item{test}{matrix containing the test statistic and the p-value for each item. 
	  "noGuess" is equal to 1 if the guessing parameter of a 3PL model was set
	  to a fixed value.}
	\item{eqmet}{the equating method used.}
	\item{DIFtype}{character indicating which parameters were tested for DIF.}
	\item{reference}{the reference group.}
	\item{focal}{the focal groups.}
	\item{names}{names of the groups.}
	\item{purification}{logical. If TRUE thre purification procedure was applied.}
	\item{signif.level}{significance level used in the purification process}
	\item{equatings}{list containing the output of function \code{\link{direc}}.}
	\item{coef_trasf}{list containing the item parameters of each group 
	  transformed to the scale of the reference group.}
	\item{var_trasf}{list containing the covariance matrix of 
	  item parameters of each group transformed to the scale of the reference group.}
	\item{items.dif}{names of the items for which the null hypothesis of the test is rejected.}
	\item{anchor}{names of the items used as anchors.}
	\item{niter}{number of iterations.}
}
\references{
	Battauz, M. (2019). On Wald tests for differential item functioning detection.
	\emph{Statistical Methods and Applications}. \bold{28}, 103-118.
	
	Candell,  G.L., Drasgow, F. (1988). An iterative procedure for linking 
	metrics and assessing item bias in item response theory. 
	\emph{Applid Psychological Measurement}, \bold{12}, 253-260.
}
\author{
	Michela Battauz
}
\examples{
# load the data
data(dataDIF)
head(dataDIF)
# estimate a 2PL model for each group using the R package mirt
library(mirt)
data1 <- dataDIF[dataDIF$group == 1, 1:20]
data2 <- dataDIF[dataDIF$group == 2, 1:20]
data3 <- dataDIF[dataDIF$group == 3, 1:20]
mod1 <- mirt(data1, SE = TRUE)
mod2 <- mirt(data2, SE = TRUE)
mod3 <- mirt(data3, SE = TRUE)

# since package version 2.5.0 it is possible to skip function import.mirt
# as follows

# perform the test for DIF on two groups
res_diftest2 <- dif.test(est.mods = list(mod1, mod2))
res_diftest2

# perform the test for DIF on three groups
res_diftest3 <- dif.test(est.mods = list(mod1, mod2, mod3))
res_diftest3

# perform the test for DIF on three groups
# reference group: 2
# equating method: Haebara
# purification applied
res_diftest3 <- dif.test(est.mods = list(mod1, mod2, mod3), 
  reference = 2, method = "Haebara", purification = TRUE)
res_diftest3
}

