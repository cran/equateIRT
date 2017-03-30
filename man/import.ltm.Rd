\name{import.ltm}
\alias{import.ltm}
\alias{import.mirt}
\alias{import.irtpro}
\alias{import.flexmirt}
\title{
	Import Item Parameters Estimates and Covariance Matrices from IRT Software
}
\description{
	Import estimated item parameters and covariance matrix from the R packages ltm
	and mirt, and from external software IRTPRO and flexMIRT.
}
\usage{
import.ltm(mod, display = TRUE, digits = 4) 

import.mirt(mod, display = TRUE, digits = 3) 

import.irtpro(fnamep, fnamev = NULL, fnameirt = NULL, display = TRUE, digits = 2) 

import.flexmirt(fnamep, fnamev = NULL, fnameirt = NULL, display = TRUE, digits = 2) 
}
\arguments{
	\item{mod}{output object from functions \code{rasch}, \code{ltm}, or \code{tpm} 
		of the \code{ltm} package or from function \code{mirt} of the \code{mirt} package.}
	\item{display}{logical; if TRUE coefficients and standard errors are printed.}
	\item{digits}{integer indicating the number of decimal places to be used
		if \code{display} is \code{TRUE}.}
	\item{fnamep}{name of the file containing the estimated item parameters. Typically, -prm.txt.}
	\item{fnamev}{name of the file containing the covariance matrix of the estimated item parameters. 
		Typically, -cov.txt.}
	\item{fnameirt}{name of the file containing additional information to link item parameters with 
		the covariance matrix. Typically, -irt.txt.}
}
\details{
	Item parameters are imported with the parameterization used by the software to
	estimate the IRT model. The usual IRT parameterization can be obtained later
	by using function \code{\link{modIRT}}.
}
\value{
	A list with components
 	\item{coef}{item parameter estimates.}
	\item{var}{covariance matrix of item parameter estimates.}
}
\references{
	Battauz, M. (2015). equateIRT: An R Package for IRT Test Equating.
	\emph{Journal of Statistical Software}, \bold{68}, 1--22.
	
	Cai L. (2013). \emph{FlexMIRT version 2: Flexible Multilevel Multidimensional Item Analysis and
	Test Scoring [Computer Software]}. Chapel Hill, NC: Vector Psychometric Group.
	
	Cai, L., du Toit, S. H. C., Thissen, D. (2011). \emph{IRTPRO: Flexible, multidimensional, multiple 
	categorical IRT modeling [Computer software]}. Chicago: Scientific Software International.
	
	Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory Package for the 
	R Environment. \emph{Journal of Statistical Software}, 
	\bold{48}, 1--29.
	
	Rizopoulos, D. (2006). ltm: an R package for latent variable modelling 
	and item response theory analyses. \emph{Journal of Statistical Software}, 
	\bold{17}, 1--25.
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link{modIRT}}
}
\examples{

#====================
# from package ltm
library(ltm)

# one-parameter logistic model
mod1pl <- rasch(LSAT)
est.mod1pl <- import.ltm(mod1pl)
est.mod1pl

# two-parameter logistic model
mod2pl <- ltm(LSAT ~ z1)
est.mod2pl <- import.ltm(mod2pl)
est.mod2pl

#====================
# from package mirt
library(mirt)

# one-parameter logistic model
data(LSAT, package = "ltm")
val <- mirt(LSAT, 1, SE = TRUE, pars = "values")
cnstr <- val[val$name == "a1",]$parnum
mod1pl.m <- mirt(LSAT, 1, SE = TRUE, SE.type = 'Richardson', constrain = list(cnstr))
est.mod1pl.m <- import.mirt(mod1pl.m, digits = 4)
est.mod1pl.m

# two-parameter logistic model
data(LSAT, package = "ltm")
mod2pl.m <- mirt(LSAT, 1, SE = TRUE, SE.type = 'Richardson')
est.mod2pl.m <- import.mirt(mod2pl.m, display = FALSE)
est.mod2pl.m

}
