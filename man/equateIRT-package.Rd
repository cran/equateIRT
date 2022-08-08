\name{equateIRT-package}
\alias{equateIRT-package}
\alias{equateIRT}
\docType{package}
\title{
	IRT Equating Methods
}
\description{
	This package computes direct, chain and average (bisector) equating coefficients 
	with standard errors using IRT methods for dichotomous items. 
	The IRT models included are the three-parameter logistic model,
	the two-parameter logistic model, the one-parameter logistic model and the
	Rasch model.
	Test scoring can be performed by true score equating and observed score equating methods.
	DIF detection can be performed using a Wald-type test.
}
\details{
\tabular{ll}{
Package: \tab equateIRT\cr
Type: \tab Package\cr
Version: \tab 2.3.0\cr
Date: \tab 2022-08-08\cr
License: \tab GPL-3\cr
}
Direct equating coefficients and their standard errors
between forms presenting common items  
can be computed using function \code{direc}.
The equating methods implemented are 
"mean-mean", "mean-geometric mean", "mean-sigma", "Haebara" and "Stocking-Lord".
Estimates of item parameters and their covariance matrix can be imported
from the R packages \code{ltm} and \code{mirt} or from the IRT programs IRTPRO and flexMIRT using
functions \code{import.ltm}, \code{import.mirt}, \code{import.irtpro} and \code{import.flexmirt}.
Item parameter estimates from other software can be imported as well
by the user.
Data should be previously organized using function \code{modIRT}.
Function \code{alldirec} computes all direct equating coefficients
(with standard errors) between all pairs of a list of forms.
Function \code{chainec} computes chain equating coefficients
(and standard errors)
given direct equating coefficients between forms directly linked.
Average equating coefficients with standard errors 
can be calculated using function \code{bisectorec}, that 
implements the (weighted) bisector method.
Once the equating coefficients are obtained, the computation of 
equated scores can be performed using function \code{score},
which implements true score equating and observed score equating.
Standard errors of equated scores are also provided.
Three simulated datasets are available for illustrative purposes.
These datasets contain item parameter coefficients and their 
covariance matrix. In particular, \code{est3pl} concerns a 
three-parameter logistic model, \code{est2pl} regards a 
two-parameter logistic model, and \code{estrasch} refers to
a Rasch model.
The estimates included in \code{est2pl} are obtained from the
dataset \code{data2pl}, also contained in the package.
Function \code{dif.test} performs a Wald-type test for the 
detection of DIF (Battauz, 2018).
The stability of the equating transformations can be assessed using
function \code{sd.test}, which tests if the equating coefficients
that link the same two forms are equal, and function \code{id.test},
which tests if the equating coefficients that link a form to itself
through a chain of forms return the identity equating.
}
\author{
	Michela Battauz
	
	Maintainer: Michela Battauz <michela.battauz@uniud.it>
}
\references{
	Battauz, M. (2013). IRT Test Equating in Complex Linkage Plans. 
	\emph{Psychometrika}, \bold{78}, 464--480.
	
	Battauz, M. (2015). equateIRT: An R Package for IRT Test Equating.
	\emph{Journal of Statistical Software}, \bold{68}, 1--22.
	
	Battauz, M. (2019). On Wald tests for differential item functioning detection.
	\emph{Statistical Methods and Applications}. \bold{28}, 103--118.
	
	Battauz, M. (2022). Testing for differences in chain equating. 
	\emph{Statistica Neerlandica}, 1--12.
	
	Cai L. (2013). \emph{FlexMIRT version 2: Flexible Multilevel Multidimensional Item Analysis and
	Test Scoring [Computer Software]}. Chapel Hill, NC: Vector Psychometric Group.
	
	Cai, L., du Toit, S. H. C., Thissen, D. (2011). \emph{IRTPRO: Flexible, multidimensional, multiple 
	categorical IRT modeling [Computer software]}. Chicago: Scientific Software International.

	Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory Package for the 
	R Environment. \emph{Journal of Statistical Software}, 
	\bold{48}, 1--29.
	
	Holland, P.W. and Strawderman, W.E. (2011). How to average equating 
	functions if you must. In A.A. von Davier (Ed.), \emph{Statistical models for 
	test equating, scaling, and linking} (pp. 89--107). New York: Springer.

	Kolen, M.J. and Brennan, R.L. (2014). \emph{Test equating, scaling, 
	and linking: methods and practices}, 3rd ed., New York: Springer.

	Ogasawara, H. (2000). Asymptotic standard errors of IRT equating coefficients 
	using moments. \emph{Economic Review (Otaru University of Commerce)}, 
	\bold{51}, 1--23.

	Ogasawara, H. (2001a). Item response theory true score equatings and
	their standard errors. \emph{Journal of Educational and Behavioral Statistics}, 
	\bold{26}, 31--50.

	Ogasawara, H. (2001b). Standard Errors of Item Response Theory Equating/Linking
	by Response Function Methods. \emph{Applied Psychological Measurement}, 
	\bold{25}, 53--67.

	Ogasawara, H. (2003). Asymptotic standard errors of IRT 
	observed-score equating methods. \emph{Psychometrika}, 
	\bold{68}, 193--211.

	Rizopoulos, D. (2006). ltm: an R package for latent variable modelling 
	and item response theory analyses. \emph{Journal of Statistical Software}, 
	\bold{17}, 1--25.
}
