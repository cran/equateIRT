\name{score}
\alias{score}
\alias{score.eqc}
\alias{score.eqclist}
\alias{score.ceqc}
\alias{score.ceqclist}
\alias{score.meqc}

\title{
	Scoring
}
\description{
	Relates number-correct scores on two forms.
}
\usage{
score(obj, \dots)

\method{score}{eqc}(obj, method = "TSE", D = 1, scores = NULL, se = TRUE, nq = 30, 
  w = 0.5, theta = NULL, weights = NULL, anchor.type = "internal", \dots)

\method{score}{eqclist}(obj, link = NULL, method = "TSE", D = 1, scores = NULL, se = TRUE, nq = 30, 
  w = 0.5, theta = NULL, weights = NULL, anchor.type = "internal", \dots)

\method{score}{ceqc}(obj, method = "TSE", D = 1, scores = NULL, se = TRUE, nq = 30, 
  w = 0.5, theta = NULL, weights = NULL, anchor.type = "internal", \dots)

\method{score}{ceqclist}(obj, path = NULL, method = "TSE", D = 1, scores = NULL, se = TRUE, nq = 30, 
  w = 0.5, theta = NULL, weights = NULL, anchor.type = "internal", \dots)

\method{score}{meqc}(obj, link = NULL, method = "TSE", D = 1, scores = NULL, se = TRUE, 
  bistype = NULL, nq = 30, w = 0.5, theta = NULL, weights = NULL, 
  anchor.type = "internal", \dots)
}
\arguments{
	\item{obj}{output object from functions \code{\link{direc}}, 
		\code{\link{alldirec}}, \code{\link{chainec}} or \code{\link{bisectorec}}.}
	\item{link}{a character string with the names of the two forms 
		being linked separated by a dot (e.g. "\code{test1.test3}"). 
		Necessary if \code{obj} is output of functions \code{alldirec} or \code{bisectorec}.}
	\item{path}{a character string with the names of the forms that constitute 
		the path separated by a dot (e.g. "\code{test1.test2.test3}").
		Necessary if \code{obj} is output of function \code{chainec}.}
	\item{method}{the scoring method to be used. This should be one of "\code{TSE}" 
		(the default) for true score equating or "\code{OSE}" for observed score equating.}
	\item{D}{constant D of the IRT model used to estimate item parameters.}
	\item{scores}{integer values to be converted.}
	\item{se}{logical; is \code{TRUE} standard errors of equated scores are computed.}
	\item{bistype}{bisector type. It should be specified when \code{obj} is an output
		object from function \code{\link{bisectorec}}.
		It should be \code{weighted} or \code{unweighted}. 
		It can be \code{NULL} if only one type was computed by function \code{\link{bisectorec}}.}
	\item{nq}{number of quadrature points used to approximate integrals with observed score equating. 
		Used only if arguments \code{theta} and \code{weights} are \code{NULL}.}
	\item{w}{synthetic weight for population 1. It should be a number between 0 and 1.}
	\item{theta}{vector of ability values used to approximate integrals with observed score equating.}
	\item{weights}{vector of weights used to approximate integrals with observed score equating.}
	\item{anchor.type}{The anchor type, which can be "\code{internal}" or "\code{external}".
	External anchors are not implemented for bisector equating coefficients. See details.}
	\item{\dots}{further arguments passed to or from other methods.}
}
\details{
	If \code{anchor.type} is \code{internal}, the common items are used for 
	computing the scores; if it is \code{external}, the common items are used
	for computing the equating coefficients but not for computing the scores.
	External anchors are not implemented for bisector equating coefficients 
	because it is possible that all items are used as achors in one of the paths
	considered.

	This function computes standard error of equated scores with the observed score equating method
	using \eqn{(\alpha_X',\alpha_{V1}',A,B)'} for the description of \eqn{\beta_{X2,x}(At_m+B)}, and
	\eqn{(\alpha_Y',\alpha_{V2}',A,B)'} for \eqn{\beta_{Y1,y}(t_m)}
	(see Ogasawara, 2003).
}
\value{
	A data frame containing theta values (only for true score equating),
	scores in the form chosen as base,
	equated scores,
	and standard errors of equated scores.
}
\references{
	Kolen, M.J. and Brennan, R.L. (2014). \emph{Test equating, scaling, 
	and linking: methods and practices}, 3nd ed., New York: Springer.

	Ogasawara, H. (2001). Item response theory true score equatings and
	their standard errors. \emph{Journal of Educational and Behavioral Statistics}, 
	\bold{26}, 31--50.

	Ogasawara, H. (2003). Asymptotic standard errors of IRT 
	observed-score equating methods. \emph{Psychometrika}, 
	\bold{68}, 193--211.
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link{alldirec}}, \code{\link{bisectorec}}, \code{\link{chainec}}, \code{\link{direc}}
}
\examples{
# two-parameter logistic model
data(est2pl)
test<-paste("test", 1:5, sep = "")
mod2pl <- modIRT(coef = est2pl$coef, var = est2pl$var, names = test, display = FALSE)
# direct equating coefficients between forms 1 and 2 using the Haebara method
l12 <- direc(mods = mod2pl, which = c(1,2), method = "Haebara")
# scoring using direct equating coefficients and the true score equating method
score(l12)
# scoring using direct equating coefficients and the observed score equating method
score(l12, method = "OSE")
# specify only scores from 10 to 15
score(l12, method = "OSE", scores = 10:15)
# external anchors and true score equating method
score(l12, anchor.type = "external")
# external anchors and observed score equating method
score(l12, method = "OSE", anchor.type = "external")

# all direct equating coefficients using the Haebara method
direclist2pl <- alldirec(mods = mod2pl, method = "Haebara")
# scoring using direct equating coefficients and the true score equating method
score(direclist2pl, link = "test1.test2")

# compute chain equating coefficients for path 1,2,3,4
pth1 <- paste("test", 1:4, sep = "")
chainec1 <- chainec(direclist = direclist2pl, pths = pth1)
# compute chain equating coefficients for path 1,5,4
pth2 <- paste("test", c(1,5,4), sep = "")
chainec2 <- chainec(direclist = direclist2pl, pths = pth2)
# scoring using chain equating coefficients and the true score equating method for score 12
score(chainec1, path = "test1.test2.test3.test4", scores = 12)
score(chainec2, path = "test1.test5.test4", scores = 12)
# scoring using chain equating coefficients, true score equating and external anchors
score(chainec1, path = "test1.test2.test3.test4", anchor.type = "external")
# scoring using chain equating coefficients, observed score equating and external anchors
score(chainec1, path = "test1.test2.test3.test4", method = "OSE", anchor.type = "external")


# create a list of objects of class ceqc
ecall <- c(chainec1, chainec2)
# compute bisector and weighted bisector coefficients
allec <- bisectorec(ecall = ecall, weighted = TRUE, unweighted = TRUE)
# scoring using weighted bisector equating coefficients and the true score 
# equating method for score 12
score(allec, bistype = "weighted", link = "test1.test4", scores = 12)
}

