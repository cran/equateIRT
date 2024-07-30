\name{itm}
\alias{itm}
\alias{itm.eqc}
\alias{itm.eqclist}
\alias{itm.ceqc}
\alias{itm.ceqclist}
\alias{itm.meqc}
\title{
	Extract Item Parameters
}
\description{
	\code{itm} is a generic function which extracts a data frame containing the
	item parameters of two forms being equated in the original scale and 
	item parameters of the first form converted in the scale of the second form.
}
\usage{
itm(x, \dots)

\method{itm}{eqc}(x, \dots)

\method{itm}{eqclist}(x, link, \dots)

\method{itm}{ceqc}(x, \dots)

\method{itm}{ceqclist}(x, path, \dots)

\method{itm}{meqc}(x, link, bistype = NULL, \dots)
}
\arguments{
	\item{x}{object of the class \code{eqc} returned by function \code{\link{direc}} 
		or of the class \code{eqclist} returned by function \code{\link{alldirec}} 
		or of the class \code{ceqc} and \code{ceqclist} returned by function \code{\link{chainec}}
		or of the class \code{meqc} return by function \code{\link{bisectorec}}.}
	\item{link}{a character string with the names of the two forms being linked separated by a dot (e.g. "\code{test1.test3}").}
	\item{path}{a character string with the names of the forms that constitute the path separated by a dot (e.g. "\code{test1.test2.test3}").}
	\item{bistype}{bisector type. It should be \code{weighted} or \code{unweighted}. 
		It could be \code{NULL} if only one type was computed by function \code{\link{bisectorec}}.}
	\item{\dots}{further arguments passed to or from other methods.}
}
\value{
	A data frame containing item names (\code{Item}), 
	item parameters of the first form (e.g. \code{test1}), 
	item parameters of the last form (e.g. \code{test3}), 
	and item parameters of the first form converted 
	in the scale of the last form (e.g. \code{test1.as.test3}).	
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link{convert}}, \code{\link{alldirec}}, \code{\link{bisectorec}}, \code{\link{chainec}}, \code{\link{direc}}, \code{\link{modIRT}}
}
\examples{
# the following code uses item parameter estimates previously obtained
# for examples that start from item response data see function modIRT

# two-parameter logistic model
data(est2pl)
test<-paste("test", 1:5, sep = "")
mod2pl <- modIRT(coef = est2pl$coef, var = est2pl$var, names = test, display = FALSE)
# direct equating coefficients between forms 1 and 2 using the Haebara method
l12 <- direc(mods = mod2pl, which = c(1,2), method = "Haebara")
# all direct equating coefficients using the Haebara method
direclist2pl <- alldirec(mods = mod2pl, method = "Haebara")
# compute all chain equating coefficients of length 3
chainec3 <- chainec(r = 3, direclist = direclist2pl)
# compute chain equating coefficients for path 1,2,3,4
pth1 <- paste("test", 1:4, sep = "")
chainec1 <- chainec(direclist = direclist2pl, pths = pth1)
# compute chain equating coefficients for path 1,5,4
pth2 <- c(paste("test", c(1,5,4), sep = ""))
chainec2 <- chainec(direclist = direclist2pl, pths = pth2)
# create a list of objects of class ceqc
ecall <- c(chainec1, chainec2)
# compute bisector and weighted bisector coefficients
allec <- bisectorec(ecall = ecall, weighted = TRUE, unweighted = TRUE)

# extract item parameters
itm(l12)
itm(direclist2pl, link = "test1.test2")
itm(chainec3, path = "test1.test2.test3")
itm(allec, bistype = "weighted")
}
