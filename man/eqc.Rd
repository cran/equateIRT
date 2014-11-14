\name{eqc}
\alias{eqc}
\alias{eqc.eqc}
\alias{eqc.eqclist}
\alias{eqc.ceqc}
\alias{eqc.ceqclist}
\alias{eqc.meqc}
\title{
	Extract Equating Coefficients
}
\description{
	\code{eqc} is a generic function which extracts the equating coefficients.
}
\usage{
eqc(x, \dots)

\method{eqc}{eqc}(x, \dots)

\method{eqc}{eqclist}(x, link = NULL, \dots)

\method{eqc}{ceqc}(x, \dots)

\method{eqc}{ceqclist}(x, link = NULL, path = NULL, \dots)

\method{eqc}{meqc}(x, link = NULL, path = NULL, \dots)
}
\arguments{
	\item{x}{object of the class \code{eqc} returned by function \code{\link{direc}} 
		or of the class \code{eqclist} returned by function \code{\link{alldirec}} 
		or of the class \code{ceqc} and \code{ceqclist} returned by function \code{\link{chainec}} 
		or of the class \code{meqc} returned by function \code{\link{bisectorec}}.}
	\item{link}{a character string with the names of the two forms being linked separated by a dot (e.g. "\code{test1.test3}").}
	\item{path}{a character string with the names of the forms that constitute the path separated by a dot (e.g. "\code{test1.test2.test3}").}
	\item{\dots}{further arguments passed to or from other methods.}
}
\value{
	A data frame containing the equating coefficients for every
	link and path.
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
l12 <- direc(mod1 = mod2pl[1], mod2 = mod2pl[2], method = "Haebara")
# all direct equating coefficients using the Haebara method
direclist2pl <- alldirec(mods = mod2pl, method = "Haebara")
# compute all chain equating coefficients of length 3
chainec3 <- chainec(r = 3, direclist = direclist2pl)
# compute chain equating coefficients for path 1,2,3,4
pth1 <- paste("test", 1:4, sep = "")
pth1 <- data.frame(t(pth1), stringsAsFactors = FALSE)
chainec1 <- chainec(direclist = direclist2pl, pths = pth1)
# compute chain equating coefficients for path 1,5,4
pth2 <- c(paste("test", c(1,5,4), sep = ""))
pth2 <- data.frame(t(pth2), stringsAsFactors = FALSE)
chainec2 <- chainec(direclist = direclist2pl, pths = pth2)
# create a list of objects of class ceqc
ecall <- c(chainec1, chainec2)
# compute bisector and weighted bisector coefficients
allec <- bisectorec(ecall = ecall, mods = mod2pl, weighted = TRUE, unweighted = TRUE)

# extract equating coefficients
eqc(l12)
eqc(direclist2pl)
eqc(direclist2pl,link="test1.test2")
eqc(chainec3)
eqc(chainec3,link="test1.test3")
eqc(allec)
eqc(allec,path="bisector")
}


