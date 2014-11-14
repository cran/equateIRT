\name{summary.meqc}
\alias{summary.meqc}
\alias{print.summary.meqc}
\title{
	Summarizing Bisector Equating Coefficients
}
\description{
	\code{summary} method for class \code{meqc}.
}
\usage{
\method{summary}{meqc}(object, ...)

\method{print}{summary.meqc}(x, ...)
}
\arguments{
	\item{object}{
		an object of the class \code{meqc} returned by function
		\code{\link{bisectorec}}.
	}
	\item{x}{
		an object of class \code{summary.meqc}, a result of a call to \code{summary.meqc}.
	}
	\item{\dots}{
		further arguments passed to or from other methods.
	}
}
\value{
	The function \code{summary.meqc} returns a list with components
	\item{link}{character vector with names of equated forms.}
	\item{method}{the equating method used.}
	\item{coefficients}{list of data frames containing \code{Path},
		\code{Estimate} and \code{StdErr} of direct, chain 
		and bisector equating coefficients.
	}
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link{bisectorec}}
}
\examples{
# three-parameter logistic model
# direct equating coefficients using the "Stocking-Lord" method
data(est3pl)
test <- paste("test", 1:5, sep = "")
mod3pl <- modIRT(coef = est3pl$coef, var = est3pl$var, names = test, display = FALSE)
direclist3pl <- alldirec(mods = mod3pl, method = "Stocking-Lord")
# compute chain equating coefficients for path 1,2,3,4
pth1 <- paste("test", 1:4, sep = "")
pth1 <- data.frame(t(pth1), stringsAsFactors = FALSE)
chainec1 <- chainec(direclist = direclist3pl, pths = pth1)
# compute chain equating coefficients for path 1,5,4
pth2 <- c(paste("test", c(1,5,4), sep = ""))
pth2 <- data.frame(t(pth2), stringsAsFactors = FALSE)
chainec2 <- chainec(direclist = direclist3pl, pths = pth2)
# compute chain equating coefficients for path 1,2,3,4,5
pth3 <- paste("test", 1:5, sep = "")
pth3 <- data.frame(t(pth3), stringsAsFactors = FALSE)
chainec3 <- chainec(direclist = direclist3pl, pths = pth3)
# create a list of objects of class eqc or ceqc
ecall <- c(chainec1, chainec2, chainec3, direclist3pl["test1.test5"])
# compute bisector and weighted bisector coefficients
allec <- bisectorec(ecall = ecall, mods = mod3pl, weighted = TRUE, unweighted = TRUE)
summary(allec)
}
