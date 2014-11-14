\name{summary.ceqclist}
\alias{summary.ceqclist}
\alias{print.summary.ceqclist}
\title{
	Summarizing a List of Estimated Chain Equating Coefficients
}
\description{
	\code{summary} method for class \code{ceqclist}.
}
\usage{
\method{summary}{ceqclist}(object, path, ...)

\method{print}{summary.ceqclist}(x, ...)
}
\arguments{
	\item{object}{
		an object of the class \code{ceqclist} returned by function
		\code{\link{chainec}}.
	}
	\item{path}{
		a vector of character strings with the names of the 
		forms that constitute 
		the path separated by a dot (e.g. "\code{test1.test2.test3}").
	}
	\item{x}{
		an object of class \code{summary.ceqclist}, a result of a call to \code{summary.ceqclist}.
	}
	\item{\dots}{
		further arguments passed to or from other methods.
	}
}
\value{
	The function \code{summary.ceqclist} returns a list containing the output
	of function \code{\link{summary.ceqc}} for each path
	contained in \code{object}.
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link{chainec}}, \code{\link{summary.eqc}}
}
\examples{
# two-parameter logistic model
# direct equating coefficients using the "Haebara" method
data(est2pl)
test<-paste("test", 1:5, sep = "")
mod2pl <- modIRT(coef = est2pl$coef, var = est2pl$var, names = test, display = FALSE)
direclist2pl <- alldirec(mods = mod2pl, method = "Haebara")
# compute all chain equating coefficients of length 3
chainec3 <- chainec(r = 3, direclist = direclist2pl)
summary(chainec3)	
summary(chainec3, path = "test1.test2.test3")	
summary(chainec3, path = c("test1.test2.test3","test1.test5.test4"))	
}



