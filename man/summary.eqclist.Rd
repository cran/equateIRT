\name{summary.eqclist}
\alias{summary.eqclist}
\alias{print.summary.eqclist}
\title{
	Summarizing a List of Estimated Direct Equating Coefficients
}
\description{
	\code{summary} method for class \code{eqclist}.
}
\usage{
\method{summary}{eqclist}(object, link, ...)

\method{print}{summary.eqclist}(x, ...)
}
\arguments{
	\item{object}{
		an object of the class \code{eqclist} returned by function
		\code{\link{alldirec}}.
	}
	\item{link}{
		a vector of character strings with the names of the two 
		forms being linked separated by a dot 
		(e.g. "\code{test1.test2}").
	}
	\item{x}{
		an object of class \code{summary.eqclist}, a result of a call to \code{summary.eqclist}.
	}
	\item{\dots}{
		further arguments passed to or from other methods.
	}
}
\value{
	The function \code{summary.eqclist} returns a list containing the output
	of function \code{\link{summary.eqc}} for each link
	contained in \code{object}.
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link{alldirec}}, \code{\link{modIRT}}, \code{\link{summary.eqc}}
}
\examples{
# the following code uses item parameter estimates previously obtained
# for examples that start from item response data see function modIRT

# Rasch model
# direct equating coefficients using the "mean-mean" method
data(estrasch)
test<-paste("test", 1:5, sep = "")
modrasch <- modIRT(coef = estrasch$coef, var = estrasch$var, names = test, 
	display = FALSE)
direclistrasch <- alldirec(mods = modrasch, method = "mean-mean", all = TRUE)
summary(direclistrasch)
summary(direclistrasch, link = "test1.test2")
summary(direclistrasch, link = c("test1.test2", "test1.test5"))
}
