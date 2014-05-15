\name{summary.eqclist}
\alias{summary.eqclist}
\alias{print.summary.eqclist}
\title{
	Summarizing a list of estimated equating coefficients
}
\description{
	\code{summary} method for class \code{eqclist}.
}
\usage{
\method{summary}{eqclist}(object, ...)

\method{print}{summary.eqclist}(x, ...)
}
\arguments{
	\item{object}{
		an object of the class \code{eqclist} returned by functions
		\code{\link{alldirec}}.
	}
	\item{x}{
		an object of class \code{summary.eqclist}, a result of a call to summary.eqclist.
	}
	\item{\dots}{
		further arguments passed to or from other methods.
	}
}
\value{
	The function \code{summary.eqclist} returns a list containing the output
	of function \code{\link{summary.eqc}} for each equating coefficient
	given in \code{object}.
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link{alldirec}}, \code{\link{summary.eqc}}
}
\examples{
# Rasch model
# direct equating coefficients using the "mean-mean" method
data(estrasch)
test<-paste("test", 1:5, sep = "")
modrasch <- modIRT(coef = estrasch$coef, var = estrasch$var, names = test, 
	display = FALSE)
direclistrasch <- alldirec(mods = modrasch, method = "mean-mean", all = TRUE)
summary(direclistrasch)
}
