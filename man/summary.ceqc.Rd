\name{summary.ceqc}
\alias{summary.ceqc}
\alias{print.summary.ceqc}
\title{
	Summarizing Estimated Chain Equating Coefficients
}
\description{
	\code{summary} method for class \code{ceqc}.
}
\usage{
\method{summary}{ceqc}(object, ...)

\method{print}{summary.ceqc}(x, ...)
}
\arguments{
	\item{object}{
		an object of the class \code{ceqc} returned by function
		\code{\link{chainec}}.
	}
	\item{x}{
		an object of class \code{summary.ceqc}, a result of a call to \code{summary.ceqc}.
	}
	\item{\dots}{
		further arguments passed to or from other methods.
	}
}
\value{
	The function \code{summary.ceqc} returns a list with components
	\item{forms}{names of equated forms.}
	\item{method}{the equating method used.}
	\item{coefficients}{a \eqn{2 \times 2}{2X2} matrix with columns for the estimated 
	coefficients A and B and standard errors. }
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link{chainec}}
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
summary(chainec3$test1.test5.test4)	
}


