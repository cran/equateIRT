\name{summary.eqc}
\alias{summary.eqc}
\alias{print.summary.eqc}
\title{
	Summarizing estimated equating coefficients
}
\description{
	\code{summary} method for class \code{eqc}.
}
\usage{
\method{summary}{eqc}(object, ...)

\method{print}{summary.eqc}(x, ...)
}
\arguments{
	\item{object}{
		an object of the class \code{eqc} returned by functions 
		\code{\link{direc}} or \code{\link{chainec}}.
	}
	\item{x}{
		an object of class \code{summary.eqc}, a result of a call to summary.eqc.
	}
	\item{\dots}{
		further arguments passed to or from other methods.
	}
}
\value{
	The function \code{summary.eqc} returns a list with components
	\item{forms}{names of equated forms.}
	\item{method}{the equating method used.}
	\item{coefficients}{a \eqn{2 \times 2}{2X2} matrix with columns for the estimated 
	coefficients A and B and standard errors. }
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link{chainec}}, \code{\link{direc}}
}
\examples{
# three-parameter logistic model
# direct equating coefficients between forms 1 and 2 using the Stocking-Lord method
data(est3pl)
test <- paste("test", 1:5, sep = "")
mod3pl <- modIRT(coef = est3pl$coef, var = est3pl$var, names = test, display = FALSE)
l12 <- direc(mod1 = mod3pl[1], mod2 = mod3pl[2], method = "Stocking-Lord")
summary(l12)

# two-parameter logistic model
# direct equating coefficients using the "Haebara" method
data(est2pl)
test<-paste("test", 1:5, sep = "")
mod2pl <- modIRT(coef = est2pl$coef, var = est2pl$var, names = test, display = FALSE)
direclist2pl <- alldirec(mods = mod2pl, method = "Haebara")
summary(direclist2pl$test1.test2)
# compute chain equating coefficients of a given path
pth <- c(paste("test", c(1,5,4), sep = ""))
pth <- data.frame(t(pth), stringsAsFactors = FALSE)
chainec154 <- chainec(r = 3, direclist = direclist2pl, pths = pth)
summary(chainec154)	
}
