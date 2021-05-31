\name{sd.test}
\alias{sd.test}
\alias{print.sdtest}

\title{
	Test on Chain Equating Coefficients
}
\description{
	Performs a statistical test to verify if the chain equating coeffients that
	link the same two forms are equal.
}
\usage{
sd.test(ecall)
}

\arguments{
	\item{ecall}{
	list of objects of class \code{eqc} or \code{ceqc} 
	returned by functions
	\code{\link{direc}} and \code{\link{chainec}}.
	}
}
\value{
	An object of class \code{idtest} with components
	\item{link}{}
	\item{paths}{}
	\item{AB}{}
	\item{statistic}{}
	\item{df}{}
	\item{p.value}{}
}
\author{
	Michela Battauz
}
\examples{
# three-parameter logistic model
# direct equating coefficients using the "Haebara" method
data(est3pl)
test <- paste("test", 1:5, sep = "")
mod3pl <- modIRT(coef = est3pl$coef, var = est3pl$var, names = test, display = FALSE)
direclist3pl <- alldirec(mods = mod3pl, method = "Haebara")
# compute chain equating coefficients for path 1,2,3,4,5
pth3 <- paste("test", 1:5, sep = "")
chainec3 <- chainec(direclist = direclist3pl, pths = pth3)
# create a list of objects of class eqc or ceqc
ecall <- c(chainec3, direclist3pl["test1.test5"])
summary(chainec3)
summary(direclist3pl$test1.test5)
# test
sd.test(ecall)
}