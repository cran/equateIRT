\name{id.test}
\alias{id.test}
\alias{print.idtest}

\title{
	Identity Test
}
\description{
	Performs a statistical test to verify if the chain equating coeffients from
	one form to itself are A=1 and B=0.
}
\usage{
id.test(chain)
}

\arguments{
	\item{chain}{an object of the class \code{ceqc} returned by function
		\code{\link{chainec}}.}
}
\value{
	An object of class \code{idtest} with components
	\item{path}{}
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
# compute chain equating coefficients for path 1,2,3,4,5,1
pth3 <- paste("test", c(1:5,1), sep = "")
chainec_circle <- chainec(direclist = direclist3pl, pths = pth3)
summary(chainec_circle)
# test
id.test(chainec_circle)
}
