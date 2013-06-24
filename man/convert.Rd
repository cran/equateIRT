\name{convert}
\alias{convert}
\title{
	Item parameters conversion
}
\description{
	Converts item parameters using equating coefficients.
}
\usage{
convert(A, B, coef)
}
\arguments{
	\item{A}{
		equating coefficient A.
	}
	\item{B}{
		equating coefficient B.
	}
	\item{coef}{
		vector of item parameters return by function \code{\link{modIRT}}.
	}
}
\details{
	Difficulty parameters \eqn{b} are converted using transformation \eqn{b \cdot A+B}{b*A+B};
	discrimination parameters \eqn{a} are converted using transformation \eqn{a/A};
	guessing parameters c are not transformed.
}
\value{
	vector of coefficients transformed.
}
\references{
	Kolen, M.J. and Brennan, R.L. (2004). \emph{Test equating, scaling, 
	and linking: methods and practices}, 2nd ed., New York: Springer
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link{modIRT}}
}
\examples{
# conversion using direct coefficients
# three-parameter logistic model
# direct equating coefficients between forms 1 and 2 using the Stocking-Lord method
data(est3pl)
test <- paste("test", 1:5, sep = "")
mod3pl <- modIRT(coef = est3pl$coef, var = est3pl$var, names = test, ltparam = TRUE, lparam = TRUE)
l12 <- direc(mod1 = mod3pl[1], mod2 = mod3pl[2], method = "Stocking-Lord", itmp = 3)
convert(A = l12$A, B = l12$B, coef = mod3pl$test1$coef)
# the conversion is given also in
l12$tab

# conversion using bisector coefficients
# two-parameter logistic model
# direct equating coefficients using the "Haebara" method
data(est2pl)
test <- paste("test", 1:5, sep = "")
mod2pl <- modIRT(coef = est2pl$coef, var = est2pl$var, names = test, ltparam = TRUE)
direclist2pl <- alldirec(mods = mod2pl, method = "Haebara", itmp = 2)
# compute chain equating coefficients for path 1,2,3,4
pth1 <- paste("test", 1:4, sep = "")
pth1 <- data.frame(t(pth1), stringsAsFactors = FALSE)
chainec1 <- chainec(r = 4, direclist = direclist2pl, pths = pth1)
# compute chain equating coefficients for path 1,5,4
pth2 <- c(paste("test", c(1,5,4), sep = ""))
pth2 <- data.frame(t(pth2), stringsAsFactors = FALSE)
chainec2 <- chainec(r = 3, direclist = direclist2pl, pths = pth2)
# compute chain equating coefficients for path 1,2,3,4,5
pth3 <- paste("test", 1:5, sep = "")
pth3 <- data.frame(t(pth3), stringsAsFactors = FALSE)
chainec3 <- chainec(r = 5, direclist = direclist2pl, pths = pth3)
# create a list of objects of class \code{eqc} 
ecall <- c(chainec1, chainec2, chainec3, direclist2pl["test1.test5"])
# compute bisector and weighted bisector coefficients
allec <- bisectorec(ecall = ecall, mod = mod2pl, weighted = TRUE, unweighted = TRUE)
summary(allec)
Ab<-allec$coef[allec$coef$link=="test1.test4" & allec$coef$path=="bisector", ]$A
Bb<-allec$coef[allec$coef$link=="test1.test4" & allec$coef$path=="bisector", ]$B
convert(A = Ab, B = Bb, coef = mod2pl$test1$coef)
}

