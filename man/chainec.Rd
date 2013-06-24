\name{chainec}
\alias{chainec}
\title{
	Chain equating coefficients
}
\description{
	Calculates chain (indirect) equating coefficients and standard errors using 
	IRT methods.
}
\usage{
chainec(r, direclist, f1 = NULL, f2 = NULL, pths = NULL)
}
\arguments{
	  \item{r}{
		length of the chain, that is the number of forms used for equating
		including extremes. It should be at least 3.
	}
	  \item{direclist}{
		an object of the class \code{eqclist} return by function \code{\link{alldirec}}
		containing direct equating coefficients	between pairs of forms.
	}
	  \item{f1}{
		the name of the first form of the chain.
	}
	  \item{f2}{
		the name of the last form of the chain.
	}
	  \item{pths}{
		data frame containing the path used for equating. The number of columns 
		is equal to r and the number of rows is equal to the equatings that
		have to be performed.
		If \code{NULL} all the chain equating coefficients of length \code{r} 
		will be computed.
	}
}
\details{
	Equating coefficients perform the conversion from the scale of the 
	first form to the scale of the last form of the path.
}
\value{
	An object of class \code{eqclist} consisting in a list with length 
	equal to the number of chain equating coefficients computed. 
	Each component of the list is an object of class \code{eqc} with components
	\item{tab1}{item parameters of the first form.}
	\item{tab2}{item parameters of the last form.}
	\item{tab}{data frame containing item parameters of the first form (value1),
		the last form (value2), and the first form converted in the last
		form (value12).}
	\item{varAll}{covariance matrix of item parameters of all forms used 
		in the chain.}
	\item{partial}{partial derivatives of equating coefficients A and B
		with respect to item parameters.}
	\item{A}{equating coefficient A.}
	\item{B}{equating coefficient B.}
	\item{varAB}{covariance matrix of the equating coefficients.}
	\item{commonitem}{list of length r-1 containing the names of common item 
		parameters between adjacent forms.}
	\item{ni}{vector containing number of common items between pairs of adjacent forms.}
	\item{forms}{names of equated forms.}
	\item{method}{the equating method used.}
	\item{itmp}{number of item parameters of the IRT model.}
}
\references{
	Battauz, M. (2013). IRT Test Equating in Complex Linkage Plans. 
	\emph{Psychometrika}, \bold{78}, 464--480.
	
	Kolen, M.J. and Brennan, R.L. (2004). \emph{Test equating, scaling, 
	and linking: methods and practices}, 2nd ed., New York: Springer
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link{alldirec}}
}
\examples{
# three-parameter logistic model
# direct equating coefficients using the "Stocking-Lord" method
data(est3pl)
test<-paste("test", 1:5, sep = "")
mod3pl <- modIRT(coef = est3pl$coef, var = est3pl$var, names = test, ltparam = TRUE, lparam = TRUE)
direclist3pl <- alldirec(mods = mod3pl, method = "Stocking-Lord", itmp = 3)
# compute all chain equating coefficients of length 4
chainec4 <- chainec(r = 4, direclist = direclist3pl)
summary(chainec4)
summary(chainec4$test1.test2.test3.test4)
# compute all chain equating coefficients of length 4 
# where the first form is test1
chainec4.1 <- chainec(r = 4, direclist = direclist3pl, f1 = "test1")
summary(chainec4.1)
# compute all chain equating coefficients of length 4 
# where the first form is test1 and the last form is test4
chainec4.14 <- chainec(r = 4, direclist = direclist3pl, f1 = "test1", f2 = "test4")
summary(chainec4.14)


# two-parameter logistic model
# direct equating coefficients using the "Haebara" method
data(est2pl)
test<-paste("test", 1:5, sep = "")
mod2pl <- modIRT(coef = est2pl$coef, var = est2pl$var, names = test, ltparam = TRUE)
direclist2pl <- alldirec(mods = mod2pl, method = "Haebara", itmp = 2)
# compute chain equating coefficients of a given path
pth <- c(paste("test", c(1,5,4), sep = ""))
pth <- data.frame(t(pth), stringsAsFactors = FALSE)
chainec154 <- chainec(r = 3, direclist = direclist2pl, pths = pth)
summary(chainec154)


# Rasch model
# direct equating coefficients using the "mean-mean" method
data(estrasch)
test<-paste("test", 1:5, sep = "")
modrasch <- modIRT(coef = estrasch$coef, var = estrasch$var, names = test, ltparam = TRUE)
direclistrasch <- alldirec(mods = modrasch, method = "mean-mean", all = TRUE, itmp = 1)
# compute chain equating coefficients of two given paths
pth1 <- paste("test", 1:3, sep = "")
pth1 <- data.frame(t(pth1), stringsAsFactors = FALSE)
pth2 <- c(paste("test", c(1,5,4), sep = ""))
pth2 <- data.frame(t(pth2), stringsAsFactors = FALSE)
pths <- rbind(pth1, pth2)
chainec1 <- chainec(r = 3, direclist = direclistrasch, pths = pths)
summary(chainec1)
}
