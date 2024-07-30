\name{linkp}
\alias{linkp}

\title{
	Linkage Plan
}
\description{
	Calculates the number of common items between a list of forms.
}
\usage{
linkp(mods = NULL, coef = NULL)
}
\arguments{
	\item{mods}{
  	an object of the class \code{modIRT} containing 
  	item parameter coefficients and their covariance matrix
  	of the forms to be equated.
	}
	\item{coef}{list of matrices (one for each form) containing item
		parameter estimates. The names of the rows of each matrix should 
		be the names of the items.
	}
}
\value{
	A matrix whose elements indicate the number of common items between
	the forms. On the diagonal of the matrix there are the number of items
	of each form.
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link{est2pl}}, \code{\link{est3pl}}, \code{\link{estrasch}}
}
\examples{

\dontrun{
library("mirt")
data("data2pl")

m1 <- mirt(data2pl[[1]], SE = TRUE)
m2 <- mirt(data2pl[[2]], SE = TRUE)
m3 <- mirt(data2pl[[3]], SE = TRUE)
m4 <- mirt(data2pl[[4]], SE = TRUE)
m5 <- mirt(data2pl[[5]], SE = TRUE)

mlist<- list(m1,m2,m3,m4,m5)
test <- paste("test", 1:5, sep = "")
mod2pl <- modIRT(est.mods = mlist, names = test, display = FALSE)
linkp(mod2pl)
}

# ===========================================================================
# the following code uses item parameter estimates previously obtained

data(est3pl)
linkp(coef = est3pl$coef)

}

