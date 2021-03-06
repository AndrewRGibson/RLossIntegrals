% Generated by roxygen2 (4.0.1): do not edit by hand
\name{qlnorm2}
\alias{qlnorm2}
\title{cover function for \code{\link{qlnorm}} using mean and standard deviation
of the resulting distribution rather than mean and sd of the logged variable}
\usage{
qlnorm2(p, mn, sd)
}
\arguments{
\item{p,}{probability (0,1)}

\item{mn}{the mean of the gamma distribution}

\item{sd}{the standard deviation of the gamma distribution}
}
\description{
cover function for \code{\link{qlnorm}} using mean and standard deviation
of the resulting distribution rather than mean and sd of the logged variable
}
\examples{
# Use qlnorm2 directly ....
qlnorm2(p = 0.8, mn = 100, sd = 60)

# as an alternative to this
p <- 0.8; mn <-100 ; sd <- 60
mnlog<- log(mn^2/(sd^2 + mn^2)^0.5)
sdlog <- (log(1 + ((sd/mn)^2)))^0.5
qlnorm(p,mnlog,sdlog)
}
\seealso{
Other utilities: \code{\link{qgamma2}}
}

