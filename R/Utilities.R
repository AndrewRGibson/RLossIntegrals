#' cover function for \code{\link{qgamma}} using mean and standard deviation
#' rather than shape and rate parameters
#' 
#' \code{qgamma2} calculates rate and shape paramaters from the provided mean
#' (mn) and standard deviation (sd) then calls the base R qgamma function
#' 
#' @param p, probability (0,1)
#' @param mn the mean of the gamma distribution
#' @param sd the standard deviation of the gamma distribution
#' 
#' @export
#' @family utilities
#' 
#' @examples
#' # Use qgamma2 directly ....
#' qgamma2(p = 0.8, mn = 100, sd = 60)
#' 
#' # as an alternative to this
#' p <- 0.8; mn <-100 ; sd <- 60
#' shape <- (mn/sd)^2
#' rate <- mn/(sd^2) 
#' qgamma(p, shape, rate)

qgamma2 <- function(p, mn, sd){
    shape <- (mn/sd)^2
    rate <- mn/(sd^2) 
    qgamma(p,shape, rate)
  }
  
#' cover function for \code{\link{qlnorm}} using mean and standard deviation
#' of the resulting distribution rather than mean and sd of the logged variable
#' 
#' 
#' @param p, probability (0,1)
#' @param mn the mean of the gamma distribution
#' @param sd the standard deviation of the gamma distribution
#' 
#' @export
#' @family utilities
#' 
#' @examples
#' # Use qlnorm2 directly ....
#' qlnorm2(p = 0.8, mn = 100, sd = 60)
#' 
#' # as an alternative to this
#' p <- 0.8; mn <-100 ; sd <- 60
#' mnlog<- log(mn^2/(sd^2 + mn^2)^0.5)
#' sdlog <- (log(1 + ((sd/mn)^2)))^0.5   
#' qlnorm(p,mnlog,sdlog)
qlnorm2 <- function(p,mn,sd) {
  mnlog<- log(mn^2/(sd^2 + mn^2)^0.5)
  sdlog <- (log(1 + ((sd/mn)^2)))^0.5   
  qlnorm(p,mnlog,sdlog)
}
