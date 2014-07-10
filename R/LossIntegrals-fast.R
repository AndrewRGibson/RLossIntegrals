## libraries ----
library(akima)


## fast interpolation functions ---- 

#' Calculate a FAST Loss Integral 
#' 
#' \code{LossIntegral} calculates a fast approximation of the Loss Integral
#' for the chosen distribution using interpolation from a grid 
#' 
#' @param mn the mean of the distribution
#' @param sd the standard deviation of the distribution
#' @param p the probability value for which you wish to calculate the loss integral
#' @param dist the assumed dstribution of demand, 'gamma' (the default), 'norm' or 'lnorm'

LossIntegralFast <- function(p, mn, sd, dist = 'gamma') {
  # push into a data-frame as an easy way of making sure all vectors
  # are the same length
  df <- data.frame(p = p, cov= sd/mn)
  # reference the right precalculated tables
  t <- switch(dist,
              norm = lossIntegralTables$Norm,
              gamma = lossIntegralTables$Gamma,
              lnorm = lossIntegralTables$LogNorm)
  
  # bicubic interpolation
  res <- bicubic(x=t$X, y = t$Y, z = t$Z, x0=df$cov, y0=df$p)
  # return the Z component - fill-rate
  return(res$z)
}


#' Calculate a FAST Inverse for a Loss Integral 
#' 
#' \code{InvLossIntegral} calculates a fast approximation of the Inverse Loss Integral
#' for the chosen distribution using interpolation from a grid 
#' 
#' @param mn the mean of the distribution
#' @param sd the standard deviation of the distribution
#' @param fillRate the probability value for which you wish to calculate the loss integral
#' @param dist the assumed dstribution of demand, 'gamma' (the default), 'norm' or 'lnorm'


InvLossIntegralFast <- function(fillRate, mn, sd, dist = 'gamma') {
  # push into a data-frame as an easy way of making sure all vectors
  # are the same length
  df <- data.frame(fillrate = fillRate, cov= sd/mn)
  # reference the right precalculated tables
  t <- switch(dist,
              norm = lossIntegralTables$NormInv,
              gamma = lossIntegralTables$GammaInv,
              lnorm = lossIntegralTables$LogNormInv)
  
  # bicubic interpolation
  res <- bicubic(x=t$X, y = t$Y, z = t$Z, x0=df$cov, y0=df$fillrate)
  # return the Z component - P
  return(res$z)
}



## user functions  ----

#' Calculate Loss Integral 
#' 
#' \code{LossIntegral} calculatesthe Loss Integral
#' for the chosen distribution using either interpolation from a grid (the default method,
#' faster but with a small loss of accuracy or numerical integration (substantially 
#' slower but slightly more accurate)
#' 
#' @param mn the mean of the distribution
#' @param sd the standard deviation of the distribution
#' @param p the probability value for which you wish to calculate the loss integral
#' @param dist the assumed dstribution of demand, 'gamma' (the default), 'norm' or 'lnorm'
#' @param fast TRUE for the faster interpolation method
#' @export
#' @family LossIntegrals
#' @examples

#' # calculate the fill rate associated with a 80% probability of 
#' # no stock outs given sales ~N(100,40)
#' LossIntegral(0.80, 100, 40, 'norm') 
#' # and assuming a gamma distribution (the default)
#' LossIntegral(0.80, 100, 40, 'gamma')  
#' # and assuming a lognormal distribution
#' LossIntegral(0.80, 100, 40, 'lnorm')  
#' 
#' # sensitivity analysis of safety stock against fill-rate for each distribution
#' xFillRate <- seq(.90, .999, length.out = 100)
#' mn <- 100
#' sd <- 80
#' plot(xFillRate, qnorm(InvLossIntegral(xFillRate, mn, sd, 'norm'),mn, sd)-mn ,  ty='l', col = 'red', xlab = 'Fill-Rate', ylab = "Safety Stock", ylim = c(0,500))
#' lines(xFillRate, qgamma2(InvLossIntegral(xFillRate, mn,sd, 'gamma'),mn, sd)-mn ,  ty='l', col = 'blue')
#' lines(xFillRate, qlnorm2(InvLossIntegral(xFillRate, mn,sd,'lnorm'),mn, sd)-mn ,  ty='l', col = 'green')
LossIntegral <- function(p, mn, sd, dist = 'gamma', fast = T) {
  if (fast == F) {
    LossIntegralExact(p,mn,sd,dist)
  } else {
    LossIntegralFast(p,mn,sd,dist)
  }
}


#' Calculate the Inverse for a Loss Integral 
#' 
#' \code{InvLossIntegral} calculates the Inverse of a Loss Integral
#' for the chosen distribution using either interpolation from a grid (the default
#' method, faster but with a small loss of accuracy) or numerical integration (substantially 
#' slower but slightly more accurate)
#' 
#' @param mn the mean of the distribution
#' @param sd the standard deviation of the distribution
#' @param fillRate the probability value for which you wish to calculate the loss integral
#' @param dist the assumed dstribution of demand, 'gamma' (the default), 'norm' or 'lnorm'
#' @param fast TRUE for the faster interpolation method
#' @export
#' @family LossIntegrals
#' @examples
#' # calculate the probability of no stock-outs associated with a fill-rate of
#' 97% # given sales ~ N(100,20)
#' InvLossIntegralExact(0.97, 100, 20, 'norm') 
#' # and assuming a gamma distribution (the default)
#' InvLossIntegralExact(0.97, 100, 20, 'gamma') 
#' # and assuming a lognormal distribution
#' InvLossIntegralExact(0.97, 100, 20, 'lnorm') 
#' 
InvLossIntegral <- function(fillRate, mn, sd, dist = 'gamma', fast = T) {
  if (fast == F) {
    InvLossIntegralExact(fillRate,mn,sd,dist)
  } else {
    InvLossIntegralFast(fillRate,mn,sd,dist)
  }
}
