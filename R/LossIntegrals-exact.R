#' Calculate an exact Loss Integral
#' 
#' \code{LossIntegralExact} calculates an exact Loss Integral
#' using numerical integration. The iterative nature of numerical
#' integration can make it perform relatively slowly. 
#' \code{\link{LossIntegralFast}} interpolates results from precalculated
#' grids which is faster but with a (very) small loss in accuracy.
#' 
#' @param mn the mean of the normal distribution
#' @param sd the standard deviation of the normal distribution
#' @param p the probability value above which you wish to calculate the loss integral
#' @param dist the assumed dstribution of demand, 'gamma' (the default), 'norm' or 'lnorm'

LossIntegralExact <- function(p, mn, sd, dist = 'gamma'){
  # push into a data-frame as an easy way of making sure all vectors are the
  # same length
  df <- data.frame(p = p, mn = mn, sd=sd)
  df$res <- 0
  
  # define additional parameters for gamma and lognormal distributions
  if (dist == 'gamma') {
    df$alpha <- (mn/sd)^2
    df$beta <- mn/(sd^2) 
  }
  
  if (dist == 'lnorm') {
    # calculate mean and sd of the logged variable
    df$mnlog<- log(mn^2/(sd^2 + mn^2)^0.5)
    df$sdlog <- (log(1 + ((sd/mn)^2)))^0.5    
  }
  
  # integrate only works with single valued vectors
  # so we'll wrap it in a for{} loop to get what we need
  for (i in 1:nrow(df)) {
    p <- df$p[i]
    
    # for the chosen distribution(norm, lnorm or gamma)
    # find the re-order point at p and define the loss integral
    if (dist == 'norm') {
      ROP <- qnorm(p,df$mn[i],df$sd[i])       
      ff <- function(x){(x-ROP)*dnorm(x, df$mn[i] ,df$sd[i])}
    }    
    else if (dist == 'lnorm') {
      ROP <- qlnorm(p,df$mnlog[i],df$sdlog[i])
      ff <- function(x){(x-ROP)*dlnorm(x, df$mnlog[i] ,df$sdlog[i])}      
    }
    else if (dist == 'gamma') {
      ROP <- qgamma(p,df$alpha[i],df$beta[i])
      ff <- function(x){(x-ROP)*dgamma(x, df$alpha[i] ,df$beta[i])}          
    }
    
    # numeric integration for ff, above ROP to find the E(lost sales)
    try(df$res[i] <- integrate(ff,lower=ROP, upper=Inf)$value,silent=T)
    
  } # next i
  
  # calculate and return the E(fill-rate) = 1 - [% lost sales=
  return(1 - df$res/df$mn)
}


#' Calculate an exact Inverse of the Loss Integral 
#' 
#' \code{InvLossIntegralExact} calculates an exact inverse of the Loss
#' Integral for a chosen distribution using a root-finder in conjunction with
#' \code{LossIntegralExact}. The iterative nature of the calculation can
#' make it perform relatively slowly. \code{\link{InvLossIntegral}}
#' interpolates results from a precalculated grid which is faster with a (very)
#' small loss in accuracy.
#' 
#' @param fillRate, the percentage of demand met from inventory (the result of a loss integral)
#' @param mn the mean of the normal distribution
#' @param sd the standard deviation of the normal distribution
#' @param dist the assumed dstribution of demand, 'gamma' (the default), 'norm' or 'lnorm'
#' 
#' 
InvLossIntegralExact <- function(fillRate, mn, sd, dist) {  
  #   push into a data-frame as an easy way of making sure all vectors are the
  #   same length
  df <- data.frame(fillRate = fillRate, mn = mn, sd=sd)
  df$p <- 0
  
  # iterate through the root-finder to get 'p'
  for (i in 1:nrow(df)) {
    fr  <- df$fillRate[i]
    cov <- df$sd[i]/df$mn[i]
    minP <- .000000001
    try(df$p[i] <- uniroot(f=function(p,cov,fr){LossIntegralExact(p, 1,cov, dist)-fr},interval=c(minP, 1-minP), cov , fr, tol=.0000000001)$root, silent = T)
  }
  
  # return p
  return(df$p)
}