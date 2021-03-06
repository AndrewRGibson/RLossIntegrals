% Generated by roxygen2 (4.0.1): do not edit by hand
\name{InvLossIntegralExact}
\alias{InvLossIntegralExact}
\title{Calculate an exact Inverse of the Loss Integral}
\usage{
InvLossIntegralExact(fillRate, mn, sd, dist)
}
\arguments{
\item{fillRate,}{the percentage of demand met from inventory (the result of a loss integral)}

\item{mn}{the mean of the normal distribution}

\item{sd}{the standard deviation of the normal distribution}

\item{dist}{the assumed dstribution of demand, 'gamma' (the default), 'norm' or 'lnorm'}
}
\description{
\code{InvLossIntegralExact} calculates an exact inverse of the Loss
Integral for a chosen distribution using a root-finder in conjunction with
\code{LossIntegralExact}. The iterative nature of the calculation can
make it perform relatively slowly. \code{\link{InvLossIntegral}}
interpolates results from a precalculated grid which is faster with a (very)
small loss in accuracy.
}

