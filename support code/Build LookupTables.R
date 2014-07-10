# Code builds 2 dimensional lookup table to get fast conversion betwen P and fill-rate 
# for a given CoV

# It will attempt to save the result to an Rdata file in the \data\ subdirectory of the 
# working directory.  If wd is set to the working directory of the  package (in dev)
# the package can then pick up the new lossIntegralTables value for use


# Utility function takes an x vector, y vector and calculates Z for all combinations of X an Y
# using func, returns a list (X, Y, Z) in a form than that akima interpolation (bicubic) can use
CalculateFunctionTable <- function(func, X, Y) {
  # calculate the function for all X, Y
  inp <- expand.grid(X = X,Y = Y)
  
  # build the result grid
  res <- list()
  res$X<- X
  res$Y<- Y
  res$Z <- matrix(func(inp$X,inp$Y), length(X))  
  
  # find names for attributes (where possible)
  attr(res,"calculation") <- deparse(substitute(func))
  attr(res,"X") <- deparse(substitute(X))
  attr(res, "Y") <- deparse(substitute(Y))
  
  return(res)
}


# Build tables with X= CoV and Y = fillRate and Z = (what we want to know, p )
len <- 50
CoV = seq(.5, 4, length.out = len)
fillRate <- seq(.5, .9999, length.out = len)
P <- seq(0.5, .9999, length.out = len)

# Create a stub for the output
lossIntegralTables <- list()

# Normal Loss Integral
lossIntegralTables$Norm <- CalculateFunctionTable(
  function(CoV, P){LossIntegralExact(P, 1,CoV, 'norm')}
  , CoV
  , P)

# Gamma Loss Integral
lossIntegralTables$Gamma <- CalculateFunctionTable(
  function(CoV, P){LossIntegralExact(P, 1,CoV,'gamma')}
  , CoV
  , P)

#LogNormal Loss Integral
lossIntegralTables$LogNorm <- CalculateFunctionTable(
  function(CoV, P){LossIntegralExact(P, 1,CoV,'lnorm')}
  , CoV
  , P)

# Normal Inverse 
lossIntegralTables$NormInv <- CalculateFunctionTable(
  function(CoV, fillRate){InvLossIntegralExact(fillRate, 1,CoV,'norm')}
  , CoV
  , fillRate)

# Gamma Inverse
lossIntegralTables$GammaInv <- CalculateFunctionTable(
  function(CoV, fillRate){InvLossIntegralExact(fillRate, 1,CoV,'gamma')}
  , CoV
  , fillRate)

# LogNormal Inverse 
lossIntegralTables$LogNormInv <- CalculateFunctionTable(
  function(CoV, fillRate){InvLossIntegralExact(fillRate, 1,CoV,'lnorm')}
  , CoV
  , fillRate)

# save the tables for use in the package
# NOTE - this assumes the working directory is set to that of the  package
save(list='lossIntegralTables', file='./data/lossIntegralTables.RData')

