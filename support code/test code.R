p <- 0.9
mn <- 100
sd <- 100
fast <- T

# Convert from p to fill-rate and back
fr <- LossIntegral(p, mn,sd, 'norm', fast) 
paste0('Normal Dist fill-rate: ', fr)
InvLossIntegral(fr, mn, sd, 'norm', fast )

fr <- LossIntegral(p, mn, sd, 'gamma', fast)
paste0('Gamma Dist fill-rate: ', fr)
InvLossIntegral(fr, mn, sd, 'gamma', fast)

fr <- LossIntegral(p, mn, sd, 'lnorm', fast)
paste0('LogNormal Dist fill-rate: ', fr)
InvLossIntegral(fr, mn, sd, 'lnorm', fast)

# Plot p vs fill-rates
xP <- seq(.5, .9999, length.out = 100)
plot(xP, 100* LossIntegral(xP, mn,sd,'norm') ,  ty='l', col = 'red', xlab = 'P(no shortages)', ylab = 'fill-rate')
lines(xP, 100* LossIntegral(xP, mn,sd,'gamma') , ty='l', col = 'blue')
lines(xP, 100* LossIntegral(xP, mn,sd,'lnorm') , ty='l', col = 'green')

# SS vs fill-rates
plot(LossIntegral(xP, mn,sd,'norm'), qnorm(xP, mn, sd)-mn ,  ty='l', col = 'red', xlab = 'P(no shortages)', ylab = "Safety-Stock")
lines(LossIntegral(xP, mn,sd,'gamma'),qgamma2(xP, mn, sd)-mn , ty='l', col = 'blue')
lines(LossIntegral(xP, mn,sd,'lnorm'),qlnorm2(xP, mn, sd)-mn , ty='l', col = 'green')

# control the axes better by starting with fill-rates (typically what I care about)
xFillRate <- seq(.90, .999, length.out = 100)
plot(xFillRate, qnorm(InvLossIntegral(xFillRate, mn,sd, 'norm'),mn, sd)-mn ,  ty='l', col = 'red', xlab = 'Fill-Rate', ylab = "Safety Stock", ylim = c(0,500))
lines(xFillRate, qgamma2(InvLossIntegral(xFillRate, mn,sd, 'gamma'),mn, sd)-mn ,  ty='l', col = 'blue')
lines(xFillRate, qlnorm2(InvLossIntegral(xFillRate, mn,sd,'lnorm'),mn, sd)-mn ,  ty='l', col = 'green')

# function to compare exact and estimated results
Acc <- function(nm, exact, estimate, ...){
  res <- list(name = nm)
  res$n <- length(exact)
  res$Correlation <- cor(exact, estimate)
  res$MeanErr <- mean(estimate-exact)
  res$MeanAbsErr <- mean(abs(estimate-exact))
  res$MAPE <- mean(abs(estimate-exact))/mean(exact)
  res$MaxAbsErr <- max(abs(estimate-exact))
  return(as.data.frame(res))
}

# compare exact and estimated loss integrals for estimation error
xP <- runif(n=1000, min=.5, max = .9999)
df <- Acc('norm', LossIntegral(p=xP, mn, sd, 'norm', F), LossIntegral(p=xP, mn, sd, 'norm'))
df <- rbind(df, Acc('gamma', LossIntegral(p=xP, mn, sd, 'gamma', F), LossIntegral(p=xP, mn, sd, 'gamma')))
df <- rbind(df, Acc('lnorm',LossIntegral(p=xP, mn, sd, 'lnorm', F), LossIntegral(p=xP, mn, sd,  'lnorm')))
View(df)


# and now for the inverse functions
# Note that the Exact values here will take some time to run
xFillRate <- runif(n=1000, min=.5, max = .9999)
df <- Acc('norm', InvLossIntegral(fillRate=xFillRate, mn, sd, 'norm', F), InvLossIntegral(fillRate=xFillRate, mn, sd, 'norm'))
df <- rbind(df, Acc('gamma', InvLossIntegral(fillRate=xFillRate, mn, sd, 'gamma', F), InvLossIntegral(fillRate=xFillRate, mn, sd, 'gamma')))
df <- rbind(df, Acc('lnorm', InvLossIntegral(fillRate=xFillRate, mn, sd, 'lnorm', F), InvLossIntegral(fillRate=xFillRate, mn, sd, 'lnorm')))
View(df)
