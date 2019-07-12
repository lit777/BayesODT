#' Goodness of Fit
#' 
#' Check goodness of fit of the model.
#' @param fit The fitted model object from BayesODT
#' @return Estimated.Prob The estimated probabilities for each classification
#' @return Estimated.Prob.ci The 95\% C.I.s for the estimated probabilities
#' @return Empirical.Prob The empirical probabilities from the observed data
#' @export

GoF <- function(fit){
  
  K <- dim(fit$alpha)[2]+1
  n <- dim(fit$A)[2]
  m <- dim(fit$U)[2]
  num <- fit$num
  
  alpha <- cbind(fit$alpha,Inf) 
  
  sp <- se <- array(dim=c(n, K, num))
  
  for(i in 1:num){
    sp[,,i] <- as.matrix(sapply(1:(K), function(x) (1/m)*(pnorm(alpha[i,x]-(matrix(fit$A[i,], ncol=m, nrow=n, byrow=FALSE)+outer(fit$B[i,],fit$U[i,]))))%*%(1-pnorm(fit$U[i,]))/mean(1-pnorm(fit$U[i,])) ) )
    se[,,i] <- as.matrix(sapply(1:(K), function(x) (1/m)*(pnorm(alpha[i,x]-(matrix(fit$A[i,], ncol=m, nrow=n, byrow=FALSE)+outer(fit$B[i,],fit$U[i,]))))%*%(pnorm(fit$U[i,]))/mean(pnorm(fit$U[i,])) ) )
  }
  
  good.temp <- matrix(nrow=num, ncol=K)
  good <- matrix(nrow=num, ncol=K)
  
  for(i in 1:num){
    for(k in 1:K){
      good.temp[i,k] <- mean(sp[,k,i])*mean(1-pnorm(fit$U[i,]))+mean(se[,k,i])*mean(pnorm(fit$U[i,]))
    }
    good[i,] <- c(good.temp[i,1], diff(good.temp[i,], lag = 1, differences = 1))
  }
  out <- list(Estimated.Prob=colMeans(good),
              Estimated.Prob.ci=apply(good, 2, function(x) quantile(x, c(0.025, 0.975))),
              Empirical.Prob = table(as.matrix(fit$w))/(length(as.matrix(fit$w))))
  return(out)
}
