#' Plot of Posterior Samples
#' 
#' Plot posterior means and 95\% credible intervals for a_j, b_j, u_i parameters.
#' @param fit The fitted model object from BayesODT
#' @return Posterior means and 95\% C.I.s of the parameters.
#' @export



plot_ABU <- function(fit){
  
  n <- dim(fit$A)[2]
  m <- dim(fit$U)[2]
  num <- fit$num
  
  par(mfrow=c(3,1), mar=c(5, 5.5, 2, 1))
  a.data <- fit$A-rowMeans(fit$A)
  plot(colMeans(a.data), ylim=c(min(a.data)-0.1,max(a.data)), ylab=expression(a[j]), xlab="rater", cex.lab=2)
  abline(h=0, col="gray")
  arrows( 1:n, apply(a.data, 2, function(x) quantile(x, 0.025)), 1:n, apply(a.data, 2, function(x) quantile(x, 0.975)), code=0)
  b.data <- fit$B
  plot(colMeans(b.data), ylim=c(min(b.data)-0.1,max(b.data)), ylab=expression(b[j]), xlab="rater", cex.lab=2)
  arrows( 1:n, apply(b.data, 2, function(x) quantile(x, 0.025)), 1:n, apply(b.data, 2, function(x) quantile(x, 0.975)), code=0)
  abline(h=mean(b.data), col="gray")
  u.data <- fit$U  
  plot(colMeans(u.data), ylim=c(min(u.data)-0.1,max(u.data)), ylab=expression(u[i]), xlab="subject", cex.lab=2)
  arrows( 1:m, apply(u.data, 2, function(x) quantile(x, 0.025)), 1:m, apply(u.data, 2, function(x) quantile(x, 0.975)), code=0)
  points(which(fit$D==1),colMeans(u.data)[which(fit$D==1)], pch=16, col="red" )
}

