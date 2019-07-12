#' Posterior Samples of the Parameters
#' 
#' Obtain posterior samples of the parameters, a_j, b_j, u_i and alpha_k.
#' @param fit The fitted model object from BayesODT
#' @return aj Posterior samples of a_j parameters
#' @return bj Posterior samples of b_j parameters
#' @return ui Posterior samples of u_i parameters
#' @return alphak Posterior samples of alpha[k] parameters
#' @return aj.ci 95\% C.I. of a_j parameters
#' @return bj.ci 95\% C.I. of b_j parameters
#' @return ui.ci 95\% C.I. of c_i parameters
#' @return alphak.ci 95\% C.I. of alpha_k parameters
#' @export


summary <- function(fit){
  
  OUT <- list(aj = colMeans(fit$A)-mean(fit$A),
              aj.ci = apply(fit$A-mean(fit$A), 2, function(x) quantile(x, c(0.025, 0.975))),
              bj = colMeans(fit$B),
              bj.ci = apply(fit$B, 2, function(x) quantile(x, c(0.025, 0.975))),
              ui = colMeans(fit$U),
              ui.ci = apply(fit$U, 2, function(x) quantile(x, c(0.025, 0.975))),
              alphak = colMeans(fit$alpha)-mean(fit$A),
              alphak.ci = apply(fit$alpha-mean(fit$A), 2, function(x) quantile(x, c(0.025, 0.975))))
  return(OUT)
}

