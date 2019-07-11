#' ROC curves
#'
#' Produce ROC curves
#' @param fit The fitted model object from BayesODT
#' @param method Select a method to produce an ROC curve ("Individual", "ROC1", "ROC2", "Smoothed.ROC", "Smoothed.pROC", "pROC")
#' @param ind For an individual ROC curve, indicate patient or subject ID
#' @return Plot of ROC curves.
#' @export

plot_ROC <- function(fit, method=c("Individual","ROC1","ROC2","Smoothed.ROC", "Smoothed.pROC", "pROC"), ind=NULL){

  #REQUIRED LIBRARIES
#  library(pROC)

  K <- dim(fit$alpha)[2]
  n <- dim(fit$A)[2]
  m <- dim(fit$U)[2]
  num <- fit$num

  sp <- fit$sp
  se <- fit$se
  sp.mean <- fit$sp.mean
  se.mean <- fit$se.mean

  par(mfrow=c(1,1))

  if(method == "Individual"){
    t <- ind
    plot(mean(sp[t,1,]),mean(se[t,1,]), xlim=c(0,1), ylim=c(0,1), xlab="False Positive Rate", ylab="True Positive Rate", main="Estimated Individual ROC",cex=1.5)
    for(k in 2:K){
      points(mean(sp[t,k,]),mean(se[t,k,]),cex=1.5)
    }
    lines(c(1,mean(sp[t,1,])),c(1,mean(se[t,1,])), lwd=1.5)
    for(k in 2:K){
      lines(c(mean(sp[t,k-1,]),mean(sp[t,k,])),c(mean(se[t,k-1,]),mean(se[t,k,])), lwd=1.5)
    }
    lines(c(mean(sp[t,K,]),0),c(mean(se[t,K,]),0), lwd=1.5)
    rater <- paste("Rater ",t, sep="")
    legend("bottomright", legend=rater, col=c(1),pch=c(1),lty=c(1), lwd=c(1.5))
  }

  if(method == "ROC1"){
    plot(mean(sp[,1,]),mean(se[,1,]), xlim=c(0,1), ylim=c(0,1), xlab="False Positive Rate", ylab="True Positive Rate", main="Estimated ROC1",cex=2)
    for(k in 2:K){
      points(mean((sp[,k,])),mean((se[,k,])),pch=1,col=1,cex=2)
    }
    lines(c(1,mean((sp[,1,]))),c(1,mean((se[,1,]))),col=1, lty=3, lwd=3)
    for(k in 2:K){
      lines(c(mean((sp[,k-1,])),mean((sp[,k,]))),c(mean((se[,k-1,])),mean((se[,k,]))),col=1, lty=3, lwd=3)
    }
    lines(c(mean((sp[,K,])),0),c(mean((se[,K,])),0),col=1, lty=3, lwd=3)
    legend("bottomright", legend="ROC1", col=c(1),pch=c(1),lty=c(1), lwd=c(1.5))
  }

  if(method == "ROC2"){

    plot(mean(sp.mean[,1]),mean(se.mean[,1]), xlim=c(0,1), ylim=c(0,1), xlab="False Positive Rate", ylab="True Positive Rate", main="Estimated ROC2",cex=2)
    for(k in 2:K){
      points(mean(sp.mean[,k]),mean(se.mean[,k]),pch=1,col=1,cex=2)
    }
    lines(c(1,mean(sp.mean[,1])),c(1,mean(se.mean[,1])),col=1, lty=3, lwd=3)
    for(k in 2:K){
      lines(c(mean(sp.mean[,k-1]),mean(sp.mean[,k])),c(mean(se.mean[,k-1]),mean(se.mean[,k])),col=1, lty=3, lwd=3)
    }
    lines(c(mean(sp.mean[,K]),0),c(mean(se.mean[,K]),0),col=1, lty=3, lwd=3)
    legend("bottomright", legend="ROC2", col=c(1),pch=c(1),lty=c(1), lwd=c(1.5))
  }

  if(method == "Smoothed.ROC"){

    SP <- function(h){
      temp <- mean(rowSums((pnorm(h-(mean(fit$A)+mean(fit$B)*fit$U)))*(1-pnorm(fit$U)))/rowSums(1-pnorm(fit$U)))
    }
    SE <- function(h){
      temp <- mean(rowSums((1-pnorm(h-(mean(fit$A)+mean(fit$B)*fit$U)))*(pnorm(fit$U)))/rowSums(pnorm(fit$U)))
    }

    x <- sapply(seq(-5,5,length=100), function(x) 1-SP(x))
    y <- sapply(seq(-5,5,length=100), function(x) SE(x))
    lo <- loess(y~x)

    plot(x,y, col=1, xlab="False Positive Rate", ylab="True Positive Rate",  type="l", main="Smoothed ROC", xlim=c(0,1), lwd=2)
    legend("bottomright", legend="Smoothed ROC", col=c(1),pch=c(1),lty=c(1), lwd=c(1.5))
  }

  if(method == "pROC"){
    roc <- roc(rep(D, n) ~ c(as.matrix(w)), smooth=F)
    plot(1-roc$specificities, roc$sensitivities, type="l", lwd=2, col=1,lty=3)
    legend("bottomright", legend="pROC", col=c(1),pch=c(1),lty=c(1), lwd=c(1.5))
  }

  if(method == "Smoothed.pROC"){
    roc <- roc(rep(D, n) ~ c(as.matrix(w)), smooth=T)
    plot(1-roc$specificities, roc$sensitivities, type="l", lwd=2, col=1,lty=3)
    legend("bottomright", legend="Smoothed pROC", col=c(1),pch=c(1),lty=c(1), lwd=c(1.5))
  }
}


