#' get the fdr and power for the two-stage method
#'
#' This function allows you to get the fdr and power for the two-stage method
#' @param n Sample size
#' @param p Variable dimension.
#' @param eta Desired FDR level.
#' @param maxbeta Signal size
#' @param alpha Tuning parameter for stage 1


teststatlist <- function(n,p,eta,maxbeta, alpha)
{
  mu <- rep(0,p)
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) -
                    (1:p - 1))
  cov <- 0.5^exponent
  l_jk = c()
  l0_jk = c()
  l1_jk = c()
  h_0 = 0
  X <- mvrnorm(n=n,mu=mu, Sigma=cov)
  l3 <- c(0,maxbeta)
  for (i in 1:(p-1)){
    x1 <- X[,i]
    for (k in (i+1):p){
      x2 <- X[,k]
      beta0 = -5
      beta1 = sample(l3,1)
      beta2 = sample(l3,1)
      beta3 = sample(l3,1)
      if (beta3 == 0){h_0 = h_0+1}
      j <- sample((i+1):p,1)
      x3 <- X[,j]
      y = beta0 + beta1*x1 + beta2*x2 + beta3*x1*x2  + x3 + rnorm(1)
      fit1 <- lm(y~x1)
      t_j = abs(fit1$coefficients[2]/sqrt(vcovHC(fit1)[2,2]))
      fit2 <- lm(y~x2)
      t_k = abs(fit2$coefficients[2]/sqrt(vcovHC(fit2)[2,2]))
      if( (t_j >= alpha)&(t_k >= alpha) ){
        fit3 <- lm(y~x1+x2+x1*x2)
        t = fit3$coefficients[4]/sqrt(vcovHC(fit3)[4,4])
        l_jk <- append(l_jk,abs(t))
        if(beta3==0){
          l0_jk <- append(l0_jk,abs(t))
        }else {l1_jk <- append(l1_jk,abs(t))}
      }
    }
  }
  h_1 = p*(p-1)/2 - h_0
  #twostagepr = sum(l1_jk >= t_hat)/h_1
  #fdr = sum(l0_jk >= t_hat)/max(sum(l_jk >= t_hat),1)

  return(list(l_jk,l0_jk,l1_jk,h_1))
}
