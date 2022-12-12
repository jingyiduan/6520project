#' search the cut-off point controlling FDR
#'
#' This function allows you to find the t_hat
#' @param l1 A list of test statistics.
#' @param eta Desired FDR level.
#' @param p Variable dimension.


searcht <- function(l1,eta,p){
  l <- sort(l1)
  sum1 = length(l)
  sum0 = sum1
  l <- append(l,0,0)
  for(i in 1:(length(l)-1)){
    t_pre = l[i]
    t = l[i+1]
    Gt_pre=2*(1-pnorm(t_pre))
    Gt=2*(1-pnorm(t))
    if (Gt/max(sum1,1)*sum0 <= eta){
      tail_p = eta*max(sum1,1)/sum0
      t_hat = qnorm((2-tail_p)/2.0)
      if (t_hat >= t_pre){
        return(t_hat)
      }else{
        return(t)
      }
    }
    sum1 = sum1-1
  }
  tail_p = eta/sum0
  t_hat = qnorm((2-tail_p)/2.0)
  if(t_hat < sqrt(2)*log(p)){
    return(t_hat)
  }else{
    return(sqrt(2)*log(p))
  }
}
