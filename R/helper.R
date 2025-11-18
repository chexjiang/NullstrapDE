library(parallel)


binary_search <- function(coef_real, coef_syn, q, conservative = FALSE) {
  # coef_syn is beta_hat from synthetic null data
  # coef_real is beta_hat from real data
  left = 0
  right = max(coef_real)
  while (abs(left - right) > 1e-8) {
    mid = ((left + right) / 2)
    if(conservative){
      num_negative = length(which(coef_syn >= mid)) + 1
    }else{
      num_negative = length(which(coef_syn >= mid)) + q/2
    }
    num_positive = max(length(which(coef_real >= mid)), 1)
    FDP = num_negative / num_positive
    if (FDP > q) {
      left = mid
    } else {
      right = mid
    }
  }
  return(right)
}



