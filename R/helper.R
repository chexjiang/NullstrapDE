#' Binary search to find FDP-controlled threshold
#'
#' @param coef_real Numeric vector of statistics from real data
#' @param coef_syn Numeric vector of statistics from synthetic null
#' @param q Target FDR level
#' @param conservative Logical; whether to use conservative mode
#'
#' @return Numeric threshold
#' @keywords internal
binary_search <- function(coef_real, coef_syn, q, conservative = FALSE) {

  # -------------------------
  # Input validation
  # -------------------------
  if (length(coef_real) == 0 || length(coef_syn) == 0) {
    stop("coef_real and coef_syn must be non-empty numeric vectors.")
  }
  if (!is.numeric(coef_real) || !is.numeric(coef_syn)) {
    stop("coef_real and coef_syn must be numeric.")
  }
  if (!is.numeric(q) || q <= 0 || q >= 1) {
    stop("q must be a number between 0 and 1.")
  }

  # Remove NA / Inf
  coef_real <- coef_real[is.finite(coef_real)]
  coef_syn  <- coef_syn[is.finite(coef_syn)]

  if (length(coef_real) == 0 || length(coef_syn) == 0) {
    stop("coef_real or coef_syn contains no finite values.")
  }

  # -------------------------
  # Binary search
  # -------------------------
  left  <- 0
  right <- max(coef_real)

  # Guard: all values same
  if (right == 0) return(0)

  while (abs(left - right) > 1e-8) {

    mid <- (left + right) / 2

    # Count positives and negatives
    num_negative <- sum(coef_syn >= mid)
    if (conservative) {
      num_negative <- num_negative + 1
    } else {
      num_negative <- num_negative + q/2
    }

    num_positive <- max(sum(coef_real >= mid), 1)

    FDP <- num_negative / num_positive

    if (FDP > q) {
      left <- mid
    } else {
      right <- mid
    }
  }

  return(right)
}
