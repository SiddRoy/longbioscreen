#' Stable log-sum-exp
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param x Numeric vector or parameter value.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
logsum <- function(x) {
  if (length(x) == 1) {
    return(x)
  }
  w <- which.max(x)
  if (is.infinite(x[w])) {
    return(log(sum(exp(x))))
  }
  x[w] + log1p(sum(exp(x[-w] - x[w])))
}

#' Stable log one-plus-exp
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param x Numeric vector or parameter value.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
log1p_exp <- function(x) -plogis(-x, log.p = TRUE)

#' Stable log one-minus-exp
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param x Numeric vector or parameter value.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
log_neg_expm1 <- function(x) pexp(-x, log.p = TRUE)
