#' Create Wald confidence intervals
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param est Estimate vector.
#' @param se Longitudinal residual standard deviation.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
ci_fit <- function(est, se) {
  tibble::tibble(est = est, se = se, lci = est - qnorm(0.975) * se, rci = est +
    qnorm(0.975) * se)
}

#' Compute two-sided normal p-values
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param t_ Input used by this helper.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
pfun_t <- function(t_) 2 * pnorm(-abs(t_))

#' Compute p-values from estimates and standard errors
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param b Subject-specific random intercept and slope values.
#' @param se Longitudinal residual standard deviation.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
pfun_bse <- function(b, se) pfun_t(-abs(b / se))

#' Compute a Wald chi-square p-value
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param x Numeric vector or parameter value.
#' @param V Variance-covariance matrix.
#' @param ... Additional arguments passed to downstream functions.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
pfun_chi2 <- function(x, V, ...) {
  R <- chol(V)
  pchisq(sum(backsolve(R, x, transpose = TRUE)^2), length(x),
    lower.tail = FALSE,
    ...
  )
}

#' Compute log two-sided normal p-values
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param t_ Input used by this helper.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
pfun_t_log <- function(t_) log(2) + pnorm(-abs(t_), log.p = TRUE)

#' Compute log p-values from estimates and standard errors
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param b Subject-specific random intercept and slope values.
#' @param se Longitudinal residual standard deviation.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
pfun_bse_log <- function(b, se) pfun_t_log(-abs(b / se))
