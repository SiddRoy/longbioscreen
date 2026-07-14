#' Build a random-effect Cholesky factor
#'
#' Advanced simulation helper for correlated random intercept/slope effects.
#' @family simulation generators
#' @param sb Random-effect standard deviations.
#' @param rho Random-effect correlation parameter.
#' @return A 2 x 2 matrix used to transform independent normal draws.
#' @export
chol_IS2 <- function(sb, rho) {
  matrix(c(1, 0, rho, sqrt(1 - rho^2)), 2) %*% diag(sb)
}

#' Draw intercept, slope, and changepoint effects
#'
#' Advanced simulation helper for subject-specific intercept, slope, and
#' changepoint effects.
#' @family simulation generators
#' @param N Number of subjects to simulate.
#' @param sb Random-effect standard deviations.
#' @param rho Random-effect correlation parameter.
#' @param sc Random changepoint slope standard deviation.
#' @param ru Function that draws changepoints.
#' @return A list with `m` random intercept/slope effects and `mC`
#'   changepoint effects.
#' @export
get_bIS <- function(N, sb, rho, sc = 0, ru = function(N) runif(N, 6, 10)) {
  b <- matrix(rnorm(N * 2), N) %*% chol_IS2(sb, rho)
  bC <- cbind(rnorm(N) * sc, ru(N))
  list(m = b, mC = bC)
}

#' Create a subject trajectory function
#'
#' Advanced simulation helper that returns a subject-specific biomarker
#' trajectory function.
#' @family simulation generators
#' @param B Fixed longitudinal intercept and slope coefficients.
#' @param Btau Fixed changepoint or acute-effect coefficient.
#' @param b Subject-specific random intercept and slope values.
#' @param bC Subject-specific random changepoint values.
#' @return A function of visit time `tim`.
#' @export
mt_IScp_i <- function(B = c(0, 0), Btau = 0, b, bC) {
  function(tim) {
    B[1] + b[1] + (B[2] + b[2]) * tim + (Btau + bC[1]) * pmax(0, tim - bC[2])
  }
}
