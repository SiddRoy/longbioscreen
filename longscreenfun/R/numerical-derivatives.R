#' Differentiate with Richardson extrapolation
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param x Numeric vector or parameter value.
#' @param f Function to differentiate.
#' @param h0 Base finite-difference step size.
#' @param r Richardson depth or residual vector, depending on context.
#' @param h Finite-difference step size.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
rich_g <- function(x, f, h0 = 1e-05, r = 4, h = NULL) {
  if (is.null(h)) {
    h <- h0 * (1 + abs(x))
  }
  if (r == 1) {
    return((f(x + h) - f(x - h)) / (2 * h))
  }
  g_r <- sapply(1:r, function(rj) {
    hj <- h / 2^(rj - 1)
    (f(x + hj) - f(x - hj)) / (2 * hj)
  })
  g_r <- matrix(g_r, ncol = r)
  r0 <- 1
  while (r0 < r) {
    r0_col <- r - r0 + 1
    g_r[, 1:(r0_col - 1)] <- (g_r[, 2:r0_col, drop = FALSE] * (4^r0) - g_r[,
      1:(r0_col - 1),
      drop = FALSE
    ]) / (4^r0 - 1)
    r0 <- r0 + 1
  }
  g_r[, 1]
}

#' Compute a numerical Jacobian
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param x Numeric vector or parameter value.
#' @param f Function to differentiate.
#' @param h0 Base finite-difference step size.
#' @param verb Whether to print finite-difference progress.
#' @param h Finite-difference step size.
#' @param r Richardson depth or residual vector, depending on context.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
g_xl <- function(x, f, h0 = 1e-05, verb = FALSE, h = NULL, r = 4) {
  p <- length(x)
  if (is.null(h) | (length(h) != length(x))) {
    h <- h0 * (1 + abs(x))
  }
  sfun <- function(j) {
    if (verb) {
      print(j)
    }
    fj <- function(x_) {
      x0 <- x
      x0[j] <- x_
      f(x0)
    }
    rich_g(x[j], fj, h0, r, h = h[j])
  }
  sapply(1:length(x), sfun)
}

#' Compute a cross-Hessian block numerically
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param x Numeric vector or parameter value.
#' @param y Second argument or response vector.
#' @param f Function to differentiate.
#' @param h0 Base finite-difference step size.
#' @param verb Whether to print finite-difference progress.
#' @param r Richardson depth or residual vector, depending on context.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
H_ij <- function(x, y, f, h0 = 1e-04, verb = FALSE, r = 4) {
  nY <- length(y)
  h_y <- h0 * (1 + abs(y))
  h_x <- h0 * (1 + abs(x))
  sy_fun <- function(m) {
    if (verb) {
      print(m)
    }
    g_m <- function(y_m) {
      y0 <- y
      y0[m] <- y_m
      g_xl(x, function(x) f(x, y0), h0, verb = verb, r = r, h = h_x)
    }
    rich_g(y[m], g_m, h0, r, h = h_y[m])
  }
  sapply(1:nY, sy_fun)
}

#' Compute a Hessian numerically
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param x Numeric vector or parameter value.
#' @param f Function to differentiate.
#' @param h0 Base finite-difference step size.
#' @param verb Whether to print finite-difference progress.
#' @param r Richardson depth or residual vector, depending on context.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
H_ii <- function(x, f, h0 = 1e-04, verb = FALSE, r = 4) {
  h_x <- h0 * (1 + abs(x))
  sy_fun <- function(m) {
    if (verb) {
      print(m)
    }
    g_m <- function(y_m) {
      x0 <- x
      x0[m] <- y_m
      g_xl(x0, f, h0, verb = verb, r = r, h = h_x)
    }
    rich_g(x[m], g_m, h0, r, h = h_x[m])
  }
  sapply(1:length(x), sy_fun)
}
