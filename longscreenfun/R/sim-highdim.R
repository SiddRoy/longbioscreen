#' Draw one subject for the high-dimensional simulator
#'
#' Advanced helper used by [jmCox_draw_WD()] to simulate one subject's
#' high-dimensional joint-model and Cox biomarkers.
#' @family high-dimensional simulation
#' @param Z Log cumulative hazard draw used to determine the latent event time.
#' @param nE Number of associated biomarkers per mechanism.
#' @param drawC Function that draws censoring times.
#' @param t_max Maximum longitudinal visit time.
#' @param B Fixed longitudinal intercept and slope coefficients.
#' @param alp_JM Joint-model association vector.
#' @param alp_Cox Cox association vector.
#' @param sb Random-effect standard deviations.
#' @param b_JM_Cox Subject-specific high-dimensional random effects.
#' @param se Longitudinal residual standard deviation.
#' @param lh_kl Weibull shape and scale for the high-dimensional simulator.
#' @return A list with longitudinal rows `X` and one survival row `S` for one subject.
#' @export
jmCox_draw_W <- function(Z = log(rexp(1)), nE = 3, drawC = function(N) {
                           pmin(rexp(
                             N,
                             0.2
                           ), 5.9)
                         }, t_max = 100, B = c(0, 0), alp_JM = 0, alp_Cox = 0, sb = c(
                           0.8, 1,
                           1.2
                         ), b_JM_Cox = NULL, se = 0.5, lh_kl = c(2, 10)) {
  if (is.null(b_JM_Cox)) {
    b_JM_Cox <- rnorm(2 * nE, sd = c(sb, sb))
  }
  tim <- 0:t_max
  r_cox <- matrix(rnorm(length(tim) * nE, sd = se), length(tim))
  mt0 <- sum((b_JM_Cox[1:3] + B[1]) * alp_JM) + sum((b_JM_Cox[4:6] + B[1]) * alp_Cox)
  a <- alp_Cox
  if (length(alp_Cox) == 1) {
    mt0 <- rowSums(r_cox) * alp_Cox + mt0
  } else {
    mt0 <- as.numeric(r_cox %*% alp_Cox) + mt0
  }
  Si <- surv_Cox_W(Z, Ci = 100, tim, mt0)
  Ci <- drawC(1)
  Sout <- tibble(
    Ti = min(Si[1], Ci), del = as.integer(Si[1] <= Ci), Si = Si[1],
    Ci = Ci
  )
  w <- which(tim < Sout$Ti)
  yCox <- sweep(r_cox[w, , drop = FALSE], 2, B[1] + b_JM_Cox[4:6], "+") %>%
    as.data.frame() %>%
    as_tibble() %>%
    setNames(paste0("y", 1:nE))
  yJM <- (tcrossprod(rep(1, length(w)), B[1] + b_JM_Cox[1:3]) + rnorm(length(w) *
    nE, sd = se)) %>%
    as.data.frame() %>%
    as_tibble() %>%
    setNames(paste0(
      "y",
      1:nE + nE
    ))
  Xdat <- bind_cols(tibble(tim = tim[w]), yCox, yJM)
  list(X = Xdat, S = Sout)
}

#' Draw a high-dimensional joint/Cox dataset
#'
#' Advanced helper used by [drawY100()] to generate the joint-model and Cox
#' biomarker block before conditional-model and null biomarkers are added.
#' @family high-dimensional simulation
#' @param N Number of subjects to simulate.
#' @param nE Number of associated biomarkers per mechanism.
#' @param seed Optional random seed.
#' @param drawC Function that draws censoring times.
#' @param t_max Maximum longitudinal visit time.
#' @param B Fixed longitudinal intercept and slope coefficients.
#' @param alp_JM Joint-model association vector.
#' @param alp_Cox Cox association vector.
#' @param sb Random-effect standard deviations.
#' @param rho Random-effect correlation parameter.
#' @param se Longitudinal residual standard deviation.
#' @param lh_kl Weibull shape and scale for the high-dimensional simulator.
#' @return A list with longitudinal rows `X` and survival rows `S`.
#' @export
jmCox_draw_WD <- function(N = 600, nE = 3, seed = NULL, drawC = function(N) {
                            pmin(rexp(
                              N,
                              0.2
                            ), 5.9)
                          }, t_max = 100, B = c(0, 0), alp_JM = 0, alp_Cox = 0, sb = c(
                            0.8, 1,
                            1.2
                          ), rho = 0, se = 0.5, lh_kl = c(2, 10)) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  Z <- log(rexp(N))
  M <- matrix(1, nE * 2, nE * 2)
  R <- chol(rho^(abs(row(M) - col(M)))) %*% diag(c(sb, sb))
  b <- crossprod(R, matrix(rnorm(N * 2 * nE), 2 * nE))
  simJC <- lapply(1:N, function(i) {
    jmCox_draw_W(Z[i],
      nE = nE, drawC = drawC, t_max = t_max, B = B, alp_JM = alp_JM,
      alp_Cox = alp_Cox, sb = sb, b_JM_Cox = b[, i], se = se, lh_kl = lh_kl
    )
  })
  Sdat <- bind_rows(lapply(1:N, function(i) simJC[[i]]$S)) %>% mutate(i = 1:n())
  Xdat <- bind_rows(lapply(1:N, function(i) {
    mutate(simJC[[i]]$X, i = i)
  }))
  list(X = Xdat, S = Sdat)
}

#' Default effect sizes for the 100-biomarker simulator
#'
#' Returns the default joint-model, Cox, and conditional-model effects used by
#' [drawY100()]. This is exported so simulation variants can start from the
#' original effect configuration and modify only selected components.
#' @family high-dimensional simulation
#' @param m Biomarker standard deviation multipliers.
#' @return A list with entries `JM`, `Cox`, and `CM`.
#' @export
default_drawY100_effects <- function(m = c(0.9, 1, 1.1)) {
  list(
    JM = 0.62 / m,
    Cox = 0.5 / m,
    CM = list(L = c(-0.06, 0), A = c(0, 0.9), B = c(-0.03, 0.45))
  )
}

#' Create a block correlation matrix
#'
#' Builds a correlation matrix with constant within-block correlation and zero
#' cross-block correlation.
#' @family high-dimensional simulation
#' @param p Matrix dimension.
#' @param rho Within-block correlation.
#' @param K Block size.
#' @return A `p` by `p` correlation matrix.
#' @export
rho_fun_blockK <- function(p, rho, K = 5) {
  if (length(p) != 1 || !is.finite(p) || p < 1 || p != as.integer(p)) {
    stop("`p` must be a positive integer.", call. = FALSE)
  }
  if (length(K) != 1 || !is.finite(K) || K < 1 || K != as.integer(K)) {
    stop("`K` must be a positive integer.", call. = FALSE)
  }
  if (length(rho) != 1 || !is.finite(rho)) {
    stop("`rho` must be a finite numeric scalar.", call. = FALSE)
  }

  p <- as.integer(p)
  K <- as.integer(K)
  block_id <- ceiling(seq_len(p) / K)
  out <- matrix(0, p, p)
  out[outer(block_id, block_id, "==")] <- rho
  diag(out) <- 1
  out
}

#' Create a banded correlation matrix
#'
#' Builds a correlation matrix with linearly decreasing off-diagonal band
#' correlations and zero correlation outside the requested bandwidth.
#' @family high-dimensional simulation
#' @param p Matrix dimension.
#' @param rho Correlation on the first off-diagonal band.
#' @param K Number of off-diagonal bands to fill.
#' @return A `p` by `p` correlation matrix.
#' @export
rho_fun_band <- function(p, rho, K = 5) {
  if (length(p) != 1 || !is.finite(p) || p < 1 || p != as.integer(p)) {
    stop("`p` must be a positive integer.", call. = FALSE)
  }
  if (length(K) != 1 || !is.finite(K) || K < 1 || K != as.integer(K)) {
    stop("`K` must be a positive integer.", call. = FALSE)
  }
  if (length(rho) != 1 || !is.finite(rho)) {
    stop("`rho` must be a finite numeric scalar.", call. = FALSE)
  }

  p <- as.integer(p)
  K <- as.integer(K)
  out <- diag(p)
  band_values <- seq(rho, 0, length.out = K)
  for (band in seq_len(min(K, p - 1))) {
    idx <- seq_len(p - band)
    out[cbind(idx, idx + band)] <- band_values[band]
    out[cbind(idx + band, idx)] <- band_values[band]
  }
  out
}

#' Draw the 100-biomarker screening dataset
#'
#' Primary high-dimensional simulator used by `35_HDsim/HDsim.R`.
#' @family high-dimensional simulation
#' @param seed Optional random seed.
#' @param th Scalar multiplier applied to the default effect sizes.
#' @param th0 Optional list of JM, Cox, and CM effect sizes. If `NULL`, uses
#'   [default_drawY100_effects()].
#' @param N Number of subjects to simulate.
#' @param p Number of biomarkers to return.
#' @param m Biomarker standard deviation multipliers.
#' @param rho_fun Optional function used to generate the null-biomarker
#'   correlation matrix from `lsb` and `rho`.
#' @param rho Optional correlation parameter passed to `rho_fun`.
#' @return A list with `X`, `S`, and biomarker matrix `y`. Rows of `y` align to
#'   rows of `X`; columns are candidate biomarkers.
#' @export
drawY100 <- function(seed = NULL, th = 1, th0 = NULL, N = 600, p = 100, m = c(
                       0.9,
                       1, 1.1
                     ), rho_fun = NULL, rho = NULL) {
  if (is.null(th0)) {
    th0 <- default_drawY100_effects(m)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  outS <- jmCox_draw_WD(N, alp_JM = th0$JM * th, alp_Cox = th0$Cox * th, sb = m)
  y_CM_j <- lapply(th0$CM, function(cm_alp) {
    CM_givT_IS(
      N = N, alp = cm_alp[1] * th, Btau = cm_alp[2] * th, Ti = outS$S$Si,
      Ci = outS$S$Ci
    )
  })
  lsb <- p - 9
  sb0 <- rep(c(1, 0.9, 1.1), ceiling(lsb / length(m)))[1:lsb]
  b0_raw <- matrix(rnorm(N * lsb), N)
  if (!is.null(rho_fun)) {
    if (is.null(rho)) {
      stop("`rho` must be provided when `rho_fun` is provided.", call. = FALSE)
    }
    Cor_b0 <- rho_fun(lsb, rho)
    if (!is.matrix(Cor_b0) || !is.numeric(Cor_b0) ||
        !all(dim(Cor_b0) == c(lsb, lsb)) || any(!is.finite(Cor_b0))) {
      stop("`rho_fun(lsb, rho)` must return a finite numeric `lsb` by `lsb` matrix.", call. = FALSE)
    }
    if (!isTRUE(all.equal(Cor_b0, t(Cor_b0)))) {
      stop("`rho_fun(lsb, rho)` must return a symmetric matrix.", call. = FALSE)
    }
    if (!isTRUE(all.equal(diag(Cor_b0), rep(1, lsb)))) {
      stop("`rho_fun(lsb, rho)` must return a matrix with diagonal entries equal to 1.", call. = FALSE)
    }
    R_b0 <- tryCatch(
      chol(Cor_b0),
      error = function(e) {
        stop("`rho_fun(lsb, rho)` must return a positive definite matrix.", call. = FALSE)
      }
    )
    b0_raw <- b0_raw %*% R_b0
  }
  b0 <- sweep(b0_raw, 2, sb0, "*")
  B1 <- rnorm(lsb, -0.25, sd = 0.1)
  y0 <- b0[outS$X$i, ] + tcrossprod(outS$X$tim, B1) + rnorm(nrow(outS$X) * lsb,
    sd = 0.5
  )
  y <- cbind(as.matrix(dplyr::select(outS$X, starts_with("y"))), do.call(
    "cbind",
    lapply(y_CM_j, function(cmy) cmy$y)
  ), y0)
  dimnames(y) <- NULL
  list(X = rename(dplyr::select(outS$X, i, tim, y1), y = y1) %>% mutate(i = factor(
    i,
    unique(i)
  )), S = dplyr::select(outS$S, i, Ti, del) %>% mutate(i = factor(
    i,
    unique(i)
  )), y = y)
}

#' Draw the 100-biomarker log-time screening dataset
#'
#' Log-time variant of [drawY100()] where the conditional-model biomarker block
#' is generated with `log(Ti)` as the event-time term. The joint-model, Cox, and
#' null biomarker blocks follow [drawY100()].
#' @family high-dimensional simulation
#' @param seed Optional random seed.
#' @param th Scalar multiplier applied to the default effect sizes.
#' @param th0 Optional list of JM, Cox, and CM effect sizes. If `NULL`, uses
#'   [default_drawY100_effects()].
#' @param N Number of subjects to simulate.
#' @param p Number of biomarkers to return.
#' @param m Biomarker standard deviation multipliers.
#' @param rho_fun Optional function used to generate the null-biomarker
#'   correlation matrix from `lsb` and `rho`.
#' @param rho Optional correlation parameter passed to `rho_fun`.
#' @return A list with `X`, `S`, and biomarker matrix `y`. Rows of `y` align to
#'   rows of `X`; columns are candidate biomarkers.
#' @export
drawY100_logT <- function(seed = NULL, th = 1, th0 = NULL, N = 600, p = 100, m = c(
                           0.9,
                           1, 1.1
                         ), rho_fun = NULL, rho = NULL) {
  if (is.null(th0)) {
    th0 <- default_drawY100_effects(m)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  outS <- jmCox_draw_WD(N, alp_JM = th0$JM * th, alp_Cox = th0$Cox * th, sb = m)
  y_CM_j <- lapply(th0$CM, function(cm_alp) {
    CM_givT_IS_logT(
      N = N, alp = cm_alp[1] * th, Btau = cm_alp[2] * th, Ti = outS$S$Si,
      Ci = outS$S$Ci
    )
  })
  lsb <- p - 9
  sb0 <- rep(c(1, 0.9, 1.1), ceiling(lsb / length(m)))[1:lsb]
  b0_raw <- matrix(rnorm(N * lsb), N)
  if (!is.null(rho_fun)) {
    if (is.null(rho)) {
      stop("`rho` must be provided when `rho_fun` is provided.", call. = FALSE)
    }
    Cor_b0 <- rho_fun(lsb, rho)
    if (!is.matrix(Cor_b0) || !is.numeric(Cor_b0) ||
        !all(dim(Cor_b0) == c(lsb, lsb)) || any(!is.finite(Cor_b0))) {
      stop("`rho_fun(lsb, rho)` must return a finite numeric `lsb` by `lsb` matrix.", call. = FALSE)
    }
    if (!isTRUE(all.equal(Cor_b0, t(Cor_b0)))) {
      stop("`rho_fun(lsb, rho)` must return a symmetric matrix.", call. = FALSE)
    }
    if (!isTRUE(all.equal(diag(Cor_b0), rep(1, lsb)))) {
      stop("`rho_fun(lsb, rho)` must return a matrix with diagonal entries equal to 1.", call. = FALSE)
    }
    R_b0 <- tryCatch(
      chol(Cor_b0),
      error = function(e) {
        stop("`rho_fun(lsb, rho)` must return a positive definite matrix.", call. = FALSE)
      }
    )
    b0_raw <- b0_raw %*% R_b0
  }
  b0 <- sweep(b0_raw, 2, sb0, "*")
  B1 <- rnorm(lsb, -0.25, sd = 0.1)
  y0 <- b0[outS$X$i, ] + tcrossprod(outS$X$tim, B1) + rnorm(nrow(outS$X) * lsb,
    sd = 0.5
  )
  y <- cbind(as.matrix(dplyr::select(outS$X, starts_with("y"))), do.call(
    "cbind",
    lapply(y_CM_j, function(cmy) cmy$y)
  ), y0)
  dimnames(y) <- NULL
  list(X = rename(dplyr::select(outS$X, i, tim, y1), y = y1) %>% mutate(i = factor(
    i,
    unique(i)
  )), S = dplyr::select(outS$S, i, Ti, del) %>% mutate(i = factor(
    i,
    unique(i)
  )), y = y)
}
