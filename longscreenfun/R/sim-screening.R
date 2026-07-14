#' Simulate a joint-model screening dataset
#'
#' Primary low-dimensional simulator for the change-point joint-model workflow.
#' @family simulation generators
#' @param N Number of subjects to simulate.
#' @param seed Optional random seed.
#' @param drawC Function that draws censoring times.
#' @param alp Association parameter.
#' @param t_max Maximum longitudinal visit time.
#' @param B Fixed longitudinal intercept and slope coefficients.
#' @param Btau Fixed changepoint or acute-effect coefficient.
#' @param sb Random-effect standard deviations.
#' @param rho Random-effect correlation parameter.
#' @param sc Random changepoint slope standard deviation.
#' @param se Longitudinal residual standard deviation.
#' @param ru Function that draws changepoints.
#' @param lh Log baseline hazard function.
#' @return A list with longitudinal rows `X` and survival rows `S`.
#' @export
JM_sim_ISCP <- function(
  N, seed = NULL, drawC = function(N) pmin(rexp(N, 0.2), 5.9),
  alp = 0, t_max = 5, B = c(0.3, -0.25), Btau = 0, sb = c(1, 0), rho = 0, sc = 0,
  se = 0.5, ru = function(N) runif(N, 6, 10), lh = function(x) lh_w(x, 2, 10)
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  Ci <- drawC(N)
  Zi <- log(rexp(N))
  b <- get_bIS(N, sb, rho, sc, ru)
  if (Btau == 0) {
    b$mC[, 2] <- Inf
  }
  Xdat <- lapply(1:N, function(j) {
    my_fun <- mt_IScp_i(B, Btau, b$m[j, ], b$mC[j, ])
    tibble(i = j, tim = 0:floor(t_max), y = my_fun(tim))
  }) %>%
    bind_rows() %>%
    mutate(y = y + rnorm(n(), sd = se))
  Si <- sapply(1:N, function(j) {
    my_fun <- mt_IScp_i(B, Btau, b$m[j, ], b$mC[j, ])
    ms_fun <- function(x) my_fun(x) * alp
    survY_jm(Ci[j], Lam_mt_IScp, ms_fun, taui = b$mC[j, 2], Zi = Zi[j], lh = lh)
  })
  Sout <- tibble(i = 1:N, Ti = Si[1, ], del = Si[2, ])
  Xdat <- Xdat %>%
    left_join(Sout, by = "i") %>%
    filter(tim <= Ti)
  list(X = mutate(Xdat, i = factor(i, sort(unique(i)))), S = mutate(Sout, i = factor(
    i,
    sort(unique(i))
  )))
}

#' Simulate a Cox screening dataset
#'
#' Low-dimensional simulator where the observed biomarker enters a stepwise Cox
#' hazard.
#' @family simulation generators
#' @param N Number of subjects to simulate.
#' @param seed Optional random seed.
#' @param drawC Function that draws censoring times.
#' @param alp Association parameter.
#' @param t_max Maximum longitudinal visit time.
#' @param B Fixed longitudinal intercept and slope coefficients.
#' @param Btau Fixed changepoint or acute-effect coefficient.
#' @param sb Random-effect standard deviations.
#' @param rho Random-effect correlation parameter.
#' @param sc Random changepoint slope standard deviation.
#' @param se Longitudinal residual standard deviation.
#' @param ru Function that draws changepoints.
#' @param lH Log baseline cumulative hazard function.
#' @return A list with longitudinal rows `X` and survival rows `S`.
#' @export
Cox_sim_ISCP <- function(
  N, seed = NULL, drawC = function(N) pmin(rexp(N, 0.2), 5.9),
  alp = 0, t_max = 5, B = c(0.3, -0.25), Btau = 0, sb = c(1, 0), rho = 0, sc = 0,
  se = 0.5, ru = function(N) runif(N, 6, 10), lH = function(x) lHw(x, 2, 10)
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  Ci <- drawC(N)
  Zi <- log(rexp(N))
  b <- get_bIS(N, sb, rho, sc, ru)
  if (Btau == 0) {
    b$mC[, 2] <- Inf
  }
  Xdat <- lapply(1:N, function(j) {
    my_fun <- mt_IScp_i(B, Btau, b$m[j, ], b$mC[j, ])
    tibble(i = j, tim = 0:floor(t_max), y = my_fun(tim))
  }) %>%
    bind_rows() %>%
    mutate(y = y + rnorm(n(), sd = se))
  Si <- sapply(1:N, function(j) {
    Xi <- filter(Xdat, i == j) %>% mutate(y = y * alp)
    surv_Cox(Zi[j], Ci[j], Xi$tim, Xi$y, lH)
  })
  Sout <- tibble(i = 1:N, Ti = Si[1, ], del = Si[2, ])
  Xdat <- Xdat %>%
    left_join(Sout, by = "i") %>%
    filter(tim <= Ti)
  list(X = mutate(Xdat, i = factor(i, sort(unique(i)))), S = mutate(Sout, i = factor(
    i,
    sort(unique(i))
  )))
}

#' Evaluate the conditional-model acute kernel
#'
#' Acute-effect kernel used by the conditional-model simulator.
#' @family simulation generators
#' @param b1i Acute-effect coefficient.
#' @param xi Conditional-model kernel decay parameter.
#' @param S Survival data frame or event time.
#' @param tx Longitudinal visit time.
#' @return Numeric acute-effect contribution at each visit time.
#' @export
gfun0 <- function(b1i, xi, S, tx) {
  b1i * exp(xi * (tx - S))
}

#' Simulate a conditional-model screening dataset
#'
#' Low-dimensional simulator where biomarker values depend conditionally on
#' event time.
#' @family simulation generators
#' @param N Number of subjects to simulate.
#' @param seed Optional random seed.
#' @param drawC Function that draws censoring times.
#' @param alp Association parameter.
#' @param t_max Maximum longitudinal visit time.
#' @param B Fixed longitudinal intercept and slope coefficients.
#' @param Btau Fixed changepoint or acute-effect coefficient.
#' @param xi Conditional-model kernel decay parameter.
#' @param sb Random-effect standard deviations.
#' @param rho Random-effect correlation parameter.
#' @param sc Random changepoint slope standard deviation.
#' @param se Longitudinal residual standard deviation.
#' @param gf Acute-effect kernel function.
#' @param lHi Inverse log cumulative baseline hazard function.
#' @param Tifun Transformation applied to event time.
#' @param addcp Whether to include random changepoints.
#' @return A list with longitudinal rows `X` and survival rows `S`.
#' @export
CM_sim_RI <- function(
  N, seed = NULL, drawC = function(N) pmin(rexp(N, 0.2), 5.9),
  alp = 0, t_max = 5, B = c(0.3, -0.25), Btau = 0, xi = 1, sb = c(1, 0), rho = 0,
  sc = 0, se = 0.5, gf = gfun0, lHi = function(x) lHwi(x, 2, 10), Tifun = identity,
  addcp = FALSE
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  Ci <- drawC(N)
  Zi <- rexp(N)
  Si <- lHi(Zi)
  Ti_ <- pmin(Si, Ci)
  del <- as.integer(Si <= Ci)
  b <- get_bIS(N, sb, rho, sc)
  if (!addcp) {
    b$mC[, 2] <- Inf
  }
  Xdat <- lapply(1:N, function(j) {
    my_fun <- mt_IScp_i(B, Btau, b$m[j, ], b$mC[j, ])
    tibble(i = j, tim = 0:floor(t_max), y = my_fun(tim))
  }) %>%
    bind_rows() %>%
    mutate(y = y + rnorm(n(), sd = se))
  Sout <- tibble(i = 1:N, Ti = Si, del = del)
  Xdat <- Xdat %>%
    left_join(Sout, by = "i") %>%
    mutate(y = y + alp * Tifun(Ti) +
      gf(Btau + b$mC[i, 1], xi, Ti, tim)) %>%
    mutate(Ti = Ti_[i], i = factor(
      i,
      sort(unique(i))
    )) %>%
    filter(tim <= Ti)
  list(X = Xdat, S = mutate(Sout, Ti = Ti_, i = factor(i, sort(unique(i)))))
}

#' Simulate a log-time conditional-model screening dataset
#'
#' Convenience wrapper around [CM_sim_RI()] using `Tifun = log`, so the event
#' time coefficient multiplies `log(Ti)`.
#' @family simulation generators
#' @param ... Arguments passed to [CM_sim_RI()].
#' @return A list with longitudinal rows `X` and survival rows `S`.
#' @export
CM_sim_RI_logT <- function(...) {
  CM_sim_RI(..., Tifun = log)
}

#' Generate conditional-model biomarkers given event times
#'
#' Advanced helper used by [drawY100()] to generate conditional-model biomarker
#' rows after event and censoring times are already available.
#' @family high-dimensional simulation
#' @param alp Association parameter.
#' @param B Fixed longitudinal intercept and slope coefficients.
#' @param Btau Fixed changepoint or acute-effect coefficient.
#' @param se Longitudinal residual standard deviation.
#' @param xi Conditional-model kernel decay parameter.
#' @param b Subject-specific random intercept and slope values.
#' @param Ti Observed or latent event times.
#' @param Ci Censoring time.
#' @param t_max Maximum longitudinal visit time.
#' @param gf Acute-effect kernel function.
#' @param Tifun Transformation applied to event time.
#' @param N Number of subjects to simulate.
#' @param sb Random-effect standard deviations.
#' @param rho Random-effect correlation parameter.
#' @return Longitudinal data frame with subject id, visit time, event time, and
#'   simulated biomarker value `y`.
#' @export
CM_givT_IS <- function(
  alp = 0, B = c(0.3, -0.25), Btau = 0, se = 0.5, xi = 1, b = NULL,
  Ti = NULL, Ci = NULL, t_max = 5, gf = gfun0, Tifun = identity, N = 600, sb = c(
    1,
    0
  ), rho = 0
) {
  if (is.null(b)) {
    b <- matrix(rnorm(N * 2), N) %*% chol_IS2(sb, rho)
  }
  N <- nrow(b)
  if (is.null(Ti)) {
    Ti <- rexp(N)
  }
  if (is.null(Ci)) {
    Ci <- Ti
  }
  Tobsi <- pmin(Ci, Ti)
  Xdat <- expand_grid(i = 1:N, tim = 0:floor(t_max)) %>%
    mutate(
      Ti = Ti[i], Tobs = Tobsi[i],
      y = (b[i, 1] + B[1]) + (b[i, 2] + B[2]) * tim + Tifun(Ti) * alp + gf(
        Btau,
        xi, Ti, tim
      ) + rnorm(n(), sd = se)
    ) %>%
    filter(tim <= Tobs)
  Xdat
}

#' Generate log-time conditional-model biomarkers given event times
#'
#' Convenience wrapper around [CM_givT_IS()] using `Tifun = log`, so the event
#' time coefficient multiplies `log(Ti)`.
#' @family high-dimensional simulation
#' @param ... Arguments passed to [CM_givT_IS()].
#' @return Longitudinal data frame with subject id, visit time, event time, and
#'   simulated biomarker value `y`.
#' @export
CM_givT_IS_logT <- function(...) {
  CM_givT_IS(..., Tifun = log)
}
