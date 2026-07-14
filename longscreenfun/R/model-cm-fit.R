#' Fit the Weibull survival model
#'
#' Fits the survival-only Weibull model used before conditional-model
#' screening.
#' @family model fitting helpers
#' @param x0 Initial parameter vector.
#' @param Ti Observed or latent event times.
#' @param del Event indicators.
#' @return A list with parameter estimates, subject score contributions,
#'   numerical Hessian, and Weibull distribution helper functions.
#' @export
fitWeibS <- function(x0 = c(log(2), log(10)), Ti, del) {
  fit_w <- optim(x0, function(x) surv_dev(x, Ti, del), method = "BFGS", hessian = FALSE)
  g_w <- g_surv_dev(fit_w$par, Ti, del)
  g_f <- function(th) colSums(g_surv_dev(th, Ti, del))
  H_ij <- g_xl(fit_w$par, g_f)
  list(est = fit_w$par, si = g_w, h = H_ij, Sfun = dfun_kl(fit_w$par))
}

#' Prepare conditional-model complete and censored data
#'
#' Splits longitudinal rows into complete-case and censored-subject pieces used
#' by the conditional-model likelihood.
#' @family model fitting helpers
#' @param X Longitudinal data frame.
#' @param S Survival data frame.
#' @param lxi0 Fixed log xi value.
#' @param addii Whether to add row indices for biomarker matrix lookup.
#' @param id_ Subject id column name.
#' @return A prepared data list with complete rows, censored rows, subject
#'   indexes, censored-row list, and survival-score match indexes.
#' @export
xifix_subX <- function(X, S, lxi0 = 0, addii = FALSE, id_ = "i") {
  if (!all(colnames(S) %in% colnames(X))) {
    X <- left_join(X, S, by = id_)
  }
  if (addii) {
    X <- mutate(X, ii = 1:n())
  }
  X_cc0 <- filter(X, del == 1) %>% mutate(tx = tim, g = gfun0(
    1, exp(lxi0), Ti,
    tim
  ))
  X_lc0 <- filter(X, del == 0) %>% rename(tx = tim)
  wcc <- wj_Xcc(X_cc0)
  X_lcl <- Xcens_jdat(mutate(X_lc0, i = as.integer(as.character(i))))
  i_s0_cc <- unique(as.integer(as.character(X_cc0$i)))
  i_s0_cens <- as.integer((as.character(unique(X_lc0$i))))
  i_s0 <- c(i_s0_cc, i_s0_cens)
  i_gw <- as.integer(as.character(S$i))
  m_s0_gw <- match(i_s0, i_gw)
  list(Xcc = X_cc0, Xrc = X_lc0, wcc = wcc, X_lcl = X_lcl, m_YS = m_s0_gw)
}

#' Build conditional-model deviance functions with fixed xi
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @family internal numerical helpers
#' @param lxi0 Fixed log xi value.
#' @param X_cc0 Complete-case longitudinal rows at fixed xi.
#' @param X_lcl List of censored-subject longitudinal rows.
#' @param wcc Complete-case row index list.
#' @param Sfun List of Weibull distribution helper functions.
#' @param lc Censored-subject integration function.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
cM_devfun_xifix <- function(lxi0 = 0, X_cc0, X_lcl, wcc, Sfun, lc = lc_fun()) {
  f_devXi <- function(th) {
    -l_full_wj(c(lxi0, th), X_cc0, X_lcl, wcc, Sfun = Sfun, lc = lc)
  }
  f_dev_lXi <- function(th) {
    -l_full_wj(c(lxi0, th), X_cc0, X_lcl, wcc, Sfun = Sfun, lc = lc, retal = TRUE)
  }
  S_i <- function(th, ...) {
    g_xl(th, f_dev_lXi, ...)
  }
  f_dev_dXi <- function(th, log_kl) {
    Sfunk <- dfun_kl(log_kl)
    -l_full_wj(c(lxi0, th), X_cc0, X_lcl, wcc, Sfun = Sfunk, lc = lc)
  }
  Lprim_ij <- function(th, thS, ...) {
    H_ij(th, thS, f_dev_dXi, ...)
  }
  list(dev = f_devXi, S_ifun = S_i, Jijfun = Lprim_ij)
}

#' Build log-time conditional-model deviance functions with fixed xi
#'
#' Log-time variant of [cM_devfun_xifix()] where the third fixed-effect
#' coefficient multiplies `log(Ti)` instead of `Ti`.
#' @family internal numerical helpers
#' @param lxi0 Fixed log xi value.
#' @param X_cc0 Complete-case longitudinal rows at fixed xi.
#' @param X_lcl List of censored-subject longitudinal rows.
#' @param wcc Complete-case row index list.
#' @param Sfun List of Weibull distribution helper functions.
#' @param lc Censored-subject integration function.
#' @return A list of functions used by [fixxi_fitcm0()] and
#'   [fixxi_fitcm_noB()].
#' @export
cM_devfun_xifix_logT <- function(lxi0 = 0, X_cc0, X_lcl, wcc, Sfun, lc = lc_fun_logT()) {
  f_devXi <- function(th) {
    -l_full_wj_logT(c(lxi0, th), X_cc0, X_lcl, wcc, Sfun = Sfun, lc = lc)
  }
  f_dev_lXi <- function(th) {
    -l_full_wj_logT(c(lxi0, th), X_cc0, X_lcl, wcc, Sfun = Sfun, lc = lc, retal = TRUE)
  }
  S_i <- function(th, ...) {
    g_xl(th, f_dev_lXi, ...)
  }
  f_dev_dXi <- function(th, log_kl) {
    Sfunk <- dfun_kl(log_kl)
    -l_full_wj_logT(c(lxi0, th), X_cc0, X_lcl, wcc, Sfun = Sfunk, lc = lc)
  }
  Lprim_ij <- function(th, thS, ...) {
    H_ij(th, thS, f_dev_dXi, ...)
  }
  list(dev = f_devXi, S_ifun = S_i, Jijfun = Lprim_ij)
}

#' Build conditional-model deviance functions with estimated xi
#'
#' @family internal numerical helpers
#' @param X_cc0 Complete-case longitudinal rows. Any precomputed `g` column is
#'   ignored so the likelihood recomputes it from the current `xi`.
#' @param X_lcl List of censored-subject longitudinal rows.
#' @param wcc Complete-case row index list.
#' @param Sfun List of Weibull distribution helper functions.
#' @param lc Censored-subject integration function.
#' @return A list of deviance, score, and cross-derivative functions.
#' @keywords internal
cM_devfun_xifit <- function(X_cc0, X_lcl, wcc, Sfun, lc = lc_fun()) {
  X_cc_fit <- X_cc0[, setdiff(names(X_cc0), "g"), drop = FALSE]
  f_devXi <- function(th) {
    -l_full_wj(th, X_cc_fit, X_lcl, wcc, Sfun = Sfun, lc = lc)
  }
  f_dev_lXi <- function(th) {
    -l_full_wj(th, X_cc_fit, X_lcl, wcc, Sfun = Sfun, lc = lc, retal = TRUE)
  }
  S_i <- function(th, ...) {
    g_xl(th, f_dev_lXi, ...)
  }
  f_dev_dXi <- function(th, log_kl) {
    Sfunk <- dfun_kl(log_kl)
    -l_full_wj(th, X_cc_fit, X_lcl, wcc, Sfun = Sfunk, lc = lc)
  }
  Lprim_ij <- function(th, thS, ...) {
    H_ij(th, thS, f_dev_dXi, ...)
  }
  list(dev = f_devXi, S_ifun = S_i, Jijfun = Lprim_ij)
}

#' Build log-time conditional-model deviance functions with estimated xi
#'
#' Log-time variant of `cM_devfun_xifit()` where the third fixed-effect
#' coefficient multiplies `log(Ti)` instead of `Ti`.
#' @family internal numerical helpers
#' @param X_cc0 Complete-case longitudinal rows. Any precomputed `g` column is
#'   ignored so the likelihood recomputes it from the current `xi`.
#' @param X_lcl List of censored-subject longitudinal rows.
#' @param wcc Complete-case row index list.
#' @param Sfun List of Weibull distribution helper functions.
#' @param lc Censored-subject integration function.
#' @return A list of deviance, score, and cross-derivative functions.
#' @keywords internal
cM_devfun_xifit_logT <- function(X_cc0, X_lcl, wcc, Sfun, lc = lc_fun_logT()) {
  X_cc_fit <- X_cc0[, setdiff(names(X_cc0), "g"), drop = FALSE]
  f_devXi <- function(th) {
    -l_full_wj_logT(th, X_cc_fit, X_lcl, wcc, Sfun = Sfun, lc = lc)
  }
  f_dev_lXi <- function(th) {
    -l_full_wj_logT(th, X_cc_fit, X_lcl, wcc, Sfun = Sfun, lc = lc, retal = TRUE)
  }
  S_i <- function(th, ...) {
    g_xl(th, f_dev_lXi, ...)
  }
  f_dev_dXi <- function(th, log_kl) {
    Sfunk <- dfun_kl(log_kl)
    -l_full_wj_logT(th, X_cc_fit, X_lcl, wcc, Sfun = Sfunk, lc = lc)
  }
  Lprim_ij <- function(th, thS, ...) {
    H_ij(th, thS, f_dev_dXi, ...)
  }
  list(dev = f_devXi, S_ifun = S_i, Jijfun = Lprim_ij)
}

#' Build fixed-xi deviance functions from prepared data
#'
#' Builds conditional-model deviance, score, and cross-derivative functions from
#' a prepared data list.
#' @family model fitting helpers
#' @param lxi0 Fixed log xi value.
#' @param datl Prepared conditional-model data list.
#' @param Sfun List of Weibull distribution helper functions.
#' @param lc Censored-subject integration function.
#' @return A list of functions used by [fixxi_fitcm0()] and
#'   [fixxi_fitcm_noB()].
#' @export
cM_devfun_xifix_datl <- function(lxi0 = 0, datl, Sfun, lc = lc_fun()) {
  cM_devfun_xifix(lxi0, datl$Xcc, datl$X_lcl, datl$wcc, Sfun = Sfun, lc = lc)
}

#' Build fixed-xi log-time deviance functions from prepared data
#'
#' Log-time variant of [cM_devfun_xifix_datl()].
#' @family model fitting helpers
#' @param lxi0 Fixed log xi value.
#' @param datl Prepared conditional-model data list.
#' @param Sfun List of Weibull distribution helper functions.
#' @param lc Censored-subject integration function.
#' @return A list of functions used by [fixxi_fitcm0()] and
#'   [fixxi_fitcm_noB()].
#' @export
cM_devfun_xifix_datl_logT <- function(lxi0 = 0, datl, Sfun, lc = lc_fun_logT()) {
  cM_devfun_xifix_logT(lxi0, datl$Xcc, datl$X_lcl, datl$wcc, Sfun = Sfun, lc = lc)
}

#' Build estimated-xi deviance functions from prepared data
#'
#' Builds conditional-model deviance, score, and cross-derivative functions from
#' a prepared data list while estimating `xi`.
#' @family internal numerical helpers
#' @param datl Prepared conditional-model data list.
#' @param Sfun List of Weibull distribution helper functions.
#' @param lc Censored-subject integration function.
#' @return A list of deviance, score, and cross-derivative functions.
#' @keywords internal
cM_devfun_xifit_datl <- function(datl, Sfun, lc = lc_fun()) {
  cM_devfun_xifit(datl$Xcc, datl$X_lcl, datl$wcc, Sfun = Sfun, lc = lc)
}

#' Build estimated-xi log-time deviance functions from prepared data
#'
#' Log-time variant of `cM_devfun_xifit_datl()`.
#' @family internal numerical helpers
#' @param datl Prepared conditional-model data list.
#' @param Sfun List of Weibull distribution helper functions.
#' @param lc Censored-subject integration function.
#' @return A list of deviance, score, and cross-derivative functions.
#' @keywords internal
cM_devfun_xifit_datl_logT <- function(datl, Sfun, lc = lc_fun_logT()) {
  cM_devfun_xifit_logT(datl$Xcc, datl$X_lcl, datl$wcc, Sfun = Sfun, lc = lc)
}

#' Initialize the fixed-xi conditional-model fit
#'
#' Builds an initial parameter vector for the fixed-xi conditional model.
#' @family model fitting helpers
#' @param datl Prepared conditional-model data list.
#' @param lxi0 Fixed log xi value.
#' @return Numeric initial parameter vector.
#' @export
init_fixxi <- function(datl, lxi0 = 0) {
  dat0 <- bind_rows(mutate(datl$Xcc, g = gfun(1, exp(lxi0), Ti, tim)), mutate(datl$Xrc,
    g = 0, Ti = 0
  ))
  fit0_cmy <- lmer(y ~ tim + (1 | i) + Ti + g, data = dat0)
  as.numeric(c(fixef(fit0_cmy), log(sigma(fit0_cmy)) + c(
    log(getME(fit0_cmy, "theta")),
    0
  )))
}

#' Initialize the fixed-xi log-time conditional-model fit
#'
#' Log-time variant of [init_fixxi()] where the third fixed-effect coefficient
#' is initialized from `log(Ti)` for complete cases. Censored starting rows use
#' `Ti_log = 0`, matching the current initialization convention for unobserved
#' event-time contributions.
#' @family model fitting helpers
#' @param datl Prepared conditional-model data list.
#' @param lxi0 Fixed log xi value.
#' @return Numeric initial parameter vector.
#' @export
init_fixxi_logT <- function(datl, lxi0 = 0) {
  dat0 <- bind_rows(
    mutate(datl$Xcc, g = gfun(1, exp(lxi0), Ti, tim), Ti_log = log(Ti)),
    mutate(datl$Xrc, g = 0, Ti_log = 0)
  )
  fit0_cmy <- lmer(y ~ tim + (1 | i) + Ti_log + g, data = dat0)
  as.numeric(c(fixef(fit0_cmy), log(sigma(fit0_cmy)) + c(
    log(getME(fit0_cmy, "theta")),
    0
  )))
}

#' Fit the fixed-xi conditional model
#'
#' Optimizes the fixed-xi conditional-model deviance and computes derivative
#' quantities needed for inference.
#' @family model fitting helpers
#' @param devF Conditional-model deviance function bundle.
#' @param fit_w Fitted Weibull survival object.
#' @param b0 Initial conditional-model parameter vector.
#' @param r Richardson depth or residual vector, depending on context.
#' @return An `optim()` fit object augmented with Hessian, cross-derivative, and
#'   score components.
#' @export
fixxi_fitcm0 <- function(devF, fit_w, b0, r = 2) {
  fit1 <- optim(b0, devF$dev, method = "BFGS", hessian = FALSE)
  fit1$hessian <- H_ii(fit1$par, devF$dev, r = r)
  fit1$Jij <- devF$Jijfun(fit1$par, fit_w$est, r = r)
  fit1$s0 <- devF$S_ifun(fit1$par, r = r)
  fit1
}

#' Compute an inverse-Hessian covariance matrix
#'
#' @param hessian Numerical Hessian matrix.
#' @return Inverse Hessian computed with Cholesky backsolves.
#' @keywords internal
cm_inverse_hessian_cov <- function(hessian) {
  RV <- chol(hessian)
  backsolve(RV, t(backsolve(RV, diag(nrow(RV)))))
}

#' Compute a sandwich covariance matrix from subject-level scores
#'
#' @param hessian Numerical Hessian matrix.
#' @param score_i Subject-level score matrix, one subject per row.
#' @return Sandwich covariance matrix.
#' @keywords internal
cm_sandwich_cov <- function(hessian, score_i) {
  score_i <- as.matrix(score_i)
  if (ncol(score_i) != nrow(hessian) && nrow(score_i) == nrow(hessian)) {
    score_i <- t(score_i)
  }
  if (ncol(score_i) != nrow(hessian)) {
    stop("`score_i` must have one column per Hessian parameter.", call. = FALSE)
  }
  RV <- chol(hessian)
  Hinv_score <- backsolve(RV, backsolve(RV, t(score_i), transpose = TRUE))
  tcrossprod(Hinv_score)
}

#' Format conditional-model association inference
#'
#' @param par Optimized conditional-model parameter vector.
#' @param Vyth Variance-covariance matrix.
#' @param assoc_idx Indexes of the two association parameters.
#' @return Numeric vector with association estimates, standard errors, and log
#'   p-values.
#' @keywords internal
cm_assoc_inference <- function(par, Vyth, assoc_idx) {
  se <- sqrt(diag(Vyth)[assoc_idx])
  p <- c(pfun_bse_log(par[assoc_idx], se), pfun_chi2(par[assoc_idx], Vyth[assoc_idx, assoc_idx],
    log.p = TRUE
  ))
  c(par[assoc_idx], se, p)
}

#' Ensure conditional-model derivative quantities are available
#'
#' @param fit1 Optional first-stage conditional-model fit.
#' @param devF Conditional-model deviance function bundle.
#' @param fit_w Fitted Weibull survival object.
#' @param b0 Initial conditional-model parameter vector.
#' @param r Richardson depth.
#' @return An `optim()` fit object augmented with Hessian, cross-derivative, and
#'   score components.
#' @keywords internal
cm_ensure_fit_derivatives <- function(fit1, devF, fit_w, b0, r = 2) {
  if (is.null(fit1)) {
    fit1 <- optim(b0, devF$dev, method = "BFGS", hessian = FALSE)
  }
  if (is.null(fit1$hessian)) {
    fit1$hessian <- H_ii(fit1$par, devF$dev, r = r)
  }
  if (is.null(fit1$Jij)) {
    fit1$Jij <- devF$Jijfun(fit1$par, fit_w$est, r = r)
  }
  if (is.null(fit1$s0)) {
    fit1$s0 <- devF$S_ifun(fit1$par, r = r)
  }
  fit1
}

#' Compute full-data adjusted conditional-model scores
#'
#' @param fit1 First-stage conditional-model fit with derivative components.
#' @param fit_w Fitted Weibull survival object.
#' @param m_YS Match index between longitudinal and survival score rows.
#' @return Subject-level adjusted score matrix.
#' @keywords internal
cm_full_adjusted_scores <- function(fit1, fit_w, m_YS) {
  RVs <- chol(fit_w$h)
  fit1$s0 -
    t(fit1$Jij %*% backsolve(
      RVs, backsolve(RVs, t(fit_w$si[m_YS, , drop = FALSE]),
        transpose = TRUE
      )
    ))
}

#' Compute fixed-xi conditional-model inference
#'
#' Computes standard errors and p-values for the fixed-xi conditional-model
#' association parameters.
#' @family model fitting helpers
#' @param devF Conditional-model deviance function bundle.
#' @param fit_w Fitted Weibull survival object.
#' @param b0 Initial conditional-model parameter vector.
#' @param m_YS Match index between longitudinal and survival score rows.
#' @param fit1 Optional first-stage conditional-model fit.
#' @param r Richardson depth or residual vector, depending on context.
#' @return Numeric vector with association estimates, standard errors, and log
#'   p-values.
#' @export
fixxi_fitcm_noB <- function(devF, fit_w, b0, m_YS, fit1 = NULL, r = 2) {
  if (is.null(fit1)) {
    fit1 <- optim(b0, devF$dev, method = "BFGS", hessian = FALSE)
    fit1$hessian <- H_ii(fit1$par, devF$dev, r = r)
    fit1$Jij <- devF$Jijfun(fit1$par, fit_w$est, r = r)
  }
  RV <- chol(fit1$hessian)
  RVs <- chol(fit_w$h)
  B <- diag(nrow(RV)) + crossprod(backsolve(RVs, t(backsolve(RV, fit1$Jij, transpose = TRUE)),
    transpose = TRUE
  ))
  Vyth <- backsolve(RV, t(backsolve(RV, B)))
  se <- sqrt(diag(Vyth)[3:4])
  p <- c(pfun_bse_log(fit1$par[3:4], se), pfun_chi2(fit1$par[3:4], Vyth[3:4, 3:4],
    log.p = TRUE
  ))
  c(fit1$par[3:4], se, p)
}

#' Compute fixed-xi sandwich conditional-model inference
#'
#' Computes standard errors and p-values for the fixed-xi conditional-model
#' association parameters using subject-level adjusted scores.
#' @family internal numerical helpers
#' @param devF Conditional-model deviance function bundle.
#' @param fit_w Fitted Weibull survival object.
#' @param b0 Initial conditional-model parameter vector.
#' @param m_YS Match index between longitudinal and survival score rows.
#' @param fit1 Optional first-stage conditional-model fit.
#' @param r Richardson depth.
#' @return Numeric vector with association estimates, standard errors, and log
#'   p-values.
#' @keywords internal
fixxi_fitcm_sandwich <- function(devF, fit_w, b0, m_YS, fit1 = NULL, r = 2) {
  fit1 <- cm_ensure_fit_derivatives(fit1, devF, fit_w, b0, r = r)
  gy_i <- cm_full_adjusted_scores(fit1, fit_w, m_YS)
  Vyth <- cm_sandwich_cov(fit1$hessian, gy_i)
  cm_assoc_inference(fit1$par, Vyth, 3:4)
}

#' Fit the estimated-xi conditional model
#'
#' @family internal numerical helpers
#' @param devF Conditional-model deviance function bundle.
#' @param fit_w Fitted Weibull survival object.
#' @param b0 Initial conditional-model parameter vector.
#' @param r Richardson depth or residual vector, depending on context.
#' @return An `optim()` fit object augmented with Hessian, cross-derivative, and
#'   score components.
#' @keywords internal
fitxi_fitcm0 <- function(devF, fit_w, b0, r = 2) {
  fit1 <- optim(b0, devF$dev, method = "BFGS", hessian = FALSE)
  fit1$hessian <- H_ii(fit1$par, devF$dev, r = r)
  fit1$Jij <- devF$Jijfun(fit1$par, fit_w$est, r = r)
  fit1$s0 <- devF$S_ifun(fit1$par, r = r)
  fit1
}

#' Compute estimated-xi conditional-model inference
#'
#' Computes standard errors and p-values for the conditional-model association
#' parameters while treating `xi` as an estimated parameter.
#' @family internal numerical helpers
#' @param devF Conditional-model deviance function bundle.
#' @param fit_w Fitted Weibull survival object.
#' @param b0 Initial conditional-model parameter vector.
#' @param m_YS Match index between longitudinal and survival score rows.
#' @param fit1 Optional first-stage conditional-model fit.
#' @param r Richardson depth or residual vector, depending on context.
#' @return Numeric vector with association estimates, standard errors, and log
#'   p-values.
#' @keywords internal
fitxi_fitcm_noB <- function(devF, fit_w, b0, m_YS, fit1 = NULL, r = 2) {
  if (is.null(fit1)) {
    fit1 <- optim(b0, devF$dev, method = "BFGS", hessian = FALSE)
    fit1$hessian <- H_ii(fit1$par, devF$dev, r = r)
    fit1$Jij <- devF$Jijfun(fit1$par, fit_w$est, r = r)
  }
  RV <- chol(fit1$hessian)
  RVs <- chol(fit_w$h)
  B <- diag(nrow(RV)) + crossprod(backsolve(RVs, t(backsolve(RV, fit1$Jij, transpose = TRUE)),
    transpose = TRUE
  ))
  Vyth <- backsolve(RV, t(backsolve(RV, B)))
  assoc_idx <- 4:5
  se <- sqrt(diag(Vyth)[assoc_idx])
  p <- c(pfun_bse_log(fit1$par[assoc_idx], se), pfun_chi2(fit1$par[assoc_idx], Vyth[assoc_idx, assoc_idx],
    log.p = TRUE
  ))
  c(fit1$par[assoc_idx], se, p)
}

#' Compute estimated-xi sandwich conditional-model inference
#'
#' Computes standard errors and p-values for the conditional-model association
#' parameters while treating `xi` as an estimated parameter and using
#' subject-level adjusted scores.
#' @family internal numerical helpers
#' @param devF Conditional-model deviance function bundle.
#' @param fit_w Fitted Weibull survival object.
#' @param b0 Initial conditional-model parameter vector.
#' @param m_YS Match index between longitudinal and survival score rows.
#' @param fit1 Optional first-stage conditional-model fit.
#' @param r Richardson depth.
#' @return Numeric vector with association estimates, standard errors, and log
#'   p-values.
#' @keywords internal
fitxi_fitcm_sandwich <- function(devF, fit_w, b0, m_YS, fit1 = NULL, r = 2) {
  fit1 <- cm_ensure_fit_derivatives(fit1, devF, fit_w, b0, r = r)
  gy_i <- cm_full_adjusted_scores(fit1, fit_w, m_YS)
  Vyth <- cm_sandwich_cov(fit1$hessian, gy_i)
  cm_assoc_inference(fit1$par, Vyth, 4:5)
}
