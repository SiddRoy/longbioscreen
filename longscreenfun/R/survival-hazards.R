#' Evaluate a Weibull log hazard
#'
#' Advanced survival helper used by the joint-model simulators.
#' @family simulation generators
#' @param x Evaluation time.
#' @param k Weibull shape parameter.
#' @param l Weibull scale parameter.
#' @return Numeric log hazard values.
#' @export
lh_w <- function(x, k, l) {
  (k - 1) * log(x) + log(k) - k * log(l)
}

#' Evaluate a Weibull log cumulative hazard
#'
#' Advanced survival helper used by Cox-style simulators.
#' @family simulation generators
#' @param x Evaluation time.
#' @param k Weibull shape parameter.
#' @param l Weibull scale parameter.
#' @return Numeric log cumulative hazard values.
#' @export
lHw <- function(x, k, l) k * log(x / l)

#' Invert a Weibull cumulative hazard draw
#'
#' Advanced survival helper used by conditional-model simulation.
#' @family simulation generators
#' @param x Cumulative hazard draw.
#' @param k Weibull shape parameter.
#' @param l Weibull scale parameter.
#' @return Event-time draw implied by the Weibull cumulative hazard.
#' @export
lHwi <- function(x, k, l) {
  qweibull(-x, k, l, lower.tail = FALSE, log.p = TRUE)
}

#' Integrate joint-model hazard before a changepoint
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @family internal numerical helpers
#' @param tim Observation or evaluation times.
#' @param mt Function for the time-varying log-hazard component.
#' @param lh Log baseline hazard function.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
Lam_IS_lt <- function(tim, mt, lh = function(x) lh_w(x, 2, 10)) {
  m0 <- mt(0)
  sapply(tim, function(ts) {
    ms <- mt(ts)
    M <- max(c(m0, ms))
    log(integrate(function(x) exp(lh(x) + mt(x) - M), 0, ts)$value) + M
  })
}

#' Integrate joint-model hazard after a changepoint
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @family internal numerical helpers
#' @param tim Observation or evaluation times.
#' @param mt Function for the time-varying log-hazard component.
#' @param lh Log baseline hazard function.
#' @param taui Subject-specific changepoint time.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
Lam_IScpp_ut <- function(tim, mt, lh = function(x) lh_w(x, 2, 10), taui = 8) {
  Lam_tau <- Lam_IS_lt(taui, mt, lh)
  mtau <- mt(taui)
  Lam_taut <- sapply(tim, function(ts) {
    ms <- mt(ts)
    M <- max(c(mtau, ms))
    log(integrate(function(x) exp(lh(x) + mt(x) - M), taui, ts)$value) + M
  })
  M <- pmax(Lam_tau, Lam_taut)
  m <- pmin(Lam_tau, Lam_taut)
  M + log1p_exp(m - M)
}

#' Integrate joint-model hazard across a changepoint
#'
#' Integrates a time-varying joint-model hazard before and after a subject's
#' changepoint.
#' @family simulation generators
#' @param tim Observation or evaluation times.
#' @param mt Function for the time-varying log-hazard component.
#' @param lh Log baseline hazard function.
#' @param taui Subject-specific changepoint time.
#' @return Numeric log cumulative hazard values.
#' @export
Lam_mt_IScp <- function(tim, mt, lh = function(x) lh_w(x, 2, 10), taui = 8) {
  out <- rep(NA, length(tim))
  w1 <- which(tim <= taui)
  if (length(w1) > 0) {
    out[w1] <- Lam_IS_lt(tim[w1], mt, lh)
  }
  w2 <- which(tim > taui)
  if (length(w2) > 0) {
    out[w2] <- Lam_IScpp_ut(tim[w2], mt, lh, taui)
  }
  out
}

#' Draw a survival time for a joint-model trajectory
#'
#' Draws or inverts one event time for a subject-specific joint-model
#' trajectory.
#' @family simulation generators
#' @param Zi Log cumulative hazard draw used to invert event time.
#' @param Ci Censoring time.
#' @param Lam_mt_IScp Function that evaluates integrated joint-model hazard.
#' @param mt Function for the time-varying log-hazard component.
#' @param lh Log baseline hazard function.
#' @param taui Subject-specific changepoint time.
#' @return Numeric vector with observed/event time, event indicator, log
#'   cumulative hazard at censoring, and censoring time.
#' @export
survY_jm <- function(Zi = NULL, Ci, Lam_mt_IScp, mt, lh = function(x) {
                       lh_w(
                         x, 2,
                         10
                       )
                     }, taui) {
  Lam_ft <- Lam_mt_IScp(Ci, mt, lh, taui)
  if (is.null(Zi)) {
    Zi <- log(rexp(1))
  }
  deli <- Zi <= Lam_ft
  if (deli == 0) {
    return(c(Ci, deli, Lam_ft, Ci))
  }
  L <- 0
  R <- Ci
  if (taui < Ci) {
    Lam_tau <- Lam_mt_IScp(taui, mt, lh, taui)
    if (Zi <= Lam_tau) {
      R <- taui
    } else {
      L <- taui
    }
  }
  Ti <- tryCatch(uniroot(function(x) {
    Zi - Lam_mt_IScp(x, mt, lh, taui)
  }, c(L, R))$root, error = function(e) NULL)
  if (is.null(Ti)) {
    return(c(NA, deli, Lam_ft, Ci))
  }
  c(Ti, deli, Lam_ft, Ci)
}

#' Compute Cox cumulative hazard at visit cutpoints
#'
#' Computes baseline and covariate-adjusted cumulative hazards at visit
#' cutpoints for a stepwise Cox covariate process.
#' @family simulation generators
#' @param tim Observation or evaluation times.
#' @param Ytim Time-varying covariate values on the hazard scale.
#' @param lH Log baseline cumulative hazard function.
#' @return A two-column matrix with baseline and adjusted log cumulative hazards.
#' @export
Lam_tc <- function(tim, Ytim, lH = function(x) lHw(x, 2, 10)) {
  Lam0t <- lHw(tim, 2, 10)
  if (length(Ytim) == 1) {
    return(cbind(Lam0t, Lam0t))
  }
  LamXt <- Lam0t
  LamXt[2] <- Ytim[1] + Lam0t[2]
  if (length(Ytim) > 2) {
    for (j in 3:length(Ytim)) {
      LamXt[j] <- LamXt[j - 1] + log1p_exp(Ytim[j - 1] - LamXt[j - 1] + Lam0t[j] +
        log_neg_expm1(Lam0t[j - 1] - Lam0t[j]))
    }
  }
  cbind(Lam0t, LamXt)
}

#' Compute Cox cumulative hazard at an arbitrary time
#'
#' Evaluates the stepwise Cox cumulative hazard between visit cutpoints.
#' @family simulation generators
#' @param ts Evaluation time.
#' @param tim Observation or evaluation times.
#' @param Ytim Time-varying covariate values on the hazard scale.
#' @param Lam_Xts Cumulative hazard values at visit cutpoints.
#' @param Lam0t Baseline cumulative hazard values at visit cutpoints.
#' @param lH Log baseline cumulative hazard function.
#' @return Numeric log cumulative hazard value at `ts`.
#' @export
LamCox_Xt <- function(ts, tim, Ytim, Lam_Xts, Lam0t, lH = function(x) lHw(x, 2, 10)) {
  Mt <- findInterval(ts, c(tim, Inf), left.open = TRUE, all.inside = TRUE)
  lhs <- Lam_Xts[Mt]
  LamTs <- lH(ts)
  rhs <- Ytim[Mt] + LamTs + log_neg_expm1(Lam0t[Mt] - LamTs)
  M <- pmax(lhs, rhs)
  m <- pmin(lhs, rhs)
  ifelse(ts == 0, -Inf, M + log1p_exp(m - M))
}

#' Draw a survival time for stepwise Cox covariates
#'
#' Draws or inverts one event time under a stepwise Cox covariate process.
#' @family simulation generators
#' @param Zi Log cumulative hazard draw used to invert event time.
#' @param Ci Censoring time.
#' @param tim Observation or evaluation times.
#' @param Ytim Time-varying covariate values on the hazard scale.
#' @param lH Log baseline cumulative hazard function.
#' @return Numeric vector with observed/event time, event indicator, log
#'   cumulative hazard at censoring, and censoring time.
#' @export
surv_Cox <- function(Zi = NULL, Ci, tim, Ytim, lH = function(x) lHw(x, 2, 10)) {
  Lam_cut <- Lam_tc(tim, Ytim, lH)
  Lam_ft <- LamCox_Xt(Ci, tim, Ytim, Lam_cut[, 2], Lam_cut[, 1])
  if (is.null(Zi)) {
    Zi <- log(rexp(1))
  }
  deli <- Zi <= Lam_ft
  if (deli == 0) {
    return(c(Ci, deli, Lam_ft, Ci))
  }
  Mzi <- findInterval(Zi, Lam_cut[, 2], left.open = TRUE)
  lZ <- pmax(0, tim[Mzi] - 0.01)
  uZ <- pmin(c(tim[-1], Ci)[Mzi] + 0.01, Ci)
  Ti <- tryCatch(uniroot(function(x) {
    Zi - LamCox_Xt(x, tim, Ytim, Lam_cut[, 2], Lam_cut[, 1])
  }, c(lZ, uZ))$root, error = function(e) NULL)
  if (is.null(Ti)) {
    return(c(NA, deli, Lam_ft, Ci))
  }
  c(Ti, deli, Lam_ft, Ci)
}

#' Draw survival time for high-dimensional Cox covariates
#'
#' Specialized event-time inverter used by the high-dimensional Cox/JM
#' simulator.
#' @family high-dimensional simulation
#' @param Zi Log cumulative hazard draw used to invert event time.
#' @param Ci Censoring time.
#' @param tim Observation or evaluation times.
#' @param Ytim Time-varying covariate values on the hazard scale.
#' @param lh_kl Weibull shape and scale for the high-dimensional simulator.
#' @return Numeric vector with observed/event time, event indicator, log
#'   cumulative hazard at censoring, and censoring time.
#' @export
surv_Cox_W <- function(Zi = NULL, Ci, tim, Ytim, lh_kl = c(2, 10)) {
  lH <- function(x) lHw(x, lh_kl[1], lh_kl[2])
  Lam_cut <- Lam_tc(tim, Ytim, lH)
  Lam_ft <- LamCox_Xt(Ci, tim, Ytim, Lam_cut[, 2], Lam_cut[, 1])
  if (is.null(Zi)) {
    Zi <- log(rexp(1))
  }
  deli <- Zi <= Lam_ft
  if (deli == 0) {
    return(c(Ci, deli, Lam_ft, Ci))
  }
  iz <- findInterval(Zi, Lam_cut[, 2], left.open = TRUE)
  lxx <- logsum(c(Zi - Ytim[iz] + log_neg_expm1(Lam_cut[iz, 2] - Zi), Lam_cut[
    iz,
    1
  ]))
  Ti <- lh_kl[2] * exp(lxx / lh_kl[1])
  if (is.null(Ti)) {
    return(c(NA, deli, Lam_ft, Ci))
  }
  c(Ti, deli, Lam_ft, Ci)
}
