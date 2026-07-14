#' Evaluate the conditional-model kernel
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param b1i Acute-effect coefficient.
#' @param xi Conditional-model kernel decay parameter.
#' @param S Survival data frame or event time.
#' @param tx Longitudinal visit time.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
gfun <- function(b1i, xi, S, tx) {
  b1i * exp(xi * (tx - S))
}

#' Compute Weibull survival negative log likelihood
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param th Parameter vector.
#' @param TT Event or censoring times.
#' @param deli Event indicators.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
surv_dev <- function(th, TT, deli) {
  kl <- exp(th)
  out <- sum(dweibull(TT[deli == 1], shape = kl[1], scale = kl[2], log = TRUE)) +
    sum(pweibull(TT[deli != 1], shape = kl[1], scale = kl[2], log.p = TRUE, lower.tail = FALSE))
  -out
}

#' Compute Weibull score contributions
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param th Parameter vector.
#' @param TT Event or censoring times.
#' @param deli Event indicators.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
g_surv_dev <- function(th, TT, deli) {
  kl <- exp(th)
  out <- (TT / kl[2])^(kl[1]) * cbind(-log(TT / kl[2]), kl[1] / kl[2])
  out[, 2] <- -kl[2] * (out[, 2] - (deli == 1) * kl[1] / kl[2])
  out[, 1] <- -kl[1] * (out[, 1] + (deli == 1) * (1 / kl[1] + log(TT / kl[2])))
  out
}

#' Create Weibull distribution helper functions
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param kl_log Log Weibull shape and scale parameters.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
dfun_kl <- function(kl_log) {
  k <- exp(kl_log[1])
  l <- exp(kl_log[2])
  ldkl <- function(Si) dweibull(Si, k, l, log = TRUE)
  qkl <- function(u) qweibull(u, k, l)
  pkl <- function(Si) pweibull(Si, k, l)
  lS <- function(Si) pweibull(Si, k, l, log.p = TRUE, lower.tail = FALSE)
  list(ldS = ldkl, qS = qkl, pS = pkl, lS = lS)
}

#' Fifteen-point Gauss-Kronrod quadrature rule
#'
#' Returns the same nodes and weights used by the unexported
#' `JM:::gaussKronrod(15)` helper. These constants were copied from
#' `JM:::gaussKronrod()` on June 24, 2026, to avoid depending on an unexported
#' function from another package.
#' @return A list with `sk` nodes and `wk` weights for the 15-point rule.
#' @noRd
gk15_rule <- function() {
  list(
    sk = c(
      -0x1.e5f178e7c622bp-1, -0x1.7ba9f9be3a1d6p-1, -0x1.9f95df119fd65p-2,
      0x0p+0, 0x1.9f95df119fd65p-2, 0x1.7ba9f9be3a1d6p-1,
      0x1.e5f178e7c622bp-1, -0x1.fba009d4d09b5p-1, -0x1.bacf827b9bb3ep-1,
      -0x1.2c13a049dfa26p-1, -0x1.a98b2892e0c78p-3, 0x1.a98b2892e0c78p-3,
      0x1.2c13a049dfa26p-1, 0x1.bacf827b9bb3ep-1, 0x1.fba009d4d09b5p-1
    ),
    wk = c(
      0x1.026cdaa7b61c4p-4, 0x1.200ed0f46e8c1p-3, 0x1.85d6861c80eb2p-3,
      0x1.ad04f90870911p-3, 0x1.85d6861c80eb2p-3, 0x1.200ed0f46e8c1p-3,
      0x1.026cdaa7b61c4p-4, 0x1.77c5b67d57472p-6, 0x1.ad384a34814c7p-4,
      0x1.5a1f266e47d5dp-3, 0x1.a2adbcbec9cd9p-3, 0x1.a2adbcbec9cd9p-3,
      0x1.5a1f266e47d5dp-3, 0x1.ad384a34814c7p-4, 0x1.77c5b67d57472p-6
    )
  )
}

#' Create a censored-subject integration function
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
lc_fun <- function() {
  rule <- gk15_rule()
  GK15 <- list(lwk = log(rule$wk), sk_u = (rule$sk + 1) / 2)
  function(th0, Xd0, Sfun, l_x = lxx_b0) {
    l_cens_i2(th0, Xd0, Sfun, l_x = lxx_b0, GK15)
  }
}

#' Cached censored-subject integration function
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param th0 Transformed conditional-model parameter list.
#' @param Xd0 One censored subject's longitudinal rows.
#' @param Sfun List of Weibull distribution helper functions.
#' @param l_x Longitudinal log-likelihood function.
#' @return A function used for censored-subject likelihood integration.
#' @keywords internal
lc2 <- lc_fun()

#' Compute centered sum of squares
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param x Numeric vector or parameter value.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
sse <- function(x) sum((x - mean(x))^2)

#' Evaluate the random-intercept longitudinal log likelihood
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param r Richardson depth or residual vector, depending on context.
#' @param sb Random-effect standard deviations.
#' @param sig Input used by this helper.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
lxx_b0 <- function(r, sb, sig) {
  n <- length(r)
  d0 <- -n * log(sig)
  tau <- (sb / sig)^2
  ssr <- -0.5 / sig^2 * (sum(r * r) + n * tau * sse(r)) / (1 + n * tau)
  d0 - 0.5 * log1p(n * tau) + ssr
}

#' Integrate a log-density with Gauss-Kronrod weights
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param lf Log-density function to integrate.
#' @param l Lower integration bound.
#' @param u Upper integration bound.
#' @param GK15 Gauss-Kronrod weights and nodes.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
Int <- function(lf, l, u, GK15) {
  si <- GK15$sk_u * (u - l) + l
  logsum(lf(si) + GK15$lwk + log((u - l) / 2))
}

#' Evaluate one censored-subject contribution
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param th0 Transformed conditional-model parameter list.
#' @param Xd0 One censored subject's longitudinal rows.
#' @param Sfun List of Weibull distribution helper functions.
#' @param l_x Longitudinal log-likelihood function.
#' @param GK15 Gauss-Kronrod weights and nodes.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
l_cens_i2 <- function(th0, Xd0, Sfun, l_x = lxx_b0, GK15) {
  r <- Xd0$y - as.numeric(cbind(1, Xd0$tx) %*% th0$B[1:2])
  lIt <- function(x) {
    sapply(Sfun$qS(x), function(y) {
      r <- r - y * th0$B[3] - gfun(1, th0$xi, y, Xd0$tx) * th0$B[4]
      l_x(r, th0$sb, th0$sig)
    })
  }
  Ti <- Sfun$pS(Xd0$Ti[1] + c(0, 2))
  a <- Int(lIt, Ti[1], Ti[2], GK15)
  b <- Int(lIt, Ti[2], 1, GK15)
  logsum(c(a, b))
}

#' Evaluate one censored-subject contribution with a log-time term
#'
#' Log-time variant of [l_cens_i2()] where the third fixed-effect coefficient
#' multiplies `log(Ti)` instead of `Ti`.
#' @param th0 Transformed conditional-model parameter list.
#' @param Xd0 One censored subject's longitudinal rows.
#' @param Sfun List of Weibull distribution helper functions.
#' @param l_x Longitudinal log-likelihood function.
#' @param GK15 Gauss-Kronrod weights and nodes.
#' @return A value used by the log-time conditional-model pipeline.
#' @keywords internal
l_cens_i2_logT <- function(th0, Xd0, Sfun, l_x = lxx_b0, GK15) {
  r <- Xd0$y - as.numeric(cbind(1, Xd0$tx) %*% th0$B[1:2])
  lIt <- function(x) {
    sapply(Sfun$qS(x), function(y) {
      r <- r - log(y) * th0$B[3] - gfun(1, th0$xi, y, Xd0$tx) * th0$B[4]
      l_x(r, th0$sb, th0$sig)
    })
  }
  Ti <- Sfun$pS(Xd0$Ti[1] + c(0, 2))
  a <- Int(lIt, Ti[1], Ti[2], GK15)
  b <- Int(lIt, Ti[2], 1, GK15)
  logsum(c(a, b))
}

#' Create a log-time censored-subject integration function
#'
#' @return A function used for censored-subject likelihood integration in the
#'   log-time conditional-model variant.
#' @keywords internal
lc_fun_logT <- function() {
  rule <- gk15_rule()
  GK15 <- list(lwk = log(rule$wk), sk_u = (rule$sk + 1) / 2)
  function(th0, Xd0, Sfun, l_x = lxx_b0) {
    l_cens_i2_logT(th0, Xd0, Sfun, l_x = lxx_b0, GK15)
  }
}

#' Cached log-time censored-subject integration function
#'
#' @keywords internal
lc2_logT <- lc_fun_logT()

#' Transform conditional-model parameters
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param th Parameter vector.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
th_set_cmI <- function(th) {
  list(xi = exp(th[1]), B = th[2:5], sb = exp(th[6]), sig = exp(th[7]))
}

#' Index complete-case rows by subject
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param X_cc Complete-case longitudinal rows.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
wj_Xcc <- function(X_cc) {
  j <- as.integer(X_cc$i)
  lapply(unique(j), function(ji) which(j == ji))
}

#' Split censored rows by subject
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param X_cens Input used by this helper.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
Xcens_jdat <- function(X_cens) {
  lapply(unique(as.integer(X_cens$i)), function(ji) {
    X_cens[which(X_cens$i == ji), ]
  })
}

#' Evaluate censored conditional-model likelihood terms
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param th_all Full conditional-model parameter vector.
#' @param X_censj List of censored-subject longitudinal rows.
#' @param l_x Longitudinal log-likelihood function.
#' @param th_fun Input used by this helper.
#' @param Sfun List of Weibull distribution helper functions.
#' @param lc Censored-subject integration function.
#' @param retal Whether to return subject-level likelihood contributions.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
L_rc_wj <- function(
  th_all, X_censj, l_x = lxx_b0, th_fun = th_set_cmI, Sfun, lc,
  retal = FALSE
) {
  N2 <- length(X_censj)
  l <- rep(NA, N2)
  th_l <- th_fun(th_all)
  for (k in 1:N2) {
    l[k] <- tryCatch(lc(th_l, X_censj[[k]], Sfun, l_x), error = function(e) -Inf)
    if (is.nan(l[k]) | is.na(l[k]) | is.infinite(l[k])) {
      return(-Inf)
    }
  }
  if (any(is.nan(l) | is.na(l) | is.infinite(l))) {
    return(-Inf)
  }
  if (retal) {
    return(l)
  }
  sum(l)
}

#' Evaluate complete-case conditional-model likelihood terms
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param th_all Full conditional-model parameter vector.
#' @param X_cc Complete-case longitudinal rows.
#' @param wj Subject row index list.
#' @param l_x Longitudinal log-likelihood function.
#' @param retal Whether to return subject-level likelihood contributions.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
L_cc_wj <- function(th_all, X_cc, wj = wj_Xcc(Xcc), l_x = lxx_b0, retal = FALSE) {
  xi <- exp(th_all[1])
  X <- cbind(1, X_cc$tx, X_cc$Ti, X_cc$g)
  B_X <- th_all[2:5]
  r <- X_cc$y - as.numeric(X %*% B_X)
  sb <- exp(th_all[6])
  s <- exp(th_all[7])
  li <- sapply(wj, function(ji) {
    lxx_b0(r[ji], sb, s)
  })
  if (any(is.nan(li) | is.na(li) | is.infinite(li))) {
    return(-Inf)
  }
  if (retal) {
    return(li)
  }
  sum(li)
}

#' Evaluate complete-case log-time conditional-model likelihood terms
#'
#' Log-time variant of [L_cc_wj()] where the third fixed-effect coefficient
#' multiplies `log(Ti)` instead of `Ti`.
#' @param th_all Full conditional-model parameter vector.
#' @param X_cc Complete-case longitudinal rows.
#' @param wj Subject row index list.
#' @param l_x Longitudinal log-likelihood function.
#' @param retal Whether to return subject-level likelihood contributions.
#' @return A value used by the log-time conditional-model pipeline.
#' @keywords internal
L_cc_wj_logT <- function(th_all, X_cc, wj = wj_Xcc(Xcc), l_x = lxx_b0, retal = FALSE) {
  xi <- exp(th_all[1])
  X <- cbind(1, X_cc$tx, log(X_cc$Ti), X_cc$g)
  B_X <- th_all[2:5]
  r <- X_cc$y - as.numeric(X %*% B_X)
  sb <- exp(th_all[6])
  s <- exp(th_all[7])
  li <- sapply(wj, function(ji) {
    lxx_b0(r[ji], sb, s)
  })
  if (any(is.nan(li) | is.na(li) | is.infinite(li))) {
    return(-Inf)
  }
  if (retal) {
    return(li)
  }
  sum(li)
}

#' Evaluate the full conditional-model likelihood
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param th_all Full conditional-model parameter vector.
#' @param Xcc Complete-case longitudinal rows.
#' @param Xrc_l Input used by this helper.
#' @param wjcc Complete-case row index list.
#' @param l_x Longitudinal log-likelihood function.
#' @param Sfun List of Weibull distribution helper functions.
#' @param lc Censored-subject integration function.
#' @param retal Whether to return subject-level likelihood contributions.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
l_full_wj <- function(th_all, Xcc, Xrc_l, wjcc, l_x = lxx_b0, Sfun, lc, retal = FALSE) {
  if (!("g" %in% names(Xcc))) {
    xi <- exp(th_all[1])
    Xcc <- mutate(Xcc, g = gfun(1, xi, Ti, tx))
  }
  if (!retal) {
    return(L_cc_wj(th_all, Xcc, wjcc, l_x, retal = FALSE) + L_rc_wj(th_all, Xrc_l,
      l_x,
      Sfun = Sfun, lc = lc, retal = FALSE
    ))
  }
  c(L_cc_wj(th_all, Xcc, wjcc, l_x, retal = TRUE), L_rc_wj(th_all, Xrc_l, l_x,
    Sfun = Sfun, lc = lc, retal = TRUE
  ))
}

#' Evaluate the full log-time conditional-model likelihood
#'
#' Log-time variant of [l_full_wj()] where the third fixed-effect coefficient
#' multiplies `log(Ti)` instead of `Ti`.
#' @param th_all Full conditional-model parameter vector.
#' @param Xcc Complete-case longitudinal rows.
#' @param Xrc_l Input used by this helper.
#' @param wjcc Complete-case row index list.
#' @param l_x Longitudinal log-likelihood function.
#' @param Sfun List of Weibull distribution helper functions.
#' @param lc Censored-subject integration function.
#' @param retal Whether to return subject-level likelihood contributions.
#' @return A value used by the log-time conditional-model pipeline.
#' @keywords internal
l_full_wj_logT <- function(th_all, Xcc, Xrc_l, wjcc, l_x = lxx_b0, Sfun, lc, retal = FALSE) {
  if (!("g" %in% names(Xcc))) {
    xi <- exp(th_all[1])
    Xcc <- mutate(Xcc, g = gfun(1, xi, Ti, tx))
  }
  if (!retal) {
    return(L_cc_wj_logT(th_all, Xcc, wjcc, l_x, retal = FALSE) + L_rc_wj(th_all, Xrc_l,
      l_x,
      Sfun = Sfun, lc = lc, retal = FALSE
    ))
  }
  c(L_cc_wj_logT(th_all, Xcc, wjcc, l_x, retal = TRUE), L_rc_wj(th_all, Xrc_l, l_x,
    Sfun = Sfun, lc = lc, retal = TRUE
  ))
}
