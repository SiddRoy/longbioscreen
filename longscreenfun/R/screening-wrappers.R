#' Fit the Cox screen to one simulated dataset
#'
#' Advanced helper used by [fit_screen_simulations()] for low-dimensional
#' simulations.
#' @family screening workflows
#' @param dat Simulated dataset list with `X` and `S` components.
#' @return A numeric vector containing the Cox estimate, standard error, and p-value.
#' @export
cox_screen_pvalue <- function(dat) {
  Xcox <- dat$X %>%
    dplyr::group_by(i) %>%
    dplyr::mutate(
      tim2 = c(tim[-1], Ti[1]),
      del = c(tim[-1] * 0, del[1])
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(i, tim, tim2, y, Ti, del)

  summary(survival::coxph(
    survival::Surv(tim, tim2, del) ~ y,
    data = Xcox,
    control = survival::coxph.control(timefix = FALSE)
  ))$coef[c(1, 3, 5)]
}

#' Fit the joint-model screen to one simulated dataset
#'
#' Advanced helper used by [fit_screen_simulations()] for low-dimensional
#' simulations.
#' @family screening workflows
#' @param dat Simulated dataset list with `X` and `S` components.
#' @param label Optional progress label printed every ten fits.
#' @param jm_random_specs Named character vector of random-effect formulas used
#'   for the longitudinal submodel before the joint-model fit.
#' @param return_timing Whether to return elapsed time in seconds for each
#'   joint-model random-effect specification.
#' @param ... Additional arguments passed to [JM::jointModel()].
#' @return A numeric vector containing the JM event-association estimate,
#'   standard error, and p-value for each requested random-effect specification.
#'   If `return_timing = TRUE`, returns a list with `values` and `time_sec`.
#' @export
jm_screen_pvalue <- function(
    dat, label = NULL, jm_random_specs = jm_default_random_specs(),
    return_timing = FALSE, ...) {
  Sjm <- dplyr::mutate(dat$S, i = factor(i, unique(i))) %>%
    dplyr::rename(tim = Ti)
  Xjm <- dat$X %>% dplyr::select(i, tim, y)
  fitSURV_coxph <- survival::coxph(
    survival::Surv(tim, del) ~ 1,
    data = Sjm, x = TRUE
  )

  if (!is.null(label) && (label %% 10) == 0) print(label)
  jm_fit_random_specs(
    Xjm, fitSURV_coxph,
    jm_random_specs = jm_random_specs,
    return_timing = return_timing,
    ...
  )
}

#' Fit the legacy all-R conditional-model screen to one dataset
#'
#' Legacy all-R implementation retained for manual reference. Package screening
#' workflows use [cm_screen()], the C++-backed implementation.
#' @family screening workflows
#' @param dat Simulated dataset list with `X` and `S` components.
#' @param time Event-time effect scale: `"T"` for `Ti` or `"logT"` for
#'   `log(Ti)`.
#' @param fit_xi Whether to estimate `xi`.
#' @param xi Fixed positive `xi` value used when `fit_xi = FALSE`.
#' @param lxi0 Initial log xi value used when `fit_xi = TRUE`.
#' @param complete_case Whether to use only event subjects and fit the
#'   complete-case likelihood.
#' @param sandwich Whether to compute conditional-model standard errors and
#'   p-values with a subject-level sandwich covariance matrix.
#' @param label Optional progress label printed before the fit.
#' @param lc Optional censored-subject integration function for full-data fits.
#' @param return_timing Whether to return elapsed time in seconds for the
#'   initialization step.
#' @return Seven conditional-model statistics used by the simulation output. If
#'   `return_timing = TRUE`, returns a list with `values` and `init_time_sec`.
#' @export
cm_screenR <- function(
    dat, time = c("T", "logT"), fit_xi = FALSE, xi = 1, lxi0 = log(xi),
    complete_case = FALSE, sandwich = FALSE, label = NULL, lc = NULL,
    return_timing = FALSE) {
  time <- match.arg(time)
  if (!is.null(label)) print(label)
  if (!isTRUE(fit_xi) && (!is.numeric(xi) || length(xi) != 1 || !is.finite(xi) || xi <= 0)) {
    stop("`xi` must be a single positive finite value.", call. = FALSE)
  }
  if (isTRUE(fit_xi) && (!is.numeric(lxi0) || length(lxi0) != 1 || !is.finite(lxi0))) {
    stop("`lxi0` must be a single finite value.", call. = FALSE)
  }
  if (!is.logical(sandwich) || length(sandwich) != 1 || is.na(sandwich)) {
    stop("`sandwich` must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.null(lc)) {
    lc <- if (time == "logT") lc2_logT else lc2
  }

  init_fun <- if (time == "logT") init_fixxi_logT else init_fixxi
  cc_lik_fun <- if (time == "logT") L_cc_wj_logT else L_cc_wj
  lxi <- if (isTRUE(fit_xi)) lxi0 else log(xi)

  if (isTRUE(complete_case)) {
    S_cc <- dplyr::filter(dat$S, .data$del == 1)
    X_cc <- dplyr::semi_join(dat$X, dplyr::select(S_cc, "i"), by = "i")
    dat_set <- xifix_subX(X_cc, S_cc, lxi0 = lxi, addii = TRUE)
    dat_set$Xrc <- dat_set$Xrc[0, , drop = FALSE]
    dat_set$X_lcl <- list()

    timed_init <- system.time({
      b0 <- init_fun(dat_set, lxi0 = lxi)
      if (isTRUE(fit_xi)) {
        b0 <- c(lxi0, b0)
      }
    })

    Xcc_fit <- dat_set$Xcc
    if (isTRUE(fit_xi)) {
      Xcc_fit <- Xcc_fit[, setdiff(names(Xcc_fit), "g"), drop = FALSE]
    }
    assoc_idx <- if (isTRUE(fit_xi)) 4:5 else 3:4
    fit_values <- function() {
      subject_dev <- function(th) {
        Xcc_i <- Xcc_fit
        th_all <- th
        if (isTRUE(fit_xi)) {
          Xcc_i <- dplyr::mutate(Xcc_i, g = gfun(1, exp(th[1]), .data$Ti, .data$tx))
        } else {
          th_all <- c(lxi, th)
        }
        -cc_lik_fun(th_all, Xcc_i, dat_set$wcc, retal = TRUE)
      }
      dev <- function(th) {
        sum(subject_dev(th))
      }
      fit_jj <- optim(b0, dev, method = "BFGS", hessian = FALSE)
      fit_jj$hessian <- H_ii(fit_jj$par, dev, r = 2)
      Vyth <- if (isTRUE(sandwich)) {
        cm_sandwich_cov(fit_jj$hessian, g_xl(fit_jj$par, subject_dev, r = 2))
      } else {
        cm_inverse_hessian_cov(fit_jj$hessian)
      }
      cm_assoc_inference(fit_jj$par, Vyth, assoc_idx)
    }
    if (return_timing) {
      values <- tryCatch(fit_values(), error = function(e) rep(NA_real_, 7))
      return(list(values = values, init_time_sec = unname(timed_init[["elapsed"]])))
    }
    return(fit_values())
  }

  fit_w <- fitWeibS(Ti = dat$S$Ti, del = dat$S$del)
  dat_set <- xifix_subX(dat$X, dat$S, lxi0 = lxi, addii = TRUE)
  if (isTRUE(fit_xi)) {
    dat_devF <- if (time == "logT") {
      cM_devfun_xifit_datl_logT(dat_set, fit_w$Sfun, lc = lc)
    } else {
      cM_devfun_xifit_datl(dat_set, fit_w$Sfun, lc = lc)
    }
    timed_init <- system.time(b0 <- c(lxi0, init_fun(dat_set, lxi0 = lxi0)))
    fit_values <- function() {
      fit_jj <- fitxi_fitcm0(dat_devF, fit_w, b0)
      outcmj <- if (isTRUE(sandwich)) {
        fitxi_fitcm_sandwich(dat_devF, fit_w, b0, dat_set$m_YS, fit_jj)
      } else {
        fitxi_fitcm_noB(dat_devF, fit_w, b0, dat_set$m_YS, fit_jj)
      }
      outcmj <- outcmj[-(1:2)]
      c(fit_jj$par[4:5], unlist(outcmj))
    }
  } else {
    dat_devF <- if (time == "logT") {
      cM_devfun_xifix_datl_logT(lxi, dat_set, fit_w$Sfun, lc = lc)
    } else {
      cM_devfun_xifix_datl(lxi, dat_set, fit_w$Sfun, lc = lc)
    }
    timed_init <- system.time(b0 <- init_fun(dat_set, lxi0 = lxi))
    fit_values <- function() {
      fit_jj <- fixxi_fitcm0(dat_devF, fit_w, b0)
      outcmj <- if (isTRUE(sandwich)) {
        fixxi_fitcm_sandwich(dat_devF, fit_w, b0, dat_set$m_YS, fit_jj)
      } else {
        fixxi_fitcm_noB(dat_devF, fit_w, b0, dat_set$m_YS, fit_jj)
      }
      outcmj <- outcmj[-(1:2)]
      c(fit_jj$par[3:4], unlist(outcmj))
    }
  }

  if (return_timing) {
    values <- tryCatch(fit_values(), error = function(e) rep(NA_real_, 7))
    return(list(values = values, init_time_sec = unname(timed_init[["elapsed"]])))
  }
  fit_values()
}

#' Fit a configurable C++ conditional-model screen to one dataset
#'
#' Fits the conditional-model screen for one simulated dataset using the
#' C++-backed random-intercept LMM implementation, with options for the
#' event-time scale, fixed or estimated `xi`, and full-data or complete-case
#' likelihood.
#' @family screening workflows
#' @param dat Simulated dataset list with `X` and `S` components.
#' @param time Event-time effect scale: `"T"` for `Ti` or `"logT"` for
#'   `log(Ti)`.
#' @param fit_xi Whether to estimate `xi`.
#' @param xi Fixed positive `xi` value used when `fit_xi = FALSE`.
#' @param lxi0 Initial log xi value used when `fit_xi = TRUE`.
#' @param complete_case Whether to use only event subjects and fit the
#'   complete-case likelihood.
#' @param sandwich Whether to compute conditional-model standard errors and
#'   p-values with a subject-level sandwich covariance matrix.
#' @param label Optional progress label printed before the fit.
#' @param whiten Whether to use the whitened optimizer. Defaults to `TRUE`.
#' @param return_fit Whether to return fitted objects and covariance details
#'   along with the seven screening statistics.
#' @return Seven named conditional-model statistics used by the simulation
#'   output, or a list with fit details when `return_fit = TRUE`.
#' @export
cm_screen <- function(
    dat, time = c("T", "logT"), fit_xi = FALSE, xi = 1, lxi0 = log(xi),
    complete_case = FALSE, sandwich = FALSE, label = NULL, whiten = TRUE,
    return_fit = FALSE) {
  time <- match.arg(time)
  if (!is.logical(fit_xi) || length(fit_xi) != 1L || is.na(fit_xi)) {
    stop("`fit_xi` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(complete_case) || length(complete_case) != 1L ||
      is.na(complete_case)) {
    stop("`complete_case` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(sandwich) || length(sandwich) != 1L || is.na(sandwich)) {
    stop("`sandwich` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(whiten) || length(whiten) != 1L || is.na(whiten)) {
    stop("`whiten` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(return_fit) || length(return_fit) != 1L || is.na(return_fit)) {
    stop("`return_fit` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!isTRUE(fit_xi) &&
      (!is.numeric(xi) || length(xi) != 1L || !is.finite(xi) || xi <= 0)) {
    stop("`xi` must be a single positive finite value.", call. = FALSE)
  }
  if (isTRUE(fit_xi) &&
      (!is.numeric(lxi0) || length(lxi0) != 1L || !is.finite(lxi0))) {
    stop("`lxi0` must be a single finite value.", call. = FALSE)
  }

  if (!is.null(label)) {
    print(label)
  }
  if (!cm_lmm_is_longscreen_data(dat)) {
    stop("`dat` must be a list with `X` and `S` tables.", call. = FALSE)
  }
  if (is.list(dat) && length(dat) == 1L && is.list(dat[[1L]]) &&
      all(c("X", "S") %in% names(dat[[1L]]))) {
    dat <- dat[[1L]]
  }

  if (isTRUE(complete_case)) {
    X <- as.data.frame(dat$X)
    S <- as.data.frame(dat$S)
    required_X <- c("i", "tim", "y")
    required_S <- c("i", "Ti", "del")
    stopifnot(all(required_X %in% names(X)), all(required_S %in% names(S)))
    S <- S[S$del == 1L, , drop = FALSE]
    if (nrow(S) == 0L) {
      stop("`complete_case` leaves no event subjects.", call. = FALSE)
    }
    X <- X[as.character(X$i) %in% as.character(S$i), , drop = FALSE]
    dat <- list(X = X, S = S)
  }

  control <- cm_control(
    ftrans_fn = if (time == "logT") log else cm_lmm_ftrans
  )
  xi_start <- if (isTRUE(fit_xi)) exp(lxi0) else xi
  xi_mode <- if (isTRUE(fit_xi)) "fitxi" else "fixxi"
  covariance_mode <- if (isTRUE(return_fit)) {
    "both"
  } else if (isTRUE(sandwich)) {
    "sandwich"
  } else {
    "model"
  }
  truth <- cm_lmm_truth(xi = xi_start)
  res <- cm_lmm_fit_two_stage(
    sim = dat,
    xi_mode = xi_mode,
    xi_fix = xi,
    xi_start = xi_start,
    whiten = whiten,
    compute_covariance = TRUE,
    covariance = covariance_mode,
    longscreen_truth = truth,
    control = control
  )

  if (isTRUE(fit_xi)) {
    theta_hat <- res$fit_y$par
    cov <- res$cov_y
  } else {
    theta_hat <- res$fit_fixxi$par
    cov <- res$cov_fixxi
  }
  V <- if (isTRUE(sandwich)) cov$V1 else cov$V2
  values <- cm_lmm_assoc_inference(theta_hat, V)

  if (!isTRUE(return_fit)) {
    return(values)
  }

  list(
    values = values,
    fit = res,
    vcov = V,
    covariance = cov,
    settings = list(
      time = time,
      fit_xi = fit_xi,
      xi = xi,
      lxi0 = lxi0,
      complete_case = complete_case,
      sandwich = sandwich,
      whiten = whiten
    )
  )
}

#' Fit the conditional-model screen to one simulated dataset
#'
#' Advanced helper used by [fit_screen_simulations()] for low-dimensional
#' simulations.
#' @family screening workflows
#' @param dat Simulated dataset list with `X` and `S` components.
#' @param label Optional progress label printed before the fit.
#' @param xi Fixed positive xi value.
#' @param complete_case Whether to use only event subjects and fit the
#'   complete-case likelihood.
#' @param sandwich Whether to compute standard errors and p-values with a
#'   subject-level sandwich covariance matrix.
#' @param whiten Whether to use the whitened optimizer. Defaults to `TRUE`.
#' @return Seven conditional-model statistics used by the simulation output.
#' @export
cm_screen_stats <- function(
    dat, label = NULL, xi = 1, complete_case = FALSE, sandwich = FALSE,
    whiten = TRUE) {
  cm_screen(dat,
    time = "T", fit_xi = FALSE, xi = xi, complete_case = complete_case,
    sandwich = sandwich, label = label, whiten = whiten
  )
}

#' Fit the estimated-xi conditional-model screen to one simulated dataset
#'
#' Advanced helper used to fit the conditional model while estimating the kernel
#' decay parameter `xi`.
#' @family screening workflows
#' @param dat Simulated dataset list with `X` and `S` components.
#' @param label Optional progress label printed before the fit.
#' @param lxi0 Initial log xi value.
#' @param complete_case Whether to use only event subjects and fit the
#'   complete-case likelihood.
#' @param sandwich Whether to compute standard errors and p-values with a
#'   subject-level sandwich covariance matrix.
#' @param whiten Whether to use the whitened optimizer. Defaults to `TRUE`.
#' @return Seven conditional-model statistics used by the simulation output.
#' @export
cm_screen_fitxi <- function(
    dat, label = NULL, lxi0 = 0, complete_case = FALSE, sandwich = FALSE,
    whiten = TRUE) {
  cm_screen(dat,
    time = "T", fit_xi = TRUE, lxi0 = lxi0, complete_case = complete_case,
    sandwich = sandwich, label = label, whiten = whiten
  )
}

screen_clean_numeric <- function(x, n) {
  out <- rep(NA_real_, n)
  x <- as.numeric(x)
  out[seq_len(min(length(x), n))] <- x[seq_len(min(length(x), n))]
  out
}

screen_clean_numeric_named <- function(x, row_names) {
  out <- setNames(rep(NA_real_, length(row_names)), row_names)
  x <- as.numeric(x)
  out[seq_len(min(length(x), length(out)))] <- x[seq_len(min(length(x), length(out)))]
  out
}

screen_lowdim_row_indices <- function(jm_rows, jm_time_rows) {
  cox_rows <- seq_len(3)
  jm_stat_rows <- seq.int(max(cox_rows) + 1L, length.out = length(jm_rows))
  cm_rows <- seq.int(max(jm_stat_rows) + 1L, length.out = 7L)
  time_cox_row <- max(cm_rows) + 1L
  time_jm_rows <- seq.int(time_cox_row + 1L, length.out = length(jm_time_rows))
  time_cm_row <- max(time_jm_rows) + 1L
  time_cm_init_row <- time_cm_row + 1L

  list(
    cox_rows = cox_rows,
    jm_stat_rows = jm_stat_rows,
    cm_rows = cm_rows,
    time_cox_row = time_cox_row,
    time_jm_rows = time_jm_rows,
    time_cm_row = time_cm_row,
    time_cm_init_row = time_cm_init_row,
    n_legacy = max(cm_rows),
    n_current_no_cm_init = time_cm_row,
    n_current = time_cm_init_row
  )
}

screen_lowdim_timing_out <- function(tCox, tJM, tCM, tCM_init) {
  rbind(
    time_cox_sec = tCox,
    tJM,
    time_cm_sec = tCM,
    time_cm_init_sec = tCM_init
  )
}

screen_validate_saved_row_names <- function(
    saved, idx, jm_rows, timing_rows = NULL, expected_timing_rows = NULL) {
  saved_row_names <- rownames(saved)
  if (is.null(saved_row_names)) {
    return(invisible(NULL))
  }

  saved_jm_rows <- saved_row_names[idx$jm_stat_rows]
  if (any(nzchar(saved_jm_rows)) && !identical(saved_jm_rows, jm_rows)) {
    stop(
      "Existing `out_file` has JM rows that do not match `jm_random_specs`.",
      call. = FALSE
    )
  }

  if (is.null(timing_rows)) {
    return(invisible(NULL))
  }

  saved_timing_rows <- saved_row_names[timing_rows]
  if (any(nzchar(saved_timing_rows)) && !identical(saved_timing_rows, expected_timing_rows)) {
    stop(
      "Existing `out_file` has timing rows that do not match the current output layout.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

screen_lowdim_resume_state <- function(out_file, n_sim, jm_rows, jm_time_rows) {
  idx <- screen_lowdim_row_indices(jm_rows, jm_time_rows)
  pCox <- matrix(NA_real_, 3, n_sim)
  pjm <- matrix(NA_real_, length(jm_rows), n_sim, dimnames = list(jm_rows, NULL))
  p_out <- matrix(NA_real_, 7, n_sim)
  tCox <- rep(NA_real_, n_sim)
  tJM <- matrix(NA_real_, length(jm_time_rows), n_sim, dimnames = list(jm_time_rows, NULL))
  tCM <- rep(NA_real_, n_sim)
  tCM_init <- rep(NA_real_, n_sim)

  state <- list(
    pCox = pCox,
    pjm = pjm,
    p_out = p_out,
    tCox = tCox,
    tJM = tJM,
    tCM = tCM,
    tCM_init = tCM_init,
    cox_todo = seq_len(n_sim),
    jm_todo = seq_len(n_sim),
    cm_todo = seq_len(n_sim),
    loaded_format = "none"
  )

  if (!file.exists(out_file)) {
    return(state)
  }

  saved <- readRDS(out_file)
  if (!is.matrix(saved) || !is.numeric(saved)) {
    stop("Existing `out_file` must contain a numeric matrix.", call. = FALSE)
  }
  if (ncol(saved) != n_sim) {
    stop(
      sprintf(
        "Existing `out_file` has %d columns, but `sim_data` has length %d.",
        ncol(saved), n_sim
      ),
      call. = FALSE
    )
  }

  loaded_format <- "current"
  if (nrow(saved) == idx$n_current_no_cm_init) {
    saved_row_names <- rownames(saved)
    saved <- rbind(saved, rep(NA_real_, ncol(saved)))
    if (is.null(saved_row_names)) {
      rownames(saved) <- NULL
    } else {
      rownames(saved) <- c(saved_row_names, "time_cm_init_sec")
    }
    loaded_format <- "current_no_cm_init"
  }

  if (nrow(saved) == idx$n_current) {
    screen_validate_saved_row_names(
      saved, idx, jm_rows,
      timing_rows = c(idx$time_cox_row, idx$time_jm_rows, idx$time_cm_row, idx$time_cm_init_row),
      expected_timing_rows = c("time_cox_sec", jm_time_rows, "time_cm_sec", "time_cm_init_sec")
    )
    state$pCox[,] <- saved[idx$cox_rows, , drop = FALSE]
    state$pjm[,] <- saved[idx$jm_stat_rows, , drop = FALSE]
    state$p_out[,] <- saved[idx$cm_rows, , drop = FALSE]
    state$tCox[] <- saved[idx$time_cox_row, ]
    state$tJM[,] <- saved[idx$time_jm_rows, , drop = FALSE]
    state$tCM[] <- saved[idx$time_cm_row, ]
    state$tCM_init[] <- saved[idx$time_cm_init_row, ]
    state$cox_todo <- integer(0)
    state$jm_todo <- integer(0)
    state$cm_todo <- which(is.na(state$tCM))
    state$loaded_format <- loaded_format
    return(state)
  }

  if (nrow(saved) == idx$n_legacy) {
    screen_validate_saved_row_names(saved, idx, jm_rows)
    state$pCox[,] <- saved[idx$cox_rows, , drop = FALSE]
    state$pjm[,] <- saved[idx$jm_stat_rows, , drop = FALSE]
    state$p_out[,] <- saved[idx$cm_rows, , drop = FALSE]
    state$cox_todo <- integer(0)
    state$jm_todo <- integer(0)
    state$cm_todo <- which(colSums(!is.na(state$p_out)) == 0L)
    state$loaded_format <- "legacy"
    return(state)
  }

  stop(
    sprintf(
      paste(
        "Existing `out_file` has %d rows; expected %d for current output,",
        "%d for current output without CM-init timing, or %d for legacy no-timing output."
      ),
      nrow(saved), idx$n_current, idx$n_current_no_cm_init, idx$n_legacy
    ),
    call. = FALSE
  )
}

fit_screen_simulations_impl <- function(
    sim_data, out_file, save_each, jm_random_specs, cm_screen_fun,
    cm_whiten) {
  j_i <- seq_along(sim_data)
  jm_random_specs <- normalize_jm_random_specs(jm_random_specs)
  jm_rows <- jm_stat_row_names(jm_random_specs)
  jm_time_rows <- jm_time_row_names(jm_random_specs, legacy_single = TRUE)
  state <- screen_lowdim_resume_state(out_file, length(j_i), jm_rows, jm_time_rows)

  pCox <- state$pCox
  tCox <- state$tCox
  for (i in state$cox_todo) {
    timed <- system.time(res <- tryCatch(
      cox_screen_pvalue(sim_data[[i]]),
      error = function(e) rep(NA_real_, 3)
    ))
    pCox[, i] <- screen_clean_numeric(res, 3)
    tCox[i] <- unname(timed["elapsed"])
  }

  pjm <- state$pjm
  tJM <- state$tJM
  for (i in state$jm_todo) {
    res <- tryCatch(
      jm_screen_pvalue(
        sim_data[[i]], i,
        jm_random_specs = jm_random_specs,
        return_timing = TRUE
      ),
      error = function(e) list(
        values = setNames(rep(NA_real_, length(jm_rows)), jm_rows),
        time_sec = setNames(rep(NA_real_, length(jm_time_rows)), jm_time_rows)
      )
    )
    pjm[, i] <- screen_clean_numeric_named(res$values, jm_rows)
    tJM[, i] <- screen_clean_numeric_named(res$time_sec, jm_time_rows)
  }

  coxjm_out <- rbind(pCox, pjm)
  p_out <- state$p_out
  tCM <- state$tCM
  tCM_init <- state$tCM_init
  timing_out <- screen_lowdim_timing_out(tCox, tJM, tCM, tCM_init)

  if (length(state$cm_todo) == 0L) {
    out <- rbind(coxjm_out, p_out, timing_out)
    if (!identical(state$loaded_format, "current")) {
      saveRDS(out, out_file)
    }
    return(invisible(out))
  }

  for (i in state$cm_todo) {
    timed <- system.time(p_i <- tryCatch(
      cm_screen_fun(sim_data[[i]], i, whiten = cm_whiten),
      error = function(e) rep(NA_real_, 7)
    ))
    p_out[, i] <- screen_clean_numeric(p_i, 7)
    tCM[i] <- unname(timed["elapsed"])
    tCM_init[i] <- NA_real_
    timing_out["time_cm_sec", i] <- tCM[i]
    timing_out["time_cm_init_sec", i] <- tCM_init[i]
    if (save_each) saveRDS(rbind(coxjm_out, p_out, timing_out), out_file)
  }

  timing_out["time_cox_sec", ] <- tCox
  timing_out[jm_time_rows, ] <- tJM
  timing_out["time_cm_sec", ] <- tCM
  timing_out["time_cm_init_sec", ] <- tCM_init
  out <- rbind(coxjm_out, p_out, timing_out)
  saveRDS(out, out_file)
  invisible(out)
}

#' Fit all three screens to replicated simulated datasets
#'
#' Primary screening wrapper used by `31_CP/CPsim.R`.
#' @family screening workflows
#' @param sim_data List of simulated datasets with `X` and `S` components.
#' @param out_file Output RDS path. If this file already exists, completed
#'   current-format Cox/JM work and CM columns with recorded CM timing are
#'   reused before continuing.
#' @param save_each Whether to save partial conditional-model progress.
#' @param jm_random_specs Named character vector of random-effect formulas used
#'   for the longitudinal submodel before the joint-model fit. The default keeps
#'   the historical random-intercept-only JM screen; pass
#'   `c(b0 = "~1 | i", b01 = "~ tim | i")` to add the random intercept+slope
#'   fit.
#' @param cm_whiten Whether to use the whitened conditional-model optimizer.
#'   Defaults to `TRUE`.
#' @return Invisibly returns the 17-row matrix written to `out_file`: three Cox
#'   rows, three JM rows, seven conditional-model rows, and three elapsed-time
#'   rows for Cox, JM, and CM fits, plus one elapsed-time row for the
#'   conditional-model initialization step.
#' @export
fit_screen_simulations <- function(
    sim_data, out_file, save_each = TRUE,
    jm_random_specs = jm_default_random_specs(), cm_whiten = TRUE) {
  fit_screen_simulations_impl(
    sim_data = sim_data,
    out_file = out_file,
    save_each = save_each,
    jm_random_specs = jm_random_specs,
    cm_screen_fun = cm_screen_stats,
    cm_whiten = cm_whiten
  )
}

#' Prepare reusable objects for high-dimensional biomarker screening
#'
#' Builds the Cox counting-process data, survival submodel, Weibull survival
#' fit, and conditional-model data split used repeatedly by the biomarker
#' screeners.
#' @family screening workflows
#' @param sim_data Dataset returned by [drawY100()].
#' @return A list containing `sim_data`, prepared Cox/JM data, the fitted
#'   survival-only Cox model, the fitted Weibull survival helper, and the
#'   conditional-model data split.
#' @export
prepare_biomarker_screen <- function(sim_data) {
  Xcox <- dplyr::left_join(sim_data$X, sim_data$S, by = "i") %>%
    dplyr::group_by(i) %>%
    dplyr::mutate(
      tim2 = c(tim[-1], Ti[1]),
      del = c(tim[-1] * 0, del[1])
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(i, tim, tim2, y, Ti, del)

  Sjm <- dplyr::mutate(sim_data$S, i = factor(i, unique(i))) %>%
    dplyr::rename(tim = Ti)
  fitSURV_coxph <- survival::coxph(
    survival::Surv(tim, del) ~ 1,
    data = Sjm, x = TRUE
  )

  list(
    sim_data = sim_data,
    Xcox = Xcox,
    Sjm = Sjm,
    Xjm = sim_data$X,
    fitSURV_coxph = fitSURV_coxph,
    fit_w = fitWeibS(Ti = sim_data$S$Ti, del = sim_data$S$del),
    dat_set0 = xifix_subX(sim_data$X, sim_data$S, addii = TRUE),
    cm_cpp_T = cm_lmm_prepare_biomarker_cpp(sim_data, time = "T")
  )
}

#' Fit Cox screens for high-dimensional biomarkers
#'
#' @family screening workflows
#' @param prepared Prepared screening object returned by [prepare_biomarker_screen()].
#' @param biomarker_ids Column ids in `sim_data$y` to screen.
#' @param return_timing Whether to return elapsed time in seconds for each
#'   biomarker fit.
#' @return Numeric vector of Cox p-values, one per biomarker id. If
#'   `return_timing = TRUE`, returns a list with `values` and `time_sec`.
#' @export
cox_screen_biomarkers <- function(prepared, biomarker_ids, return_timing = FALSE) {
  Xcox <- prepared$Xcox
  sim_data <- prepared$sim_data
  p_values <- rep(NA_real_, length(biomarker_ids))
  time_sec <- rep(NA_real_, length(biomarker_ids))

  for (pos in seq_along(biomarker_ids)) {
    j <- biomarker_ids[pos]
    timed <- system.time(p_values[pos] <- tryCatch({
      Xcox$y <- sim_data$y[, j]
      summary(survival::coxph(
        survival::Surv(tim, tim2, del) ~ y,
        data = Xcox,
        control = survival::coxph.control(timefix = FALSE)
      ))$coef[5]
    }, error = function(e) NA_real_))
    time_sec[pos] <- unname(timed[["elapsed"]])
  }

  if (return_timing) {
    return(list(values = p_values, time_sec = time_sec))
  }
  p_values
}

#' Fit joint-model screens for high-dimensional biomarkers
#'
#' @family screening workflows
#' @param prepared Prepared screening object returned by [prepare_biomarker_screen()].
#' @param biomarker_ids Column ids in `sim_data$y` to screen.
#' @param return_timing Whether to return elapsed time in seconds for each
#'   biomarker fit.
#' @param jm_random_specs Named character vector of random-effect formulas used
#'   for the longitudinal submodel before the joint-model fit.
#' @param jm_stats Character vector of joint-model event-association statistics
#'   to retain. The default keeps the historical p-value-only high-dimensional
#'   output; pass `c("est", "se", "p")` to retain estimates and standard errors.
#' @return Numeric vector of joint-model p-values for the historical default, or
#'   a matrix of requested joint-model statistics for expanded settings. If
#'   `return_timing = TRUE`, returns a list with `values` and `time_sec`.
#' @export
jm_screen_biomarkers <- function(
    prepared, biomarker_ids, return_timing = FALSE,
    jm_random_specs = jm_default_random_specs(), jm_stats = "p") {
  jm_random_specs <- normalize_jm_random_specs(jm_random_specs)
  jm_stats <- normalize_jm_stats(jm_stats)
  jm_rows <- jm_stat_row_names(jm_random_specs, jm_stats)
  jm_time_rows <- jm_time_row_names(jm_random_specs)
  legacy_p_only <- length(jm_random_specs) == 1 &&
    names(jm_random_specs)[1] == "b0" &&
    identical(jm_stats, "p")
  Xjm <- prepared$Xjm
  sim_data <- prepared$sim_data
  values <- matrix(NA_real_, length(jm_rows), length(biomarker_ids), dimnames = list(jm_rows, NULL))
  time_sec <- matrix(NA_real_, length(jm_time_rows), length(biomarker_ids), dimnames = list(jm_time_rows, NULL))

  for (pos in seq_along(biomarker_ids)) {
    j <- biomarker_ids[pos]
    res <- tryCatch({
      if ((j %% 10) == 0) print(j)
      Xjm$y <- sim_data$y[, j]
      jm_fit_random_specs(
        Xjm, prepared$fitSURV_coxph,
        jm_random_specs = jm_random_specs,
        jm_stats = jm_stats,
        return_timing = TRUE
      )
    }, error = function(e) list(
      values = setNames(rep(NA_real_, length(jm_rows)), jm_rows),
      time_sec = setNames(rep(NA_real_, length(jm_time_rows)), jm_time_rows)
    ))
    values[, pos] <- as.numeric(res$values[jm_rows])
    time_sec[, pos] <- as.numeric(res$time_sec[jm_time_rows])
  }

  if (legacy_p_only) {
    values <- values[1, ]
    time_sec <- time_sec[1, ]
    names(values) <- NULL
    names(time_sec) <- NULL
  }

  if (return_timing) {
    return(list(values = values, time_sec = time_sec))
  }
  values
}

screen_print_biomarker_progress <- function(biomarker_id, pos) {
  if (is.numeric(biomarker_id) &&
      length(biomarker_id) == 1L &&
      !is.na(biomarker_id)) {
    if ((biomarker_id %% 10) == 0) print(biomarker_id)
    return(invisible(NULL))
  }
  if ((pos %% 10) == 0) print(biomarker_id)
  invisible(NULL)
}

#' Fit one C++ conditional-model screen for a high-dimensional biomarker
#'
#' @family screening workflows
#' @param prepared Prepared screening object returned by [prepare_biomarker_screen()].
#' @param biomarker_id Column id in `sim_data$y` to screen.
#' @param whiten Whether to use the whitened optimizer. Defaults to `TRUE`.
#' @param return_timing Whether to return elapsed time in seconds for the
#'   C++ CM starting-value step.
#' @return Seven conditional-model statistics for one biomarker. If
#'   `return_timing = TRUE`, returns a list with `values` and `init_time_sec`.
#' @export
cm_screen_biomarker <- function(prepared, biomarker_id, whiten = TRUE, return_timing = FALSE) {
  prepared <- cm_lmm_ensure_biomarker_cpp(
    prepared,
    time = "T",
    biomarker_ids = biomarker_id
  )
  cm_lmm_fit_biomarker_cpp(
    prepared[[cm_lmm_biomarker_cpp_field("T")]],
    biomarker_id,
    whiten = whiten,
    return_timing = return_timing
  )
}

#' Fit one legacy all-R conditional-model screen for a high-dimensional biomarker
#'
#' Legacy all-R implementation retained for manual reference. Package screening
#' workflows use [cm_screen_biomarker()], the C++-backed implementation.
#' @family screening workflows
#' @param prepared Prepared screening object returned by [prepare_biomarker_screen()].
#' @param biomarker_id Column id in `sim_data$y` to screen.
#' @param lc Censored-subject integration function.
#' @param return_timing Whether to return elapsed time in seconds for the
#'   `init_fixxi()` initialization step.
#' @return Seven conditional-model statistics for one biomarker. If
#'   `return_timing = TRUE`, returns a list with `values` and `init_time_sec`.
#' @export
cm_screen_biomarkerR <- function(prepared, biomarker_id, lc = lc2, return_timing = FALSE) {
  sim_data <- prepared$sim_data
  dat_set <- prepared$dat_set0
  dat_set$Xcc$y <- sim_data$y[dat_set$Xcc$ii, biomarker_id]
  dat_set$Xrc$y <- sim_data$y[dat_set$Xrc$ii, biomarker_id]
  dat_set$X_lcl <- Xcens_jdat(dplyr::mutate(
    dat_set$Xrc,
    i = as.integer(as.character(i))
  ))
  dat_devF <- cM_devfun_xifix_datl(0, dat_set, prepared$fit_w$Sfun, lc = lc)
  timed_init <- system.time(b0 <- init_fixxi(dat_set))
  fit_values <- function() {
    fit_jj <- fixxi_fitcm0(dat_devF, prepared$fit_w, b0)
    outcmj <- fixxi_fitcm_noB(dat_devF, prepared$fit_w, b0, dat_set$m_YS, fit_jj)[-(1:2)]
    c(fit_jj$par[3:4], unlist(outcmj))
  }
  if (return_timing) {
    values <- tryCatch(fit_values(), error = function(e) rep(NA_real_, 7))
    return(list(values = values, init_time_sec = unname(timed_init[["elapsed"]])))
  }
  fit_values()
}

#' Fit conditional-model screens for high-dimensional biomarkers
#'
#' @family screening workflows
#' @param prepared Prepared screening object returned by [prepare_biomarker_screen()].
#' @param biomarker_ids Column ids in `sim_data$y` to screen.
#' @param out_file Optional output RDS path for partial progress.
#' @param coxjm_out Optional Cox/JM result block to prepend when saving progress.
#' @param save_each Whether to save partial conditional-model progress.
#' @param whiten Whether to use the whitened optimizer. Defaults to `TRUE`.
#' @param return_timing Whether to return elapsed time in seconds for each
#'   biomarker fit.
#' @param timing_out Optional timing block to append when saving progress.
#' @return A 7 x `length(biomarker_ids)` matrix of conditional-model statistics.
#'   If `return_timing = TRUE`, returns a list with `values`, `time_sec`, and
#'   `init_time_sec`.
#' @export
cm_screen_biomarkers <- function(
    prepared, biomarker_ids, out_file = NULL, coxjm_out = NULL,
    save_each = FALSE, whiten = TRUE, return_timing = FALSE, timing_out = NULL) {
  prepared <- cm_lmm_ensure_biomarker_cpp(
    prepared,
    time = "T",
    biomarker_ids = biomarker_ids
  )
  p_out <- matrix(
    NA_real_,
    length(cm_lmm_assoc_stat_names),
    length(biomarker_ids),
    dimnames = list(cm_lmm_assoc_stat_names, NULL)
  )
  time_sec <- rep(NA_real_, length(biomarker_ids))
  init_time_sec <- rep(NA_real_, length(biomarker_ids))
  for (pos in seq_along(biomarker_ids)) {
    j <- biomarker_ids[pos]
    screen_print_biomarker_progress(j, pos)
    timed <- system.time(p_i <- tryCatch(
      cm_screen_biomarker(prepared, j, whiten = whiten, return_timing = TRUE),
      error = function(e) list(
        values = setNames(rep(NA_real_, length(cm_lmm_assoc_stat_names)), cm_lmm_assoc_stat_names),
        init_time_sec = NA_real_
      )
    ))
    p_out[, pos] <- p_i$values
    time_sec[pos] <- unname(timed[["elapsed"]])
    init_time_sec[pos] <- p_i$init_time_sec
    if (save_each && !is.null(out_file) && !is.null(coxjm_out)) {
      if (!is.null(timing_out)) {
        timing_out["time_cm_sec", pos] <- time_sec[pos]
        timing_out["time_cm_init_sec", pos] <- init_time_sec[pos]
        saveRDS(rbind(coxjm_out, p_out, timing_out), out_file)
      } else {
        saveRDS(rbind(coxjm_out, p_out), out_file)
      }
    }
  }
  if (return_timing) {
    return(list(
      values = p_out,
      time_sec = time_sec,
      init_time_sec = init_time_sec
    ))
  }
  p_out
}

#' Fit legacy all-R conditional-model screens for high-dimensional biomarkers
#'
#' Legacy all-R implementation retained for manual reference. Package screening
#' workflows use [cm_screen_biomarkers()], the C++-backed implementation.
#' @family screening workflows
#' @param prepared Prepared screening object returned by [prepare_biomarker_screen()].
#' @param biomarker_ids Column ids in `sim_data$y` to screen.
#' @param out_file Optional output RDS path for partial progress.
#' @param coxjm_out Optional Cox/JM result block to prepend when saving progress.
#' @param save_each Whether to save partial conditional-model progress.
#' @param lc Censored-subject integration function.
#' @param return_timing Whether to return elapsed time in seconds for each
#'   biomarker fit.
#' @param timing_out Optional timing block to append when saving progress.
#' @return A 7 x `length(biomarker_ids)` matrix of conditional-model statistics.
#'   If `return_timing = TRUE`, returns a list with `values`, `time_sec`, and
#'   `init_time_sec`.
#' @export
cm_screen_biomarkersR <- function(
    prepared, biomarker_ids, out_file = NULL, coxjm_out = NULL,
    save_each = FALSE, lc = lc2, return_timing = FALSE, timing_out = NULL) {
  p_out <- matrix(NA_real_, 7, length(biomarker_ids))
  time_sec <- rep(NA_real_, length(biomarker_ids))
  init_time_sec <- rep(NA_real_, length(biomarker_ids))
  for (pos in seq_along(biomarker_ids)) {
    j <- biomarker_ids[pos]
    screen_print_biomarker_progress(j, pos)
    timed <- system.time(p_i <- tryCatch(
      cm_screen_biomarkerR(prepared, j, lc = lc, return_timing = TRUE),
      error = function(e) list(values = rep(NA_real_, 7), init_time_sec = NA_real_)
    ))
    p_out[, pos] <- p_i$values
    time_sec[pos] <- unname(timed[["elapsed"]])
    init_time_sec[pos] <- p_i$init_time_sec
    if (save_each && !is.null(out_file) && !is.null(coxjm_out)) {
      if (!is.null(timing_out)) {
        timing_out["time_cm_sec", pos] <- time_sec[pos]
        timing_out["time_cm_init_sec", pos] <- init_time_sec[pos]
        saveRDS(rbind(coxjm_out, p_out, timing_out), out_file)
      } else {
        saveRDS(rbind(coxjm_out, p_out), out_file)
      }
    }
  }
  if (return_timing) {
    return(list(
      values = p_out,
      time_sec = time_sec,
      init_time_sec = init_time_sec
    ))
  }
  p_out
}

#' Fit all three screens to a high-dimensional biomarker block
#'
#' Primary screening wrapper used by `35_HDsim/HDsim.R`.
#' @family screening workflows
#' @param sim_data Dataset returned by [drawY100()].
#' @param biomarker_ids Column ids in `sim_data$y` to screen.
#' @param out_file Output RDS path.
#' @param save_each Whether to save partial conditional-model progress.
#' @param jm_random_specs Named character vector of random-effect formulas used
#'   for the longitudinal submodel before the joint-model fit.
#' @param jm_stats Character vector of joint-model event-association statistics
#'   to retain. The default keeps the historical p-value-only high-dimensional
#'   output; pass `c("est", "se", "p")` to retain estimates and standard errors.
#' @param cm_whiten Whether to use the whitened optimizer for conditional-model
#'   screens. Defaults to `TRUE`.
#' @return Invisibly returns the 13-row matrix written to `out_file`: one Cox
#'   row, one JM row, seven conditional-model rows, and three elapsed-time rows
#'   for Cox, JM, and CM fits, plus one elapsed-time row for the
#'   conditional-model initialization step.
#' @export
fit_screen_biomarkers <- function(
    sim_data, biomarker_ids, out_file, save_each = TRUE,
    jm_random_specs = jm_default_random_specs(), jm_stats = "p",
    cm_whiten = TRUE) {
  prepared <- prepare_biomarker_screen(sim_data)
  pCox <- cox_screen_biomarkers(prepared, biomarker_ids, return_timing = TRUE)
  pjm <- jm_screen_biomarkers(
    prepared, biomarker_ids,
    return_timing = TRUE,
    jm_random_specs = jm_random_specs,
    jm_stats = jm_stats
  )
  if (is.matrix(pjm$values)) {
    coxjm_out <- rbind(pCox = pCox$values, pjm$values)
  } else {
    coxjm_out <- rbind(pCox = pCox$values, pjm = pjm$values)
  }
  jm_timing <- if (is.matrix(pjm$time_sec)) {
    pjm$time_sec
  } else {
    rbind(time_jm_sec = pjm$time_sec)
  }
  timing_out <- rbind(
    time_cox_sec = pCox$time_sec,
    jm_timing,
    time_cm_sec = rep(NA_real_, length(biomarker_ids)),
    time_cm_init_sec = rep(NA_real_, length(biomarker_ids))
  )
  p_out <- cm_screen_biomarkers(
    prepared, biomarker_ids,
    out_file = out_file, coxjm_out = coxjm_out, timing_out = timing_out,
    save_each = save_each, whiten = cm_whiten, return_timing = TRUE
  )
  timing_out["time_cm_sec", ] <- p_out$time_sec
  timing_out["time_cm_init_sec", ] <- p_out$init_time_sec
  out <- rbind(coxjm_out, p_out$values, timing_out)
  saveRDS(out, out_file)
  invisible(out)
}

#' Fit the log-time conditional-model screen to one simulated dataset
#'
#' Log-time variant of [cm_screen_stats()] where the third fixed-effect
#' coefficient multiplies `log(Ti)` instead of `Ti`.
#' @family screening workflows
#' @param dat Simulated dataset list with `X` and `S` components.
#' @param label Optional progress label printed before the fit.
#' @param xi Fixed positive xi value.
#' @param complete_case Whether to use only event subjects and fit the
#'   complete-case likelihood.
#' @param sandwich Whether to compute standard errors and p-values with a
#'   subject-level sandwich covariance matrix.
#' @param whiten Whether to use the whitened optimizer. Defaults to `TRUE`.
#' @return Seven conditional-model statistics used by the simulation output.
#' @export
cm_screen_stats_logT <- function(
    dat, label = NULL, xi = 1, complete_case = FALSE, sandwich = FALSE,
    whiten = TRUE) {
  cm_screen(dat,
    time = "logT", fit_xi = FALSE, xi = xi, complete_case = complete_case,
    sandwich = sandwich, label = label, whiten = whiten
  )
}

#' Fit the estimated-xi log-time conditional-model screen
#'
#' Log-time variant of `cm_screen_fitxi()` where the third fixed-effect
#' coefficient multiplies `log(Ti)` instead of `Ti`.
#' @family screening workflows
#' @param dat Simulated dataset list with `X` and `S` components.
#' @param label Optional progress label printed before the fit.
#' @param lxi0 Initial log xi value.
#' @param complete_case Whether to use only event subjects and fit the
#'   complete-case likelihood.
#' @param sandwich Whether to compute standard errors and p-values with a
#'   subject-level sandwich covariance matrix.
#' @param whiten Whether to use the whitened optimizer. Defaults to `TRUE`.
#' @return Seven conditional-model statistics used by the simulation output.
#' @export
cm_screen_fitxi_logT <- function(
    dat, label = NULL, lxi0 = 0, complete_case = FALSE, sandwich = FALSE,
    whiten = TRUE) {
  cm_screen(dat,
    time = "logT", fit_xi = TRUE, lxi0 = lxi0, complete_case = complete_case,
    sandwich = sandwich, label = label, whiten = whiten
  )
}

#' Fit all three screens with the log-time conditional model
#'
#' Log-time variant of [fit_screen_simulations()]. Cox and joint-model rows are
#' unchanged; conditional-model rows use `log(Ti)` for the event-time term.
#' @family screening workflows
#' @param sim_data List of simulated datasets with `X` and `S` components.
#' @param out_file Output RDS path. If this file already exists, completed
#'   current-format Cox/JM work and CM columns with recorded CM timing are
#'   reused before continuing.
#' @param save_each Whether to save partial conditional-model progress.
#' @param jm_random_specs Named character vector of random-effect formulas used
#'   for the longitudinal submodel before the joint-model fit. The default keeps
#'   the historical random-intercept-only JM screen; pass
#'   `c(b0 = "~1 | i", b01 = "~ tim | i")` to add the random intercept+slope
#'   fit.
#' @param cm_whiten Whether to use the whitened conditional-model optimizer.
#'   Defaults to `TRUE`.
#' @return Invisibly returns the 17-row matrix written to `out_file`.
#' @export
fit_screen_simulations_logT <- function(
    sim_data, out_file, save_each = TRUE,
    jm_random_specs = jm_default_random_specs(), cm_whiten = TRUE) {
  fit_screen_simulations_impl(
    sim_data = sim_data,
    out_file = out_file,
    save_each = save_each,
    jm_random_specs = jm_random_specs,
    cm_screen_fun = cm_screen_stats_logT,
    cm_whiten = cm_whiten
  )
}

#' Fit one C++ log-time conditional-model screen for a high-dimensional biomarker
#'
#' Log-time variant of [cm_screen_biomarker()] where the third fixed-effect
#' coefficient multiplies `log(Ti)` instead of `Ti`.
#' @family screening workflows
#' @param prepared Prepared screening object returned by [prepare_biomarker_screen()].
#' @param biomarker_id Column id in `sim_data$y` to screen.
#' @param whiten Whether to use the whitened optimizer. Defaults to `TRUE`.
#' @param return_timing Whether to return elapsed time in seconds for the
#'   C++ CM starting-value step.
#' @return Seven conditional-model statistics for one biomarker. If
#'   `return_timing = TRUE`, returns a list with `values` and `init_time_sec`.
#' @export
cm_screen_biomarker_logT <- function(prepared, biomarker_id, whiten = TRUE, return_timing = FALSE) {
  prepared <- cm_lmm_ensure_biomarker_cpp(
    prepared,
    time = "logT",
    biomarker_ids = biomarker_id
  )
  cm_lmm_fit_biomarker_cpp(
    prepared[[cm_lmm_biomarker_cpp_field("logT")]],
    biomarker_id,
    whiten = whiten,
    return_timing = return_timing
  )
}

#' Fit one legacy all-R log-time CM screen for a high-dimensional biomarker
#'
#' Legacy all-R implementation retained for manual reference. Package screening
#' workflows use [cm_screen_biomarker_logT()], the C++-backed implementation.
#' @family screening workflows
#' @param prepared Prepared screening object returned by [prepare_biomarker_screen()].
#' @param biomarker_id Column id in `sim_data$y` to screen.
#' @param lc Censored-subject integration function.
#' @param return_timing Whether to return elapsed time in seconds for the
#'   `init_fixxi_logT()` initialization step.
#' @return Seven conditional-model statistics for one biomarker. If
#'   `return_timing = TRUE`, returns a list with `values` and `init_time_sec`.
#' @export
cm_screen_biomarker_logTR <- function(prepared, biomarker_id, lc = lc2_logT, return_timing = FALSE) {
  sim_data <- prepared$sim_data
  dat_set <- prepared$dat_set0
  dat_set$Xcc$y <- sim_data$y[dat_set$Xcc$ii, biomarker_id]
  dat_set$Xrc$y <- sim_data$y[dat_set$Xrc$ii, biomarker_id]
  dat_set$X_lcl <- Xcens_jdat(dplyr::mutate(
    dat_set$Xrc,
    i = as.integer(as.character(i))
  ))
  dat_devF <- cM_devfun_xifix_datl_logT(0, dat_set, prepared$fit_w$Sfun, lc = lc)
  timed_init <- system.time(b0 <- init_fixxi_logT(dat_set))
  fit_values <- function() {
    fit_jj <- fixxi_fitcm0(dat_devF, prepared$fit_w, b0)
    outcmj <- fixxi_fitcm_noB(dat_devF, prepared$fit_w, b0, dat_set$m_YS, fit_jj)[-(1:2)]
    c(fit_jj$par[3:4], unlist(outcmj))
  }
  if (return_timing) {
    values <- tryCatch(fit_values(), error = function(e) rep(NA_real_, 7))
    return(list(values = values, init_time_sec = unname(timed_init[["elapsed"]])))
  }
  fit_values()
}

#' Fit log-time conditional-model screens for high-dimensional biomarkers
#'
#' Log-time variant of [cm_screen_biomarkers()].
#' @family screening workflows
#' @param prepared Prepared screening object returned by [prepare_biomarker_screen()].
#' @param biomarker_ids Column ids in `sim_data$y` to screen.
#' @param out_file Optional output RDS path for partial progress.
#' @param coxjm_out Optional Cox/JM result block to prepend when saving progress.
#' @param save_each Whether to save partial conditional-model progress.
#' @param whiten Whether to use the whitened optimizer. Defaults to `TRUE`.
#' @param return_timing Whether to return elapsed time in seconds for each
#'   biomarker fit.
#' @param timing_out Optional timing block to append when saving progress.
#' @return A 7 x `length(biomarker_ids)` matrix of conditional-model statistics.
#'   If `return_timing = TRUE`, returns a list with `values`, `time_sec`, and
#'   `init_time_sec`.
#' @export
cm_screen_biomarkers_logT <- function(
    prepared, biomarker_ids, out_file = NULL, coxjm_out = NULL,
    save_each = FALSE, whiten = TRUE, return_timing = FALSE, timing_out = NULL) {
  prepared <- cm_lmm_ensure_biomarker_cpp(
    prepared,
    time = "logT",
    biomarker_ids = biomarker_ids
  )
  p_out <- matrix(
    NA_real_,
    length(cm_lmm_assoc_stat_names),
    length(biomarker_ids),
    dimnames = list(cm_lmm_assoc_stat_names, NULL)
  )
  time_sec <- rep(NA_real_, length(biomarker_ids))
  init_time_sec <- rep(NA_real_, length(biomarker_ids))
  for (pos in seq_along(biomarker_ids)) {
    j <- biomarker_ids[pos]
    screen_print_biomarker_progress(j, pos)
    timed <- system.time(p_i <- tryCatch(
      cm_screen_biomarker_logT(prepared, j, whiten = whiten, return_timing = TRUE),
      error = function(e) list(
        values = setNames(rep(NA_real_, length(cm_lmm_assoc_stat_names)), cm_lmm_assoc_stat_names),
        init_time_sec = NA_real_
      )
    ))
    p_out[, pos] <- p_i$values
    time_sec[pos] <- unname(timed[["elapsed"]])
    init_time_sec[pos] <- p_i$init_time_sec
    if (save_each && !is.null(out_file) && !is.null(coxjm_out)) {
      if (!is.null(timing_out)) {
        timing_out["time_cm_sec", pos] <- time_sec[pos]
        timing_out["time_cm_init_sec", pos] <- init_time_sec[pos]
        saveRDS(rbind(coxjm_out, p_out, timing_out), out_file)
      } else {
        saveRDS(rbind(coxjm_out, p_out), out_file)
      }
    }
  }
  if (return_timing) {
    return(list(
      values = p_out,
      time_sec = time_sec,
      init_time_sec = init_time_sec
    ))
  }
  p_out
}

#' Fit legacy all-R log-time CM screens for high-dimensional biomarkers
#'
#' Legacy all-R implementation retained for manual reference. Package screening
#' workflows use [cm_screen_biomarkers_logT()], the C++-backed implementation.
#' @family screening workflows
#' @param prepared Prepared screening object returned by [prepare_biomarker_screen()].
#' @param biomarker_ids Column ids in `sim_data$y` to screen.
#' @param out_file Optional output RDS path for partial progress.
#' @param coxjm_out Optional Cox/JM result block to prepend when saving progress.
#' @param save_each Whether to save partial conditional-model progress.
#' @param lc Censored-subject integration function.
#' @param return_timing Whether to return elapsed time in seconds for each
#'   biomarker fit.
#' @param timing_out Optional timing block to append when saving progress.
#' @return A 7 x `length(biomarker_ids)` matrix of conditional-model statistics.
#'   If `return_timing = TRUE`, returns a list with `values`, `time_sec`, and
#'   `init_time_sec`.
#' @export
cm_screen_biomarkers_logTR <- function(
    prepared, biomarker_ids, out_file = NULL, coxjm_out = NULL,
    save_each = FALSE, lc = lc2_logT, return_timing = FALSE, timing_out = NULL) {
  p_out <- matrix(NA_real_, 7, length(biomarker_ids))
  time_sec <- rep(NA_real_, length(biomarker_ids))
  init_time_sec <- rep(NA_real_, length(biomarker_ids))
  for (pos in seq_along(biomarker_ids)) {
    j <- biomarker_ids[pos]
    screen_print_biomarker_progress(j, pos)
    timed <- system.time(p_i <- tryCatch(
      cm_screen_biomarker_logTR(prepared, j, lc = lc, return_timing = TRUE),
      error = function(e) list(values = rep(NA_real_, 7), init_time_sec = NA_real_)
    ))
    p_out[, pos] <- p_i$values
    time_sec[pos] <- unname(timed[["elapsed"]])
    init_time_sec[pos] <- p_i$init_time_sec
    if (save_each && !is.null(out_file) && !is.null(coxjm_out)) {
      if (!is.null(timing_out)) {
        timing_out["time_cm_sec", pos] <- time_sec[pos]
        timing_out["time_cm_init_sec", pos] <- init_time_sec[pos]
        saveRDS(rbind(coxjm_out, p_out, timing_out), out_file)
      } else {
        saveRDS(rbind(coxjm_out, p_out), out_file)
      }
    }
  }
  if (return_timing) {
    return(list(
      values = p_out,
      time_sec = time_sec,
      init_time_sec = init_time_sec
    ))
  }
  p_out
}

#' Fit all three screens to a high-dimensional biomarker block with log-time CM
#'
#' Log-time variant of [fit_screen_biomarkers()]. Cox and joint-model rows are
#' unchanged; conditional-model rows use `log(Ti)` for the event-time term.
#' @family screening workflows
#' @param sim_data Dataset returned by [drawY100()] or `drawY100_logT()`.
#' @param biomarker_ids Column ids in `sim_data$y` to screen.
#' @param out_file Output RDS path.
#' @param save_each Whether to save partial conditional-model progress.
#' @param jm_random_specs Named character vector of random-effect formulas used
#'   for the longitudinal submodel before the joint-model fit.
#' @param jm_stats Character vector of joint-model event-association statistics
#'   to retain. The default keeps the historical p-value-only high-dimensional
#'   output; pass `c("est", "se", "p")` to retain estimates and standard errors.
#' @param cm_whiten Whether to use the whitened optimizer for conditional-model
#'   screens. Defaults to `TRUE`.
#' @return Invisibly returns the 13-row matrix written to `out_file`.
#' @export
fit_screen_biomarkers_logT <- function(
    sim_data, biomarker_ids, out_file, save_each = TRUE,
    jm_random_specs = jm_default_random_specs(), jm_stats = "p",
    cm_whiten = TRUE) {
  prepared <- prepare_biomarker_screen(sim_data)
  pCox <- cox_screen_biomarkers(prepared, biomarker_ids, return_timing = TRUE)
  pjm <- jm_screen_biomarkers(
    prepared, biomarker_ids,
    return_timing = TRUE,
    jm_random_specs = jm_random_specs,
    jm_stats = jm_stats
  )
  if (is.matrix(pjm$values)) {
    coxjm_out <- rbind(pCox = pCox$values, pjm$values)
  } else {
    coxjm_out <- rbind(pCox = pCox$values, pjm = pjm$values)
  }
  jm_timing <- if (is.matrix(pjm$time_sec)) {
    pjm$time_sec
  } else {
    rbind(time_jm_sec = pjm$time_sec)
  }
  timing_out <- rbind(
    time_cox_sec = pCox$time_sec,
    jm_timing,
    time_cm_sec = rep(NA_real_, length(biomarker_ids)),
    time_cm_init_sec = rep(NA_real_, length(biomarker_ids))
  )
  p_out <- cm_screen_biomarkers_logT(
    prepared, biomarker_ids,
    out_file = out_file, coxjm_out = coxjm_out, timing_out = timing_out,
    save_each = save_each, whiten = cm_whiten, return_timing = TRUE
  )
  timing_out["time_cm_sec", ] <- p_out$time_sec
  timing_out["time_cm_init_sec", ] <- p_out$init_time_sec
  out <- rbind(coxjm_out, p_out$values, timing_out)
  saveRDS(out, out_file)
  invisible(out)
}
