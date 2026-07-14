#' Extract the joint-model association estimate
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @family internal numerical helpers
#' @param jmF Fitted joint model.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
ci_fit_jm1 <- function(jmF) {
  sfJ <- summary(jmF)
  nam <- row.names(sfJ$`CoefTable-Event`)
  w <- which(str_detect(nam, "Assoc"))
  ci_fit(sfJ$`CoefTable-Event`[w, 1], sfJ$`CoefTable-Event`[w, 2]) %>% mutate(
    nam = "alp",
    mod = "jm"
  )
}

#' Fit a joint model and return association inference
#'
#' Advanced helper that wraps [JM::jointModel()] and extracts association
#' inference for the screening workflows.
#' @family model fitting helpers
#' @param lmeF Fitted longitudinal mixed model.
#' @param fitS Fitted survival model.
#' @param mod_ Suffix identifying the joint-model random-effect structure.
#' @param tV Time variable name passed to the joint-model fit.
#' @param ... Additional arguments passed to [JM::jointModel()].
#' @return A tibble with association estimate, standard error, Wald interval,
#'   association name, and model label; returns an empty tibble if fitting fails.
#' @export
jm_fit_wmods <- function(lmeF, fitS, mod_ = "", tV = "tim", ...) {
  fitJOINT <- tryCatch(JM::jointModel(lmeF, fitS, timeVar = tV, ...), error = function(e) NULL)
  if (is.null(fitJOINT)) {
    return(tibble())
  }
  out <- tryCatch(ci_fit_jm1(fitJOINT), error = function(e) NULL)
  if (!is.null(out)) {
    out$mod <- paste0("jmcox-", mod_)
  }
  return(out)
}

#' Default joint-model random-effect specifications
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @return A named character vector of random-effect formulas.
#' @keywords internal
jm_default_random_specs <- function() c(b0 = "~1 | i")

#' Validate joint-model random-effect specifications
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param jm_random_specs Named character vector of random-effect formulas.
#' @return A named character vector of random-effect formulas.
#' @keywords internal
normalize_jm_random_specs <- function(jm_random_specs = jm_default_random_specs()) {
  if (!is.character(jm_random_specs) || length(jm_random_specs) == 0) {
    stop("`jm_random_specs` must be a non-empty named character vector.", call. = FALSE)
  }
  spec_names <- names(jm_random_specs)
  if (is.null(spec_names) || any(spec_names == "") || any(is.na(spec_names))) {
    stop("`jm_random_specs` must name each random-effect formula.", call. = FALSE)
  }
  if (anyDuplicated(spec_names)) {
    stop("`jm_random_specs` names must be unique.", call. = FALSE)
  }
  jm_random_specs
}

#' Validate requested joint-model statistics
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param jm_stats Character vector of statistics to retain.
#' @return A character vector containing any of `est`, `se`, and `p`.
#' @keywords internal
normalize_jm_stats <- function(jm_stats = c("est", "se", "p")) {
  unique(match.arg(jm_stats, c("est", "se", "p"), several.ok = TRUE))
}

#' Build joint-model statistic row names
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param jm_random_specs Named character vector of random-effect formulas.
#' @param jm_stats Character vector of statistics to retain.
#' @return Row names for joint-model statistic output.
#' @keywords internal
jm_stat_row_names <- function(
    jm_random_specs = jm_default_random_specs(),
    jm_stats = c("est", "se", "p")) {
  jm_random_specs <- normalize_jm_random_specs(jm_random_specs)
  jm_stats <- normalize_jm_stats(jm_stats)
  unlist(lapply(names(jm_random_specs), function(spec_name) {
    paste0("jm_", spec_name, "_", jm_stats)
  }), use.names = FALSE)
}

#' Build joint-model timing row names
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param jm_random_specs Named character vector of random-effect formulas.
#' @param legacy_single Whether to keep the historical single-JM timing row name
#'   for the default random-intercept-only fit.
#' @return Row names for joint-model timing output.
#' @keywords internal
jm_time_row_names <- function(
    jm_random_specs = jm_default_random_specs(),
    legacy_single = FALSE) {
  jm_random_specs <- normalize_jm_random_specs(jm_random_specs)
  if (legacy_single && length(jm_random_specs) == 1 && names(jm_random_specs)[1] == "b0") {
    return("time_jm_sec")
  }
  paste0("time_jm_", names(jm_random_specs), "_sec")
}

#' Fit joint models over random-effect specifications
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @param Xjm Longitudinal data frame with `i`, `tim`, and `y`.
#' @param fitS Fitted survival model.
#' @param jm_random_specs Named character vector of random-effect formulas.
#' @param jm_stats Character vector of statistics to retain.
#' @param return_timing Whether to return elapsed time in seconds for each
#'   joint-model random-effect specification.
#' @param ... Additional arguments passed to [JM::jointModel()].
#' @return A named numeric vector of requested statistics, or a list with
#'   `values` and `time_sec` if `return_timing = TRUE`.
#' @keywords internal
jm_fit_random_specs <- function(
    Xjm, fitS, jm_random_specs = jm_default_random_specs(),
    jm_stats = c("est", "se", "p"), return_timing = FALSE, ...) {
  jm_random_specs <- normalize_jm_random_specs(jm_random_specs)
  jm_stats <- normalize_jm_stats(jm_stats)
  value_names <- jm_stat_row_names(jm_random_specs, jm_stats)
  time_names <- jm_time_row_names(jm_random_specs)
  values <- setNames(rep(NA_real_, length(value_names)), value_names)
  time_sec <- setNames(rep(NA_real_, length(time_names)), time_names)

  for (spec_pos in seq_along(jm_random_specs)) {
    spec_name <- names(jm_random_specs)[spec_pos]
    spec_rows <- paste0("jm_", spec_name, "_", jm_stats)
    timed <- system.time(res <- tryCatch({
      fitLM <- lmm_fit_jm0(Xjm, jm_random_specs[[spec_pos]])
      fitJM <- jm_fit_wmods(fitLM, fitS, spec_name, ...)
      if (is.null(fitJM) || nrow(fitJM) == 0) {
        structure(setNames(rep(NA_real_, length(jm_stats)), jm_stats), fit_ok = FALSE)
      } else {
        est <- fitJM$est[1]
        se <- fitJM$se[1]
        structure(c(est = est, se = se, p = pfun_bse(est, se))[jm_stats], fit_ok = TRUE)
      }
    }, error = function(e) {
      structure(setNames(rep(NA_real_, length(jm_stats)), jm_stats), fit_ok = FALSE)
    }))

    values[spec_rows] <- as.numeric(res[jm_stats])
    if (isTRUE(attr(res, "fit_ok"))) {
      time_sec[spec_pos] <- unname(timed["elapsed"])
    }
  }

  if (return_timing) {
    return(list(values = values, time_sec = time_sec))
  }
  values
}

#' Fit the longitudinal submodel for the joint model
#'
#' Fits the longitudinal mixed model used before the JM screen, retrying with
#' the alternate ML/REML method if the first fit fails.
#' @family model fitting helpers
#' @param dat Longitudinal data frame or prepared data list.
#' @param rand_chr Random-effect formula as a character string.
#' @param defm Mixed-model fitting method.
#' @return An `nlme::lme` fit or `NULL` if both fitting attempts fail.
#' @export
lmm_fit_jm0 <- function(dat, rand_chr = "~1 | i", defm = "REML") {
  lm_fit <- tryCatch(lmm_fit_jm0_m(dat, rand_chr, defm), error = function(e) NULL)
  if (is.null(lm_fit)) {
    if (defm == "REML") {
      defm2 <- "ML"
    }
    if (defm == "ML") {
      defm2 <- "REML"
    }
    lm_fit <- tryCatch(lmm_fit_jm0_m(dat, rand_chr, defm2), error = function(e) NULL)
  }
  lm_fit
}

#' Run one longitudinal mixed-model fit
#'
#' Helper retained for the longbioscreen simulation pipeline.
#' @family internal numerical helpers
#' @param dat Longitudinal data frame or prepared data list.
#' @param rand_chr Random-effect formula as a character string.
#' @param defm Mixed-model fitting method.
#' @return A value used by the longbioscreen simulation pipeline.
#' @keywords internal
lmm_fit_jm0_m <- function(dat, rand_chr = "~1 | i", defm = "REML") {
  r_f <- as.formula(rand_chr)
  nlme::lme(y ~ tim,
    random = r_f, dat, control = nlme::lmeControl(maxIter = 500),
    method = defm
  )
}
