test_that("new cm_screen exposes the C++ screen API", {
  cm_args <- names(formals(cm_screen))
  expect_true("whiten" %in% cm_args)
  expect_identical(formals(cm_screen)$whiten, TRUE)
  expect_false("lc" %in% cm_args)
  expect_false("return_timing" %in% cm_args)

  legacy_args <- names(formals(cm_screenR))
  expect_true("lc" %in% legacy_args)
  expect_true("return_timing" %in% legacy_args)
})

test_that("low-dimensional wrappers expose whitening without legacy controls", {
  wrappers <- list(
    cm_screen_stats,
    cm_screen_fitxi,
    cm_screen_stats_logT,
    cm_screen_fitxi_logT
  )
  for (fun in wrappers) {
    args <- names(formals(fun))
    expect_true("whiten" %in% args)
    expect_identical(formals(fun)$whiten, TRUE)
    expect_false("lc" %in% args)
    expect_false("return_timing" %in% args)
  }

  expect_true("cm_whiten" %in% names(formals(fit_screen_simulations)))
  expect_identical(formals(fit_screen_simulations)$cm_whiten, TRUE)
  expect_true("cm_whiten" %in% names(formals(fit_screen_simulations_logT)))
  expect_identical(formals(fit_screen_simulations_logT)$cm_whiten, TRUE)
})

test_that("package workflows do not call legacy cm_screenR", {
  workflow_names <- c(
    "cm_screen_stats",
    "cm_screen_fitxi",
    "cm_screen_stats_logT",
    "cm_screen_fitxi_logT",
    "fit_screen_simulations",
    "fit_screen_simulations_logT"
  )
  workflow_bodies <- vapply(
    workflow_names,
    function(name) paste(deparse(body(get(name, envir = asNamespace("longscreenfun")))), collapse = "\n"),
    character(1)
  )
  expect_false(any(grepl("cm_screenR", workflow_bodies, fixed = TRUE)))
})

test_that("C++ covariance modes compute only requested covariance blocks", {
  fit_two_stage <- getFromNamespace("cm_lmm_fit_two_stage", "longscreenfun")
  sim <- CM_sim_RI(N = 30, seed = 3)

  model_fit <- fit_two_stage(
    sim,
    xi_mode = "fixxi",
    xi_fix = 1,
    xi_start = 1,
    whiten = FALSE,
    covariance = "model"
  )
  model_cov <- model_fit$cov_fixxi
  expect_identical(model_cov$covariance, "model")
  expect_true(is.matrix(model_cov$V2))
  expect_null(model_cov$V1)
  expect_null(model_cov$U)
  expect_null(model_cov$meat)

  sandwich_fit <- fit_two_stage(
    sim,
    xi_mode = "fixxi",
    xi_fix = 1,
    xi_start = 1,
    whiten = FALSE,
    covariance = "sandwich"
  )
  sandwich_cov <- sandwich_fit$cov_fixxi
  expect_identical(sandwich_cov$covariance, "sandwich")
  expect_true(is.matrix(sandwich_cov$V1))
  expect_null(sandwich_cov$V2)
  expect_true(is.matrix(sandwich_cov$U))
  expect_true(is.matrix(sandwich_cov$meat))
})

test_that("cm_screen return_fit preserves both covariance blocks", {
  sim <- CM_sim_RI(N = 30, seed = 4)
  out <- cm_screen(sim, whiten = FALSE, return_fit = TRUE)
  expect_identical(out$covariance$covariance, "both")
  expect_true(is.matrix(out$covariance$V1))
  expect_true(is.matrix(out$covariance$V2))
  expect_true(is.matrix(out$vcov))
  expect_identical(out$vcov, out$covariance$V2)
})

test_that("high-dimensional CM wrappers expose C++ API and legacy references", {
  wrappers <- list(
    cm_screen_biomarker,
    cm_screen_biomarkers,
    cm_screen_biomarker_logT,
    cm_screen_biomarkers_logT
  )
  for (fun in wrappers) {
    args <- names(formals(fun))
    expect_true("whiten" %in% args)
    expect_identical(formals(fun)$whiten, TRUE)
    expect_false("lc" %in% args)
  }

  legacy_wrappers <- list(
    cm_screen_biomarkerR,
    cm_screen_biomarkersR,
    cm_screen_biomarker_logTR,
    cm_screen_biomarkers_logTR
  )
  for (fun in legacy_wrappers) {
    expect_true("lc" %in% names(formals(fun)))
  }

  expect_true("cm_whiten" %in% names(formals(fit_screen_biomarkers)))
  expect_identical(formals(fit_screen_biomarkers)$cm_whiten, TRUE)
  expect_true("cm_whiten" %in% names(formals(fit_screen_biomarkers_logT)))
  expect_identical(formals(fit_screen_biomarkers_logT)$cm_whiten, TRUE)
})

test_that("high-dimensional package workflows do not call legacy CM wrappers", {
  workflow_names <- c(
    "cm_screen_biomarker",
    "cm_screen_biomarkers",
    "fit_screen_biomarkers",
    "cm_screen_biomarker_logT",
    "cm_screen_biomarkers_logT",
    "fit_screen_biomarkers_logT"
  )
  workflow_bodies <- vapply(
    workflow_names,
    function(name) paste(deparse(body(get(name, envir = asNamespace("longscreenfun")))), collapse = "\n"),
    character(1)
  )
  expect_false(any(grepl("cm_screen_biomarkerR", workflow_bodies, fixed = TRUE)))
  expect_false(any(grepl("cm_screen_biomarkersR", workflow_bodies, fixed = TRUE)))
  expect_false(any(grepl("cm_screen_biomarker_logTR", workflow_bodies, fixed = TRUE)))
  expect_false(any(grepl("cm_screen_biomarkers_logTR", workflow_bodies, fixed = TRUE)))
})

test_that("high-dimensional C++ prep aligns y_all after splitting", {
  sim <- drawY100(seed = 1, N = 30, p = 10)
  prepared <- prepare_biomarker_screen(sim)
  cm <- prepared$cm_cpp_T

  expect_equal(
    cm$split$datX1$y_all[, 2],
    sim$y[cm$split$datX1$row_index, 2],
    ignore_attr = TRUE
  )
  expect_equal(
    cm$split$datX0$y_all[, 2],
    sim$y[cm$split$datX0$row_index, 2],
    ignore_attr = TRUE
  )

  ensure_cpp <- getFromNamespace("cm_lmm_ensure_biomarker_cpp", "longscreenfun")
  prepared_logT <- ensure_cpp(prepared, time = "logT", biomarker_ids = 2)
  cm_logT <- prepared_logT$cm_cpp_logT
  expect_equal(
    cm_logT$split$datX1$y_all[, 2],
    sim$y[cm_logT$split$datX1$row_index, 2],
    ignore_attr = TRUE
  )
  expect_equal(
    cm_logT$split$datX0$y_all[, 2],
    sim$y[cm_logT$split$datX0$row_index, 2],
    ignore_attr = TRUE
  )
})

test_that("high-dimensional C++ CM returns seven named statistics", {
  cm_names <- c(
    "beta_T", "beta_g", "se_beta_T", "se_beta_g",
    "logp_beta_T", "logp_beta_g", "logp_joint"
  )
  sim <- drawY100(seed = 1, N = 30, p = 10)
  prepared <- prepare_biomarker_screen(sim)

  out <- cm_screen_biomarker(prepared, 2, whiten = FALSE)
  expect_named(out, cm_names)
  expect_length(out, 7)
  expect_true(all(is.finite(out)))

  out_mat <- cm_screen_biomarkers(prepared, 2, whiten = FALSE)
  expect_equal(dim(out_mat), c(7L, 1L))
  expect_identical(rownames(out_mat), cm_names)
})

test_that("high-dimensional CM progress printing is sparse and quiet", {
  sim <- drawY100(seed = 1, N = 30, p = 10)
  prepared <- prepare_biomarker_screen(sim)

  one_out <- NULL
  one_printed <- capture.output(
    one_out <- cm_screen_biomarker(prepared, 2, whiten = FALSE)
  )
  expect_length(one_out, 7)
  expect_identical(one_printed, character(0))

  batch_out <- NULL
  batch_printed <- capture.output(
    batch_out <- cm_screen_biomarkers(prepared, c(9, 10), whiten = FALSE)
  )
  expect_equal(dim(batch_out), c(7L, 2L))
  expect_identical(batch_printed, "[1] 10")
  expect_false(any(grepl("user|system|elapsed", batch_printed)))

  sim_logT <- drawY100_logT(seed = 2, N = 30, p = 10)
  prepared_logT <- prepare_biomarker_screen(sim_logT)
  logT_out <- NULL
  logT_printed <- capture.output(
    logT_out <- cm_screen_biomarker_logT(prepared_logT, 2, whiten = FALSE)
  )
  expect_length(logT_out, 7)
  expect_identical(logT_printed, character(0))

  progress <- getFromNamespace("screen_print_biomarker_progress", "longscreenfun")
  char_printed <- capture.output(progress("marker_10", 10))
  expect_identical(char_printed, "[1] \"marker_10\"")
  char_quiet <- capture.output(progress("marker_9", 9))
  expect_identical(char_quiet, character(0))
})

test_that("high-dimensional fit wrappers preserve output layout", {
  sim <- drawY100(seed = 1, N = 30, p = 10)
  out_file <- tempfile(fileext = ".rds")
  out <- fit_screen_biomarkers(
    sim,
    biomarker_ids = 2,
    out_file = out_file,
    save_each = FALSE,
    cm_whiten = TRUE
  )
  expect_equal(nrow(out), 13L)
  expect_true(file.exists(out_file))
  expect_identical(readRDS(out_file), out)

  sim_logT <- drawY100_logT(seed = 2, N = 30, p = 10)
  out_file_logT <- tempfile(fileext = ".rds")
  out_logT <- fit_screen_biomarkers_logT(
    sim_logT,
    biomarker_ids = 2,
    out_file = out_file_logT,
    save_each = FALSE,
    cm_whiten = TRUE
  )
  expect_equal(nrow(out_logT), 13L)
  expect_true(file.exists(out_file_logT))
  expect_identical(readRDS(out_file_logT), out_logT)
})
