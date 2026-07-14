# Conditional-model random-intercept LMM implementation used by the C++ screen.

cm_lmm_theta_names <- c(
  "xi", "beta0", "beta1", "beta_T", "beta_g", "log_sigma_b", "log_sigma_e"
)

cm_lmm_assoc_stat_names <- c(
  "beta_T", "beta_g", "se_beta_T", "se_beta_g",
  "logp_beta_T", "logp_beta_g", "logp_joint"
)

cm_lmm_ftrans <- function(Ti) {
  Ti
}

cm_lmm_gfun <- function(xi, S, tx) {
  gfun(1, xi, S, tx)
}

cm_control <- function(ftrans_fn = cm_lmm_ftrans,
                       candidate_times = 0:5) {
  stopifnot(is.function(ftrans_fn))
  structure(
    list(
      ftrans_fn = ftrans_fn,
      candidate_times = as.double(candidate_times)
    ),
    class = "lmm_cm_control"
  )
}

cm_normalize_control <- function(control = NULL, candidate_times = NULL) {
  if (is.null(control)) {
    control <- cm_control()
  }
  if (!is.list(control)) {
    stop("control must be a cm_control() list", call. = FALSE)
  }
  if (is.null(control$ftrans_fn)) {
    control$ftrans_fn <- cm_lmm_ftrans
  }
  if (is.null(control$candidate_times)) {
    control$candidate_times <- 0:5
  }
  if (!is.null(candidate_times)) {
    control$candidate_times <- candidate_times
  }
  stopifnot(is.function(control$ftrans_fn))
  control$candidate_times <- as.double(control$candidate_times)
  class(control) <- unique(c("lmm_cm_control", class(control)))
  control
}

cm_get_control <- function(sim = NULL, control = NULL, candidate_times = NULL) {
  if (is.null(control) && !is.null(sim) && !is.null(sim$control)) {
    control <- sim$control
  }
  cm_normalize_control(control, candidate_times = candidate_times)
}

cm_get_ftrans <- function(sim = NULL, control = NULL) {
  cm_get_control(sim = sim, control = control)$ftrans_fn
}

cm_lmm_truth <- function(k = 2,
                         lam = 10,
                         xi = 1,
                         beta0 = 0.3,
                         beta1 = -0.25,
                         beta_T = -0.1,
                         beta_g = 1.2,
                         sigma_b = 1,
                         sigma_e = 0.5,
                         B = NULL,
                         Btau = NULL,
                         alp = NULL,
                         sb = NULL,
                         se = NULL) {
  if (!is.null(B)) {
    stopifnot(length(B) >= 2)
    beta0 <- B[1]
    beta1 <- B[2]
  }
  if (!is.null(Btau)) {
    beta_g <- Btau
  }
  if (!is.null(alp)) {
    beta_T <- alp
  }
  if (!is.null(sb)) {
    stopifnot(length(sb) >= 1)
    sigma_b <- sb[1]
  }
  if (!is.null(se)) {
    sigma_e <- se
  }

  list(
    k = k,
    lam = lam,
    log_klam = c(log_k = log(k), log_lam = log(lam)),
    xi = xi,
    beta0 = beta0,
    beta1 = beta1,
    beta_T = beta_T,
    beta_t = beta_T,
    beta_g = beta_g,
    sigma_b = sigma_b,
    sigma_e = sigma_e,
    theta_y = c(
      xi = xi,
      beta0 = beta0,
      beta1 = beta1,
      beta_T = beta_T,
      beta_g = beta_g,
      log_sigma_b = log(sigma_b),
      log_sigma_e = log(sigma_e)
    )
  )
}

cm_lmm_is_longscreen_data <- function(x) {
  if (is.list(x) && length(x) == 1L && is.list(x[[1L]]) &&
      all(c("X", "S") %in% names(x[[1L]]))) {
    return(TRUE)
  }
  is.list(x) && all(c("X", "S") %in% names(x))
}

cm_lmm_as_sim <- function(x,
                          truth = cm_lmm_truth(),
                          candidate_times = NULL,
                          control = cm_control()) {
  if (!cm_lmm_is_longscreen_data(x)) {
    return(x)
  }
  if (is.list(x) && length(x) == 1L && is.list(x[[1L]]) &&
      all(c("X", "S") %in% names(x[[1L]]))) {
    x <- x[[1L]]
  }

  X <- as.data.frame(x$X)
  S <- as.data.frame(x$S)
  X$.cm_row_index <- seq_len(nrow(X))
  required_X <- c("i", "tim", "y")
  required_S <- c("i", "Ti", "del")
  stopifnot(all(required_X %in% names(X)), all(required_S %in% names(S)))

  subject_id <- as.character(S$i)
  if (anyDuplicated(subject_id)) {
    stop("S$i must contain one row per subject", call. = FALSE)
  }
  X$id_int <- match(as.character(X$i), subject_id)
  if (anyNA(X$id_int)) {
    stop("all X$i values must appear in S$i", call. = FALSE)
  }

  X <- X[order(X$id_int, X$tim), , drop = FALSE]
  J <- nrow(S)
  n_i <- as.integer(tabulate(X$id_int, nbins = J))
  if (any(n_i == 0L)) {
    stop("all subjects must have at least one longitudinal row", call. = FALSE)
  }

  Ti_til <- as.double(S$Ti)
  del <- as.integer(S$del)
  id <- as.integer(X$id_int)
  tim <- as.double(X$tim)
  y <- as.double(X$y)
  row_index <- as.integer(X$.cm_row_index)
  Ti_row <- Ti_til[id]

  cluster_end <- cumsum(n_i)
  cluster_start <- cluster_end - n_i + 1L
  if (is.null(candidate_times)) {
    candidate_times <- sort(unique(tim))
  }
  control <- cm_normalize_control(control, candidate_times = candidate_times)
  ftrans_fn <- control$ftrans_fn

  data <- data.frame(
    id = factor(id, levels = seq_len(J)),
    id_int = id,
    row_index = row_index,
    y = y,
    tim = tim,
    t = tim,
    Ti = as.double(Ti_row),
    fTi = as.double(ftrans_fn(Ti_row)),
    Ti_til = as.double(Ti_til[id]),
    del = as.integer(del[id]),
    fT = as.double(ftrans_fn(Ti_row)),
    g = as.double(cm_lmm_gfun(truth$xi, Ti_row, tim))
  )

  individual_data <- data.frame(
    id = factor(seq_len(J), levels = seq_len(J)),
    id_int = seq_len(J),
    Ti = as.double(Ti_til),
    fTi = as.double(ftrans_fn(Ti_til)),
    C = as.double(Ti_til),
    Ti_til = as.double(Ti_til),
    del = as.integer(del),
    n_i = as.integer(n_i)
  )

  list(
    y = y,
    tim = tim,
    t = tim,
    Ti = as.double(Ti_til),
    fTi = as.double(ftrans_fn(Ti_til)),
    C = as.double(Ti_til),
    Ti_til = as.double(Ti_til),
    del = as.integer(del),
    id = id,
    row_index = row_index,
    n_i = as.integer(n_i),
    cluster_start = as.integer(cluster_start),
    cluster_end = as.integer(cluster_end),
    candidate_times = as.double(candidate_times),
    control = control,
    data = data,
    individual_data = individual_data,
    truth = truth
  )
}

cm_lmm_make_GK15 <- function(rule = gk15_rule()) {
  list(lwk = log(rule$wk), sk_u = (rule$sk + 1) / 2)
}

cm_lmm_mod_pts <- function(l, u, GK15 = cm_lmm_make_GK15()) {
  u_Ti <- GK15$sk_u * (u - l) + l
  lw <- GK15$lwk + log((u - l) / 2)
  list(Ti = u_Ti, lw = lw)
}

cm_lmm_set_quadpts_Ti <- function(Ti_til, GK15, pS, qS) {
  Ti_cuts <- pS(Ti_til + c(0, 2))
  pI1 <- cm_lmm_mod_pts(Ti_cuts[1], Ti_cuts[2], GK15)
  pI2 <- cm_lmm_mod_pts(Ti_cuts[2], 1, GK15)
  list(
    Ti = qS(c(pI1$Ti, pI2$Ti)),
    lw = c(pI1$lw, pI2$lw)
  )
}

cm_lmm_set_quadpts_Ti_weib <- function(Ti_til,
                                       GK15 = cm_lmm_make_GK15(),
                                       log_klam) {
  stopifnot(length(log_klam) == 2, all(is.finite(log_klam)))

  n_nodes <- 2L * length(GK15$sk_u)
  if (length(Ti_til) == 0) {
    return(list(
      Ti_mat = matrix(numeric(0), nrow = n_nodes, ncol = 0),
      lw_mat = matrix(numeric(0), nrow = n_nodes, ncol = 0)
    ))
  }

  k <- exp(log_klam[1])
  lam <- exp(log_klam[2])
  pS <- function(ti) pweibull(ti, shape = k, scale = lam)
  qS <- function(ui) qweibull(ui, shape = k, scale = lam)
  all_pts <- lapply(Ti_til, cm_lmm_set_quadpts_Ti, GK15 = GK15, pS = pS, qS = qS)

  list(
    Ti_mat = vapply(all_pts, function(x) x$Ti, numeric(n_nodes)),
    lw_mat = vapply(all_pts, function(x) x$lw, numeric(n_nodes))
  )
}

cm_lmm_make_group <- function(sim, ids, include_Ti,
                              ftrans_fn = cm_get_ftrans(sim)) {
  ids <- as.integer(ids)
  rows <- if (length(ids) == 0) {
    integer(0)
  } else {
    unlist(lapply(ids, function(i) sim$cluster_start[i]:sim$cluster_end[i]),
      use.names = FALSE
    )
  }

  n_i <- as.integer(sim$n_i[ids])
  cluster_end <- cumsum(n_i)
  cluster_start <- cluster_end - n_i + 1L

  group <- list(
    J = length(ids),
    id = ids,
    row_id = as.integer(sim$id[rows]),
    row_index = as.integer(sim$row_index[rows]),
    y = as.double(sim$y[rows]),
    tim = as.double(sim$tim[rows]),
    t = as.double(sim$tim[rows]),
    cluster_start = as.integer(cluster_start),
    cluster_end = as.integer(cluster_end),
    n_i = n_i,
    Ti_til = as.double(sim$Ti_til[ids]),
    del = as.integer(sim$del[ids])
  )

  group$data <- data.frame(
    id = factor(group$row_id, levels = seq_along(sim$Ti)),
    id_int = group$row_id,
    row_index = group$row_index,
    y = group$y,
    tim = group$tim,
    t = group$t,
    Ti_til = as.double(sim$Ti_til[group$row_id]),
    del = as.integer(sim$del[group$row_id])
  )

  if (include_Ti) {
    group$Ti <- as.double(sim$Ti[ids])
    group$fTi <- as.double(ftrans_fn(group$Ti))
    group$Ti_row <- as.double(sim$Ti[group$row_id])
    group$fTi_row <- as.double(ftrans_fn(group$Ti_row))
    group$data$Ti <- group$Ti_row
    group$data$fTi <- group$fTi_row
  }

  group
}

cm_lmm_check_quad_alignment <- function(datX0) {
  stopifnot(!is.null(datX0$Ti_mat), !is.null(datX0$lw_mat), !is.null(datX0$fTi_mat))
  stopifnot(
    ncol(datX0$Ti_mat) == datX0$J,
    ncol(datX0$lw_mat) == datX0$J,
    ncol(datX0$fTi_mat) == datX0$J,
    nrow(datX0$Ti_mat) == nrow(datX0$lw_mat),
    nrow(datX0$Ti_mat) == nrow(datX0$fTi_mat)
  )

  if (datX0$J > 0) {
    ok <- vapply(seq_len(datX0$J), function(j) {
      all(datX0$Ti_mat[, j] > datX0$Ti_til[j])
    }, logical(1))
    if (!all(ok)) {
      stop("quadrature node/truncation alignment failed", call. = FALSE)
    }
  }

  TRUE
}

cm_lmm_split_data <- function(sim,
                              log_klam = NULL,
                              GK15 = cm_lmm_make_GK15(),
                              attach_quad = TRUE,
                              ftrans_fn = cm_get_ftrans(sim)) {
  if (is.null(log_klam)) {
    log_klam <- sim$truth$log_klam
  }

  ids1 <- which(sim$del == 1L)
  ids0 <- which(sim$del == 0L)

  datX1 <- cm_lmm_make_group(sim, ids1, include_Ti = TRUE, ftrans_fn = ftrans_fn)
  datX0 <- cm_lmm_make_group(sim, ids0, include_Ti = FALSE, ftrans_fn = ftrans_fn)

  if (attach_quad) {
    quad <- cm_lmm_set_quadpts_Ti_weib(datX0$Ti_til, GK15 = GK15, log_klam = log_klam)
    datX0$Ti_mat <- quad$Ti_mat
    datX0$lw_mat <- quad$lw_mat
    datX0$fTi_mat <- ftrans_fn(quad$Ti_mat)
    datX0$quad_log_klam <- log_klam
    cm_lmm_check_quad_alignment(datX0)
  }

  list(datX1 = datX1, datX0 = datX0, log_klam = log_klam, GK15 = GK15)
}

cm_lmm_surv_ll_vec <- function(phi, TT, deli) {
  k <- exp(phi[1])
  lam <- exp(phi[2])
  ev <- dweibull(TT, shape = k, scale = lam, log = TRUE)
  cn <- pweibull(TT, shape = k, scale = lam, lower.tail = FALSE, log.p = TRUE)
  ifelse(deli == 1L, ev, cn)
}

cm_lmm_surv_score_mat <- function(phi, TT, deli) {
  k <- exp(phi[1])
  lam <- exp(phi[2])
  u <- log(TT) - log(lam)
  q <- exp(k * u)

  score_log_k <- ifelse(
    deli == 1L,
    1 + k * u * (1 - q),
    -q * k * u
  )
  score_log_lam <- ifelse(
    deli == 1L,
    k * (q - 1),
    k * q
  )

  cbind(log_k = score_log_k, log_lam = score_log_lam)
}

cm_lmm_fit_stage1_weibull <- function(Ti_til, del,
                                      start = log(c(k = 2, lam = 10))) {
  dev <- function(phi) -sum(cm_lmm_surv_ll_vec(phi, Ti_til, del))
  gr <- function(phi) -colSums(cm_lmm_surv_score_mat(phi, Ti_til, del))

  fit <- optim(
    start,
    dev,
    gr,
    method = "L-BFGS-B",
    lower = log(c(0.05, 0.5)),
    upper = log(c(20, 80))
  )

  phi_hat <- fit$par
  s1 <- cm_lmm_surv_score_mat(phi_hat, Ti_til, del)
  H <- numDeriv::jacobian(
    function(phi) colSums(cm_lmm_surv_score_mat(phi, Ti_til, del)),
    phi_hat
  )
  H <- (H + t(H)) / 2

  list(
    fit = fit,
    phi = phi_hat,
    estimate = c(k = exp(phi_hat[1]), lam = exp(phi_hat[2])),
    score = s1,
    H_phiphi = H,
    logLik = sum(cm_lmm_surv_ll_vec(phi_hat, Ti_til, del))
  )
}

cm_lmm_build_cpp_handles <- function(split) {
  datX1 <- split$datX1
  datX0 <- split$datX0

  long <- setdata_03_lmm_CM(
    as.double(datX1$y),
    as.double(datX1$tim),
    as.double(datX1$Ti),
    as.double(datX1$fTi),
    as.integer(datX1$cluster_start),
    as.integer(datX1$cluster_end),
    as.double(datX0$y),
    as.double(datX0$tim),
    as.double(datX0$Ti_til),
    as.integer(datX0$cluster_start),
    as.integer(datX0$cluster_end)
  )

  quad <- setquad_03_lmm_CM(
    datX0$Ti_mat,
    datX0$lw_mat,
    datX0$fTi_mat,
    as.double(datX0$Ti_til)
  )

  list(long = long, quad = quad)
}

cm_lmm_stage2_obj <- function(theta, data) {
  ll_R_cpp(as.double(theta), data$long, data$quad)
}

cm_lmm_stage2_vec <- function(theta, data) {
  ll_R_cpp_vec(as.double(theta), data$long, data$quad)
}

cm_lmm_theta_from_eta <- function(eta) {
  theta <- eta
  theta[1] <- exp(eta[1])
  names(theta) <- cm_lmm_theta_names
  theta
}

cm_lmm_eta_from_theta <- function(theta) {
  eta <- theta
  eta[1] <- log(theta[1])
  names(eta) <- c("log_xi", cm_lmm_theta_names[-1])
  eta
}

cm_lmm_stage2_obj_eta <- function(eta, data) {
  theta <- cm_lmm_theta_from_eta(eta)
  out <- cm_lmm_stage2_obj(theta, data)
  grad <- out$gradient
  grad[1] <- grad[1] * theta[1]
  names(grad) <- names(eta)
  list(value = out$value, gradient = grad)
}

cm_lmm_stage2_obj_fixxi <- function(theta_free, data, xi_fix = 1) {
  out <- cm_lmm_stage2_obj(c(xi_fix, theta_free), data)
  list(value = out$value, gradient = out$gradient[-1])
}

cm_lmm_stage2_start_Ti <- function(Ti_til, del) {
  Ti_til <- as.double(Ti_til)
  del <- as.integer(del)
  stopifnot(length(Ti_til) == length(del), length(Ti_til) > 0)
  ifelse(del == 1L, Ti_til, max(Ti_til) + 2)
}

cm_lmm_stage2_starting_theta <- function(sim, xi_start = 1,
                                         ftrans_fn = cm_get_ftrans(sim)) {
  df <- sim$data
  Ti_init <- cm_lmm_stage2_start_Ti(df$Ti_til, df$del)
  df$fT_init <- ftrans_fn(Ti_init)
  df$g_init <- cm_lmm_gfun(xi_start, Ti_init, df$t)

  fit_lm <- stats::lm(y ~ 1 + t + fT_init + g_init, data = df)
  beta <- stats::coef(fit_lm)[c("(Intercept)", "t", "fT_init", "g_init")]
  sigma_b <- 0.8
  sigma_e <- sigma(fit_lm)

  if (requireNamespace("lme4", quietly = TRUE)) {
    fit_lmer <- try(
      lme4::lmer(
        y ~ 1 + t + fT_init + g_init + (1 | id),
        data = df,
        REML = FALSE
      ),
      silent = TRUE
    )
    if (!inherits(fit_lmer, "try-error")) {
      beta <- lme4::fixef(fit_lmer)[c("(Intercept)", "t", "fT_init", "g_init")]
      vc <- as.data.frame(lme4::VarCorr(fit_lmer))
      sigma_b <- vc$sdcor[vc$grp == "id" & vc$var1 == "(Intercept)"]
      sigma_e <- sigma(fit_lmer)
    }
  }

  theta <- c(
    xi = xi_start,
    beta0 = beta[[1]],
    beta1 = beta[[2]],
    beta_T = beta[[3]],
    beta_g = beta[[4]],
    log_sigma_b = log(max(sigma_b, 0.05)),
    log_sigma_e = log(max(sigma_e, 0.05))
  )
  names(theta) <- cm_lmm_theta_names
  theta
}

cm_lmm_dev_4opt <- function(obj_fn, data, maximize = TRUE) {
  last_th <- NULL
  l_b <- NULL
  sign <- if (maximize) -1 else 1

  refresh <- function(th) {
    if (is.null(last_th) || any(th != last_th)) {
      l_b <<- obj_fn(th, data)
      stopifnot(
        "obj_fn must return a list" = is.list(l_b),
        "obj_fn's result must have a 'value' element" = !is.null(l_b$value),
        "obj_fn's result must have a 'gradient' element" = !is.null(l_b$gradient)
      )
      last_th <<- th
    }
  }

  list(
    fn = function(th) {
      refresh(th)
      sign * l_b$value
    },
    gr = function(th) {
      refresh(th)
      sign * l_b$gradient
    }
  )
}

cm_lmm_make_stage2_whitener <- function(obj_fn,
                                        data,
                                        par_start,
                                        ridge = 1e-6) {
  H <- numDeriv::jacobian(
    function(par) obj_fn(par, data)$gradient,
    par_start,
    method = "simple"
  )
  info <- -(H + t(H)) / 2
  eig <- eigen(info, symmetric = TRUE)
  floor_val <- max(max(eig$values), 1) * ridge
  eig_values_pd <- pmax(eig$values, floor_val)
  info_pd <- eig$vectors %*% diag(eig_values_pd, nrow = length(eig_values_pd)) %*%
    t(eig$vectors)
  R <- chol(info_pd)

  from_z <- function(z) {
    par <- par_start + as.numeric(backsolve(R, z))
    names(par) <- names(par_start)
    par
  }
  to_z <- function(par) {
    z <- as.numeric(R %*% (par - par_start))
    names(z) <- paste0("z_", names(par_start))
    z
  }

  list(
    par_start = par_start,
    H_start = H,
    info_start = info,
    info_pd = info_pd,
    eig_values = eig$values,
    eig_values_pd = eig_values_pd,
    R = R,
    from_z = from_z,
    to_z = to_z,
    z_start = to_z(par_start)
  )
}

cm_lmm_make_stage2_whitened_obj <- function(obj_fn,
                                            data,
                                            par_start,
                                            ridge = 1e-6,
                                            lower = rep(-Inf, length(par_start)),
                                            upper = rep(Inf, length(par_start)),
                                            boundary_penalty = 1e8) {
  whitener <- cm_lmm_make_stage2_whitener(obj_fn, data, par_start, ridge = ridge)

  obj <- function(z, ignored_data = NULL) {
    par <- whitener$from_z(z)
    below <- pmax(lower - par, 0)
    above <- pmax(par - upper, 0)
    if (any(!is.finite(par)) || any(below > 0) || any(above > 0)) {
      grad_par <- 2 * boundary_penalty * below - 2 * boundary_penalty * above
      grad_z <- as.numeric(backsolve(whitener$R, grad_par, transpose = TRUE))
      names(grad_z) <- names(z)
      return(list(
        value = -boundary_penalty * (1 + sum(below^2) + sum(above^2)),
        gradient = grad_z
      ))
    }

    out <- obj_fn(par, data)
    grad_z <- as.numeric(backsolve(whitener$R, out$gradient, transpose = TRUE))
    names(grad_z) <- names(z)
    list(value = out$value, gradient = grad_z)
  }

  list(obj = obj, whitener = whitener, z_start = whitener$z_start)
}

cm_lmm_fit_stage2_fitxi <- function(cpp_data,
                                    theta_start,
                                    whiten = TRUE,
                                    whiten_ridge = 1e-6,
                                    fallback = TRUE) {
  eta0 <- cm_lmm_eta_from_theta(theta_start)
  lower <- rep(-Inf, length(eta0))
  upper <- rep(Inf, length(eta0))
  lower[1] <- log(1e-3)
  upper[1] <- log(10)
  lower[6:7] <- -20

  if (isTRUE(whiten)) {
    fit_whitened <- try({
      w <- cm_lmm_make_stage2_whitened_obj(
        cm_lmm_stage2_obj_eta,
        cpp_data,
        eta0,
        ridge = whiten_ridge,
        lower = lower,
        upper = upper
      )
      dev <- cm_lmm_dev_4opt(w$obj, NULL)
      fit <- optim(
        w$z_start,
        dev$fn,
        dev$gr,
        method = "BFGS",
        control = list(maxit = 500, reltol = 1e-9)
      )
      eta_hat <- w$whitener$from_z(fit$par)
      theta_hat <- cm_lmm_theta_from_eta(eta_hat)
      fit$z <- fit$par
      fit$eta <- eta_hat
      fit$par <- theta_hat
      fit$whitened <- TRUE
      fit$whitener <- w$whitener
      fit
    }, silent = TRUE)

    if (!inherits(fit_whitened, "try-error") &&
        identical(fit_whitened$convergence, 0L) &&
        all(is.finite(fit_whitened$par))) {
      return(fit_whitened)
    }
    if (!isTRUE(fallback)) {
      stop("whitened fit-xi optimization failed", call. = FALSE)
    }
    whitened_attempt <- fit_whitened
  } else {
    whitened_attempt <- NULL
  }

  dev <- cm_lmm_dev_4opt(cm_lmm_stage2_obj_eta, cpp_data)
  fit <- optim(
    eta0,
    dev$fn,
    dev$gr,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 500, pgtol = 1e-8, factr = 1e5)
  )

  theta_hat <- cm_lmm_theta_from_eta(fit$par)
  fit$eta <- fit$par
  fit$par <- theta_hat
  fit$whitened <- FALSE
  if (!is.null(whitened_attempt)) {
    fit$whitened_attempt <- whitened_attempt
  }
  fit
}

cm_lmm_fit_stage2_fixxi <- function(cpp_data,
                                    theta_start,
                                    xi_fix = 1,
                                    whiten = TRUE,
                                    whiten_ridge = 1e-6,
                                    fallback = TRUE) {
  free0 <- theta_start[-1]
  lower <- rep(-Inf, length(free0))
  lower[5:6] <- -20

  obj <- function(theta_free, data) {
    cm_lmm_stage2_obj_fixxi(theta_free, data, xi_fix)
  }

  if (isTRUE(whiten)) {
    fit_whitened <- try({
      w <- cm_lmm_make_stage2_whitened_obj(
        obj,
        cpp_data,
        free0,
        ridge = whiten_ridge,
        lower = lower
      )
      dev <- cm_lmm_dev_4opt(w$obj, NULL)
      fit <- optim(
        w$z_start,
        dev$fn,
        dev$gr,
        method = "BFGS",
        control = list(maxit = 500, reltol = 1e-9)
      )
      free_hat <- w$whitener$from_z(fit$par)
      fit$z <- fit$par
      fit$par_free <- free_hat
      fit$par <- c(xi = xi_fix, fit$par_free)
      names(fit$par) <- cm_lmm_theta_names
      fit$whitened <- TRUE
      fit$whitener <- w$whitener
      fit
    }, silent = TRUE)

    if (!inherits(fit_whitened, "try-error") &&
        identical(fit_whitened$convergence, 0L) &&
        all(is.finite(fit_whitened$par))) {
      return(fit_whitened)
    }
    if (!isTRUE(fallback)) {
      stop("whitened fixed-xi optimization failed", call. = FALSE)
    }
    whitened_attempt <- fit_whitened
  } else {
    whitened_attempt <- NULL
  }

  dev <- cm_lmm_dev_4opt(obj, cpp_data)
  fit <- optim(
    free0,
    dev$fn,
    dev$gr,
    method = "L-BFGS-B",
    lower = lower,
    control = list(maxit = 500, pgtol = 1e-8, factr = 1e5)
  )

  fit$par_free <- fit$par
  fit$par <- c(xi = xi_fix, fit$par_free)
  names(fit$par) <- cm_lmm_theta_names
  fit$whitened <- FALSE
  if (!is.null(whitened_attempt)) {
    fit$whitened_attempt <- whitened_attempt
  }
  fit
}

cm_lmm_align_stage2_scores <- function(score_compact, split, J) {
  p <- ncol(score_compact)
  U <- matrix(0, nrow = J, ncol = p)
  colnames(U) <- colnames(score_compact)

  J1 <- split$datX1$J
  if (J1 > 0) {
    U[split$datX1$id, ] <- score_compact[seq_len(J1), , drop = FALSE]
  }
  if (split$datX0$J > 0) {
    idx0 <- J1 + seq_len(split$datX0$J)
    U[split$datX0$id, ] <- score_compact[idx0, , drop = FALSE]
  }

  U
}

cm_lmm_stage2_gradient_at_phi <- function(theta, phi, sim, split_template,
                                          long_ptr,
                                          ftrans_fn = cm_get_ftrans(sim)) {
  quad <- cm_lmm_set_quadpts_Ti_weib(
    split_template$datX0$Ti_til,
    GK15 = split_template$GK15,
    log_klam = phi
  )
  quad_ptr <- setquad_03_lmm_CM(
    quad$Ti_mat,
    quad$lw_mat,
    ftrans_fn(quad$Ti_mat),
    as.double(split_template$datX0$Ti_til)
  )
  cm_lmm_stage2_obj(theta, list(long = long_ptr, quad = quad_ptr))$gradient
}

cm_lmm_stage2_covariance <- function(theta_hat,
                                     phi_hat,
                                     sim,
                                     split,
                                     cpp_data,
                                     stage1,
                                     fixed_xi = FALSE,
                                     s1_override = NULL,
                                     covariance = c("both", "model", "sandwich")) {
  covariance <- match.arg(covariance)
  need_model <- covariance %in% c("both", "model")
  need_sandwich <- covariance %in% c("both", "sandwich")

  A <- numDeriv::jacobian(
    function(theta) cm_lmm_stage2_obj(theta, cpp_data)$gradient,
    theta_hat
  )
  A <- (A + t(A)) / 2
  rownames(A) <- colnames(A) <- cm_lmm_theta_names

  H_thphi <- numDeriv::jacobian(
    function(phi) {
      cm_lmm_stage2_gradient_at_phi(
        theta = theta_hat,
        phi = phi,
        sim = sim,
        split_template = split,
        long_ptr = cpp_data$long,
        ftrans_fn = cm_get_ftrans(sim)
      )
    },
    phi_hat
  )
  rownames(H_thphi) <- cm_lmm_theta_names
  colnames(H_thphi) <- c("log_k", "log_lam")

  H_phiphi <- stage1$H_phiphi

  idx <- if (fixed_xi) 2:7 else 1:7
  A_use <- A[idx, idx, drop = FALSE]
  H_use <- H_thphi[idx, , drop = FALSE]

  bread <- solve(-A_use)
  V_naive_model <- bread

  out <- list(
    A_thth = A,
    H_phiphi = H_phiphi,
    H_thphi = H_thphi,
    bread = bread,
    V_naive_model = V_naive_model,
    se_naive_model = sqrt(diag(V_naive_model)),
    fixed_xi = fixed_xi,
    free_index = idx,
    covariance = covariance
  )

  if (need_model) {
    V2 <- bread + bread %*% H_use %*% solve(-H_phiphi) %*% t(H_use) %*% bread
    out$V2 <- V2
    out$se_V2 <- sqrt(diag(V2))
  }

  if (need_sandwich) {
    score_compact <- numDeriv::jacobian(
      function(theta) cm_lmm_stage2_vec(theta, cpp_data),
      theta_hat
    )
    colnames(score_compact) <- cm_lmm_theta_names
    U <- cm_lmm_align_stage2_scores(score_compact, split, length(sim$Ti))
    U_use <- U[, idx, drop = FALSE]
    s1 <- if (is.null(s1_override)) stage1$score else s1_override
    adjust <- t(H_use %*% solve(H_phiphi, t(s1)))
    Utilde <- U_use - adjust
    meat <- crossprod(Utilde)
    V1 <- bread %*% meat %*% bread

    out$U <- U
    out$s1 <- s1
    out$Utilde <- Utilde
    out$meat <- meat
    out$V1 <- V1
    out$se_V1 <- sqrt(diag(V1))

    if (covariance == "both") {
      naive_meat <- crossprod(U_use)
      V_naive_sandwich <- bread %*% naive_meat %*% bread
      out$naive_meat <- naive_meat
      out$V_naive_sandwich <- V_naive_sandwich
      out$naive <- V_naive_sandwich
      out$se_naive_sandwich <- sqrt(diag(V_naive_sandwich))
    }
  }

  out
}

cm_lmm_assoc_inference <- function(theta,
                                   V,
                                   assoc_names = c("beta_T", "beta_g")) {
  theta <- unlist(theta)
  storage.mode(theta) <- "double"
  if (is.null(names(theta))) {
    stop("theta must be named", call. = FALSE)
  }
  if (is.null(rownames(V)) || is.null(colnames(V))) {
    stop("V must have row and column names", call. = FALSE)
  }
  if (!all(assoc_names %in% names(theta))) {
    stop("theta is missing association parameters", call. = FALSE)
  }
  if (!all(assoc_names %in% rownames(V)) ||
      !all(assoc_names %in% colnames(V))) {
    stop("V is missing association parameters", call. = FALSE)
  }

  est <- theta[assoc_names]
  V_assoc <- V[assoc_names, assoc_names, drop = FALSE]
  se <- sqrt(diag(V_assoc))
  logp <- log(2) + pnorm(-abs(est / se), log.p = TRUE)
  R <- chol(V_assoc)
  q <- sum(backsolve(R, est, transpose = TRUE)^2)
  logp_joint <- pchisq(
    q,
    df = length(assoc_names),
    lower.tail = FALSE,
    log.p = TRUE
  )

  out <- c(est, se, logp, logp_joint)
  names(out) <- c(
    assoc_names,
    paste0("se_", assoc_names),
    paste0("logp_", assoc_names),
    "logp_joint"
  )
  out
}

cm_lmm_fit_two_stage <- function(sim,
                                 xi_mode = c("fitxi", "fixxi"),
                                 xi_fix = 1,
                                 xi_start = 1,
                                 whiten = TRUE,
                                 whiten_ridge = 1e-6,
                                 compute_covariance = TRUE,
                                 covariance = c("both", "model", "sandwich"),
                                 longscreen_truth = cm_lmm_truth(),
                                 control = cm_control()) {
  xi_mode <- match.arg(xi_mode)
  covariance <- match.arg(covariance)
  sim <- cm_lmm_as_sim(sim, truth = longscreen_truth, control = control)

  stage1_start <- if (!is.null(sim$truth$log_klam)) {
    sim$truth$log_klam
  } else {
    log(c(k = 2, lam = 10))
  }
  stage1 <- cm_lmm_fit_stage1_weibull(sim$Ti_til, sim$del, start = stage1_start)
  split <- cm_lmm_split_data(sim, log_klam = stage1$phi)
  cpp_data <- cm_lmm_build_cpp_handles(split)
  theta_start <- cm_lmm_stage2_starting_theta(sim, xi_start = xi_start)

  res <- list(
    sim = sim,
    xi_mode = xi_mode,
    xi_fix = xi_fix,
    xi_start = xi_start,
    stage1 = stage1,
    split = split,
    cpp_data = cpp_data,
    theta_start = theta_start
  )

  if (xi_mode == "fitxi") {
    fit_y <- cm_lmm_fit_stage2_fitxi(
      cpp_data,
      theta_start,
      whiten = whiten,
      whiten_ridge = whiten_ridge
    )
    res$fit_y <- fit_y
    if (isTRUE(compute_covariance)) {
      res$cov_y <- cm_lmm_stage2_covariance(
        theta_hat = fit_y$par,
        phi_hat = stage1$phi,
        sim = sim,
        split = split,
        cpp_data = cpp_data,
        stage1 = stage1,
        covariance = covariance
      )
    }
  }

  if (xi_mode == "fixxi") {
    fit_fixxi <- cm_lmm_fit_stage2_fixxi(
      cpp_data,
      theta_start,
      xi_fix = xi_fix,
      whiten = whiten,
      whiten_ridge = whiten_ridge
    )
    res$fit_fixxi <- fit_fixxi
    if (isTRUE(compute_covariance)) {
      res$cov_fixxi <- cm_lmm_stage2_covariance(
        theta_hat = fit_fixxi$par,
        phi_hat = stage1$phi,
        sim = sim,
        split = split,
        cpp_data = cpp_data,
        stage1 = stage1,
        fixed_xi = TRUE,
        covariance = covariance
      )
    }
  }

  res
}

cm_lmm_biomarker_cpp_field <- function(time = c("T", "logT")) {
  time <- match.arg(time)
  if (time == "logT") "cm_cpp_logT" else "cm_cpp_T"
}

cm_lmm_biomarker_y_matrix <- function(sim_data) {
  if (is.null(sim_data$y)) {
    stop("`sim_data$y` is required for high-dimensional CM screening.", call. = FALSE)
  }
  y_all <- as.matrix(sim_data$y)
  storage.mode(y_all) <- "double"
  if (nrow(y_all) != nrow(sim_data$X)) {
    stop("`sim_data$y` rows must align with `sim_data$X` rows.", call. = FALSE)
  }
  y_all
}

cm_lmm_normalize_biomarker_ids <- function(y_all, biomarker_ids = NULL) {
  if (is.null(biomarker_ids)) {
    biomarker_ids <- seq_len(ncol(y_all))
  }

  if (is.character(biomarker_ids) && !is.null(colnames(y_all))) {
    idx <- match(biomarker_ids, colnames(y_all))
    if (anyNA(idx)) {
      stop("all `biomarker_ids` must be column names in `sim_data$y`.", call. = FALSE)
    }
    return(idx)
  }

  if (!is.numeric(biomarker_ids) && !is.integer(biomarker_ids)) {
    stop("`biomarker_ids` must be numeric column indices.", call. = FALSE)
  }
  idx <- as.integer(biomarker_ids)
  if (length(idx) != length(biomarker_ids) ||
      anyNA(idx) ||
      any(idx < 1L | idx > ncol(y_all)) ||
      any(idx != biomarker_ids)) {
    stop("`biomarker_ids` must be valid integer column indices in `sim_data$y`.", call. = FALSE)
  }
  idx
}

cm_lmm_prepare_biomarker_cpp <- function(sim_data,
                                         time = c("T", "logT"),
                                         xi = 1) {
  time <- match.arg(time)
  if (!is.numeric(xi) || length(xi) != 1L || !is.finite(xi) || xi <= 0) {
    stop("`xi` must be a single positive finite value.", call. = FALSE)
  }

  timed <- system.time({
    y_all <- cm_lmm_biomarker_y_matrix(sim_data)
    biomarker_ids <- seq_len(ncol(y_all))
    control <- cm_control(
      ftrans_fn = if (time == "logT") log else cm_lmm_ftrans
    )
    truth <- cm_lmm_truth(xi = xi)
    sim <- cm_lmm_as_sim(sim_data, truth = truth, control = control)
    if (is.null(sim$row_index)) {
      stop("C++ CM data conversion did not preserve row indices.", call. = FALSE)
    }
    if (anyNA(sim$row_index) || any(sim$row_index < 1L | sim$row_index > nrow(y_all))) {
      stop("C++ CM row indices are not aligned with `sim_data$y`.", call. = FALSE)
    }

    stage1_start <- if (!is.null(sim$truth$log_klam)) {
      sim$truth$log_klam
    } else {
      log(c(k = 2, lam = 10))
    }
    stage1 <- cm_lmm_fit_stage1_weibull(
      sim$Ti_til,
      sim$del,
      start = stage1_start
    )
    split <- cm_lmm_split_data(sim, log_klam = stage1$phi)
    split$datX1$y_all <- y_all[split$datX1$row_index, biomarker_ids, drop = FALSE]
    split$datX0$y_all <- y_all[split$datX0$row_index, biomarker_ids, drop = FALSE]

    out <- list(
      time = time,
      xi = xi,
      control = control,
      truth = truth,
      sim = sim,
      stage1 = stage1,
      split = split,
      y_all_ordered = y_all[sim$row_index, biomarker_ids, drop = FALSE],
      biomarker_ids = biomarker_ids,
      biomarker_names = colnames(y_all)
    )
  })

  out$setup_time_sec <- unname(timed[["elapsed"]])
  out
}

cm_lmm_ensure_biomarker_cpp <- function(prepared,
                                        time = c("T", "logT"),
                                        biomarker_ids = NULL) {
  time <- match.arg(time)
  field <- cm_lmm_biomarker_cpp_field(time)
  y_all <- cm_lmm_biomarker_y_matrix(prepared$sim_data)
  needed <- cm_lmm_normalize_biomarker_ids(y_all, biomarker_ids)
  existing <- prepared[[field]]

  if (is.null(existing) ||
      is.null(existing$biomarker_ids) ||
      !all(needed %in% existing$biomarker_ids)) {
    prepared[[field]] <- cm_lmm_prepare_biomarker_cpp(
      prepared$sim_data,
      time = time,
      xi = 1
    )
  }

  prepared
}

cm_lmm_assign_biomarker_y <- function(cm_data, biomarker_id) {
  if (is.character(biomarker_id) && !is.null(cm_data$biomarker_names)) {
    pos <- match(biomarker_id, cm_data$biomarker_names)
  } else {
    pos <- match(as.character(biomarker_id), as.character(cm_data$biomarker_ids))
  }
  if (is.na(pos)) {
    stop("`biomarker_id` was not prepared for high-dimensional CM screening.",
      call. = FALSE
    )
  }

  split <- cm_data$split
  y1 <- as.double(split$datX1$y_all[, pos])
  y0 <- as.double(split$datX0$y_all[, pos])
  split$datX1$y <- y1
  split$datX0$y <- y0
  split$datX1$data$y <- y1
  split$datX0$data$y <- y0

  sim <- cm_data$sim
  sim_y <- as.double(cm_data$y_all_ordered[, pos])
  sim$y <- sim_y
  sim$data$y <- sim_y

  list(sim = sim, split = split, pos = pos)
}

cm_lmm_fit_biomarker_cpp <- function(cm_data,
                                     biomarker_id,
                                     whiten = TRUE,
                                     whiten_ridge = 1e-6,
                                     sandwich = FALSE,
                                     return_timing = FALSE) {
  if (!is.logical(whiten) || length(whiten) != 1L || is.na(whiten)) {
    stop("`whiten` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(sandwich) || length(sandwich) != 1L || is.na(sandwich)) {
    stop("`sandwich` must be TRUE or FALSE.", call. = FALSE)
  }

  current <- cm_lmm_assign_biomarker_y(cm_data, biomarker_id)
  cpp_data <- cm_lmm_build_cpp_handles(current$split)
  timed_init <- system.time({
    theta_start <- cm_lmm_stage2_starting_theta(
      current$sim,
      xi_start = cm_data$xi
    )
  })

  fit_values <- function() {
    fit_fixxi <- cm_lmm_fit_stage2_fixxi(
      cpp_data,
      theta_start,
      xi_fix = cm_data$xi,
      whiten = whiten,
      whiten_ridge = whiten_ridge
    )
    cov_fixxi <- cm_lmm_stage2_covariance(
      theta_hat = fit_fixxi$par,
      phi_hat = cm_data$stage1$phi,
      sim = current$sim,
      split = current$split,
      cpp_data = cpp_data,
      stage1 = cm_data$stage1,
      fixed_xi = TRUE,
      covariance = if (isTRUE(sandwich)) "sandwich" else "model"
    )
    V <- if (isTRUE(sandwich)) cov_fixxi$V1 else cov_fixxi$V2
    cm_lmm_assoc_inference(fit_fixxi$par, V)
  }

  if (return_timing) {
    return(list(
      values = fit_values(),
      init_time_sec = unname(timed_init[["elapsed"]])
    ))
  }

  fit_values()
}
