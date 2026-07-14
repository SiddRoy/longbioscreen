#' longscreenfun: longitudinal biomarker screening simulation helpers
#'
#' Simulation, model-fitting, and screening helpers for longitudinal biomarker
#' studies. The package supports the batch workflows in `31_CP/CPsim.R` and
#' `35_HDsim/HDsim.R`, plus lower-level helpers for trying edits to individual
#' simulation and screening components.
#'
#' Low-dimensional simulators such as [JM_sim_ISCP()] return `list(X, S)`.
#' `X` contains longitudinal rows with subject id `i`, visit time `tim`, and
#' biomarker value `y`; `S` contains one survival row per subject with `i`,
#' observed time `Ti`, and event indicator `del`.
#'
#' The high-dimensional simulator [drawY100()] returns `list(X, S, y)`. `X`
#' and `S` use the same longitudinal/survival contract, while `y` is a matrix
#' whose columns are candidate biomarkers aligned to the rows of `X`.
#'
#' @importFrom dplyr bind_cols bind_rows filter group_by left_join mutate n rename select starts_with ungroup
#' @importFrom JM jointModel
#' @importFrom lme4 fixef getME lmer
#' @importFrom magrittr %>%
#' @importFrom nlme lme lmeControl
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom rlang .data
#' @importFrom stats as.formula dweibull integrate optim pchisq pexp plogis pnorm pweibull qnorm qweibull rexp rnorm runif setNames sigma uniroot
#' @importFrom stringr str_detect
#' @importFrom survival Surv coxph coxph.control
#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr expand_grid
#' @importFrom utils globalVariables
#' @useDynLib longscreenfun, .registration = TRUE
"_PACKAGE"
