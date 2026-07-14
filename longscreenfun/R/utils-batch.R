#' Parse batch command line arguments
#'
#' Parse `--key=value` trailing command-line arguments used by the batch scripts.
#' @family batch helpers
#' @param defaults Default command line string used when no trailing args are supplied.
#' @return A named list of numeric argument values.
#' @export
parse_args <- function(defaults = "--aid=1 --taskid=1") {
  raw <- commandArgs(trailingOnly = TRUE)
  if (length(raw) == 0) {
    raw <- strsplit(trimws(defaults), "\\s+")[[1]]
  }
  result <- list()
  for (a in raw) {
    parts <- strsplit(a, "=")[[1]]
    key <- sub("^--", "", parts[1])
    value <- as.numeric(parts[2])
    result[[key]] <- value
  }
  result
}

#' Draw reproducible integer seeds
#'
#' Draw integer seeds used to fan out reproducible simulation jobs.
#' @family batch helpers
#' @param n Number of values to draw.
#' @param s0 Optional seed used before drawing new seeds.
#' @return Numeric vector of length `n`.
#' @export
get_seed <- function(n = 1, s0 = NULL) {
  if (!is.null(s0)) {
    set.seed(s0)
  }
  round(runif(n, -1, 1) * .Machine$integer.max)
}

#' Left-pad values for file names
#'
#' Left-pad ids for stable batch output file names.
#' @family batch helpers
#' @param x Values to pad.
#' @param l Target string width.
#' @param char Padding character.
#' @return Character vector with padded values.
#' @export
padLeft <- function(x, l, char = "0") {
  n <- pmax(0, l - nchar(x))
  pL <- sapply(n, function(N) {
    paste0(rep(char, N), collapse = "")
  })
  paste0(pL, x)
}

#' Set up a batch simulation run
#'
#' Creates input/output path helpers, ensures the output directory exists, and
#' parses `--aid` and `--taskid` command-line arguments.
#' @family batch helpers
#' @param wd Directory containing the method-specific script.
#' @param defaults Default argument string used when no trailing args are supplied.
#' @return A list with `io_f` path helpers and parsed `args`.
#' @export
batch_setup <- function(wd = ".", defaults = "--aid=1 --taskid=1") {
  io_f <- list(
    i = function(f) file.path(wd, f),
    o = function(f) file.path(wd, "out", f)
  )
  if (!dir.exists(io_f$o(""))) dir.create(io_f$o(""), recursive = TRUE)
  list(io_f = io_f, args = parse_args(defaults))
}

#' Add scheduler ids to a simulation grid
#'
#' Adds `seed`, `a_id`, and `task_id` columns in the same two-task layout used
#' by the existing batch scripts, then filters to the requested task.
#' @family batch helpers
#' @param grid Data frame with one row per simulation file to produce.
#' @param args Parsed batch arguments from [parse_args()].
#' @param seed0 Optional seed used to draw one seed per grid row.
#' @return The single grid row assigned to the requested array/task id.
#' @export
batch_task_grid <- function(grid, args, seed0 = NULL) {
  n_grid <- nrow(grid)
  if (n_grid %% 2 != 0) stop("batch_task_grid() requires an even number of rows")
  if (!is.null(seed0)) grid$seed <- get_seed(n_grid, seed0)
  grid$a_id <- rep(seq_len(n_grid / 2), each = 2)
  grid$task_id <- rep(1:2, n_grid / 2)
  dplyr::filter(grid, .data$a_id == args$aid, .data$task_id == args$taskid)
}

#' Build a stable replicate-level seed grid
#'
#' Creates one row per simulated replicate. Seeds are generated at the
#' replicate level, then assigned to output files by `n_per`. Keep
#' `reference_n_per` fixed to preserve replicate seeds when changing `n_per`.
#' @family batch helpers
#' @param param_grid Data frame with one row per simulation scenario.
#' @param n_reps Number of simulated replicates per scenario.
#' @param n_per Number of simulated replicates to bundle into each output file.
#' @param seed0 Master seed used to draw stable replicate seeds.
#' @param reference_n_per Historical or reference bundle size used to derive
#'   replicate seeds. Leave this fixed when changing `n_per`.
#' @return A tibble with scenario columns plus `sim_rep`, `seed`, `file_rep`,
#'   and `within_file`.
#' @export
replicate_seed_grid <- function(param_grid, n_reps, n_per, seed0,
                                reference_n_per = n_per) {
  if (!is.data.frame(param_grid)) {
    stop("param_grid must be a data frame")
  }
  if (nrow(param_grid) < 1) {
    stop("param_grid must contain at least one scenario row")
  }
  if (length(n_reps) != 1 || n_reps < 1 || n_reps != floor(n_reps)) {
    stop("n_reps must be a positive whole number")
  }
  if (length(n_per) != 1 || n_per < 1 || n_per != floor(n_per)) {
    stop("n_per must be a positive whole number")
  }
  if (length(reference_n_per) != 1 ||
      reference_n_per < 1 ||
      reference_n_per != floor(reference_n_per)) {
    stop("reference_n_per must be a positive whole number")
  }
  if (n_reps %% n_per != 0) {
    stop("n_reps must be divisible by n_per")
  }
  if (n_reps %% reference_n_per != 0) {
    stop("n_reps must be divisible by reference_n_per")
  }

  n_reps <- as.integer(n_reps)
  n_per <- as.integer(n_per)
  reference_n_per <- as.integer(reference_n_per)
  param_grid <- tibble::as_tibble(param_grid)
  n_scenarios <- nrow(param_grid)
  ref_files_per_scenario <- n_reps / reference_n_per
  file_seeds <- get_seed(n_scenarios * ref_files_per_scenario, seed0)
  replicate_seeds <- unlist(
    lapply(file_seeds, function(s) get_seed(reference_n_per, s)),
    use.names = FALSE
  )

  out <- param_grid[rep(seq_len(n_scenarios), each = n_reps), , drop = FALSE]
  out$sim_rep <- rep(seq_len(n_reps), times = n_scenarios)
  out$seed <- replicate_seeds
  out$file_rep <- as.integer(((out$sim_rep - 1L) %/% n_per) + 1L)
  out$within_file <- as.integer(((out$sim_rep - 1L) %% n_per) + 1L)
  out
}

#' Simulate one batch file from replicated seeds
#'
#' Draws `n_per` child seeds from a file-level seed and calls a simulator once
#' per child seed.
#' @family batch helpers
#' @param n_per Number of simulated datasets in the output file.
#' @param seed File-level seed.
#' @param sim_fun Simulation function accepting a `seed` argument.
#' @param ... Additional arguments passed to `sim_fun`.
#' @return A list of simulated datasets.
#' @export
simulate_seed_reps <- function(n_per, seed, sim_fun, ...) {
  simulate_seed_vec(get_seed(n_per, seed), sim_fun, ...)
}

#' Simulate replicated datasets from explicit seeds
#'
#' Calls a simulator once per supplied seed, passing each value through the
#' simulator's `seed` argument.
#' @family batch helpers
#' @param seeds Numeric vector of replicate-level seeds.
#' @param sim_fun Simulation function accepting a `seed` argument.
#' @param ... Additional arguments passed to `sim_fun`.
#' @return A list of simulated datasets.
#' @export
simulate_seed_vec <- function(seeds, sim_fun, ...) {
  lapply(seeds, function(s) sim_fun(seed = s, ...))
}
