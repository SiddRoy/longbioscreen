#' Conditional-model association simulation batch
#'
#' Output is a 16 x `n_per` matrix: Cox estimates, standard errors, and
#' p-values in rows 1-3; JM estimates, standard errors, and p-values in rows
#' 4-6; seven conditional-model statistics in rows 7-13; and model timing
#' rows in rows 14-16.
#'
#' This simulation generates form xi = 2, to provide results
#' for mismatched xi = 1. Later, I will expand this.

wd <- "."
if (interactive()) wd <- "longbioscreen/36_CM_xi2"
library(longscreenfun)

batch <- batch_setup(wd)
io_f <- batch$io_f
args <- batch$args

Nrep0 <- 2e3
n_per <- 400
seed_ref_n_per <- 50

seed_grid <- replicate_seed_grid(
  tidyr::expand_grid(BG = 0:9, BS = 0:9),
  n_reps = Nrep0,
  n_per = n_per,
  seed0 = 6105523,
  reference_n_per = seed_ref_n_per
) |>
  dplyr::group_by(BG, BS, file_rep) |>
  dplyr::summarise(
    seeds = list(seed),
    sim_reps = list(sim_rep),
    fout = paste0(BG[1], "_", BS[1], "_", padLeft(file_rep[1], 3), ".rda"),
    .groups = "drop"
  )

facs <- seed_grid |>
  batch_task_grid(args = args)

sim_data <- simulate_seed_vec(
  facs$seeds[[1]], CM_sim_RI,
  N = 600,
  Btau = 1.2 * facs$BG / 10,
  alp = -0.06 * facs$BS / 10,
  xi = 2
)

fit_screen_simulations(sim_data, io_f$o(facs$fout))
