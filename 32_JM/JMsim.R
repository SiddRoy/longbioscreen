#' Joint-model association simulation batch
#'
#' Output is a 17 x `n_per` matrix: Cox estimates, standard errors, and
#' p-values in rows 1-3; JM estimates, standard errors, and p-values in rows
#' 4-6; seven conditional-model statistics in rows 7-13; and model timing
#' rows in rows 14-17.

wd <- "."
if (interactive()) wd <- "longbioscreen/32_JM"
library(longscreenfun)

batch <- batch_setup(wd)
io_f <- batch$io_f
args <- batch$args

Nrep0 <- 2e3
n_per <- 400
seed_ref_n_per <- 5

seed_grid <- replicate_seed_grid(
  tidyr::expand_grid(aC = 0:9),
  n_reps = Nrep0,
  n_per = n_per,
  seed0 = -1527155019,
  reference_n_per = seed_ref_n_per
) |>
  dplyr::group_by(aC, file_rep) |>
  dplyr::summarise(
    seeds = list(seed),
    sim_reps = list(sim_rep),
    fout = paste0(aC[1], "_", padLeft(file_rep[1], 3), ".rda"),
    .groups = "drop"
  )

facs <- seed_grid |>
  batch_task_grid(args = args)

sim_data <- simulate_seed_vec(
  facs$seeds[[1]], JM_sim_ISCP,
  N = 600,
  alp = facs$aC / 10 / 3
)

fit_screen_simulations(sim_data, io_f$o(facs$fout))
