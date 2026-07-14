#' High-dimensional biomarker simulation batch
#'
#' Output is a 13 x 25 matrix: Cox and JM p-values in rows 1-2,
#' conditional-model estimates, standard errors, and p-values in rows 3-9,
#' followed by per-biomarker elapsed times in rows 10-13.
#'
#' This runs

wd <- "."
if (interactive()) wd <- "longbioscreen/38_HDsim_rho"
library(longscreenfun)
batch <- batch_setup(wd)
io_f <- batch$io_f
args <- batch$args

Nrep0 <- 1e3
j_grp <- 100
j_ <- 1:(100 / j_grp)
facs <- tidyr::expand_grid(
  rhoS = 1:2,
  rhoF = 1:2,
  rep = seq_len(Nrep0)
) |>
  dplyr::mutate(seed = get_seed(dplyr::n(), 1557711060)) |>
  tidyr::expand_grid(j = j_) |>
  dplyr::mutate(
    fout = paste0(rhoS, "_", rhoF, "_", padLeft(rep, 4), "_", j, ".rda")) |>
  batch_task_grid(args = args)

set.seed(facs$seed)
rho_b0 <- c(0.3, 0.5)[facs$rhoS]
rho_f <- list(rho_fun_band, rho_fun_blockK)[[facs$rhoF]]
sim_data <- drawY100(th = 1, rho_fun = rho_f, rho = rho_b0)
biomarker_ids <- (facs$j - 1) * j_grp + seq_len(j_grp)

fit_screen_biomarkers(sim_data, biomarker_ids, io_f$o(facs$fout))
