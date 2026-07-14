#' High-dimensional biomarker simulation batch
#'
#' Output is a 13 x 25 matrix: Cox and JM p-values in rows 1-2,
#' conditional-model estimates, standard errors, and p-values in rows 3-9,
#' followed by per-biomarker elapsed times in rows 10-13.

wd <- "."
if (interactive()) wd <- "longbioscreen/37_HDsim_bigN"
library(longscreenfun)

batch <- batch_setup(wd)
io_f <- batch$io_f
args <- batch$args

Nrep0 <- 1e3
j_grp <- 50
j_ <- 1:(100 / j_grp)
facs <- tibble::tibble(
  rep = seq_len(Nrep0),
  seed = get_seed(Nrep0, -1725655327)
) |>
  tidyr::expand_grid(j = j_) |>
  dplyr::mutate(fout = paste0(padLeft(rep, 4), "_", j, ".rda")) |>
  batch_task_grid(args = args)

set.seed(facs$seed)
sim_data <- drawY100(th = 1, N = 2400)
biomarker_ids <- (facs$j - 1) * j_grp + seq_len(j_grp)

fit_screen_biomarkers(sim_data, biomarker_ids, io_f$o(facs$fout))
