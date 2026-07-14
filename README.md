# Comparison of Testing Procedures for Detecting Short and Long Term Effects of Biomarkers for Cancer Screening

The local R package `longscreenfun` contains all functions 
used in simulations. 

To replicate the results, please use the following steps: 

1. Clone the GitHub repository via
   ```bash
   git clone https://github.com/SiddRoy/longbioscreen.git
   cd longbioscreen
   ```
2. Install the local R package in `longscreenfun/` 
   ```bash
   Rscript -e 'devtools::install("longscreenfun", dependencies = TRUE)'
   ```
   - All exported functions are documented. 
3. Run the worked examples below
4. Then run R scripts inside the subfolder to generate the 
   results. As specified in supplement 2 of the simulation,
   fitting all three models takes less than 10 seconds per biomarker
   at the $N=600$ sample size. 
   - Each subsection from section 3 has a corresponding 
     subfolder to generate results. Folders 36-39 correspond 
     to supplement sections and extra results referred to in
     the paper. 
   - The full set of repetitions were split and implemented 
     using SLURM. Variables Nrep0, n_per, determine the
     total number of repetitions and the number of 
     repetitions done per split. 

## Running the simulations

Each simulation script determines its own directory, so it can be run from
the repository root regardless of the repository folder name. For each `aid`,
run both `taskid=1` and `taskid=2`. The following replicates the conditional
model simulations from Section 3.4.
 
```bash
for i in {1..250}
do
  Rscript 34_CM/CMsim.R --aid=$i --taskid=1 &
  Rscript 34_CM/CMsim.R --aid=$i --taskid=2 &
  wait
done   
```

The maximum array IDs and conservative runtime estimates per
split used for each simulation are shown below. 

| Simulation folder | Script | Maximum `aid` | Total simulation blocks | Runtime per block |
|:------------------|:-------|--------------:|------------------------:|:------------------|
| `31_CP`            | `CPsim.R`  | 10    | 20    | < 70 min |
| `32_JM`            | `JMsim.R`  | 25    | 50    | < 70 min |
| `33_Cox`           | `Coxsim.R` | 25    | 50    | < 70 min |
| `34_CM`            | `CMsim.R`  | 250   | 500   | < 70 min |
| `35_HDsim`         | `HDsim.R`  | 500   | 1,000 | < 20 min |
| `36_CM_xi2`        | `CMsim.R`  | 250   | 500   | < 70 min |
| `37_HDsim_bigN`    | `HDsim.R`  | 1,000 | 2,000 | < 20 min ($N=2400$) |
| `38_HDsim_rho`     | `HDsim.R`  | 2,000 | 4,000 | < 20 min |
| `39_HDsim_CS`      | `HDsim.R`  | 1,500 | 3,000 | < 20 min |
   
Outputs go to each folder's out/ directory. 

For the low-dimensional simulations, the 17 rows are:

- rows 1-3: Cox estimate, standard error, and p-value
- rows 4-6: joint-model estimate, standard error, and p-value
- rows 7-13: conditional-model statistics
- row 14: Cox elapsed time in seconds
- row 15: joint-model elapsed time in seconds
- row 16: conditional-model elapsed time in seconds
- row 17: conditional-model initialization elapsed time in seconds

For the high-dimensional simulations, the 13 rows are:

- row 1: Cox p-values
- row 2: joint-model p-values
- rows 3-9: conditional-model estimates, standard errors, and p-values/statistics
- row 10: Cox elapsed time in seconds
- row 11: joint-model elapsed time in seconds
- row 12: conditional-model elapsed time in seconds
- row 13: conditional-model initialization elapsed time in seconds

## Worked examples

Generate and fit one low-dimensional simulated dataset:

```r
library(longscreenfun)

set.seed(1)
dat <- JM_sim_ISCP(N = 600)

cox_screen_pvalue(dat)
jm_screen_pvalue(dat)
jm_screen_pvalue(dat, jm_random_specs = c(b0 = "~1 | i", b01 = "~ tim | i"))
cm_screen_stats(dat)
```

Generate and analyze one high-dimensional dataset:

```r
library(longscreenfun)

set.seed(1)
dat <- drawY100(th = 1)
prepared <- prepare_biomarker_screen(dat)

cox_screen_biomarkers(prepared, 1:5)
jm_screen_biomarkers(prepared, 1:5)
jm_screen_biomarkers(
  prepared, 1:5,
  jm_random_specs = c(b0 = "~1 | i", b01 = "~ tim | i"),
  jm_stats = c("est", "se", "p")
)
cm_screen_biomarkers(prepared, 1:5)
```
