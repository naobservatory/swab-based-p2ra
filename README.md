# Estimating the sensitivity of wastewater metagenomic sequencing using nasal swabs

Simon Grimm, Dan Rice, Mike McLaren

Code and documentation supporting the NAO blog post [Estimating the sensitivity of wastewater metagenomic sequencing using nasal swabs](https://naobservatory.org/blog/swab-based-p2ra/)
and its [Appendix](https://naobservatory.github.io/swab-based-p2ra/).

## Reproducing the statistical analysis

To reproduce the model fits and diagnostics:

1. Clone this git repository.
2. From the R console, install the necessary R packages with: `install.packages(c("tidyverse", "here", "rstan"))`.
   RStan may require additional setup; see the [RStan docs](https://mc-stan.org/rstan/).
3. Run the model-fitting script with: `Rscript statistics/scripts/fit_model.R`.
   (This script should take 5--10 minutes to run.)
   
The model-fitting script outputs four files (all of which are already checked into the repository):

- `tables/model_diagnostics.tsv` contains Stan diagnostic information,
  including the number of divergent transitions and the maximum rhat over all parameters.
- `tables/posteriors.tsv` contains posterior samples for all viruses.
  Note that re-running the fitting script may shuffle the order of the posterior samples in this file.
- `tables/{prevalence,ra01}_posterior_summary.tsv` contain posterior 15th,
  50th, and 85th percentiles for prevalence and $RA_p(1\%)$ for all viruses.

## Publishing the appendix (NAO-internal only)

To publish an update to the appendix, make any changes, then run:
```
quarto publish gh-pages statistics/appendix/index.qmd
```
from the root of the repository.
