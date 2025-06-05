library(readr)      # read_tsv(), write_tsv()
library(dplyr)      # mutate(), select(), bind_rows(), etc.
library(tidyr)      # crossing(), nest(), unnest(), unnest_wider()
library(purrr)      # map(), map2(), compose()
library(lubridate)  # ymd()
library(stringr)    # str_interp()
library(tibble)     # as_tibble()
library(rstan)
options(mc.cores = parallel::detectCores())

# Run this from the root of the project as
# Rscript statistics/scripts/fit_model.R

# Load data

swab_metadata <- read_tsv("tables/swab-sample-metadata.tsv") %>%
  mutate(date = ymd(paste0("20", date)))

# Use the per-day swab data, which disregards treatment
swab_reads <- read_tsv("tables/swabs-ra-summary.tsv") %>%
  mutate(date = ymd(paste0("20", date))) %>%
  # Redundant with the metadata table
  select(-all_reads)

wastewater_metadata <- read_tsv("tables/wastewater-sample-metadata.tsv") %>%
  mutate(date = ymd(paste0("20", date)))

wastewater_reads <- read_tsv("tables/ww-ra-summary.tsv") %>%
  mutate(date = ymd(paste0("20", date))) %>%
  # Redundant with the metadata table
  select(-all_reads)

# Aggregate to species level and fill in zero counts for species unobserved in
# each sample
prepare_data <- function(all_species, metadata, read_data){
  all_species_and_samples <- crossing(all_species, metadata)
  # The taxonomic assignments go to below species level, but we want to analyze
  # at the species-level assignments
  read_data_species <- summarise(
    read_data,
    .by = c(date, location, species, group),
    viral_reads = sum(dedup_hv)
  )
  left_join(
    all_species_and_samples,
    read_data_species,
    by = join_by(date, location, species, group),
  ) %>%
    arrange(group, species) %>%
    mutate(viral_reads = replace_na(viral_reads, 0))
}

all_species <- bind_rows(swab_reads, wastewater_reads) %>%
  select(species, group) %>%
  distinct %>%
  arrange(group, species)

swab_reads_complete <- prepare_data(
  all_species,
  swab_metadata,
  swab_reads
  )
wastewater_reads_complete <- prepare_data(
  all_species,
  wastewater_metadata,
  wastewater_reads
  )

# Define model and priors

## Logit-normal prior on prevalence
mean_logit_p <- -6
sd_logit_p <- 2

## Logit-normal prior on mu_swab and mu_ww
mean_logit_mu_swab <- -10
sd_logit_mu_swab <- 6
mean_logit_mu_ww <- -14
sd_logit_mu_ww <- 6

## Gamma prior on phi_swab and phi_ww
shape_phi_swab <- 2
rate_phi_swab <- 7
shape_phi_ww <- 1
rate_phi_ww <- 1

model_code = "
data {
  // Pooled swab sequence data
  int<lower=1> N_pools;                   // number of pools
  int<lower=1> n_swabs[N_pools];          // number of swabs in each pool
  int<lower=0> r_total_pool[N_pools];     // total read count
  int<lower=0> r_viral_pool[N_pools];     // viral read count
  
  // Wastewater MGS data
  int<lower=0> N_ww;                      // number of wastewater samples
  int<lower=0> r_total_ww[N_ww];          // total read count
  int<lower=0> r_viral_ww[N_ww];          // viral read count
}
  
parameters {
  real logit_p;                   // prevalence
  
  // Swab parameters (using logit transformation)
  real logit_mu_swab;             // logit of expected viral fraction per infected swab
  real<lower=0> phi_swab;         // inverse-dispersion in a pool containing one infected swab
  
  // Wastewater MGS parameters (using logit transformation)
  real logit_mu_ww;               // logit of expected viral fraction in wastewater per prevalence
  real<lower=0> phi_ww;           // inverse-dispersion in wastewater MGS
}

// Back-transform logistic parametres to [0, 1] scale
transformed parameters {
  real<lower=0, upper=1> p = inv_logit(logit_p);
  real<lower=0, upper=1> mu_swab = inv_logit(logit_mu_swab);
  real<lower=0, upper=1> mu_ww = inv_logit(logit_mu_ww);
}
  
model {
  // Priors
  logit_p ~ normal(${mean_logit_p}, ${sd_logit_p});
  logit_mu_swab ~ normal(${mean_logit_mu_swab}, ${sd_logit_mu_swab});
  phi_swab ~ gamma(${shape_phi_swab}, ${rate_phi_swab});
  logit_mu_ww ~ normal(${mean_logit_mu_ww}, ${sd_logit_mu_ww});
  phi_ww ~ gamma(${shape_phi_ww}, ${rate_phi_ww});
  
  // Marginalized likelihood for pooled testing data
  for (i in 1:N_pools) {
    vector[n_swabs[i] + 1] lpmf_sum;  // to store log probabilities for each possible n_pos value
    for (n_pos in 0:n_swabs[i]) {
      // Probability of n_pos positive samples out of n_swabs[i]
      real binomial_lpmf_val = binomial_lpmf(n_pos | n_swabs[i], p);
      
      // Probability of observing r_viral_pool[i] given n_pos positive samples
      real nb_lpmf_val; 
      if (n_pos == 0){
        // If no positive samples in the pool, we should see zero viral reads
        nb_lpmf_val = r_viral_pool[i] == 0 ? 0 : negative_infinity(); 
      } else {
        nb_lpmf_val = neg_binomial_2_lpmf(r_viral_pool[i] | r_total_pool[i] * n_pos * mu_swab / n_swabs[i], n_pos * phi_swab);
      }
      
      // Store the sum of log probabilities
      lpmf_sum[n_pos + 1] = binomial_lpmf_val + nb_lpmf_val;
    }
    
    // Marginalize out n_pos[i] using log-sum-exp for numerical stability
    target += log_sum_exp(lpmf_sum);
  }
  
  // Wastewater MGS
  for (j in 1:N_ww) {
    real mu = r_total_ww[j] * p * mu_ww;
    r_viral_ww[j] ~ neg_binomial_2(mu, phi_ww);
  }
}
"
stan_model <- stan_model(model_code = str_interp(model_code))

# Fit model

species_data <- full_join(
    swab_reads_complete %>% nest(.by = c(species, group)),
    wastewater_reads_complete %>% nest(.by = c(species, group)),
    by = join_by(species, group),
    suffix = c(".swab", ".wastewater")
)

create_stan_data <- function(df.sw, df.ww) {
  list(
    N_pools = nrow(df.sw),
    n_swabs = df.sw$pool_size,
    r_total_pool = df.sw$all_reads,
    r_viral_pool = df.sw$viral_reads,
    N_ww = nrow(df.ww),
    r_total_ww = df.ww$all_reads,
    r_viral_ww = df.ww$viral_reads
  )
}

fit_model <- function(d) {
  sampling(stan_model, d, iter = 4000, warmup = 1000, chains = 4, cores = 4, seed = 381928)
}

fits <- species_data %>%
  mutate(
    stan_data = map2(data.swab, data.wastewater, create_stan_data),
    fit = map(stan_data, fit_model)
  )

# Output diagnostics and posteriors

get_diagnostics <- function(f) {
  check <- summary(f)$summary
  sampler_params <- get_sampler_params(f, inc_warmup = FALSE)
  
  n_divergent <- sum(do.call(rbind, sampler_params)[, "divergent__"])
  max_treedepth <- max(do.call(rbind, sampler_params)[, "treedepth__"])
  # Stan default max_treedepth = 10
  n_max_treedepth <- sum(do.call(rbind, sampler_params)[, "treedepth__"] == 10)
  # Get the maximum Rhat across all the parameters
  # Efficient sampling is < 1.01
  rhat_max <- max(check[, "Rhat"], na.rm = TRUE)
  
  list(
    n_divergent = n_divergent,
    max_treedepth = max_treedepth,
    n_max_treedepth = n_max_treedepth,
    rhat_max = rhat_max
  )
}

fits %>%
  mutate(diagnostics = map(fit, get_diagnostics)) %>%
  select(species, group, diagnostics) %>%
  unnest_wider(diagnostics) %>%
  write_tsv("tables/model_diagnostics.tsv")

posteriors <- fits %>%
  mutate(posterior = map(fit, compose(as_tibble, extract))) %>%
  select(species, group, posterior) %>%
  unnest(cols = posterior) %>%
  mutate(across(-c(species, group), as.double))

write_tsv(posteriors, "tables/posteriors.tsv")

prev_posteriors <- left_join(
  swab_reads_complete %>%
    summarise(
      swab_samples_with_viral_reads = sum(viral_reads > 0),
      total_swab_viral_reads = sum(viral_reads),
      .by = c(species, group),
    ),
  posteriors %>%
    summarise(
      q15 = quantile(p, .15),
      median = median(p),
      q85 = quantile(p, .85),
      ratio_q85_q15 = q85 / q15,
      .by = c(species, group),
    )
  )

write_tsv(prev_posteriors, "tables/prevalence_posterior_summary.tsv")

ra_posteriors <- posteriors %>%
  mutate(ra01 = mu_ww * 0.01) %>%
  summarise(
    q15 = quantile(ra01, .15),
    median = median(ra01),
    q85 = quantile(ra01, .85),
    ratio_q85_q15 = q85 / q15,
    .by = c(species, group),
  )

write_tsv(ra_posteriors, "tables/ra01_posterior_summary.tsv")
