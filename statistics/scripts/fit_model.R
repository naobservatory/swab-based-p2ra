library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(here)
library(rstan)
options(mc.cores = parallel::detectCores())

# Fit models to viruses and output diagnostics, posteriors, and summaries
# Run script from anywhere within the project directory

source(here("statistics/scripts/prepare_data.R"))
source(here("statistics/scripts/define_model.R"))

# Compile model
stan_model <- stan_model(model_code = model_code)

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
  write_tsv(here(tables_dir, "model_diagnostics.tsv"))

posteriors <- fits %>%
  mutate(posterior = map(fit, compose(as_tibble, extract))) %>%
  select(species, group, posterior) %>%
  unnest(cols = posterior) %>%
  mutate(across(-c(species, group), as.double))

write_tsv(posteriors, here(tables_dir, "posteriors.tsv"))

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

write_tsv(prev_posteriors, here(tables_dir, "prevalence_posterior_summary.tsv"))

ra_posteriors <- posteriors %>%
  mutate(ra01 = mu_ww * 0.01) %>%
  summarise(
    q15 = quantile(ra01, .15),
    median = median(ra01),
    q85 = quantile(ra01, .85),
    ratio_q85_q15 = q85 / q15,
    .by = c(species, group),
  )

write_tsv(ra_posteriors, here(tables_dir, "ra01_posterior_summary.tsv"))
