library(tidyverse)
library(rstan)

options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)

stan_model <- stan_model(model_code = "
data {
  // Pooled testing data
  int<lower=0> N_pools;              // number of pools
  int<lower=0> x[N_pools];           // test results (0 = negative, 1 = positive)
  int<lower=1> n[N_pools];           // size of each pool
  
  // MGS data
  int<lower=0> N_mgs;                // number of MGS samples
  int<lower=0> viral_reads[N_mgs];   // observed viral reads in each MGS sample
  int<lower=0> total_reads[N_mgs];   // total reads in each MGS sample
  
  // Probability a positive individual makes a pool positive
  real<lower=0, upper=1> sensitivity;
}

parameters {
  real<lower=0, upper=1> p;          // prevalence (fraction of positive individuals)
  real b;                            // log10(p2ra factor)
  real<lower=0> phi;                 // dispersion parameter for negative binomial
}

model {
  // Priors
  p ~ beta(1, 1);                    // Uniform prior on prevalence
  b ~ normal(-6, 4);                 // Derived from eyeballing P2RA report and picking a wide range
  phi ~ exponential(1);              // Prior for dispersion parameter
  
  // Likelihood for pooled testing data
  for (i in 1:N_pools) {
    real prob_neg = pow(1 - p * sensitivity, n[i]); // Probability of negative pool
    if (x[i] == 0)
      target += log(prob_neg);
    else
      target += log(1 - prob_neg);    
  }
  
  // Likelihood for MGS data
  for (j in 1:N_mgs) {
    real mu = total_reads[j] * p * pow(10, b);
    viral_reads[j] ~ neg_binomial_2(mu, phi);
  }
}
")

viruses <- c("SARS-CoV-2", "Coronaviruses (seasonal)", "Mononegavirales", "Rhinoviruses")

pool_positivity <- read_csv("data/swabs.csv") %>%
  select(sample, pool_size, all_of(viruses)) %>%
  pivot_longer(all_of(viruses), names_to = "virus", values_to = "positivity") %>%
  mutate(pool_size = as.integer(pool_size), positivity = as.integer(positivity)) %>%
  nest(.by = virus)

mgs_counts <- read_csv("data/ww_mgs.csv") %>%
  select(sample:n_reads) %>%
  rename(virus = pathogen) %>%
  filter(virus %in% viruses) %>%
  nest(.by = virus)

fits <- left_join(pool_positivity, mgs_counts, by = "virus", suffix = c(".pool", ".mgs")) %>%
  mutate(
    pool_data = map(data.pool, \(df) list(N_pools = nrow(df), x = df$positivity, n = df$pool_size)),
    mgs_data = map(data.mgs, \(df) list(N_mgs = nrow(df), total_reads = df$total_reads, viral_reads = df$n_reads)),
    fit = map2(pool_data, mgs_data, \(x, y)
               sampling(stan_model, c(x, y, sensitivity = 1.0), iter = 2000, warmup = 1000, chains = 4, cores = 4, seed = 381928)
               ),
  )

posteriors <- fits %>%
  mutate(posterior = map(fit, compose(as_tibble, extract))) %>%
  select(virus, posterior) %>%
  unnest(cols = c(posterior)) %>%
  mutate(
    across(-virus, as.double),
    virus = factor(virus, rev(viruses)),
    ) %>%
  glimpse

ggplot(posteriors, mapping=aes(p, virus)) +
  geom_violin(draw_quantiles = c(0.05, 0.5, 0.95)) +
  scale_x_continuous(
    trans = "log10", 
    limits = c(1e-4, 2e-1),
    name = "Prevalence",
    ) 

ggplot(posteriors, mapping=aes(b - 2, virus)) +
  geom_violin(draw_quantiles = c(0.05, 0.5, 0.95)) +
  scale_x_continuous(
    labels = scales::math_format(10^.x),
    name = "RA_p(1%)",
    limits = c(-9, -5),
  )
  