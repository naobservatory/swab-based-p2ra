library(tidyverse)
library(rstan)
library(knitr)
theme_set(theme_minimal())

# Colors borrowed from https://github.com/JLSteenwyk/ggpubfigs
wong_eight <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
options(ggplot2.discrete.colour = function() scale_colour_manual(values = wong_eight))
options(ggplot2.discrete.fill = function() scale_fill_manual(values = wong_eight))

options(mc.cores = parallel::detectCores())

# Model

## Hyperparameters

### prevalence
alpha_p <- 0.5
beta_p <- 10
prior_p <- function(x) dbeta(x, alpha_p, beta_p)

### swab: logit(mu_1)
mean_logit_mu_1 <- -7  # Prior mean on logit scale
sd_logit_mu_1 <- 3     # Prior standard deviation on logit scale
prior_logit_mu_1 <- function(x) dnorm(x, mean_logit_mu_1, sd_logit_mu_1)

### swab: phi_1
shape_phi_1 <- 2
rate_phi_1 <- 2
prior_phi_1 <- function(x) dgamma(x, shape_phi_1, rate_phi_1)

### wastewater: logit(mu_ww)
mean_logit_mu_ww <- -11  # Prior mean on logit scale
sd_logit_mu_ww <- 3      # Prior standard deviation on logit scale
prior_logit_mu_ww <- function(x) dnorm(x, mean_logit_mu_ww, sd_logit_mu_ww)

### wastewater: phi_ww
shape_phi_ww <- 2
rate_phi_ww <- 2
prior_phi_ww <- function(x) dgamma(x, shape_phi_ww, rate_phi_ww)

## Plot priors (in transformed space)

crossing(
  datatype = factor(c("swab", "wastewater")),
  logit_mu = seq(-18, 2, length.out=1000),
) %>%
  mutate(
    mu = plogis(logit_mu),
    prior = if_else(datatype == "swab", 
                    prior_logit_mu_1(logit_mu), 
                    prior_logit_mu_ww(logit_mu)),
  ) %>%
  ggplot(aes(logit_mu, prior, color=datatype)) +
  geom_line() +
  labs(title = "Prior distributions on logit(mu)")

## Plot priors (in original space)
crossing(
  datatype = factor(c("swab", "wastewater")),
  logit_mu = seq(-20, 4, length.out=1000),
) %>%
  mutate(
    mu = plogis(logit_mu),
    prior = if_else(datatype == "swab", 
                          prior_logit_mu_1(logit_mu), 
                          prior_logit_mu_ww(logit_mu)),
    prior_log10 = prior * mu * log(10),
  ) %>%
  ggplot(aes(mu, prior, color=datatype)) +
  geom_line() +
  scale_x_log10(breaks = 10^seq(-8, 0, 1)) +
  labs(title = "Implied prior distributions on mu")

## Model code

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
  real<lower=0, upper=1> p;       // prevalence
  
  // Swab parameters (using logit transformation)
  real logit_mu_1;                // logit of expected viral fraction in a pool containing a single positive swab
  real<lower=0> phi_1;            // inverse-dispersion in a pool containing a single positive swab
  
  // Wastewater MGS parameters (using logit transformation)
  real logit_mu_ww;               // logit of expected viral fraction in wastewater if everyone were infected
  real<lower=0> phi_ww;           // inverse-dispersion in wastewater MGS
}

transformed parameters {
  real<lower=0, upper=1> mu_1 = inv_logit(logit_mu_1);    // back-transform to [0,1] scale
  real<lower=0, upper=1> mu_ww = inv_logit(logit_mu_ww);  // back-transform to [0,1] scale
}
  
model {
  // Priors
  p ~ beta(${alpha_p}, ${beta_p});
  logit_mu_1 ~ normal(${mean_logit_mu_1}, ${sd_logit_mu_1});
  phi_1 ~ gamma(${shape_phi_1}, ${rate_phi_1});
  logit_mu_ww ~ normal(${mean_logit_mu_ww}, ${sd_logit_mu_ww});
  phi_ww ~ gamma(${shape_phi_ww}, ${rate_phi_ww});  // Updated from placeholder
  
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
        nb_lpmf_val = neg_binomial_2_lpmf(r_viral_pool[i] | r_total_pool[i] * n_pos * mu_1 / n_swabs[i], n_pos * phi_1);
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

# Swab data

swab_metadata <- read_tsv("../tables/swab-sample-metadata.tsv") %>%
  mutate(date = ymd(paste0("20", date))) 

# Use the per-day data. Caveat about the different treatments
swab_reads <- read_tsv("../tables/swabs-ra-summary.tsv") %>%
  mutate(date = ymd(paste0("20", date))) %>%
  select(-all_reads)

# Wastewater data

wastewater_metadata <- read_tsv("../tables/wastewater-sample-metadata.tsv") %>%
  mutate(date = ymd(paste0("20", date))) 

wastewater_reads <- read_tsv("../tables/ww-ra-summary.tsv") %>%
  mutate(date = ymd(paste0("20", date))) %>%
  select(-all_reads)

# Fit model to each virus at the species level

swab_species <- swab_reads %>%
  select(species, group) %>%
  distinct

wastewater_species <- wastewater_reads %>%
  select(species, group) %>%
  distinct

all_species <- full_join(wastewater_species, swab_species)

species_order <- all_species %>%
  arrange(group, species) %>%
  pull(species)

# Fill in zero counts
swab_reads_complete <- crossing(all_species, swab_metadata) %>%
  left_join(
    summarise(swab_reads, .by = c(date, location, species, group), viral_reads = sum(dedup_hv)),
    by = join_by(date, location, species, group),
  ) %>%
  arrange(group, species) %>%
  mutate(
    viral_reads= replace_na(viral_reads, 0),
    species = fct_inorder(species),
    group = fct_inorder(group)
    )

wastewater_reads_complete <- wastewater_metadata %>%
  crossing(all_species) %>%
  left_join(
    summarise(wastewater_reads, .by = c(date, location, species, group), viral_reads = sum(dedup_hv)),
    by = join_by(date, location, species, group),
  ) %>%
  arrange(group, species) %>%
  mutate(
    viral_reads= replace_na(viral_reads, 0),
    species = fct_inorder(species),
    group = fct_inorder(group)
    )

# Assemble full dataset and fit model

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

# Diagnostics

get_diagnostics <- function(f) {
  check <- summary(f)$summary
  sampler_params <- get_sampler_params(f, inc_warmup = FALSE)
  
  n_divergent <- sum(do.call(rbind, sampler_params)[, "divergent__"])
  max_treedepth <- max(do.call(rbind, sampler_params)[, "treedepth__"])
  n_max_treedepth <- sum(do.call(rbind, sampler_params)[, "treedepth__"] == 10)  # Assuming max_treedepth of 10
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
  print

posteriors <- fits %>%
  mutate(posterior = map(fit, compose(as_tibble, extract))) %>%
  select(species, group, posterior) %>%
  unnest(cols = posterior) %>%
  mutate(
    across(-c(species, group), as.double),
    )

# Compare posterior median prevalence with naive estimate from positive pools

left_join(
  swab_reads_complete %>%
    summarise(
      pos_pools_div_by_total_swabs = sum(viral_reads > 0) / sum(pool_size),
      .by = c(species, group),
      ),
  posteriors %>%
    summarise(
      posterior_median = median(p),
      .by = species,
    )
  ) %>%
  mutate(ratio = posterior_median / pos_pools_div_by_total_swabs) %>%
  print

#  Plot posteriors

prior_density <- function(x, density_function, log10x = FALSE){
  if (log10x){
    data <- tibble(
      x = x,
      density = density_function(x) * x * log(10),
    )
  } else {
    data <- tibble(
      x = x,
      density = density_function(x),
    )
  }
  geom_line(
    aes(x, density),
    data = data,
    color="grey",
    linetype="dashed",
    inherit.aes = FALSE
    )
}

posteriors %>%
  ggplot(aes(p, color=group)) +
  geom_density() +
  prior_density(10^seq(-8, 0, length.out=100), prior_p, log10x = TRUE) +
  scale_x_log10(limits = c(1e-8, 1), name = "prevalence") +
  scale_y_continuous(name = "") +
  facet_wrap(vars(species))

posteriors %>%
  ggplot(aes(logit_mu_1, color=group)) +
  geom_density() +
  prior_density(seq(-20, 4, length.out=100), prior_logit_mu_1) +
  scale_y_continuous(name = "") +
  facet_wrap(vars(species))

posteriors %>%
  ggplot(aes(phi_1, color=group)) +
  geom_density() +
  prior_density(10^seq(-3, 3, length.out=100), prior_phi_1, log10x = TRUE) +
  scale_x_log10() +
  scale_y_continuous(name = "") +
  facet_wrap(vars(species))

posteriors %>%
  ggplot(aes(logit_mu_ww, color=group)) +
  geom_density() +
  prior_density(seq(-20, 0, length.out=100), prior_logit_mu_ww) +
  scale_y_continuous(name = "") +
  facet_wrap(vars(species))

posteriors %>%
  ggplot(aes(phi_ww, color=group)) +
  geom_density() +
  prior_density(10^seq(-3, 3, length.out=100), prior_phi_ww, log10x = TRUE) +
  scale_x_log10(limits = c(1e-3, 1e3)) +
  facet_wrap(vars(species))

# Joint posteriors

## Effects on prevalence estimates

posteriors %>%
  ggplot(aes(mu_1, p)) +
  geom_point(alpha = 0.05, size = 0.1) +
  geom_smooth() +
  scale_x_log10(limits = c(1e-7, 1)) +
  scale_y_log10(limits = c(1e-6, 1)) +
  facet_wrap(~species)

posteriors %>%
  ggplot(aes(phi_1, p)) +
  geom_point(alpha = 0.05, size = 0.1) +
  geom_smooth() +
  scale_x_log10() +
  scale_y_log10(limits = c(1e-6, 1)) +
  facet_wrap(~species)

## Effects on P2RA estimates

posteriors %>%
  ggplot(aes(p, mu_ww)) +
  geom_point(alpha = 0.05, size = 0.1) +
  geom_smooth() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~species)

posteriors %>%
  ggplot(aes(mu_1, mu_ww)) +
  geom_point(alpha = 0.05, size = 0.1) +
  geom_smooth() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~species)

posteriors %>%
  ggplot(aes(phi_1, mu_ww)) +
  geom_point(alpha = 0.05, size = 0.1) +
  geom_smooth() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~species)

posteriors %>%
  ggplot(aes(phi_ww, mu_ww)) +
  geom_point(alpha = 0.05, size = 0.1) +
  geom_smooth() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~species)

# Posterior predictive distributions

## Swabs 

posterior_predictive_swab <- posteriors %>%
  mutate(rep = 1:n(), .by=species) %>%
  cross_join(swab_metadata) %>%
  mutate(
    n_pos = map2_int(pool_size, p, ~ rbinom(1, .x, .y)),
    mu = n_pos * mu_1 * all_reads / pool_size,
    phi = n_pos * phi_1,
    viral_reads = map2_int(mu, phi, ~ifelse(.x > 0, rnbinom(1, size = .y, mu = .x), 0)),
    false_negative = n_pos > 0 & viral_reads == 0,
  )

posterior_predictive_swab_summary <- posterior_predictive_swab %>%
  summarise(
    pos_pools = sum(viral_reads > 0),
    viral_reads_total = sum(viral_reads),
    viral_reads_max = max(viral_reads),
    viral_reads_sd = sd(viral_reads),
    false_negative_pools = sum(false_negative),
    across(c(p, mu_1, phi_1), first),
    .by = c(species, group, rep),
  )

species_swab_summary <- swab_reads_complete %>%
  summarise(
    pos_pools = sum(viral_reads > 0),
    viral_reads_total = sum(viral_reads),
    viral_reads_max = max(viral_reads),
    viral_reads_sd = sd(viral_reads),
    .by = c(species, group),
  ) %>%
  arrange(species) %>%
  print

posterior_predictive_swab_summary %>%
  ggplot(aes(x = pos_pools, fill = group, after_stat(density))) +
  geom_histogram(binwidth = 1, center = 0, color = "grey") +
  geom_vline(data = species_swab_summary, mapping = aes(xintercept=pos_pools)) +
  facet_wrap(vars(species)) +
  scale_x_continuous(breaks = seq(0, 12, 2), name = "positive pools")

posterior_predictive_swab_summary %>%
  ggplot(aes(x = viral_reads_total, y = species, fill = group)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(data = species_swab_summary) +
  scale_x_log10(
    breaks = 10^seq(0, 4),
    name = "total viral read count (conditional on > 0)",
    oob = scales::oob_keep,
  ) +
  scale_y_discrete(limits = rev(species_order))

posterior_predictive_swab_summary %>%
  ggplot(aes(x = viral_reads_max, y = species, fill = group)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(data = species_swab_summary) +
  scale_x_log10(
    breaks = 10^seq(0, 4),
    name = "largest viral read count (conditional on > 0)",
    oob = scales::oob_keep,
  ) +
  scale_y_discrete(limits = rev(species_order))

posterior_predictive_swab_summary %>%
  ggplot(aes(x = viral_reads_sd, y = species, fill = group)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(data = species_swab_summary) +
  scale_x_log10(
    breaks = 10^seq(0, 4),
    name = "stdev viral read count (conditional on > 0)",
    oob = scales::oob_keep,
  ) +
  scale_y_discrete(limits = rev(species_order))

posterior_predictive_swab_summary %>%
  ggplot(aes(x = false_negative_pools, fill = group, after_stat(density))) +
  geom_histogram(binwidth = 1, center = 0, color = "grey") +
  facet_wrap(vars(species)) +
  scale_x_continuous(breaks = seq(0, 10), limits = c(-1, 8)) +
  theme_minimal()

posterior_predictive_swab_summary %>%
  ggplot(aes(phi_1, pos_pools, color = group)) +
  geom_point(alpha = 0.05, size = 0.1) +
  geom_smooth() +
  geom_hline(data = species_swab_summary, mapping = aes(yintercept=pos_pools)) +
  scale_x_log10() +
  facet_wrap(vars(species))

## Wastewater

posterior_predictive_ww <- posteriors %>%
  mutate(rep = 1:n(), .by=species) %>%
  cross_join(wastewater_metadata) %>%
  mutate(
    mu = p * mu_ww * all_reads,
    viral_reads = map2_int(mu, phi_ww, ~ifelse(.x > 0, rnbinom(1, size = .y, mu = .x), 0)),
  )

posterior_predictive_ww_summary <- posterior_predictive_ww %>%
  summarise(
    pos_samples = sum(viral_reads > 0),
    viral_reads_total = sum(viral_reads),
    viral_reads_max = max(viral_reads),
    viral_reads_sd = sd(viral_reads),
    across(c(p, mu_ww, phi_ww), first),
    viral_reads_expected = p * mu_ww * sum(all_reads),
    .by = c(species, group, rep),
  )

species_ww_summary <- wastewater_reads_complete %>%
  summarise(
    pos_samples = sum(viral_reads > 0),
    viral_reads_total = sum(viral_reads),
    viral_reads_max = max(viral_reads),
    viral_reads_sd = sd(viral_reads),
    .by = c(species, group),
  ) %>%
  arrange(species) %>%
  print

posterior_predictive_ww_summary %>%
  ggplot(aes(x = pos_samples, fill = group, after_stat(density))) +
  geom_histogram(binwidth = 1, center = 0, color = "grey") +
  geom_vline(data = species_ww_summary, mapping = aes(xintercept = pos_samples)) +
  facet_wrap(vars(species)) +
  scale_x_continuous(breaks = seq(0, 12, 2), name = "positive samples")

posterior_predictive_ww_summary %>%
  ggplot(aes(x = viral_reads_total, y = species, fill = group)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(data = species_ww_summary) +
  scale_x_log10(
    breaks = 10^seq(0, 4),
    name = "total viral read count (conditional on > 0)",
    oob = scales::oob_keep,
  ) +
  scale_y_discrete(limits = rev(species_order))

posterior_predictive_ww_summary %>%
  ggplot(aes(x = viral_reads_expected, y = species, fill = group)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(data = species_ww_summary, mapping = aes(x = viral_reads_total)) +
  scale_x_log10(
    breaks = 10^seq(0, 4),
    name = "total viral reads expected (observed)",
    oob = scales::oob_keep,
  ) +
  scale_y_discrete(limits = rev(species_order))

posterior_predictive_ww_summary %>%
  ggplot(aes(x = viral_reads_max, y = species, fill = group)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(data = species_ww_summary) +
  scale_x_log10(
    breaks = 10^seq(0, 4),
    name = "largest viral read count (conditional on > 0)",
    oob = scales::oob_keep,
  ) +
  scale_y_discrete(limits = rev(species_order))

posterior_predictive_ww_summary %>%
  ggplot(aes(x = viral_reads_sd, y = species, fill = group)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(data = species_ww_summary) +
  scale_x_log10(
    breaks = 10^seq(0, 4),
    name = "stdev viral read count (conditional on > 0)",
    oob = scales::oob_keep,
  ) +
  scale_y_discrete(limits = rev(species_order))

# RA_i(1%) and related quantities

posteriors %>%
  mutate(ra01 = mu_ww * 0.01) %>%
  ggplot(aes(ra01, species, fill=group)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_x_log10(name = "wastewater RA_p(1%)") +
  scale_y_discrete(limits=rev(species_order))

posteriors %>%
  mutate(ratio = mu_1 / mu_ww) %>%
  ggplot(aes(ratio, species, fill=group)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_x_log10(name = "swab RA / wastewater RA", limits=c(1e-2,1e8), breaks=10^seq(0,6,2)) +
  scale_y_discrete(limits=rev(species_order))

posteriors %>%
  mutate(ra01 = mu_ww * 0.01) %>%
  summarise(
    q25 = quantile(ra01, 0.25),
    median = median(ra01),
    q75 = quantile(ra01, 0.75),
    .by = c(species, group)
  ) %>%
  mutate(across(-c(species, group), ~format(.x, scientific=TRUE, digits = 2))) %>%
  kable

posteriors %>%
  mutate(ratio = mu_1 / mu_ww) %>%
  summarise(
    q25 = quantile(ratio, 0.25),
    median = median(ratio),
    q75 = quantile(ratio, 0.75),
    .by = c(species, group)
  ) %>%
  mutate(across(-c(species, group), ~format(.x, scientific=TRUE, digits = 2))) %>%
  kable
