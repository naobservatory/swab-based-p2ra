library(tidyverse)
library(rstan)

options(mc.cores = parallel::detectCores())

stan_model <- stan_model(model_code = "
data {
  // Pooled testing data
  int<lower=1> N_pools;              // number of pools
  int<lower=1> n_swabs[N_pools];     // number of swabs in each pool
  int<lower=0> r_total[N_pools];     // total read count
  int<lower=0> r_viral[N_pools];     // viral read count
  }

parameters {
  real<lower=0, upper=1> p;          // prevalence
  real<lower=0, upper=1> mu_1;       // expected viral fraction in a pool containing only a single positive swab
  real<lower=0> phi_1;               // inverse-dispersion in a pool containing a single positive swab
  }

model {
  // Priors
  p ~ beta(0.5, 10.5);                 // Jeffreys left tail, informative right tail (10 pseudo-observations)
  mu_1 ~ beta(0.5, 100.5);             // Jeffreys left tail, informative right tail (100 pseudo-observations)
  phi_1 ~ exponential(1);            // PLACEHOLDER
  
  // Marginalized likelihood for pooled testing data
  for (i in 1:N_pools) {
    vector[n_swabs[i] + 1] lpmf_sum;  // to store log probabilities for each possible n_pos value
    
    for (j in 0:n_swabs[i]) {
      // Probability of j positive samples out of n_swabs[i]
      real binomial_lpmf_val = binomial_lpmf(j | n_swabs[i], p);
      
      // Probability of observing r_viral[i] given j positive samples
      real nb_lpmf_val; 
      if (j == 0){
        nb_lpmf_val = r_viral[i] == 0 ? 1 : negative_infinity();
      } else {
        nb_lpmf_val = neg_binomial_2_lpmf(r_viral[i] | r_total[i] * j * mu_1 / n_swabs[i], j * phi_1);
      }
      
      // Store the sum of log probabilities
      lpmf_sum[j + 1] = binomial_lpmf_val + nb_lpmf_val;
    }
    
    // Marginalize out n_pos[i] using log-sum-exp for numerical stability
    target += log_sum_exp(lpmf_sum);
  }
}
")

swab_metadata <- read_tsv("../tables/swab-sample-metadata.tsv") %>%
  mutate(date = ymd(paste0("20", date))) %>%
  glimpse

# Use the per-day data. Caveat about the different treatments
swab_reads <- read_tsv("../tables/swabs-ra-summary.tsv") %>%
  mutate(date = ymd(paste0("20", date))) %>%
  select(-all_reads)

all_assignments <- swab_reads %>%
  select(assignment, species, group) %>%
  distinct()

# Fill in zero counts
swab_reads_complete <- swab_metadata %>%
  select(date, location) %>%
  crossing(all_assignments) %>%
  left_join(
    swab_reads,
    by = join_by(date, location, assignment, species, group),
  ) %>%
  mutate(
    non_dedup_hv = replace_na(non_dedup_hv, 0),
    dedup_hv = replace_na(dedup_hv, 0)
  ) %>%
  glimpse

# Just look at Rhinoviruses

swab_reads %>%
  filter(group == "Rhinoviruses") %>%
  print

rhino_data <- swab_reads_complete %>%
  filter(startsWith(species, "Rhinovirus")) %>%
  # Provisionally work on the species level
  summarise(.by = c(date, location, species), viral_reads = sum(dedup_hv)) %>%
  left_join(swab_metadata)

rhino_data %>%
  pivot_wider(id_cols = date, names_from = species, values_from = viral_reads) %>%
  left_join(swab_metadata) %>%
  select(-c(sample, location))

rhino_fits <- rhino_data %>%
  nest(.by = species) %>%
  mutate(
    stan_data = map(
      data,
      \(df) list(
        N_pools = nrow(df),
        n_swabs = df$pool_size,
        r_total = df$all_reads,
        r_viral = df$viral_reads
        )
    ),
    fit = map(
      stan_data,
      \(d) sampling(stan_model, d, iter = 4000, warmup = 1000, chains = 4, cores = 4, seed = 381928)
    ),
  )

rhino_posteriors <- rhino_fits %>%
  mutate(posterior = map(fit, compose(as_tibble, extract))) %>%
  select(species, posterior) %>%
  unnest(cols = c(posterior)) %>%
  mutate(
    across(-species, as.double),
    species = factor(species),
    ) %>%
  glimpse

# Compare posterior median with naive estimate

rhino_data %>%
  summarise(
    pos_pools_div_by_total_swabs = sum(viral_reads > 0) / sum(pool_size),
    .by = species,
    ) %>%
  left_join(
    rhino_posteriors %>%
      summarise(
        posterior_median = median(p),
        .by = species
      )
  )
 
rhino_data %>%
  summarise(
    pos_pools = sum(viral_reads > 0),
    total_pools = sum(pool_size),
    .by = species,
    ) %>%
  left_join(
    rhino_posteriors %>%
      summarise(
        posterior_mean_p = mean(p),
        .by = species
      )
  ) %>%
  mutate(mean_pos_swabs = posterior_mean_p * total_pools)

#  Plot posteriors

p_prior <- tibble(
  p = seq(0, 0.1, length.out = 501),
  density = dbeta(p, 0.5, 10.5),
  ) 

mu_prior <- tibble(
  mu_1 = seq(0, 0.01, length.out = 101),
  density = dbeta(mu_1, 0.5, 100.5),
  ) 

phi_prior <- tibble(
  phi_1 = seq(0, 10, length.out = 101),
  density = dexp(phi_1, rate = 1)
  ) 

rhino_posteriors %>%
  ggplot(aes(p, color=species)) +
  geom_density() +
  geom_line(data = p_prior, mapping = aes(y = density, color="Prior"), linetype="dashed") +
  scale_x_continuous(limits=c(0, 0.05), name = "prevalence") 

rhino_posteriors %>%
  ggplot(aes(mu_1, color=species)) +
  geom_density() +
  geom_line(data = mu_prior, mapping = aes(y = density, color="Prior"), linetype="dashed") +
  scale_x_continuous(limits=c(0, 0.01))

rhino_posteriors %>%
  ggplot(aes(phi_1, color=species)) +
  geom_density() +
  geom_line(data = phi_prior, mapping = aes(y = density, color="Prior"), linetype="dashed")

rhino_posteriors %>%
  ggplot(aes(mu_1, p)) +
  geom_density_2d() +
  facet_wrap(~species, scales="free")

rhino_posteriors %>%
  ggplot(aes(phi_1, p)) +
  geom_density_2d() +
  facet_wrap(~species, scales="free")

rhino_posteriors %>%
  ggplot(aes(phi_1, mu_1)) +
  geom_density_2d() +
  facet_wrap(~species)

# Posterior predictive distributions

posterior_predictive <- rhino_posteriors %>%
  # group_by(species) %>%
  # slice_sample(n = n_samples) %>%
  # mutate(rep = 1:n_samples) %>%
  # ungroup %>%
  mutate(rep = 1:n(), .by=species) %>%
  cross_join(swab_metadata) %>%
  mutate(
    n_pos = map2_int(pool_size, p, ~ rbinom(1, .x, .y)),
    mu = n_pos * mu_1 * all_reads / pool_size,
    phi = n_pos * phi_1,
    viral_reads = map2_int(mu, phi, ~ifelse(.x > 0, rnbinom(1, size = .y, mu = .x), 0)),
    false_negative = n_pos > 0 & viral_reads == 0,
  ) %>%
  glimpse 

posterior_predictive_summary <- posterior_predictive %>%
  summarise(
    n_pos_total = sum(n_pos),
    pos_pools = sum(n_pos > 0),
    viral_reads_total = sum(viral_reads),
    viral_reads_max = max(viral_reads),
    false_negative_pools = sum(false_negative),
    .by = c(species, rep),
  ) %>%
  glimpse

species_summary <- rhino_data %>%
  summarise(
    pos_pools = sum(viral_reads > 0),
    viral_reads_total = sum(viral_reads),
    viral_reads_max = max(viral_reads),
    .by = species
  )

posterior_predictive_summary %>%
  ggplot(aes(x = pos_pools)) +
  geom_histogram(binwidth = 1, center = 0, color = "grey") +
  geom_vline(data = species_summary, mapping = aes(xintercept=pos_pools)) +
  facet_grid(rows = vars(species)) +
  scale_x_continuous(breaks = seq(0, 12)) +
  theme_minimal()

posterior_predictive_summary %>%
  ggplot(aes(x = viral_reads_total, y = species, color = species)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(data = species_summary) +
  scale_x_continuous(limits = c(0, 400)) +
  theme_minimal()

posterior_predictive_summary %>%
  ggplot(aes(x = viral_reads_total, y = species, color = species)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(data = species_summary) +
  scale_x_continuous(transform = "log10") +
  theme_minimal()

posterior_predictive_summary %>%
  ggplot(aes(x = viral_reads_max, y = species, color = species)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(data = species_summary) +
  scale_x_continuous(limits = c(0, 200)) +
  theme_minimal()

posterior_predictive_summary %>%
  ggplot(aes(x = viral_reads_max, y = species, color = species)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(data = species_summary) +
  scale_x_continuous(transform = "log10") +
  theme_minimal()

posterior_predictive_summary %>%
  ggplot(aes(x = false_negative_pools)) +
  geom_histogram(binwidth = 1, center = 0, color = "grey") +
  facet_grid(rows = vars(species)) +
  scale_x_continuous(breaks = seq(0, 10)) +
  theme_minimal()
