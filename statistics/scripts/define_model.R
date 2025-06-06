library(stringr)

# Priors

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

# Model

model_code <-  str_interp("
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
")