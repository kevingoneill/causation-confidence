//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
  real ord_beta_lpdf(real y, real mu, real phi, vector cutpoints) {  
    //auxiliary variables
    real mu_logit = logit(mu);
        
    if (y == 0) {
      return log1m_inv_logit(mu_logit - cutpoints[1]);
    } else if (y == 1) {
      return log_inv_logit(mu_logit - cutpoints[2]);
    } else {
      return log(inv_logit(mu_logit - cutpoints[1]) - inv_logit(mu_logit - cutpoints[2])) +
        beta_proportion_lpdf(y | mu, phi);
    }
  }
  
  // vectorized version over y/mu
  real ord_beta_lpdf(vector y, vector mu, vector phi, vector cutpoints) {  
    int N = size(y);
    vector[N] lp;
    for (n in 1:N)
      lp[n] = ord_beta_lpdf(y[n] | mu[n], phi[n], cutpoints);
    return sum(lp);
  }
  
  // sample from ordered beta distribution
  real ord_beta_rng(real mu, real phi, vector cutpoints) {
    //auxiliary variables
    real mu_logit = logit(mu);
    real p_0 = 1 - inv_logit(logit(mu) - cutpoints[1]);
    real p_01 = inv_logit(logit(mu) - cutpoints[1]) - inv_logit(logit(mu) - cutpoints[2]);
    real p_1 = inv_logit(logit(mu) - cutpoints[2]);
    int mix = categorical_rng([p_0, p_01, p_1]');
    
    if (mix == 1)
      return(0);
    else if (mix == 2)
      return(beta_proportion_rng(mu, phi));
    else
      return(1);
  }
  
  // calculate expected value
  real ord_beta_epred(real mu, real phi, vector cutpoints) {
    //auxiliary variables
    real p_01 = inv_logit(logit(mu) - cutpoints[1]) - inv_logit(logit(mu) - cutpoints[2]);
    real p_1 = inv_logit(logit(mu) - cutpoints[2]);
    return p_01*mu + p_1;
  }
}

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0, upper=1> prior_only; // sample from the prior?
  int<lower=1> N;                   // number of data points
  int<lower=1> K;                   // number of predictors
  matrix[N, K] X;                   // design matrix
  int<lower=1> V;                   // number of vignettes
  int<lower=1> K_V;                 // number of vignette-level predictors
  array[N] int<lower=1> v;          // vignette indicator
  matrix[N, K_V] X_V;               // vignette-level design matrix
  vector[N] y;                      // response variable

  int<lower=1> N_pred;           // number of predictions
  matrix[N_pred, K] X_pred;      // design matrix for predictions
  matrix[N_pred, K] X_V_pred;    // vignette-level design matrix for predictions
}

transformed data {
  // double the number of vignette-level predictors
  // (one set for mu, one set for phi)
  int Kv = K_V*2;
}

parameters {
  ordered[2] cutpoints;  // bounds on logit scale to force 0/1 values
  
  // population-level coefficients
  vector[K] b;      // coefficients on mu
  vector[K] b_phi;  // coefficients on phi
  
  // vignette-level effects
  vector<lower=0>[Kv] sd_v;      // vignette-level SDs
  matrix[Kv, V] z_v;             // standardized vignette-level effects
  cholesky_factor_corr[Kv] L_v; // vignette-level effect correlations
}

transformed parameters {
  // actual vignette-level effects
  array[V] vector[K_V] r_v;
  array[V] vector[K_V] r_v_phi;

  {
    // pre-index vignette-level effects for efficiency
    matrix[Kv, V] R = diag_pre_multiply(sd_v, L_v) * z_v;
    for (i in 1:V) {
      r_v[i] = R[1:K_V, i];
      r_v_phi[i] = R[(K_V+1):Kv, i];
    }
  }
}

model {
  // priors for coefficients
  b ~ std_normal();
  b_phi ~ std_normal();
  cutpoints ~ normal(0, 2.5);
  
  // priors for vignette-level effects
  sd_v ~ normal(0, 1);
  to_vector(z_v) ~ std_normal();
  L_v ~ lkj_corr_cholesky(2);
  
  if (!prior_only) {
    vector[N] mu_logit = X*b;
    vector[N] phi_log = X*b_phi;
    
    for (n in 1:N) {
      mu_logit[n] += X_V[n,]*r_v[v[n]];
      phi_log[n] += X_V[n,]*r_v_phi[v[n]];
    }
    
    target += ord_beta_lpdf(y | inv_logit(mu_logit), exp(phi_log), cutpoints);
  }
}

generated quantities {
  // vignette-level effect correlations
  matrix[Kv, Kv] Omega = L_v * transpose(L_v);
  
  // vignette-average conditional parameters
  vector[N_pred] mu_pred = inv_logit(X_pred * b);
  vector[N_pred] phi_pred = exp(X_pred * b_phi);
  vector[N_pred] e_pred; // conditional means
  vector[N_pred] y_pred; // posterior predictive distribution
  for (i in 1:N_pred) {
    e_pred[i] = ord_beta_epred(mu_pred[i], phi_pred[i], cutpoints);
    y_pred[i] = ord_beta_rng(mu_pred[i], phi_pred[i], cutpoints);
  }
  
  // conditional parameters per vignette
  matrix[N_pred, V] mu_pred_v;
  matrix[N_pred, V] phi_pred_v;
  matrix[N_pred, V] e_pred_v;
  for (j in 1:V) {
    mu_pred_v[, j] = inv_logit(X_pred*b + X_V_pred*r_v[j]);
    phi_pred_v[, j] = exp(X_pred*b_phi + X_V_pred*r_v_phi[j]);
    
    for (i in 1:N_pred) {
      e_pred_v[i, j] = ord_beta_epred(mu_pred_v[i, j], phi_pred_v[i, j], cutpoints);
    }
  }
}
