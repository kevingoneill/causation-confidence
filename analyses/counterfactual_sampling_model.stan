functions {
  // ordered beta log likelihood
  real ord_beta_lpdf(real y, real mu, real phi, vector cutpoints) {  
    //auxiliary variables
    real mu_logit = logit(mu);
        
    if (y == 0)
      return log1m_inv_logit(mu_logit - cutpoints[1]);
    else if (y == 1)
      return log_inv_logit(mu_logit - cutpoints[2]);
    else
      return log_inv_logit_diff(mu_logit - cutpoints[1], mu_logit - cutpoints[2]) + beta_proportion_lpdf(y | mu, phi);
  }

  // vectorized version
  real ord_beta_lpdf(vector y, vector mu, vector phi, array[] vector cutpoints) {  
    int N = size(y);
    vector[N] lp;
    for (n in 1:N)
      lp[n] = ord_beta_lpdf(y[n] | mu[n], phi[n], cutpoints[n]);
    return sum(lp);
  }

  // get a 3-simplex of ordered beta response probabilities [P(y=0), P(0 < y < 1), P(y=1)]
  vector ord_beta_probs(real mu, real phi, vector cutpoints) {
    real p_0 = 1 - inv_logit(logit(mu) - cutpoints[1]);
    real p_01 = inv_logit(logit(mu) - cutpoints[1]) - inv_logit(logit(mu) - cutpoints[2]);
    real p_1 = inv_logit(logit(mu) - cutpoints[2]);
    return [p_0, p_01, p_1]';
  }

  // get the expected value of the ordered beta distribution
  real ord_beta_epred(real mu, real phi, vector cutpoints) {
    vector[3] p = ord_beta_probs(mu, phi, cutpoints);
    return p[2]*mu + p[3];
  }

  // get the expected variance of the ordered beta distribution
  real ord_beta_varpred(real mu, real phi, vector cutpoints) {
    vector[3] p = ord_beta_probs(mu, phi, cutpoints);
    return p[3]*(1-p[3]) + p[2]*mu*(1-mu)/(phi+1) + p[2]*(1-p[2])*square(mu) - 2*p[2]*p[3]*mu;
  }

  // sample from the ordered beta distribution
  real ord_beta_rng(real mu, real phi, vector cutpoints) {
    // determine if response is 0, continuous, or 1
    int i = categorical_rng(ord_beta_probs(mu, phi, cutpoints));
    
    if (i == 1)
      return 0;
    else if (i == 2)
      return beta_proportion_rng(mu, phi);
    else
      return 1;
  }
    
  /*
   * Calculate causal strength given P(C), P(A), and causal structure using one of five measures:
   *   measure=1: delta P model
   *   measure=2: power PC model
   *   measure=3: crediting causality model
   *   measure=4: necessity-sufficiency model
   *   measure=5: counterfactual effect size model
   *
   *  returned value is a tuple of two matrices:
   *    the first element is a matrix of causal strength values, and
   *    the second element is a matrix of normalization constants.
   *    for all measures except 5 (CES), the constant is fixed at 1.
   *    for measure=5, the constant is the SD_C / SD_E
   */
  tuple(matrix, matrix) causal_strength(vector pC, vector pA, int structure, int measure) {
    int M = size(pC);
    int N = size(pA);
    matrix[M, N] PC = rep_matrix(pC, N);
    matrix[M, N] PA = rep_matrix(pA', M);
    matrix[M, N] lambda = rep_matrix(1.0, M, N);
    
    if (measure == 1) {             // Delta P model
      if (structure == 1)
	return (PA, lambda);
      else
	return (1-PA, lambda);
    } else if (measure == 2) {      // Power PC model
      if (structure == 1)
	return (PA, lambda);
      else
	return (rep_matrix(1-1e-3, M, N), lambda);   // subtract small constant to avoid infinite gradients
    } else if (measure == 3) {      // Crediting Causality model
      if (structure == 1)
	return ((1-PC).*PA, lambda);
      else
	return ((1-PC).*(1-PA), lambda);
    } else if (measure == 4) {      // Necessity-sufficiency model
      if (structure == 1)
	return (PC.*PA - PC + 1, lambda);
      else
	return (PC, lambda);
    } else {                        // Counterfactual effect size model
      if (structure == 1)
	return (PA, sqrt(PC.*(1-PC) ./ (PC.*PA .* (1 - PC.*PA))));
      else
	return (1-PA, sqrt(PC.*(1-PC) ./ ((PC+PA-PC.*PA) .* (1-PC-PA+PC.*PA))));
    }
  }
  
  /*
   * Calculate inverse precision based on causal strength using one of four measures
   *   measure=1: variance
   *   measure=2: standard deviation
   *   measure=3: coefficient of variation
   *   measure=4: entropy
   */
  matrix inv_precision(matrix K, int measure, matrix lambda) {
    if (measure == 1)
      return square(lambda) .* K .* (1-K);
    else if (measure == 2)
      return lambda .* sqrt(K .* (1-K));
    else if (measure == 3)
      return log(lambda) + .5*log(1-K) - .5*log(K);   // put CV on log scale to avoid infinite values
    else
      return -K.*log(K) - (1-K).*log(1-K);
  }  
  
  /*
   * apply the ordered transform (as defined in the stan manual)
   * to map an unconstrained vector v to an ordered vector o such that:
   *  o[n] := v[1] + sum_{2:N} exp(v[n])
   */
  vector order(vector v) {
    int N = size(v);
    vector[N] o = rep_vector(v[1], N);
    o[2:N] += cumulative_sum(exp(v[2:N]));
    return o;
  }
}

data {
  int<lower=0, upper=1> prior_only;          // sample from the prior?
  int<lower=1, upper=5> cause_measure;       // which measure of causal strength? (1=DP, 2=PPC, 3=SP, 4=NS, 5=CES)
  int<lower=1, upper=4> confidence_measure;  // use which measure of confidence? (1=Var, 2=SD, 3=CV, 4=Entropy)
  
  int<lower=1> N;     // number of data points
  int<lower=1> S;     // number of causal structures (2: conjunctive and disjunctive)
  int<lower=1> P;     // number of grid points (10: .1 - 1)
  
  array[N] int<lower=1, upper=S> s;    // causal structure index for each trial (1:2)
  array[N] int<lower=1, upper=P> pc;   // P(C) index for each trial (1:10)
  array[N] int<lower=1, upper=P> pa;   // P(A) index for each trial

  vector<lower=0, upper=1>[N] cause;
  vector<lower=0, upper=1>[N] confidence;
}

parameters {
  // shift + scale parameters
  real shift_cause;
  real shift_confidence;
  real<lower=0> scale_cause;
  real<upper=0> scale_confidence;
  
  // parameters for ordered beta
  vector<lower=0>[S] phi_cause;
  vector<lower=0>[S] phi_confidence;
  array[S] ordered[2] cutpoints_cause;
  array[S] ordered[2] cutpoints_confidence;
  
  // counterfactual sampling probabilities on transformed scale to ensure 
  // (a) bounds of [0,1] and (b) positive effect of probability manipulation
  array[S] vector[P] z_logit_pC;
  array[S] vector[P] z_logit_pA;
}

transformed parameters {
  array[S] vector[P] pC;
  array[S] vector[P] pA;
  array[S] matrix[P, P] mu_cause;
  array[S] matrix[P, P] mu_confidence;

  for (i in 1:S) {
    // calculate counterfactual sampling probabilities
    pC[i] = inv_logit(order(z_logit_pC[i]));
    pA[i] = inv_logit(order(z_logit_pA[i]));

    // calculate mean & precision of counterfactual effects
    matrix[P, P] K;
    matrix[P, P] lambda;
    (K, lambda) = causal_strength(pC[i], pA[i], i, cause_measure);
    matrix[P, P] K_p = inv_precision(K, confidence_measure, lambda);
    
    // calculate mean parameters for causal judgments/confidence
    mu_cause[i] = inv_logit(lambda.*K*scale_cause + shift_cause);
    mu_confidence[i] = inv_logit(K_p*scale_confidence + shift_confidence);
  }
}

model {
  // model priors
  target += normal_lpdf(shift_cause | 0, .5);
  target += normal_lpdf(shift_confidence | 0, .5);
  target += normal_lpdf(scale_cause | 0, .5);
  target += normal_lpdf(scale_confidence | 0, .5);
  target += lognormal_lpdf(phi_cause | 0, .66);
  target += lognormal_lpdf(phi_confidence | 0, .66);
  for (i in 1:S) {
    target += std_normal_lpdf(cutpoints_cause[i]);
    target += std_normal_lpdf(cutpoints_confidence[i]);

    // put priors on P(C) and P(A) that loosely encourage linearity when transformed
    target += normal_lpdf(z_logit_pC[i, 1] | -1.5, 1.5);
    target += normal_lpdf(z_logit_pA[i, 1] | -1.5, 1.5);
    target += normal_lpdf(z_logit_pC[i, 2:P] | -1.25, .75);
    target += normal_lpdf(z_logit_pA[i, 2:P] | -1.25, .75);
  }
  
  if (!prior_only) {
    // pre-index means for efficiency
    vector[N] mu_cause_n;
    vector[N] mu_confidence_n;
    for (n in 1:N) {
      mu_cause_n[n] = mu_cause[s[n], pc[n], pa[n]];
      mu_confidence_n[n] = mu_confidence[s[n], pc[n], pa[n]];
    }
    
    // evaluate model likelihood
    target += ord_beta_lpdf(cause | mu_cause_n, phi_cause[s], cutpoints_cause[s]);
    target += ord_beta_lpdf(confidence | mu_confidence_n, phi_confidence[s], cutpoints_confidence[s]);
  }
}

generated quantities {
  // calculate expected value the posterior predictive distribution (i.e., conditional means)
  // and sample from the posterior predictive distribution
  array[S] matrix[P, P] epred_cause;
  array[S] matrix[P, P] epred_confidence;
  array[S] matrix[P, P] cause_hat;
  array[S] matrix[P, P] confidence_hat;
  for (i in 1:S) {
    for (j in 1:P) {
      for (k in 1:P) {
	epred_cause[i, j, k] = ord_beta_epred(mu_cause[i, j, k], phi_cause[i], cutpoints_cause[i]);
	epred_confidence[i, j, k] = ord_beta_epred(mu_confidence[i, j, k], phi_confidence[i], cutpoints_confidence[i]);
	cause_hat[i, j, k] = ord_beta_rng(mu_cause[i, j, k], phi_cause[i], cutpoints_cause[i]);
	confidence_hat[i, j, k] = ord_beta_rng(mu_confidence[i, j, k], phi_confidence[i], cutpoints_confidence[i]);
      }
    }
  }
  
  // calculate predictive log likelihood and variance explained (per Gelman et al., 2019)
  vector[N] log_lik_cause;
  vector[N] log_lik_confidence;
  real R2_cause;
  real R2_confidence;
  {
    // temporarily calculate pointwise fit & residual variance
    vector[N] epred_cause_n;
    vector[N] varpred_cause_n;
    vector[N] epred_confidence_n;
    vector[N] varpred_confidence_n;
    for (n in 1:N) {
      log_lik_cause[n] = ord_beta_lpdf(cause[n] | mu_cause[s[n], pc[n], pa[n]], phi_cause[s[n]], cutpoints_cause[s[n]]);
      log_lik_confidence[n] = ord_beta_lpdf(confidence[n] | mu_confidence[s[n], pc[n], pa[n]], phi_confidence[s[n]], cutpoints_confidence[s[n]]);
      epred_cause_n[n] = epred_cause[s[n], pc[n], pa[n]];
      epred_confidence_n[n] = epred_confidence[s[n], pc[n], pa[n]];
      varpred_cause_n[n] = ord_beta_varpred(mu_cause[s[n], pc[n], pa[n]], phi_cause[s[n]], cutpoints_cause[s[n]]);
      varpred_confidence_n[n] = ord_beta_varpred(mu_confidence[s[n], pc[n], pa[n]], phi_confidence[s[n]], cutpoints_confidence[s[n]]);
    }

    R2_cause = variance(epred_cause_n) / (variance(epred_cause_n) + mean(varpred_cause_n));
    R2_confidence = variance(epred_confidence_n) / (variance(epred_confidence_n) + mean(varpred_confidence_n));
  }
}
