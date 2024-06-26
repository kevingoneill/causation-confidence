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

  // calculate the kernel between two observations
  real k_00(vector x_i, vector x_j, real sq_alpha, vector rho) {
    return sq_alpha * exp(-0.5 * dot_self((x_i - x_j) ./ rho));
  }

  // calculate the kernel between two gradients dy/dx1 and dy/dx2
  real k_11(vector x_i, vector x_j, real sq_alpha, vector rho, int di, int dj) {
    int del = di == dj;
    return k_00(x_i, x_j, sq_alpha, rho) *
      (del - (x_i[di] - x_j[di])*(x_i[dj] - x_j[dj])/square(rho[dj])) / square(rho[di]);
  }

  // calculate the kernel between an observation and a gradient observation
  real k_10(vector x_i, vector x_j, real sq_alpha, vector rho, int di) {
    return k_00(x_i, x_j, sq_alpha, rho) * ((x_j[di] - x_i[di]) / square(rho[di]));
  }
  
  matrix L_cov_exp_quad_ARD(array[] vector x,
			    array[] int d,
                            real alpha,
                            vector rho,
                            real delta) {
    int N = size(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    
    for (i in 1:(N-1)) {
      if(d[i] == 0) {
        K[i, i] = sq_alpha + delta;
      } else {
        K[i, i] = sq_alpha / square(rho[d[i]]) + delta;
      }
      
      for (j in (i + 1):N) {
	if (d[i] == 0 && d[j] == 0) {
	  // both function observations
	  K[i, j] = k_00(x[i], x[j], sq_alpha, rho);
	} else if (d[i] > 0 && d[j] == 0) {
	  // one derivative observation
	  K[i, j] = k_10(x[i], x[j], sq_alpha, rho, d[i]);
	} else if (d[i] == 0 && d[j] > 0) {
	  // one derivative observation
	  K[i, j] = k_10(x[j], x[i], sq_alpha, rho, d[j]);
	} else {
	  // both derivative observations
	  K[i, j] = k_11(x[i], x[j], sq_alpha, rho, d[i], d[j]);
	}
	
	K[j, i] = K[i, j];
      }      
    }
    
    if(d[N] == 0) {
      K[N, N] = sq_alpha + delta;
    } else {
      K[N, N] = sq_alpha / square(rho[d[N]]) + delta;
    }
    return cholesky_decompose(K);
  }

  
}

data {
  int<lower=0, upper=1> prior_only; // sample from the prior?
  
  int<lower=1> N;     // number of data points
  int<lower=1> G;     // number of grid points (100: 10 P(C) X 10 P(A))
  int<lower=1> Dx;    // number of X dimensions (2: P(C) and P(A))
  int<lower=1> Dy;    // number of Y dimensions (2: causal judgments & confidence ratings)
  int<lower=1> C;     // number of conditions (2: conjunctive and disjunctive)
  int<lower=1> V;     // number of vignettes (6)

  array[N] int<lower=1, upper=G> g;  // grid index for each trial
  array[N] int<lower=1, upper=C> c;  // condition index for each trial
  array[N] int<lower=1, upper=V> v;  // vignette index for each trial
    
  array[G] vector[Dx] x;     // input data (same across trials)
  matrix<lower=0, upper=1>[N, Dy] y;
}

transformed data {
  real delta = 1e-9;
  int<lower=1> DY = Dy * 2;  // number of Y dimensions * 2 (for mean/precision estimators)

  // input data for raw/gradient observations
  // first dimension for GP function, all others for partial derivatives
  int GD = G * (Dx+1);
  array[GD] vector[Dx] X;

  // dimension of derivative:
  //   0 = to function values
  //   1 = partial derivative over Dx==1
  //   2 = partial derivative over Dx==2, etc
  array[GD] int<lower=0, upper=Dx> d;

  for (dim in 0:Dx) {
    for (i in 1:G) {
      X[G*dim + i] = x[i];
      d[G*dim + i] = dim;
    }
  }
}

parameters {
  // GP (hyper)parameters
  array[C] vector<lower=0>[Dx] rho;             // GP length-scale
  array[C] vector<lower=0>[Dx] rho_tilde;       // vignette-level GP length-scale
  array[C] vector<lower=0>[DY] alpha;           // group-level GP marginal SD
  array[C, V] vector<lower=0>[DY] alpha_tilde;  // vignette-level GP marginal SD
  array[C] matrix[GD, DY] eta;                  // group-level GP mean/variance variates
  array[C, V] matrix[GD, DY] eta_tilde;         // vignette-level GP mean/variance variates
  array[C] cholesky_factor_corr[DY] L_Omega;              // group-level correlation matrix between DVs
  array[C] cholesky_factor_corr[DY] L_Omega_tilde;   // vignette-level correlation matrix between DVs
  
  // parameters for ordered beta
  array[Dy] ordered[2] cutpoints;           // bounds to force 0/1 values
}

transformed parameters {
  array[C] matrix[GD, DY] f;                 // group-level GP mean/gradient
  array[C, V] matrix[GD, DY] f_tilde;        // vignette-level GP mean/gradient
  
  for (i in 1:C) {
    matrix[GD, GD] L_x_tilde = L_cov_exp_quad_ARD(X, d, 1.0, rho_tilde[i], delta);
    
    // (fix alpha=1 to remove colinearity with L_Omega)
    f[i] = L_cov_exp_quad_ARD(X, d, 1.0, rho[i], delta) * eta[i] * diag_pre_multiply(alpha[i], L_Omega[i])';
    for (j in 1:V) {
      f_tilde[i, j] = L_x_tilde * eta_tilde[i, j] * diag_pre_multiply(alpha_tilde[i, j], L_Omega_tilde[i])';
    }
  }
}

model {
  for (i in 1:C) {
    target += inv_gamma_lpdf(rho[i] | 6, 2);
    target += std_normal_lpdf(alpha[i]);
    target += lkj_corr_cholesky_lpdf(L_Omega[i] | 3);
    target += std_normal_lpdf(to_vector(eta[i]));
    
    target += inv_gamma_lpdf(rho_tilde[i] | 6, 2);
    target += lkj_corr_cholesky_lpdf(L_Omega_tilde[i] | 3);
    for (j in 1:V) {      
      target += std_normal_lpdf(alpha_tilde[i, j]);
      target += std_normal_lpdf(to_vector(eta_tilde[i, j]));
    }
  }
  
  for (dy in 1:Dy)
    target += normal_lpdf(cutpoints[dy] | 0, 10);
  
  if (!prior_only) {
    for (n in 1:N) {
      for (dy in 1:Dy) {
	target += ord_beta_lpdf(y[n, dy] | inv_logit(f[c[n], g[n], dy] + f_tilde[c[n], v[n], g[n], dy]),
				exp(f[c[n], g[n], Dy+dy] + f_tilde[c[n], v[n], g[n], Dy+dy]),
				cutpoints[dy]);
      }
    }
  }
}

generated quantities {
  array[C] matrix[DY, DY] Omega;
  array[C] matrix[DY, DY] Omega_tilde;
  for (i in 1:C) {
    Omega[i] = L_Omega[i] * L_Omega[i]';
    Omega_tilde[i] = L_Omega_tilde[i] * L_Omega_tilde[i]';
  }  

  // beta mean parameter
  array[C] matrix[G, Dy] mu_logit = f[:, 1:G, 1:Dy];
  array[C] matrix[G, Dy] mu = inv_logit(mu_logit);
  array[C, V] matrix[G, Dy] mu_logit_tilde;
  array[C, V] matrix[G, Dy] mu_tilde;
  array[C, Dx] matrix[G, Dy] mu_logit_grad;
  array[C, V, Dx] matrix[G, Dy] mu_logit_tilde_grad;
  array[C, Dx] matrix[G, Dy] mu_grad;
  array[C, V, Dx] matrix[G, Dy] mu_tilde_grad;

  // beta precision paremter
  array[C] matrix[G, Dy] phi_log = f[:, 1:G, (Dy+1):DY];
  array[C] matrix[G, Dy] phi = exp(phi_log);
  array[C, V] matrix[G, Dy] phi_log_tilde; 
  array[C, V] matrix[G, Dy] phi_tilde;
  array[C, Dx] matrix[G, Dy] phi_log_grad;
  array[C, V, Dx] matrix[G, Dy] phi_log_tilde_grad;
  array[C, Dx] matrix[G, Dy] phi_grad;
  array[C, V, Dx] matrix[G, Dy] phi_tilde_grad;

  // mixture probabilities
  array[C] matrix[G, Dy] p_0;
  array[C] matrix[G, Dy] p_1;
  array[C] matrix[G, Dy] p_01;
  array[C, V] matrix[G, Dy] p_0_tilde;
  array[C, V] matrix[G, Dy] p_1_tilde;
  array[C, V] matrix[G, Dy] p_01_tilde;
  array[C, Dx] matrix[G, Dy] p_0_grad;
  array[C, Dx] matrix[G, Dy] p_1_grad;
  array[C, Dx] matrix[G, Dy] p_01_grad;
  array[C, V, Dx] matrix[G, Dy] p_0_tilde_grad;
  array[C, V, Dx] matrix[G, Dy] p_1_tilde_grad;
  array[C, V, Dx] matrix[G, Dy] p_01_tilde_grad;
  
  // posterior predictive mean and variance
  array[C] matrix[G, Dy] e_pred;
  array[C, Dx] matrix[G, Dy] e_pred_grad;
  array[C, V] matrix[G, Dy] e_pred_tilde;
  array[C, V, Dx] matrix[G, Dy] e_pred_tilde_grad;
  array[C] matrix[G, Dy] var_pred;
  array[C, Dx] matrix[G, Dy] var_pred_grad;
  array[C, V] matrix[G, Dy] var_pred_tilde;
  array[C, V, Dx] matrix[G, Dy] var_pred_tilde_grad;

  for (i in 1:C) {
    for (j in 1:V) {
      mu_logit_tilde[i, j] = mu_logit[i] + f_tilde[i, j, 1:G, 1:Dy];
      mu_tilde[i, j] = inv_logit(mu_logit_tilde[i, j]);
      phi_log_tilde[i, j] = phi_log[i] + f_tilde[i, j, 1:G, (Dy+1):DY];
      phi_tilde[i, j] = exp(phi_log_tilde[i, j]);
    }
      
    // mixture probabilities
    for (j in 1:G) {
      for (k in 1:Dy) {
	p_0[i, j, k] = 1 - inv_logit(mu_logit[i,j,k] - cutpoints[k,1]);
	p_1[i, j, k] = inv_logit(mu_logit[i,j,k] - cutpoints[k,2]);
	p_01[i, j, k] = inv_logit(mu_logit[i,j,k] - cutpoints[k,1]) - inv_logit(mu_logit[i,j,k] - cutpoints[k,2]);

	for (l in 1:V) {
	  p_0_tilde[i, l, j, k] = 1 - inv_logit(mu_logit_tilde[i,l,j,k] - cutpoints[k,1]);
	  p_1_tilde[i, l, j, k] = inv_logit(mu_logit_tilde[i,l,j,k] - cutpoints[k,2]);
	  p_01_tilde[i, l, j, k] = inv_logit(mu_logit_tilde[i,l,j,k] - cutpoints[k,1]) - inv_logit(mu_logit_tilde[i,l,j,k] - cutpoints[k,2]);
	}
      }
    }

    // posterior predictive mean and variance
    e_pred[i] = p_01[i] .* mu[i] + p_1[i];
    var_pred[i] = p_1[i].*(1-p_1[i]) + p_01[i].*mu[i].*(1-mu[i])./(phi[i]+1) + p_01[i].*(1-p_01[i]).*square(mu[i]) - 2.*p_1[i].*p_01[i].*mu[i];
    for (j in 1:V) {
      e_pred_tilde[i,j] = p_01_tilde[i,j] .* mu_tilde[i,j] + p_1_tilde[i,j];
      var_pred_tilde[i,j] = p_1_tilde[i,j].*(1-p_1_tilde[i,j]) + p_01_tilde[i,j].*mu_tilde[i,j].*(1-mu_tilde[i,j])./(phi_tilde[i,j]+1) + p_01_tilde[i,j].*(1-p_01_tilde[i,j]).*square(mu_tilde[i,j]) - 2.*p_1_tilde[i,j].*p_01_tilde[i,j].*mu_tilde[i,j];
    }
    
    // gradients
    for (j in 1:Dx) {
      mu_logit_grad[i,j] = f[i, (G*j+1):(G*(j+1)), 1:Dy];
      mu_grad[i,j] = mu[i] .* (1-mu[i]) .* mu_logit_grad[i,j];

      phi_log_grad[i,j] = f[i, (G*j+1):(G*(j+1)), (Dy+1):DY];
      phi_grad[i,j] = phi[i] .* phi_log_grad[i,j];

      p_0_grad[i,j] = -p_0[i] .* (1-p_0[i]) .* mu_logit_grad[i,j];
      p_1_grad[i,j] = p_1[i] .* (1-p_1[i]) .* mu_logit_grad[i,j];
      p_01_grad[i,j] = (p_0[i].*(1-p_0[i]) - p_1[i].*(1-p_1[i])) .* mu_logit_grad[i,j];

      e_pred_grad[i,j] = p_1_grad[i,j] + p_01_grad[i,j].*mu[i] + p_01[i].*mu_grad[i,j];
      var_pred_grad[i,j] = (1 - 2*p_1[i]).*p_1_grad[i,j] +
	p_01_grad[i,j].*mu[i].*(1 - mu[i])./(phi[i] + 1) +
	p_01[i].* ((1-2*mu[i]).*mu_grad[i,j].*(phi[i] + 1) - mu[i].*(1-mu[i]).*phi_grad[i,j]) ./ square(phi[i] + 1) +
	p_01_grad[i,j].*(1-2*p_01[i]).*square(mu[i]) + 2.*p_01[i].*(1-p_01[i]).*mu[i].*mu_grad[i,j] -
	2*p_1_grad[i,j].*p_01[i].*mu[i] - 2*p_1[i].*p_01_grad[i,j].*mu[i] - 2*p_1[i].*p_01[i].*mu_grad[i,j];

      for (k in 1:V) {
	mu_logit_tilde_grad[i,k,j] = mu_logit_grad[i,j] + f_tilde[i, k, (G*j+1):(G*(j+1)), 1:Dy];
	mu_tilde_grad[i,k,j] = mu_tilde[i,k] .* (1-mu_tilde[i,k]) .* mu_logit_tilde_grad[i,k,j];

	phi_log_tilde_grad[i,k,j] = phi_log_grad[i,j] + f_tilde[i, k, (G*j+1):(G*(j+1)), (Dy+1):DY];
	phi_tilde_grad[i,k,j] = phi_tilde[i,k] .* phi_log_tilde_grad[i,k,j];
	
	p_0_tilde_grad[i,k,j] = -p_0_tilde[i,k] .* (1-p_0_tilde[i,k]) .* mu_logit_tilde_grad[i,k,j];
	p_1_tilde_grad[i,k,j] = p_1_tilde[i,k] .* (1-p_1_tilde[i,k]) .* mu_logit_tilde_grad[i,k,j];
	p_01_tilde_grad[i,k,j] = (p_0_tilde[i,k].*(1-p_0_tilde[i,k]) - p_1_tilde[i,k].*(1-p_1_tilde[i,k])) .* mu_logit_tilde_grad[i,k,j];
	
	e_pred_tilde_grad[i,k,j] = p_1_tilde_grad[i,k,j] + p_01_tilde_grad[i,k,j].*mu_tilde[i,k] + p_01_tilde[i,k].*mu_tilde_grad[i,k,j];
	var_pred_tilde_grad[i,k,j] = (1 - 2*p_1_tilde[i,k]).*p_1_tilde_grad[i,k,j] +
	  p_01_tilde_grad[i,k,j].*mu_tilde[i,k].*(1 - mu_tilde[i,k])./(phi_tilde[i,k] + 1) +
	  p_01_tilde[i,k].* ((1-2*mu_tilde[i,k]).*mu_tilde_grad[i,k,j].*(phi_tilde[i,k] + 1) - mu_tilde[i,k].*(1-mu_tilde[i,k]).*phi_tilde_grad[i,k,j]) ./ square(phi_tilde[i,k] + 1) +
	  p_01_tilde_grad[i,k,j].*(1-2*p_01_tilde[i,k]).*square(mu_tilde[i,k]) + 2.*p_01_tilde[i,k].*(1-p_01_tilde[i,k]).*mu_tilde[i,k].*mu_tilde_grad[i,k,j] -
	  2*p_1_tilde_grad[i,k,j].*p_01_tilde[i,k].*mu_tilde[i,k] - 2*p_1_tilde[i,k].*p_01_tilde_grad[i,k,j].*mu_tilde[i,k] - 2*p_1_tilde[i,k].*p_01_tilde[i,k].*mu_tilde_grad[i,k,j];
      }
    }
  }

  matrix[N, Dy] log_lik;
  matrix<lower=0, upper=1>[N, Dy] y_hat;
  for (n in 1:N) {
    for (dy in 1:Dy) {
      log_lik[n, dy] = ord_beta_lpdf(y[n, dy] | inv_logit(f[c[n], g[n], dy] + f_tilde[c[n], v[n], g[n], dy]),
				     exp(f[c[n], g[n], Dy+dy] + f_tilde[c[n], v[n], g[n], Dy+dy]),
				     cutpoints[dy]);
      
      int mixture = categorical_rng([p_0_tilde[c[n], v[n], g[n], dy],
				     p_01_tilde[c[n], v[n], g[n], dy],
				     p_1_tilde[c[n], v[n], g[n], dy]]');
      if (mixture == 1)
	y_hat[n, dy] = 0;
      else if (mixture == 2)
	y_hat[n, dy] = beta_proportion_rng(mu_tilde[c[n], v[n], g[n], dy], phi_tilde[c[n], v[n], g[n], dy]);
      else
	y_hat[n, dy] = 1;      
    }
  }
}
