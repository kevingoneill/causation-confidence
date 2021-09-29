functions {
  real ord_beta_reg_lpdf(real y, real mu, real phi, vector cutpoints) {  
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
  
  matrix L_cov_exp_quad_ARD(vector[] x,
			    int[] d,
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

  int<lower=1, upper=G> g[N];  // grid index for each trial
  int<lower=1, upper=C> c[N];  // condition index for each trial
  int<lower=1, upper=V> v[N];  // vignette index for each trial
    
  vector[Dx] x[G];     // input data (same across trials)
  matrix<lower=0, upper=1>[N, Dy] y;
}

transformed data {
  real delta = 1e-9;

  // input data for raw/gradient observations
  // first dimension for GP function, all others for partial derivatives
  int GD = G * (Dx+1);
  vector[Dx] X[GD];

  // dimension of derivative:
  //   0 = to function values
  //   1 = partial derivative over Dx==1
  //   2 = partial derivative over Dx==2, etc
  int<lower=0, upper=Dx> d[GD];

  for (dim in 0:Dx) {
    for (i in 1:G) {
      X[G*dim + i] = x[i];
      d[G*dim + i] = dim;
    }
  }

  print(d);
}

parameters {
  // GP (hyper)parameters
  vector<lower=0>[Dx] rho[C];            // GP length-scale
  vector<lower=0>[Dy] alpha[C];          // group-level GP marginal SD
  cholesky_factor_corr[Dy] L_Omega[C];   // Correlation matrix between DVs
  matrix[GD, Dy] eta[C];                  // group-level GP mean variates
  
  // parameters for ordered beta
  vector<lower=0>[Dy] phi;            // beta precision
  ordered[2] cutpoints[Dy];           // bounds to force 0/1 values
}

transformed parameters {
  matrix[GD, Dy] f[C];                 // group-level GP mean/gradient
  
  for (i in 1:C) {
    // cache kernel for the covariate grid (L_x)
    // (fix alpha=1 to remove colinearity with L_Omega)
    matrix[GD, GD] L_x = L_cov_exp_quad_ARD(X, d, 1.0, rho[i], delta);
    
    f[i] =  L_x * eta[i] * diag_pre_multiply(alpha[i], L_Omega[i])';
  }
}

model {
  for (i in 1:C) {
    rho[i] ~ inv_gamma(4, 0.5);
    alpha[i] ~ normal(0, 2);
    L_Omega[i] ~ lkj_corr_cholesky(3);
    
    to_vector(eta[i]) ~ std_normal();
  }
  
  for (dy in 1:Dy)
    cutpoints[dy, 2] - cutpoints[dy, 1] ~ normal(0, 3);
  
  phi ~ std_normal();
  
  
  if (!prior_only) {
    for (n in 1:N) {
      for (dy in 1:Dy)
	target += ord_beta_reg_lpdf(y[n, dy] | inv_logit(f[c[n], g[n], dy]), phi[dy], cutpoints[dy]);
    }
  }
}

generated quantities {
  matrix[Dy, Dy] Omega[C];
  for (i in 1:C)
    Omega[i] = L_Omega[i] * L_Omega[i]';
}
