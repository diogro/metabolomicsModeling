data {
    int<lower=1> N;                           // total number of rows in data
    int<lower=1> N_ids;
    array[N] int<lower=1, upper=N_ids> id;
    array[N] real t;
    array[N] real y;
}
transformed data {
  real sigma_intercept = 0.1;
}
parameters {
  array[N_ids] real level;
  real<lower=0> sigma_residual;

  vector[N] deltas_tilde;

  // Per-batch parameters 
  real<lower=0> rho;   
  real<lower=0> alpha; 
}
transformed parameters {
  vector[N] deltas;
  {
    matrix[N, N] K = rep_matrix(0, N, N);
    matrix[N, N] L_K;

    /* Gaussian Processes for the deltas */
    K = gp_exp_quad_cov(tn, alpha, rho);
    L_K = cholesky_decompose(add_diag(K, sigma_intercept));
    deltas = L_K * deltas_tilde;
  }
}

model {
  array[N] real theta;

  /* Linear model for the mean */
  for (i in 1:N){
    theta[i] = level[id[i]] + deltas[i];
  }

  /* Likelihood */
  yn ~ normal(theta, sigma_residual);

  /* Priors */
  level ~ std_normal();
  sigma_residual  ~ std_normal();

  /* GP Priors */
  deltas_tilde ~ std_normal();
  rho ~ std_normal();
  alpha ~ std_normal();
}
