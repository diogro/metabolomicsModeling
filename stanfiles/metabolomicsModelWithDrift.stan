data {
    int<lower=1> N;                           // total number of rows in data
    int<lower=1> N_batch;
    int<lower=1> N_ids;
    array[N_batch] int batch_size;
    array[N] int<lower=1, upper=N_ids> id;
    array[N] int<lower=1> batch;
    array[N] int t;
    array[N] real y;
}
transformed data {
  real epsilon = 0.01;
}
parameters {
  real mu_0;
  array[N_batch] real b;
  array[N_ids] real level;
  real<lower=0> sigma_residual;

  vector[N] deltas_tilde;

  // Per-batch parameters 
  array[N_batch] real<lower=0> rho;   
  array[N_batch] real<lower=0> alpha; 
}
transformed parameters {
  vector[N] deltas;
  {
    matrix[N, N] K = rep_matrix(0, N, N);
    matrix[N, N] L_K;
    int pos;

    /* Gaussian Processes for the deltas */
    pos = 1;
    for (i in 1:N_batch) {
      K[pos:batch_size[i], pos:batch_size[i]] = gp_exp_quad_cov(t[pos:batch_size[i]], 
                                                                alpha[i], rho[i]);
      pos = pos + batch_size[i];
    }

    L_K = cholesky_decompose(add_diag(K, epsilon));
    deltas = L_K * deltas_tilde;
  }
}

model {
    array[N] real theta;

    /* Linear model for the mean */
    for (i in 1:N){
      theta[i] = mu_0 + b[batch[i]] + level[id[i]] + deltas[i];
    }

    /* Likelihood */
    y ~ normal(theta, sigma_residual);

    /* Priors */
    mu_0 ~ normal(0, 5);
    b ~ normal(0, 0.3);
    level ~ std_normal();
    sigma_residual ~ exponential(1);

    /* GP Priors */
    deltas_tilde ~ std_normal();
    rho ~ normal(0, 0.3);
    alpha ~ normal(0, 0.3);
}
generated quantities{
  array[N_batch] real batch_means;

  for(i in 1:N_batch){
    batch_means[i] = mu_0 + b[i];
  }
}
