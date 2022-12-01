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

  vector[N] deltas;

  // Per-batch parameters (non-centered parameterization)
  array[N_batch] real rho_tilde;   //non-centered std of length scale
  array[N_batch] real alpha_tilde; //non-centered std of signal std dev
  array[N_batch] real sigma_tilde; //non-centered std of noise std dev
   
  // Population-level parameters
  real<lower=0> rho_m;   //median of rho population distribution
  real<lower=0> rho_s;   //std of rho population distribution
  real<lower=0> alpha_m; //median of alpha population distribution
  real<lower=0> alpha_s; //std of alpha population distribution
  real<lower=0> sigma_m; //median of sigma population distribution
  real<lower=0> sigma_s; //std of sigma population distribution
}
 
transformed parameters {
  // Per-batch parameters
  array[N_batch] real<lower=0> rho;   //length scale
  array[N_batch] real<lower=0> alpha; //signal standard deviation
  array[N_batch] real<lower=0> sigma; //noise standard deviation
   
  // Non-centered parameterization of per-batch parameters
  for (s in 1:N_batch) {
    rho[s] = exp(log(rho_m) + rho_s * rho_tilde[s]);
    alpha[s] = exp(log(alpha_m) + alpha_s * alpha_tilde[s]);
    sigma[s] = exp(log(sigma_m) + sigma_s * sigma_tilde[s]);
  }
}

model {
   
    array[N] real theta;
    matrix[N, N] K = rep_matrix(0, N, N);
    matrix[N, N] L_K;
    int pos;

    /* Gaussian Processes for the deltas */
    pos = 1;
    for (i in 1:N_batch) {
      K[pos:batch_size[i], pos:batch_size[i]] = add_diag(gp_exp_quad_cov(t[pos:batch_size[i]], alpha[i], rho[i]), sigma[i]);
      pos = pos + batch_size[i];
    }

    L_K = cholesky_decompose(add_diag(K, epsilon));
    deltas ~ multi_normal_cholesky(rep_vector(0, N), L_K);

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

    rho_m ~ inv_gamma(2, 0.5);  //median of rho population distribution
    rho_s ~ normal(0, .5);      //std of rho population distribution
    alpha_m ~ normal(0, 1);     //median of alpha population distribution
    alpha_s ~ normal(0, .5);    //std of alpha population distribution
    sigma_m ~ normal(0, 1);     //median of sigma population distribution
    sigma_s ~ normal(0, .5);    //std of sigma population distribution
      
    rho_tilde ~ std_normal();
    alpha_tilde ~ std_normal();
    sigma_tilde ~ std_normal();

}
generated quantities{
  array[N_batch] real batch_means;

  for(i in 1:N_batch){
    batch_means[i] = mu_0 + b[i];
  }
}
