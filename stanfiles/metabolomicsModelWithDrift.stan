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
  real epslon = 0.01;
}
parameters {
    real mu_0;
    array[N_batch] real b;
    array[N_ids] real level;
    real<lower=0> sigma;
    vector[N] delta_tilde;
    array[N_batch] real<lower=0> eta;
    array[N_batch] real<lower=0> rho;
}
model {
    vector[N] deltas;
    array[N] real theta;
    matrix[N, N] K = rep_matrix(0, N, N);
    matrix[N, N] L_K;
    int pos;

    /* Gaussian Processes for the deltas */
    pos = 1;
    for (i in 1:N_batch) {
      K[pos:batch_size[i], pos:batch_size[i]] = gp_exp_quad_cov(t[pos:batch_size[i]], eta[i], rho[i]);
      pos = pos + batch_size[i];
    }

    L_K = cholesky_decompose(add_diag(K, epslon));
    deltas = L_K * delta_tilde;

    /* Linear model for the mean */
    for (i in 1:N){
        theta[i] = mu_0 + b[batch[i]] + level[id[i]] + deltas[i];
    }

    /* Likelihood */
    y ~ normal(theta, sigma);

    /* Priors */
    mu_0 ~ normal(0, 5);
    b ~ normal(0, 0.3);
    level ~ std_normal();
    sigma ~ exponential(1);

    /* GP Priors */
    delta_tilde ~ std_normal();
    eta ~ exponential(2);
    rho ~ exponential(1);

}
generated quantities{
  array[N_batch] real batch_means;

  for(i in 1:N_batch){
    batch_means[i] = mu_0 + b[i];
  }
}
