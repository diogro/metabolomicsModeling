data {
    int<lower=1> N;
    int<lower=1> N_batch;
    int<lower=1> N_ids;
    array[N] int<lower=1, upper=N_ids> id;
    array[N] int<lower=1> batch;

    array[N] real y;
}

parameters {
    real mu_0;
    array[N_batch] real b;
    array[N_ids] real level;
    real<lower=0> sigma;
}

model {
    array[N] real theta;

    /* Linear model for the mean */
    for (i in 1:N){
        theta[i] = mu_0 + b[batch[i]] + level[id[i]];
    }

    /* Likelihood */
    y ~ normal(theta, sigma);

    /* Priors */
    mu_0 ~ normal(0, 5);
    b ~ normal(0, 0.3);
    level ~ normal(0, 1);
    sigma ~ exponential(1);
}
generated quantities{
  array[N_batch] real batch_means;

  for(i in 1:N_batch){
    batch_means[i] = mu_0 + b[i];
  }
}
