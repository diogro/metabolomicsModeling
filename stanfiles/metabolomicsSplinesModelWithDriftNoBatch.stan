functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}

data {
  int<lower=1> N;    // total number of rows in data
  int<lower=1> N_ids;
  array[N] int<lower=1, upper=N_ids> id;
  int N_knots;            // num of knots
  vector[N_knots] knots;  // the sequence of knots
  int spline_degree;        // the degree of spline (is equal to order - 1)
  array[N] real Y;
  array[N] real X;
}

transformed data {
  int N_basis = N_knots + spline_degree - 1;      // total number of B-splines
  matrix[N_basis, N] B;                           // matrix of B-splines
  vector[spline_degree + N_knots] ext_knots_temp;
  vector[2*spline_degree + N_knots] ext_knots;    // set of extended knots
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[N_knots], spline_degree));
  for (ind in 1:N_basis)
    B[ind,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));
  B[N_knots + spline_degree - 1, N] = 1;
}

parameters {
  array[N_ids] real level;
  real<lower=0> sigma_residual;

  row_vector[N_basis] a_raw;
  real a0;  // intercept
  real<lower=0> sigma;
  real<lower=0> tau;
}

transformed parameters {
  row_vector[N_basis] a;
  vector[N] deltas;
  a[1] = a_raw[1];
  for (i in 2:N_basis)
    a[i] = a[i-1] + a_raw[i]*tau;
  deltas = a0*to_vector(X) + to_vector(a*B);
}

model {
  array[N] real theta;

  /* Linear model for the mean */
  for (i in 1:N){
    theta[i] = level[id[i]] + deltas[i];
  }

  /* Likelihood */
  Y ~ normal(theta, sigma_residual);

  /* Priors */
  level ~ std_normal();
  sigma_residual  ~ std_normal();

  // Splines Priors
  a_raw ~ normal(0, 1);
  a0 ~ normal(0, 1);
  tau ~ normal(0, 1);
  sigma ~ normal(0, 1);
}
