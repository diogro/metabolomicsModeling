# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
cmdstanr::install_cmdstan(overwrite=FALSE)

pak::pkg_install(c("rmcelreath/rethinking", "bayesplot", "posterior", "ggplot2", "cowplot"))

library(rethinking)
library(cmdstanr)
library(ggplot2)
library(cowplot)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

n_ind = 300
n_pooled = 20
level_i = matrix(rnorm(n_ind))
mean_qc = mean(level_i)
level_qc = rnorm(n_pooled, mean_qc, sd = 0.1)
dens(level_qc)
abline(v = mean_qc)

n_batch = 3
beta_batch = rnorm(n_batch, sd = 1)

mu = rnorm(1, 3)

etasq <- rexp(n_batch,2)
rhosq <- rexp(n_batch,1)
deltas_tilde = rnorm(n_ind + n_pooled)

batch_id = c(sample(1:n_batch, n_pooled + n_ind, T))

data = data.frame(batch = batch_id,
           QC = c(rep(1, n_pooled), rep(0, n_ind)),
           ID = c(rep(1, n_pooled), 2:(n_ind+1)))

Ds = vector("list", n_batch)
data$t = 0
data$deltas = rep(1, n_ind + n_pooled)
for(b in 1:n_batch){
  n_in_b = sum(data$batch==b)
  n_qc_in_B = sum(data[data$batch==b,"QC"]==1)
  t_QC = floor(seq(1, n_in_b, length.out = n_qc_in_B))
  data[data$batch==b & data$QC==1,"t"] = t_QC

  t_ind = seq(1, n_in_b)
  t_ind = t_ind[!t_ind %in% t_QC]
  data[data$batch==b & data$QC!=1,"t"] = t_ind

  ts = data[data$batch==b,"t"]
  D = outer(ts, ts, function(x, y) abs(x - y))
  Ds[[b]] = D/max(D) * 5

  K = rethinking::cov_GPL2(Ds[[b]], etasq[b], rhosq[b], 0.01)[[1]]
  data$deltas[data$batch==b] = t(chol(K)) %*% deltas_tilde[data$batch==b]
}

plot( NULL , xlim=c(0,5) , ylim=c(0,2) , xlab="Scaled time difference" , ylab="covariance" )
for ( i in 1:n_batch )
  curve( etasq[i]*exp(-rhosq[i]*x^2) , add=TRUE , lwd=4 , col=col.alpha(2,0.5) )


data = within(data, {
  y = mu + beta_batch[batch] + c(level_qc, level_i)
})
ggplot(data, aes(t, y, color = as.factor(QC))) +
  geom_point() +
  geom_point(aes(y = deltas), color  = "DarkOrange") +
  geom_line(data = dplyr::filter(data, QC == 1)) +
  facet_wrap(~batch) +
  theme_cowplot() + theme(legend.position = "none")

mod <- cmdstan_model("metabolomicsModelNoDrift.stan")

# data {
#     int<lower=1> N;
#     int<lower=1> N_batch;
#     int<lower=1> N_ids;

#     array[N] int<lower=1> id;
#     array[N] int<lower=1> batch;

#     array[N] real<lower=0> y;
# }

data_list = list(
  N = nrow(data),
  N_batch = n_batch,
  N_ids = length(unique(data$ID)),
  id = data$ID,
  batch = data$batch,
  y = data$y
)
fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)
# plot posterior kernel
fit$summary()
mcmc_recover_hist(fit$draws("mu_0"), mu)
mcmc_recover_hist(fit$draws("b"), beta_batch)
mcmc_recover_hist(fit$draws("batch_means"), mu + beta_batch)

mcmc_recover_intervals(fit$draws("level"), c(mean_qc, level_i))




