# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
cmdstanr::install_cmdstan()

pak::pkg_install(c("rmcelreath/rethinking", "bayesplot", "posterior", "ggplot2", "cowplot", "here", "patchwork", "invgamma"))

library(rethinking)
library(cmdstanr)
library(ggplot2)
library(cowplot)
library(posterior)
library(bayesplot)
library(patchwork)
library(invgamma)
color_scheme_set("brightblue")

n_ind = 200
n_pooled = 50
level_i = matrix(rnorm(n_ind))
mean_qc = mean(level_i)
level_qc = rnorm(n_pooled, mean_qc, sd = 0.1)
dens(level_qc)
abline(v = mean_qc)

n_batch = 2
beta_batch = rnorm(n_batch, sd = 1)

mu = rnorm(1, 3)

meta_eta = rinvgamal
alphasq <- abs(rnorm(2)) #rexp(n_batch,2)
rhosq <- abs(rnorm(2) #rexp(n_batch,1)
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

  K = rethinking::cov_GPL2(Ds[[b]], alphasq[b], rhosq[b], 0.01)[[1]]
  data$deltas[data$batch==b] = t(chol(K)) %*% deltas_tilde[data$batch==b]
}

png("test.png", height= 1080, width = 1080)
plot( NULL , xlim=c(0,5) , ylim=c(0,2) , xlab="Scaled time difference" , ylab="covariance" )
for ( i in 1:n_batch )
  curve( alphasq[i]*exp(-rhosq[i]*x^2) , add=TRUE , lwd=4 , col=col.alpha(2,0.5) )
dev.off()

data = within(data, {
  levels = c(level_qc, level_i)
  y = mu + beta_batch[batch] + levels + deltas * 2
})
p = ggplot(data, aes(t, y, color = as.factor(QC))) +
  geom_point() +
  geom_point(aes(y = deltas), color  = "DarkOrange") +
  geom_line(data = dplyr::filter(data, QC == 1)) +
  facet_wrap(~batch) +
  theme_cowplot() + theme(legend.position = "none", 
                          plot.background = element_rect(fill = "white"),
                          panel.background = element_rect(fill = "white"))
save_plot("test.png", p, base_height = 8)

mod <- cmdstan_model(here::here("stanfiles/metabolomicsModelWithDrift.stan"))

# data {
#   int<lower=1> N;                           // total number of rows in data
#   int<lower=1> N_batch;
#   int<lower=1> N_ids;
#   array[N_batch] int<lower=1, upper=N> batch_size;
#   array[N] int<lower=1, upper=N_ids> id;
#   array[N] int<lower=1> batch;
#   array[N] int t;
#   array[N] real y;
# }
data = dplyr::arrange(data, batch, t)
data_list = list(
  N = nrow(data),
  N_batch = n_batch,
  N_ids = length(unique(data$ID)),
  batch_size = table(data$batch),
  id = data$ID,
  batch = data$batch,
  t = data$t / max(data$t) * 6,
  y = data$y
)
fit <- mod$sample(
  data = data_list,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 200 # print update every 500 iters
)
# plot posterior kernel
batch_mean = tapply(data$y, data$batch, mean)
fit$summary()

p = mcmc_recover_hist(fit$draws("batch_means"), batch_mean) /
(mcmc_recover_intervals(fit$draws("eta"), etasq) + mcmc_recover_intervals(fit$draws("rho"), rhosq))
save_plot("test.png", p, base_height = 9)

p = mcmc_recover_intervals(fit$draws("level"), c(mean_qc, level_i)) /
mcmc_recover_intervals(fit$draws("deltas"), data$deltas) 
save_plot("test.png", p, base_height = 9)



fit$draws("deltas")
