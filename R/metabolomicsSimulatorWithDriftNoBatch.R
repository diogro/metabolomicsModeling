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
n_pooled = 200
N = n_ind + n_pooled
level_i = matrix(rnorm(n_ind))
mean_qc = mean(level_i)
level_qc = rnorm(n_pooled, mean_qc, sd = 0.1)
dens(level_qc)
abline(v = mean_qc)

mu = rnorm(1, 3)

meta_eta = rinvgamma(1, 2, 0.5)
alphasq <- 1#abs(rnorm(2)) #rexp(n_batch,2)
rhosq <- 1#abs(rnorm(2) #rexp(n_batch,1)
deltas_tilde = rnorm(n_ind + n_pooled)


data = data.frame(
           QC = c(rep(1, n_pooled), rep(0, n_ind)),
           ID = c(rep(1, n_pooled), 2:(n_ind+1)))

data$t = 0
data$deltas = rep(1, n_ind + n_pooled)
t_QC = floor(seq(1, n_ind + n_pooled, length.out = n_pooled))
data[data$QC==1,"t"] = t_QC

t_ind = seq(1, N)
t_ind = t_ind[!t_ind %in% t_QC]
data[data$QC!=1,"t"] = t_ind

s_ts = (data$t - mean(data$t))/sd(data$t)
Ds = outer(s_ts, s_ts, function(x, y) abs(x - y))

K = rethinking::cov_GPL2(Ds, alphasq, rhosq, 0.001)[[1]]
data$deltas = t(chol(K)) %*% deltas_tilde


png("test.png", height= 1080, width = 1080)
plot( NULL , xlim=c(min(D),max(D)) , ylim=c(0,2) , xlab="Scaled time difference" , ylab="covariance" )
for ( i in 1:n_batch )
  curve( alphasq[i]*exp(-rhosq[i]*x^2) , add=TRUE , lwd=4 , col=col.alpha(2,0.5) )
dev.off()

data = within(data, {
  levels = c(level_qc, level_i)
  y = mu + levels + deltas * 2
})
p = ggplot(data, aes(t, y, color = as.factor(QC))) +
  geom_point() +
  geom_point(aes(y = deltas), color  = "DarkOrange") +
  geom_line(data = dplyr::filter(data, QC == 1)) +
  theme_cowplot() + theme(legend.position = "none", 
                          plot.background = element_rect(fill = "white"),
                          panel.background = element_rect(fill = "white"))
save_plot("test.png", p, base_height = 8)

mod <- cmdstan_model(here::here("stanfiles/metabolomicsModelWithDriftNoBatch.stan"))

# data {
#     int<lower=1> N;                           // total number of rows in data
#     int<lower=1> N_ids;
#     array[N] int<lower=1, upper=N_ids> id;
#     array[N] int t;
#     array[N] real y;
# }
data_list = list(
  N = nrow(data),
  N_ids = length(unique(data$ID)),
  id = data$ID,
  t = data$t,
  y = as.vector(data$y)
)
fit <- mod$sample(
  data = data_list,
  iter_warmup = 20,
  iter_sampling = 20,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 5 # print update every 500 iters
)
# plot posterior kernel
fit$summary()

p = mcmc_recover_hist(fit$draws("mu_0"), mean(data$y)) /
(mcmc_recover_intervals(fit$draws("alpha"), alphasq) + mcmc_recover_intervals(fit$draws("rho"), rhosq))
save_plot("test.png", p, base_height = 9)

p = ggplot(data, aes(t, deltas)) + geom_point() + geom_point(col = 2, data = data.frame(deltas = fit$draws("deltas") |> colMeans() |> colMeans(), t = data$t))
save_plot("test.png", p, base_height = 9)


p = mcmc_recover_intervals(fit$draws("level"), c(mean_qc, level_i)) /
mcmc_recover_intervals(fit$draws("deltas"), data$deltas[order(data$t)]) 
save_plot("test.png", p, base_height = 9)



fit$draws("deltas")
