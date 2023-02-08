cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
cmdstanr::install_cmdstan(overwrite=FALSE)

pak::pkg_install(c("rmcelreath/rethinking", "bayesplot", "posterior", "ggplot2", "cowplot"))

library(rethinking)
library(cmdstanr)
library(ggplot2)
library(cowplot)
library(posterior)
library(bayesplot)
library(splines)
color_scheme_set("brightblue")

set.seed(123)
num_knots <- 10 # true number of knots
spline_degree <- 3
num_basis <- num_knots + spline_degree - 1
X <- seq(from=-1, to=1, by=0.01)
knots <- unname(quantile(X,probs=seq(from=0, to=1, length.out = num_knots)))
num_data <- length(X)
a0 <- 0.
a <- rnorm(num_basis, 0, 2)
B_true <- t(bs(X, df=num_basis, degree=spline_degree, intercept = TRUE))
drift <- as.vector(a0 + a%*%B_true)

n_pooled = 10
n_ind = num_data - n_pooled

level_i = matrix(rnorm(n_ind))
mean_qc = mean(level_i)
level_qc = rnorm(n_pooled, mean_qc, sd = 0.1)

data = data.frame(
           QC = c(rep(1, n_pooled), rep(0, n_ind)),
           ID = c(rep(1, n_pooled), 2:(n_ind+1)))

data$t = 0
data$deltas = rep(1, num_data)
t_QC = floor(seq(1, num_data, length.out = n_pooled))
data[data$QC==1,"t"] = t_QC

t_ind = seq(1, num_data)
t_ind = t_ind[!t_ind %in% t_QC]
data[data$QC!=1,"t"] = t_ind

data$deltas = drift[data$t]

mu = rnorm(1, 3)
data = within(data, {
  levels = c(level_qc, level_i)
  y = mu + levels + deltas
})
p = ggplot(data, aes(t, y, color = as.factor(QC))) +
  geom_point() +
  geom_point(aes(y = deltas), color  = "DarkOrange") +
  geom_line(data = dplyr::filter(data, QC == 1)) +
  theme_cowplot() + theme(legend.position = "none", 
                          plot.background = element_rect(fill = "white"),
                          panel.background = element_rect(fill = "white"))
save_plot("test.png", p, base_height = 8)

num_knots <- 20; # number of knots for fitting
num_basis <- num_knots + spline_degree - 1
knots <- unname(quantile(X,probs=seq(from=0, to=1, length.out = num_knots)))
sm <- cmdstan_model(here::here("stanfiles/metabolomicsSplinesModelWithDriftNoBatch.stan"))

# data {
#   int<lower=1> N;    // total number of rows in data
#   int<lower=1> N_ids;
#   array[N] int<lower=1, upper=N_ids> id;
#   int N_knots;            // num of knots
#   vector[N_knots] knots;  // the sequence of knots
#   int spline_degree;        // the degree of spline (is equal to order - 1)
#   array[N] real Y;
#   array[N] real X;
# }
data = data[order(data$t),]
data_list <- list(N = num_data, 
                  N_ids = n_ind + 1,
                  id = data$ID,
                  N_knots = num_knots, 
                  knots = knots,
                  spline_degree = 3,
                  Y = data$y, X = X)
fit <- sm$sample(data = data_list, chains = 4, parallel_chains = 4, adapt_delta = 0.99)
a0_e = mean(fit$draws('a0'))
a_e = colMeans(fit$draws('a', f = "draws_matrix"))
colMeans(fit$draws('a', f = "draws_matrix"))
num_basis <- num_knots + spline_degree - 1
B <- t(bs(X, df=num_basis, degree=spline_degree, intercept = TRUE))
y_e = as.vector(a0_e*X + a_e%*%B)
png("test.png")
plot(y_e ~ X, pch = 19, col = 2)
points(drift + mu ~ X, pch = 19)
dev.off()

levels_e = colMeans(fit$draws('level', f = "draws_matrix"))
levels = tapply(data$levels, data$ID, mean)

png("test.png")
plot(abs(levels_e - levels), pch = 19, col = 2)
dev.off()

png("test.png")
plot(levels_e ~ levels, pch = 19, col = 2)
abline(0, 1)
dev.off()

p = ggplot(data, aes(t, y, color = as.factor(QC))) +
  geom_point() +
  geom_point(aes(y = deltas), color  = "DarkOrange") +
  geom_line(data = dplyr::filter(data, QC == 1)) +
  geom_point(data = data.frame(deltas = colMeans(fit$draws('deltas', f = "draws_matrix")),
                               t = data$t), 
  aes(y = deltas), color  = 3) +
  theme_cowplot() + theme(legend.position = "none", 
                          plot.background = element_rect(fill = "white"),
                          panel.background = element_rect(fill = "white"))
save_plot("test.png", p, base_height = 8)
