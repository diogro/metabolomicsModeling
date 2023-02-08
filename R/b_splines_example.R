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

set.seed(1234)
num_knots <- 10 # true number of knots
spline_degree <- 3
num_basis <- num_knots + spline_degree - 1
X <- seq(from=-10, to=10, by=.1)
knots <- unname(quantile(X,probs=seq(from=0, to=1, length.out = num_knots)))
num_data <- length(X)
a0 <- 0.2
a <- rnorm(num_basis, 0, 1)
B_true <- t(bs(X, df=num_basis, degree=spline_degree, intercept = TRUE))
Y_true <- as.vector(a0*X + a%*%B_true)
Y <- Y_true + rnorm(length(X), 0, 0.2)

num_knots <- 100; # number of knots for fitting
num_basis <- num_knots + spline_degree - 1
knots <- unname(quantile(X,probs=seq(from=0, to=1, length.out = num_knots)))
sm <- cmdstan_model(here::here("stanfiles/fit_basis.stan"))

# data {
#   int num_data;             // number of data points
#   int num_knots;            // num of knots
#   vector[num_knots] knots;  // the sequence of knots
#   int spline_degree;        // the degree of spline (is equal to order - 1)
#   real Y[num_data];
#   real X[num_data];
# }
data_list <- list(num_data = num_data, 
                  num_knots = num_knots, 
                  knots = knots,
                  spline_degree = 3,
                  Y = Y, X = X)
fit <- sm$sample(data = data_list, chains = 4, parallel_chains = 4)
a0_e = mean(fit$draws('a0'))
a_e = colMeans(fit$draws('a', f = "draws_matrix"))
y_e = as.vector(a0_e*X + a_e%*%B)
png("test.png")
plot(y_e ~ X, pch = 19, col = 2)
points(Y_true ~ X, pch = 19)
dev.off()
