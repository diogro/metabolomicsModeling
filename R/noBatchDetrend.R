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
library(rio)
library(tidyverse)
color_scheme_set("brightblue")

raw_data = import(here::here("output/ordered_samples.csv"))
miss_df = import(here::here("output/missingness.csv")) %>% filter(missingness < 1e-8)

metabolite_list = miss_df$metabolite
batch_list = unique(raw_data$batch)

metab = metabolite_list[2]
b = batch_list[1]
data = raw_data %>% 
    filter(metabolite == metab, batch == b) %>%
    mutate(ID = if_else(Pooled == 1, "Pool", ID), 
           num_ID = as.numeric(factor(ID)))
head(data)
p = ggplot(data, aes(ColPos, concentration, color = as.factor(Pooled))) +
  geom_point() +
  geom_line(data = dplyr::filter(data, Pooled == 1)) +
  theme_cowplot() + theme(legend.position = "none", 
                          plot.background = element_rect(fill = "white"),
                          panel.background = element_rect(fill = "white"))
save_plot("test.png", p, base_height = 8)

num_knots <- 4; # number of knots for fitting
num_basis <- num_knots + spline_degree - 1
knots <- unname(quantile(X,probs=seq(from=0, to=1, length.out = num_knots)))
sm <- cmdstan_model(here::here("stanfiles/metabolomicsSplinesModelWithDriftNoBatchRobust.stan"))

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
data = data[order(data$ColPos),]
x = data$ColPos
data_list <- list(N = nrow(data), 
                  N_ids = length(unique(data$ID)),
                  id = data$num_ID,
                  N_knots = num_knots, 
                  knots = knots,
                  spline_degree = 3,
                  Y = log(data$concentration), X = ((x - min(x))/(max(x) - min(x)) - 0.5 ) * 2)
fit <- sm$sample(data = data_list, chains = 4, parallel_chains = 4, adapt_delta = 0.99, max_treedepth = 12)
fit
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

p = ggplot(data, aes(ColPos, log(concentration), color = as.factor(Pooled))) +
  geom_point() +
  geom_line(data = dplyr::filter(data, Pooled == 1)) +
  geom_point(data = data.frame(deltas = colMeans(fit$draws('deltas', f = "draws_matrix")),
                               ColPos = data$ColPos), 
  aes(y = deltas), color  = 3) +
  theme_cowplot() + theme(legend.position = "none", 
                          plot.background = element_rect(fill = "white"),
                          panel.background = element_rect(fill = "white"))
save_plot("test.png", p, base_height = 8)

png("test.png", height = 1000, width = 1000)
colMeans(fit$draws("level", f = "draws_matrix")) |> hist()
dev.off()

