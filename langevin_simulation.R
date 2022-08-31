
setwd("~/Code/nested_rhat/summer_2022")
.libPaths("~/Rlib/")
library(ggplot2)
library(rootSolve)
library(latex2exp)

rho_t <- function (t) {
  (1 - exp(-t)) / t
}

xi_t <- function (t) {
  (1 - exp(- 2 * t)) / (2 * t)
}

eta_t <- function (t) {
  2 * (1 - rho_t(t)) / t
}

# N > 1 case for $\hat R$.
B <- function (sigma2, sigma2_0, t) {
  (sigma2_0 - sigma2) * rho_t(t)^2 + sigma2 * eta_t(t)
}

W <- function (sigma2, sigma2_0, mu0, t) {
  
  (sigma2_0 - sigma2 + mu0^2) * (xi_t(t) - rho_t(t)^2) + sigma2 * (1 - eta_t(t))
}


sigma2 <- 1
sigma2_0 <- 0
mu0 <- 5
t_array <- c(seq(from = 0.00001, to = 1, by = 0.01), 1:100)
B_array <- rep(NA, length(t_array))
W_array <- rep(NA, length(t_array))

for (i in 1:length(t_array)) {
  B_array[i] <- B(sigma2, sigma2_0, t_array[i]) 
  W_array[i] <- W(sigma2, sigma2_0, mu0, t_array[i])
}

plot_data <- data.frame(y = c(B_array, W_array, B_array / W_array),
                        time = rep(t_array, 3),
                        label = rep(c("B", "W", "B / W"), each = length(t_array)))

p <- ggplot(data = plot_data, aes(x = time, y = y, color = label)) +
  geom_line() + scale_y_continuous(trans='log') +
  theme_bw() 
p


f <- function (sigma2, sigma2_0, mu0, t) {
  B(sigma2, sigma2_0, t) / (sigma2 + (sigma2_0 - sigma2) * xi_t(t) + mu0^2 * (xi_t(t) - rho_t(t)^2))
}

f_array <- rep(NA, length(t_array))
for (i in 1:length(t_array)) {
  f_array[i] <- f(sigma2, sigma2_0, mu0, t_array[i])
}

p <- ggplot(data = data.frame(t = t_array, f = f_array), aes(x = t, y = f)) +
  geom_line() + theme_bw()
p


###############################################################################
## Condition for reliability

sigma2_0_lower <- function (delta, delta_prime, M, mu_0, sigma2 = 1) {
  # T_star = max(0, 0.5 * log (mu_0^2 / (delta_prime * sigma2)))
  # (delta - 1 / M) * (M / (2 * M - 1)) * (exp(2 * T_star - 1)) * sigma2
  (delta - 1/ M) * (M / (2 * M - 1)) * (mu_0^2 / (delta_prime * sigma2) - 1) * sigma2
}

M = 32
delta = 1 / M + 0.01
delta_prime = 0.01
mu_0 = seq(from = 0, to = 100, by = 0.1)
sigma2_0 = rep(NA, length(mu_0))

for (i in 1:length(sigma2_0)) {
  sigma2_0[i] = sigma2_0_lower(delta, delta_prime, M, mu_0[i])
}

p <- ggplot(data = data.frame(bias = mu_0, variance = sigma2_0), 
            aes(x = bias, y = variance)) +
  geom_line() +
  theme_bw()
p

###############################################################################
## Upper bound on delta

coth <- function(t) {
  (exp(2 * t) + 1) / (exp(2 * t) - 1)
}

delta_bound <- function(t) {
  1 / (0.5 * t * coth(t / 2) - 1)
}

bound_search <- function(delta_prime, mu0, sigma2 = 1) {
  # 1) find t*
  fun <- function (t) {
    mu0 * rho_t(t)^2 / sigma2 - delta_prime
  }

  t_star = uniroot(fun, c(0.001, 1000))$root

  # 2) find bound
  delta_bound(t_star)
}

mu0 = seq(from = 1, to = 20, by = 1)
delta_prime = seq(from = 0.001, to = 0.99, by = 0.005)

bound_mat = array(NA, dim = c(length(delta_prime), length(mu0)))
for (i in 1:length(mu0)) {
  for (j in 1:length(delta_prime)) {
    bound_mat[j, i] = bound_search(delta_prime[j], mu0[i]) 
  }
}

bound_vec = c(bound_mat)

plot_data = data.frame(bound = bound_vec,
                       mu0 = as.factor(rep(mu0, each = length(delta_prime))),
                       delta_prime = rep(delta_prime, length(mu0)))
plot_data$bound_norm <- plot_data$bound / plot_data$delta_prime

p <- ggplot(data = plot_data, aes(x = delta_prime, y = bound, color = mu0)) +
  geom_line() + theme_bw() + scale_y_continuous(trans='log10') +
  # geom_line(aes(x = delta_prime, y = delta_prime, color = 'delta'), # color = "black",
  #          linetype = "dashed") +
  expand_limits(x = 0, y = 1) +
  ylab(TeX("$\\delta$")) + xlab(TeX("$\\delta'$")) +
  # scale_color_discrete(labels = c("1", "4", "9", "16", "25") ) +
  labs(col = TeX("$(\\mu - \\mu_0)^2$")) +
  # guides(title = TeX("$(\\mu - \\mu_0)^2$")) +
  theme(text = element_text(size = 20))
p

###############################################################################
## Lower bound on sigma

sigma0_bound <- function (delta, delta_prime, mu0, sigma2 = 1) {
  # 1) find t*
  fun <- function (t) {
    mu0 * rho_t(t)^2 / sigma2 - delta_prime
  }
  
  t_star = uniroot(fun, c(0.001, 1000))$root
  
  # 2) compute lower bound
  A = delta * (xi_t(t_star) - rho_t(t_star)^2) * mu0^2
  offset = eta_t(t_star) - rho_t(t_star)^2
  B = (delta * (1 + rho_t(t_star)^2 - eta_t(t_star) - xi_t(t_star)) 
         - offset) * sigma2
  C = (1 + delta) * rho_t(t_star)^2 - delta * xi_t(t_star)
  
  list(bound = (A + B) / C, offset = - offset * sigma2 / C)
}

delta = 0.01
delta_prime = 0.005
mu0 = seq(from = 0.05, to = 20, by = 0.01)
bound_vec = rep(NA, length(mu0))
offset_vec = rep(NA, length(mu0))

for (i in 1:length(mu0)) {
  temp_list = sigma0_bound(delta, delta_prime, mu0[i])
  bound_vec[i] = temp_list$bound
  offset_vec[i] = temp_list$offset
}

p <- ggplot(data = data.frame(mu0 = mu0, sigma0 = bound_vec, 
                              offset = offset_vec),
            aes(x = mu0, y = bound_vec)) +
  geom_line() + theme_bw() +
  geom_line(aes(x = mu0, y = offset_vec), linetype = "dashed") +
  geom_line(aes(x = mu0, y = bound_vec - offset_vec), linetype = "dotted")
p





