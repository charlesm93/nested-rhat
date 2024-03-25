
rm(list = ls())
gc()
set.seed(1954)

.libPaths("~/Rlib/")
library(ggplot2)
library(ggrepel)
library(scales)
library(latex2exp)

#####################################################################
## utility functions

covariance_L <- function(N, rho, sigma2 = 1, spatial = TRUE) {
  # Returns Cholesky decomposition of the covariance.
  # If dist = TRUE then there is a spacial component in the covariance;
  # else the off-diagonal terms are all equal.
  Sigma = matrix(NA, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:(i - 1)) {
      if (spatial) {
        correlation = rho^(abs(i - j))
      } else {
        correlation = rho
      }
      Sigma[i, j] = correlation * sigma2
      Sigma[j, i] = Sigma[i, j]
    }
  }
  for (i in 1:N) Sigma[i, i] = sigma2

  chol(Sigma)
}


multi_normal <- function(K, M, L) {
  N = nrow(L)
  theta = array(NA, dim = c(K, M, N))
  for (k in 1:K) {
    for (m in 1:M) {
      theta[k, m, ] = L %*% rnorm(N)
    }
  }
  theta
}

B_W <- function(K, M, N, theta) {
  # NOTE: code for this in Python is much more elegant...
  # K: number of super-chains
  # M: number of chains
  # samples: should be an array with dimension K x M x N
  #
  # Returns: B_hat / W_hat

  theta_... = mean(theta)
  theta_..K = rep(NA, K)
  for (k in 1:K) {
    theta_..K[k] = mean(theta[k, , ])
  }

  B_hat = (1 / (K - 1)) * sum((theta_..K - theta_...)^2)
  
  theta_.MK = array(NA, dim = c(K, M))
  for (k in 1:K) {
    for (m in 1:M) {
      theta_.MK[k, m] = mean(theta[k, m, ])
    }
  }
  
  sum_between_chains = array(NA, c(K, M))
  for (k in 1:K) {
    for (m in 1:M) {
      sum_between_chains[k, m] = (theta_.MK[k, m] - theta_..K[k])^2
    }
  }
  
  sum_within_chains = array(NA, c(K, M, N))
  for (k in 1:K) {
    for (m in 1:M) {
      for (n in 1:N) {
        sum_within_chains[k, m, n] = (theta[k, m, n] - theta_.MK[k, m])^2
      }
    }
  }
  
  W_hat = (1 / K) * (1 / (M - 1)) * sum(sum_between_chains)
                     
  if (N > 1) W_hat = W_hat + (1 / K) * (1 / M) * (1 / (N - 1)) * sum(sum_within_chains)

  B_hat / W_hat
}

# Monte Carlo variance assuming stationarity
mc_variance <- function(n_sim, K, M, N, rho, sigma2 = 1, fractional = FALSE) {
  sim <- rep(NA, n_sim)
  L <- covariance_L(N, rho, sigma2)
  for (n in 1:n_sim) {
    theta = multi_normal(K, M, L)
    sim[n] = B_W(K, M, N, theta)
  }

  if (fractional) {
    mean_B_W = mean(sim)^2
  } else {
    mean_B_W = 1
  }

  var(sim) / mean_B_W
}


# Variance of F distribution (analytical solution)
var_F <- function(K, M, log = FALSE) {
  o = (1 / M^2) *
    2 * K^2 * (M - 1)^2 * (K - 1 + K * (M - 1) - 2) /
    ((K - 1) * (K * (M - 1) - 2)^2 * (K * (M - 1) - 4))

  if (log) {
    mu = K * M
    o = 2 * log(M - 1) - log(mu - M) - 2 * log((mu - 3) * M - mu) - 
      log((mu - 4) * M - mu)
  }
  
  return(o)
}


# Monte Carlo variance assuming initialization from stationarity
# (but the super chains are not yet stationary)
# Code currently assume N = 1.
mc_variance_warm <- function(n_sim, K, M, rho, num_warm, sigma2 = 1, 
                             fractional = TRUE) {
  sim <- rep(NA, n_sim)
  L <- covariance_L(M, rho^(2 * num_warm), sigma2, spatial = FALSE)
  N = 1
  for (n in 1:n_sim) {
    theta = array(NA, dim = c(K, M, 1))
    theta[, , 1] = multi_normal(1, K, L)[1, , ]
    sim[n] = B_W(K, M, N, theta)
  }

  if (fractional) {
    mean_B_W = mean(sim)^2
  } else {
    mean_B_W = 1
  }

  list(variance = var(sim),
       variance_frac = var(sim) / mean_B_W)
}



#####################################################################
## simulations (varying the number of initializations and samples)

mu = 2048  # 124
K = c(2, 4, 8, 16, 32, 64)
N = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
rho = 0.5
n_sim = 1e3

mc_variance_mat = array(NA, dim = c(length(K), length(N)))
for (n in 1:length(N)) {
  for (k in 1:length(K)) {
    print(paste0("k: ", K[k]))
    M = mu / K[k]
    mc_variance_mat[k, n] = mc_variance(n_sim, K[k], M, N[n], rho)
  }
}

mc_variance_vec = c(mc_variance_mat) # / min(mc_variance_vec)

plot_data = data.frame(K = factor(rep(K, length(N)),
                                  levels = c("2", "4", "8", "16", "32", "64")),
                       var = mc_variance_vec,
                       N = rep(N, each = length(K)))

plot_data$labels = ifelse(plot_data$N == 10, paste0("K = ", as.character(plot_data$K)),
                          NA)
# plot_data$labels[plot_data$K == 2] = NA
# plot_data$labelx = NA
# plot_data$labelx[plot_data$N == 8 & plot_data$K == 2] = "K = 2"

# point plot for comparison
p <- ggplot(data = plot_data,
            aes(x = N, y = var, color = K)) +
  # geom_point(size = 1.5) + 
  geom_line() + geom_point() + theme_bw() +
  ylab(" ") + xlab("sampling length") +
  theme(text = element_text(size = 20)) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(breaks=1:10) +
  theme(legend.position = "None") +
  geom_label(aes(label = labels), nudge_x = -0.5, nudge_y = 0.0) # +
  # geom_label(aes(label = labelx), nudge_x = 0, nudge_y = -0.25)

p


#####################################################################
## Varying K, the number of initializations.
# Can reframe the problem to not assume stationarity of the chains 
# (although we do assume the initialization is drawn from the 
# stationary dist). Assume the correlation between two subsequent
# samples is rho.

set.seed(1954)

mu = 2048 # 128
K = c(2, 4, 8, 16, 32, 64)
# K = c(2, 4, 64)
rho = 0.9
n_sim = 1e3
num_warm = 1:50
fractional = FALSE  # if TRUE computed the fractional variance
                    # (i.e. standardized by the squared mean of B / W).

mc_variance_mat = array(NA, dim = c(length(K), length(num_warm), 2))

for (n in 1:length(num_warm)) {
  print(paste0("num_warm: ", num_warm[n]))
  for (k in 1:length(K)) {
    M = mu / K[k]
    list_var = mc_variance_warm(n_sim, K[k], M, rho, num_warm[n],
                                fractional = fractional)
    mc_variance_mat[k, n, 1] = list_var$variance
    mc_variance_mat[k, n, 2] = list_var$variance_frac
  }
}


mc_variance_vec = c(mc_variance_mat)

plot_data = data.frame(K = as.factor(rep(K, 2 * length(num_warm))),
                       var = mc_variance_vec,
                       num_warm = rep(rep(num_warm, each = length(K)), 2),
                       variance = factor(rep(c("variance", "fractional variance"),
                                             each = length(K) * length(num_warm)),
                       levels = c("variance", "fractional variance"))
)

# point plot for comparison
# NOTE: if fractional = FALSE, fractional variance is not computed.
p <- ggplot(data = plot_data,
            aes(x = num_warm, y = var, color = K)) +
  geom_point(size = 0.75, alpha = 0.5) +
  geom_line(linewidth = 0.5, alpha = 1) + theme_bw() + 
  ylab(" ") + xlab("warmup length") +
  theme(text = element_text(size = 20)) +
  scale_y_continuous(trans='log10') +
  facet_wrap(~ variance, scale = "free")
p

# remove the fractional variance which is misleading
plot_data_var = plot_data[plot_data$variance == "variance",]
plot_data_var$label = ifelse(plot_data_var$num_warm == 50, 
                             paste0("K = ", as.character(plot_data_var$K)), NA)
plot_data_var$label[plot_data_var$K == 2] = NA
plot_data_var$labelx = NA
plot_data_var$labelx[plot_data_var$num_warm == 50 & plot_data_var$K == 2] = "K = 2"

p <- ggplot(data = plot_data_var,
            aes(x = num_warm, y = var, color = K)) +
  geom_point(size = 0.75, alpha = 0.5) + 
  geom_line(size = 0.5, alpha = 1) + theme_bw() + 
  ylab(TeX("variance of $n \\widehat{B} / n \\widehat{W}$")) + xlab("warmup length") +
  theme(text = element_text(size = 20)) +
  scale_y_continuous(trans='log10') +
  geom_label(aes(label = label), nudge_x = -2.5) +
  geom_label(aes(label = labelx), nudge_x = -1.9, nudge_y = -0.2) +
  theme(legend.position="none")
p

# can use size = 6 for the labels


#######################################################################
## Varying the number of chains.
# Let's investigate how the variance changes with the number of chains
# We'll focus on the M = mu / 2, which is the maximum number of inits.

Mu = c(4, 32, 128, 516, 1024, 2048)
K = Mu / 2
rho = 0.9
n_sim = 1e3
num_warm = 1:50

mc_variance_mat = array(NA, dim = c(length(Mu), length(num_warm)))
for (n in 1:length(num_warm)) {
  print(paste0("num_warm: ", num_warm[n]))
  for (mu in 1:length(Mu)) {
    M = Mu[mu] / K[mu]  # should always be 2
    mc_variance_mat[mu, n] = mc_variance_warm(n_sim, K[mu], M, rho, num_warm[n])$variance
  }
}

mc_variance_vec = c(mc_variance_mat)

plot_data = data.frame(mu = as.factor(rep(Mu, length(num_warm))),
                       var = mc_variance_vec,
                       num_warm = rep(num_warm, each = length(Mu)))


plot_data_label = plot_data
plot_data_label$label = ifelse(plot_data_label$num_warm == 50, 
                             paste0("KM = ", as.character(plot_data_label$mu)), NA)

# point plot for comparison
p <- ggplot(data = plot_data_label[plot_data$mu != 4, ],
            aes(x = num_warm, y = var, color = mu)) +
  # geom_point(size = 0.75, alpha = 0.5) + 
  geom_line(size = 0.5, alpha = 1) + theme_bw() + 
  # scale_x_continuous(trans='log2') +
  ylab("variance") + xlab("warmup length") +
  theme(text = element_text(size = 20)) +
  theme(legend.position="none") +
  geom_label(aes(label = label), nudge_x = -2, label.size = NA)

p <- p + scale_y_continuous(trans='log10') + 
  geom_hline(yintercept = 1e-5, linetype = "dashed") +
  annotate(geom="text", x = 12, y = 2e-5, label="Asymptotic variance, KM = 128, K = 2")
p
