library(OptimalRerandExpDesigns)

set.seed(1984)
n = 100
p = 10
X = matrix(rnorm(n * p), nrow = n, ncol = p)
X = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})

W_base_obj = generate_W_base_and_sort(X)
W_base_obj
#Fig 1
plot(W_base_obj, title = "", xlab = "Mahalanobis Distance")


#Fig 2

rerand_res_ex = optimal_rerandomization_exact(W_base_obj, z_sim_fun = function(){rnorm(n)})
rerand_res_ex
plot(rerand_res_ex, s_min = 5)

library(rmutil)
rerand_res_ex_laplace = optimal_rerandomization_exact(W_base_obj, z_sim_fun = function(){rlaplace(n)})
rerand_res_ex_laplace
plot(rerand_res_ex_laplace, s_min = 5)

rerand_res_ex_tdist = optimal_rerandomization_exact(W_base_obj, z_sim_fun = function(){rt(n, df = 6)})
rerand_res_ex_tdist
plot(rerand_res_ex_tdist, s_min = 5)

rerand_res_norm = optimal_rerandomization_normality_assumed(W_base_obj)
rerand_res_norm
plot(rerand_res_norm, s_min = 5)


rerand_res_approx = optimal_rerandomization_tail_approx(W_base_obj)
rerand_res_approx
plot(rerand_res_approx, s_min = 5)
head(rerand_res_approx$all_data_from_run, 15)

rerand_res_approx_kurt_3 = optimal_rerandomization_tail_approx(W_base_obj, excess_kurtosis_z = 3)
rerand_res_approx_kurt_3
plot(rerand_res_approx_kurt_3, s_min = 5)
head(rerand_res_approx_kurt_3$all_data_from_run, 15)

max_designs = W_base_obj$max_designs
imbalances = rerand_res_ex$all_data_from_run$imbalance_by_w_sorted

data_rerand_res_ex = rerand_res_ex$all_data_from_run$Q_primes
data_rerand_res_ex_laplace = rerand_res_ex_laplace$all_data_from_run$Q_primes
data_rerand_res_ex_tdist = rerand_res_ex_tdist$all_data_from_run$Q_primes
data_rerand_res_norm = rerand_res_norm$all_data_from_run$Q_primes
data_rerand_res_approx = rerand_res_approx$all_data_from_run$Q_primes
data_rerand_res_approx_kurt_3 = rerand_res_approx_kurt_3$all_data_from_run$Q_primes

data_rerand_res_ex = data_rerand_res_ex / data_rerand_res_ex[max_designs]
data_rerand_res_ex_laplace = data_rerand_res_ex_laplace / data_rerand_res_ex_laplace[max_designs]
data_rerand_res_ex_tdist = data_rerand_res_ex_tdist / data_rerand_res_ex_tdist[max_designs]
data_rerand_res_norm = data_rerand_res_norm / data_rerand_res_norm[max_designs]
data_rerand_res_approx = data_rerand_res_approx / data_rerand_res_approx[max_designs]
data_rerand_res_approx_kurt_3 = data_rerand_res_approx_kurt_3 / data_rerand_res_approx_kurt_3[max_designs]


ggplot(data.frame(
  data_rerand_res_ex = data_rerand_res_ex,
  data_rerand_res_ex_laplace = data_rerand_res_ex_laplace,
  data_rerand_res_ex_tdist = data_rerand_res_ex_tdist,
  data_rerand_res_norm = data_rerand_res_norm,
  data_rerand_res_approx = data_rerand_res_approx,
  data_rerand_res_approx_kurt_3 = data_rerand_res_approx_kurt_3
)) +
  geom_line(aes(x = imbalances, y = data_rerand_res_ex), color = "green") + 
  geom_line(aes(x = imbalances, y = data_rerand_res_ex_laplace), color = "purple") + 
  geom_line(aes(x = imbalances, y = data_rerand_res_ex_tdist), color = "black") + 
  geom_line(aes(x = imbalances, y = data_rerand_res_norm), color = "red") + 
  geom_line(aes(x = imbalances, y = data_rerand_res_approx), color = "blue") + 
  geom_line(aes(x = imbalances, y = data_rerand_res_approx_kurt_3), color = "orange") + 
  ylim(0, 2) + xlim(0.04, max(imbalances))

# smoothing_degree = 1
# smoothing_span = 0.1
# smoothing_fit = loess(rerand_res_ex$all_data_from_run$Q_primes ~ imbalances, degree = smoothing_degree, span = smoothing_span)
# Q_primes_smoothed = predict(smoothing_fit, data.frame(X = imbalances))
# imbalances[which.min(Q_primes_smoothed)]
# rerand_res_ex$a_star
# 
# ggplot(data.frame(
#   imbalances = imbalances,
#   data_rerand_res_ex = data_rerand_res_ex,
#   data_rerand_res_ex_sm = Q_primes_smoothed / Q_primes_smoothed[max_designs]
# )) +
#   geom_line(aes(x = imbalances, y = data_rerand_res_ex), color = "green") + 
#   geom_line(aes(x = imbalances, y = data_rerand_res_ex_sm), color = "purple") + 
#   ylim(0, 2) + xlim(0.04, max(imbalances))


#Fig 3
set.seed(1984)
bbeta = rnorm(p)
betaT = 1
N_w = 1000
N_z = 1000
sigma_z = 3

mse_z_all = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, N_w)
  for (n_w in 1 : N_w){
    w = W_base_obj$W_base_sorted[sample(1 : 25000, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_all[n_z] = mean(mse_w)
}


mse_z_top_5000 = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, N_w)
  for (n_w in 1 : N_w){
    w = W_base_obj$W_base_sorted[sample(1 : 5000, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_top_5000[n_z] = mean(mse_w)
}


mse_z_worst_5000 = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, N_w)
  for (n_w in 1 : N_w){
    w = W_base_obj$W_base_sorted[sample(20001 : 25000, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_worst_5000[n_z] = mean(mse_w)
}


mse_z_opt = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, N_w)
  for (n_w in 1 : N_w){
    w = rerand_res_norm$W_star[sample(1 : rerand_res_norm$w_star_size, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_opt[n_z] = mean(mse_w)
}


mse_z_determ = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  w = W_base_obj$W_base_sorted[1, ]
  y = X %*% bbeta + z + betaT * w
  betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
  mse_z_determ[n_z] = (betaThatLinRegr - betaT)^2
}





ggplot(x$all_data_from_run) +
  ggtitle(title, subtitle = subtitle) +
  xlab(xlab) +
  ylab(ylab) +
  geom_line(aes(x = imbalance_by_w_sorted, y = Q_primes)) +
  geom_vline(xintercept = x$a_star, col = "green") +
  scale_x_log10(
    # breaks = scales::trans_breaks("log10", function(x) 10^x, n = 10),
    # labels = scales::trans_format("log10", scales::scientific_format(digits = 2))
    ) +
  scale_y_log10()
  # coord_trans(y = "log10")










max_designs = 25000


c_vals = seq(1, 10, by = 0.25)

res = list()

for (c_val in c_vals){
  res[[as.character(c_val)]] = optimal_rerandomization(X, max_designs, c_val = c_val)
}

#save(res, file = "res_thresholds.RData")
load("res_thresholds.RData")
unlist(lapply(res, function(l){l$Q_star}))
unlist(lapply(res, function(l){l$a_star}))
unlist(lapply(res, function(l){l$W_star_size}))



# Frob norm sqd debiasing works?
library(OptimalRerandExpDesigns)

set.seed(1984)
Nsim = 1000
n = 100
S = 1000
frob_norm_sqs_naive = array(NA, Nsim)
frob_norm_sqs_debiased = array(NA, Nsim)
for (nsim in 1 : Nsim){
  ws = complete_randomization_with_forced_balance_plus_one_min_one(n, S)
  Sigmahat = t(ws) %*% ws / S
  #naive
  frob_norm_sqs_naive[nsim] = frob_norm_sq(Sigmahat)
  frob_norm_sqs_debiased[nsim] = frob_norm_sq_debiased_times_matrix(Sigmahat, diag(n), S, n)
}
mean(frob_norm_sqs_naive)
mean(frob_norm_sqs_debiased)
sd(frob_norm_sqs_debiased) / sqrt(Nsim)


#Clinical Trial Data
library(PTE)


