source("fig2.R")
# pacman::p_load(gridExtra, ggExtra)

rerand_res_norm$W_star_size
rerand_res_approx$W_star_size
rerand_res_ex$W_star_size_smoothed
rerand_res_norm$a_star
rerand_res_approx$a_star
rerand_res_ex$a_star_smoothed

#Fig 3
set.seed(10)
bbeta = rnorm(p)
betaT = 1
N_w = 1000
N_z = 5000
S_good = 125
W_base_obj$imbalance_by_w_sorted[S_good]
sigma_z = 3

mse_z_all = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, N_w)
  for (n_w in 1 : N_w){
    w = W_base_obj$W_base_sorted[sample(1 : S, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_all[n_z] = mean(mse_w)
}


mse_z_top = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, min(N_w, S_good))
  for (n_w in 1 : min(N_w, S_good)){
    w = W_base_obj$W_base_sorted[n_w, ]
    # w = W_base_obj$W_base_sorted[sample(1 : S_good, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_top[n_z] = mean(mse_w)
}


mse_z_bad = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, min(N_w, S_good))
  for (n_w in 1 : min(N_w, S_good)){
    # w = rerand_res_norm$W_star[n_w, ]
    w = W_base_obj$W_base_sorted[S - n_w + 1, ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_bad[n_z] = mean(mse_w)
}


mse_z_opt1 = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, N_w)
  for (n_w in 1 : N_w){
    w = rerand_res_norm$W_star[sample(1 : rerand_res_norm$W_star_size, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_opt1[n_z] = mean(mse_w)
}

mse_z_opt2 = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, N_w)
  for (n_w in 1 : N_w){
    w = rerand_res_approx$W_star[sample(1 : rerand_res_approx$W_star_size, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_opt2[n_z] = mean(mse_w)
}

mse_z_opt3 = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, N_w)
  for (n_w in 1 : N_w){
    w = rerand_res_ex$W_star[sample(1 : rerand_res_ex$W_star_size, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_opt3[n_z] = mean(mse_w)
}

w_opt = W_base_obj$W_base_sorted[1, ]
mse_z_determ = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  y = X %*% bbeta + z + betaT * w_opt
  betaThatLinRegr = coef(lm(y ~ 0 + X + w_opt))[p + 1]
  mse_z_determ[n_z] = (betaThatLinRegr - betaT)^2
}


all_res = rbind(
  data.table(sqerr = mse_z_all, design = "BCRD"),
  data.table(sqerr = mse_z_top, design = "GOOD"),
  data.table(sqerr = mse_z_bad, design = "BAD"),
  data.table(sqerr = mse_z_opt1, design = "OPT-chisq"),
  data.table(sqerr = mse_z_opt2, design = "OPT-tail"),
  data.table(sqerr = mse_z_opt3, design = "OPT-exact"),
  data.table(sqerr = mse_z_determ, design = "DET")
)
all_res_summary = all_res[, .(avg = mean(sqerr), q = quantile(sqerr, .95)), by = design]
xL = quantile(mse_z_all, 0.08)
xU = quantile(mse_z_determ, 0.96)
manual_colors = c(
  "darkgrey", "black", "yellow3", "green", "blue", "red", "purple"
)
ggplot(all_res) +
  xlab("OLS Estimator Squared Error over the Z distribution") +
  scale_color_manual(values = alpha(manual_colors, 0)) +
  scale_fill_manual(values = manual_colors) +
  geom_density(aes(x = sqerr, fill = design), alpha = 0.3) +
  geom_vline(data = all_res_summary,
             aes(xintercept = avg, color = design),
             lwd = 1, alpha = 1, linetype = "solid") +
  geom_vline(data = all_res_summary,
             aes(xintercept = q, color = design),
             lwd = 1, alpha = 1, linetype = "dashed") +
  guides(color = FALSE, design = TRUE, fill = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_log10(limits = c(xL, xU))

all_res = rbind(
  data.table(sqerr = mse_z_opt1, design = "OPT-chisq"),
  data.table(sqerr = mse_z_opt2, design = "OPT-tail"),
  data.table(sqerr = mse_z_opt3, design = "OPT-exact")
)
all_res_summary = all_res[, .(avg = mean(sqerr), q = quantile(sqerr, .95)), by = design]
xL = quantile(mse_z_opt1, 0.50)
xU = quantile(mse_z_opt1, 0.96)
manual_colors = c(
  "blue", "red", "purple"
)
ggplot(all_res) +
  xlab("OLS Estimator Squared Error over the Z distribution") +
  scale_color_manual(values = alpha(manual_colors, 0)) +
  scale_fill_manual(values = manual_colors) +
  geom_density(aes(x = sqerr, fill = design), alpha = 0.3) +
  geom_vline(data = all_res_summary,
             aes(xintercept = avg, color = design),
             lwd = 1, alpha = 1, linetype = "solid") +
  geom_vline(data = all_res_summary,
             aes(xintercept = q, color = design),
             lwd = 1, alpha = 1, linetype = "dashed") +
  guides(color = FALSE, design = TRUE, fill = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_log10(limits = c(xL, xU))