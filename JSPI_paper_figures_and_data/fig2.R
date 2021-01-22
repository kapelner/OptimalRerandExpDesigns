source("fig1.R")
pacman::p_load(rmutil, data.table)

set.seed(10)
skip_search_length = 10

#Fig 2
rerand_res_ex = optimal_rerandomization_exact(W_base_obj, z_sim_fun = function(){rnorm(n)}, skip_search_length = skip_search_length)
rerand_res_ex
plot(rerand_res_ex, s_min = 5)

rerand_res_ex_laplace = optimal_rerandomization_exact(W_base_obj, z_sim_fun = function(){rlaplace(n)}, skip_search_length = skip_search_length)
rerand_res_ex_laplace
# plot(rerand_res_ex_laplace, s_min = 5)

rerand_res_ex_tdist = optimal_rerandomization_exact(W_base_obj, z_sim_fun = function(){rt(n, df = 6)}, skip_search_length = skip_search_length)
rerand_res_ex_tdist
# plot(rerand_res_ex_tdist, s_min = 5)

rerand_res_norm = optimal_rerandomization_normality_assumed(W_base_obj, skip_search_length = skip_search_length)
rerand_res_norm
# plot(rerand_res_norm, s_min = 5)


rerand_res_approx = optimal_rerandomization_tail_approx(W_base_obj, skip_search_length = skip_search_length, use_frob_norm_sq_unbiased_estimator = FALSE)
rerand_res_approx
# plot(rerand_res_approx, s_min = 5)
# head(rerand_res_approx$all_data_from_run, 15)

rerand_res_approx_kurt_3 = optimal_rerandomization_tail_approx(W_base_obj, excess_kurtosis_z = 3, skip_search_length = skip_search_length, use_frob_norm_sq_unbiased_estimator = FALSE)
rerand_res_approx_kurt_3
# plot(rerand_res_approx_kurt_3, s_min = 5)
# head(rerand_res_approx_kurt_3$all_data_from_run, 15)

data_rerand_res_ex = rerand_res_ex$all_data_from_run$Q_primes
data_rerand_res_ex_laplace = rerand_res_ex_laplace$all_data_from_run$Q_primes
data_rerand_res_ex_tdist = rerand_res_ex_tdist$all_data_from_run$Q_primes
data_rerand_res_norm = rerand_res_norm$all_data_from_run$Q_primes
data_rerand_res_approx = rerand_res_approx$all_data_from_run$Q_primes
data_rerand_res_approx_kurt_3 = rerand_res_approx_kurt_3$all_data_from_run$Q_primes

q_pin = 0.95
data_rerand_res_ex = data_rerand_res_ex / quantile(data_rerand_res_ex, q_pin, na.rm = TRUE)
data_rerand_res_ex_laplace = data_rerand_res_ex_laplace / quantile(data_rerand_res_ex_laplace, q_pin, na.rm = TRUE)
data_rerand_res_ex_tdist = data_rerand_res_ex_tdist / quantile(data_rerand_res_ex_tdist, q_pin, na.rm = TRUE)
data_rerand_res_norm = data_rerand_res_norm / quantile(data_rerand_res_norm, q_pin, na.rm = TRUE)
data_rerand_res_approx = data_rerand_res_approx / quantile(data_rerand_res_approx, q_pin, na.rm = TRUE)
data_rerand_res_approx_kurt_3 = data_rerand_res_approx_kurt_3 / quantile(data_rerand_res_approx_kurt_3, q_pin, na.rm = TRUE)

imbalance = rerand_res_ex$all_data_from_run$imbalance_by_w_sorted
design_names = c("OPT-chisq", "OPT-tail-kappa0", "OPT-exact-Normal", "OPT-exact-Laplace", "OPT-exact-T6", "OPT-tail-kappa3")
res = rbind(
  data.table(imbalance = imbalance, relativeQ = data_rerand_res_norm,            design = design_names[1]),
  data.table(imbalance = imbalance, relativeQ = data_rerand_res_approx,          design = design_names[2]),
  data.table(imbalance = imbalance, relativeQ = data_rerand_res_ex,              design = design_names[3]),
  data.table(imbalance = imbalance, relativeQ = data_rerand_res_ex_laplace,      design = design_names[4]),
  data.table(imbalance = imbalance, relativeQ = data_rerand_res_ex_tdist,        design = design_names[5]),
  data.table(imbalance = imbalance, relativeQ = data_rerand_res_approx_kurt_3,   design = design_names[6])
)
a_stars = data.table(
  imbalance = c(
    rerand_res_norm$a_star, 
    rerand_res_approx$a_star, 
    rerand_res_ex$a_star_smoothed, 
    rerand_res_ex_laplace$a_star_smoothed, 
    rerand_res_ex_tdist$a_star_smoothed, 
    rerand_res_approx_kurt_3$a_star
  ),
  design = design_names
)
manual_colors = c(
  "blue", "cyan", "red", "yellow", "purple", "darkgreen"
)
ggplot(na.omit(res)) +
  geom_line(aes(x = imbalance, y = relativeQ, color = design)) +
  scale_color_manual(values = manual_colors) +
  geom_vline(data = a_stars,
             aes(xintercept = imbalance, color = design),
             lwd = 1, alpha = 1, linetype = "solid") +
  xlab("imbalance threshold (a)") +
  ylab("relative tail criterion value (Q')") +
  scale_x_log10(limits = c(0.05, 0.8)) +
  scale_y_log10(breaks = NULL)
