if (!require("pacman")) install.packages("pacman")
pacman::p_load(OptimalRerandExpDesigns, rmutil, tidyverse, magrittr, gridExtra, ggExtra, data.table, YARF)

#### Fig 5
set.seed(10)

#load the data
X = read.csv("JSPI_paper_figures_and_data/star.csv") %>%
  select(
    gender, 
    race, 
    contains("birth"),
    contains("g3")
  ) %>%
  select(-c(g3schid, g3tchid)) %>%
  filter(g3classtype != 3) %>%
  mutate(year = birthyear + birthmonth / 12 + birthday / 365.25) %>%
  select(-contains("birth")) %>%
  mutate(g3classtype = ifelse(g3classtype == 1, 1, 0)) %>%
  mutate(g3freelunch = ifelse(g3freelunch == 2, 1, 0)) %>%
  mutate(g3surban4 = ifelse(g3surban == 4, 1, 0)) %>%
  mutate(g3surban3 = ifelse(g3surban == 3, 1, 0)) %>%
  mutate(g3surban2 = ifelse(g3surban == 2, 1, 0)) %>%
  select(-g3surban) %>%
  rename(trt = g3classtype) %>%
  rowwise %>%
  mutate(y = mean(c(g3treadss, g3tmathss, g3tlangss, g3tlistss), na.rm = TRUE))  %>%
  select(-c(g3treadss, g3tmathss, g3tlangss, g3tlistss)) %>%
  na.omit %>%
  data.table

summary(X)

#measure the ATE and treat as real
betaT = mean(X[trt == 1, y]) - mean(X[trt == 0, y])

#copy the dataset and delete the effect of the ATE
Xcopy = copy(X)
Xcopy[trt == 1, y := y - betaT]
#fit RF model and get back z's
yarf_mod = YARF(Xcopy[, !"y"], Xcopy$y, seed = 1984)
z = yarf_mod$residuals
Xcopy[, y := NULL]
Xcopy[, fhat := yarf_mod$y_oob]
rm(yarf_mod, X); gc()


n = 100
p = 8
Nsim = 1
N_z = 1000
N_w = 1000
S = 25000
skip_search_length = 10

res = data.table(
  nsim = integer(),
  mse = numeric(),
  design = character()
)

for (nsim in 1 : Nsim){
  mses_bcrd = array(NA, N_z)
  mses_opt_chisq = array(NA, N_z) 
  mses_opt_tail = array(NA, N_z)
  mses_det = array(NA, N_z) 
  
  #sample a dataset
  X_0 = rbind(
    sample_n(Xcopy[trt == 1], n / 2),
    sample_n(Xcopy[trt == 0], n / 2)
  ) %>% 
    select(-trt)
  fhat_0 = X_0$fhat
  #pull out only the covariates and standardize
  X_0$fhat = NULL
  X_0 = apply(X_0, 2, function(xj){(xj - mean(xj)) / sd(xj)})
  
  #begin with S BCRD vectors
  W_base_obj = generate_W_base_and_sort(X_0, max_designs = S)
  
  #now compute the strategies
  w_opt = W_base_obj$W_base_sorted[1, ]
  design_opt_chisq = optimal_rerandomization_normality_assumed(W_base_obj, skip_search_length = skip_search_length)
  plot(design_opt_chisq, s_min = 5)
  design_opt_tail = optimal_rerandomization_tail_approx(W_base_obj, skip_search_length = skip_search_length)
  plot(design_opt_tail, s_min = 5)
  
  for (n_z in 1 : N_z){
    z_rep = sample(z, n)
    mse_w = array(NA, N_w)
    for (n_w in 1 : N_w){
      w = W_base_obj$W_base_sorted[sample(1 : S, 1), ]
      y_0 = fhat_0 + z_rep + betaT * w
      betaThatLinRegr = coef(lm(y_0 ~ X_0 + w))[p + 2]
      mse_w[n_w] = (betaThatLinRegr - betaT)^2
    }
    mses_bcrd[n_z] = mean(mse_w)
  }
  cat("B ")
  
  for (n_z in 1 : N_z){
    z_rep = sample(z, n)
    mse_w = array(NA, N_w)
    for (n_w in 1 : N_w){
      w = design_opt_chisq$W_star[sample(1 : design_opt_chisq$W_star_size, 1), ]
      y_0 = fhat_0 + z_rep + betaT * w
      betaThatLinRegr = coef(lm(y_0 ~ X_0 + w))[p + 2]
      mse_w[n_w] = (betaThatLinRegr - betaT)^2
    }
    mses_opt_chisq[n_z] = mean(mse_w)
  }
  cat("Ox ")
  
  for (n_z in 1 : N_z){
    z_rep = sample(z, n)
    mse_w = array(NA, N_w)
    for (n_w in 1 : N_w){
      w = design_opt_tail$W_star[sample(1 : design_opt_tail$W_star_size, 1), ]
      y_0 = fhat_0 + z_rep + betaT * w
      betaThatLinRegr = coef(lm(y_0 ~ X_0 + w))[p + 2]
      mse_w[n_w] = (betaThatLinRegr - betaT)^2
    }
    mses_opt_tail[n_z] = mean(mse_w)
  }
  cat("OT ")  
  
  mses_det = array(NA, N_z)
  for (n_z in 1 : N_z){
    z_rep = sample(z, n)
    y_0 = fhat_0 + z_rep + betaT * w_opt
    betaThatLinRegr = coef(lm(y_0 ~ X_0 + w_opt))[p + 2]
    mses_det[n_z] = (betaThatLinRegr - betaT)^2
  }
  cat("D\n")
  
  #tally the results
  res = rbind(res, 
    data.table(nsim = nsim, mse = mses_bcrd, design = "BCRD"),
    data.table(nsim = nsim, mse = mses_opt_chisq, design = "OPT-chisq"),
    data.table(nsim = nsim, mse = mses_opt_tail, design = "OPT-tail"),
    data.table(nsim = nsim, mse = mses_det, design = "DET")
  )
  res_summary = res[, .(avg = mean(mse), q = quantile(mse, .95)), by = design]
  xL = quantile(mses_bcrd, 0.08)
  xU = quantile(mses_det, 0.96)
  manual_colors = c(
    "black", "yellow3", "blue", "purple"
  )
  plot(ggplot(res) +
   xlab("OLS Estimator Squared Error over the Z distribution") +
   scale_color_manual(values = alpha(manual_colors, 0)) +
   scale_fill_manual(values = manual_colors) +
   geom_density(aes(x = mse, fill = design), alpha = 0.3) +
   geom_vline(data = res_summary,
              aes(xintercept = avg, color = design),
              lwd = 1, alpha = 1, linetype = "solid") +
   geom_vline(data = res_summary,
              aes(xintercept = q, color = design),
              lwd = 1, alpha = 1, linetype = "dashed") +
   guides(color = FALSE, design = TRUE, fill = guide_legend(override.aes = list(alpha = 1))) +
   scale_x_log10(limits = c(xL, xU)))
}
