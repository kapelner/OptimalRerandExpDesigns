if (!require("pacman")) install.packages("pacman")
pacman::p_load(OptimalRerandExpDesigns)

#### Fig 4
set.seed(10)

n = 200
S = 25000
skip_search_length = 10
ps = 1 : (n - 1)
X = matrix(rnorm(n * max(ps)), nrow = n)

res = matrix(NA, nrow = length(ps), ncol = 3)
colnames(res) = c("p", "a_star", "fraction")

for (i in 1 : length(ps)){
  p = ps[i]
  cat("p = ", p, "\n")
  Xi = X[, 1 : p, drop = FALSE]
  W_base_obj = generate_W_base_and_sort(Xi, max_designs = S)
  design = optimal_rerandomization_normality_assumed(W_base_obj, skip_search_length = skip_search_length)
  res[i, ] = c(p, design$a_star, design$W_star_size / S)
  plot(ggplot(data.frame(res)) + geom_line(aes(x = p, y = fraction)))
}
