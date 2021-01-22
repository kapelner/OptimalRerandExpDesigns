if (!require("pacman")) install.packages("pacman")
pacman::p_load(OptimalRerandExpDesigns)
#R CMD INSTALL -l ~/Documents/R/win-library/3.6/ OptimalRerandExpDesigns

#Fig 1
set.seed(10)
n = 100
p = 10
X = matrix(rnorm(n * p), nrow = n, ncol = p)
X = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
S = 25000

W_base_obj = generate_W_base_and_sort(X, max_designs = S)
W_base_obj
plot(W_base_obj, title = "", xlab = "scaled Mahalanobis Distance")
