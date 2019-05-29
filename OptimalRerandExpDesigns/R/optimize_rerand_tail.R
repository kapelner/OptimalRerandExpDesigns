

#' Find the Optimal Rerandomization Design Under the Tail and Kurtosis Approximation
#' 
#' Finds the optimal rerandomization threshold based on a user-defined quantile
#' and kurtosis based on an approximation of tail standard errors
#' 
#' @param W_base_object			An object that contains the assignments to begin with sorted by imbalance.
#' @param estimator 			"linear" for the covariate-adjusted linear regression estimator (default).
#' @param q 					The tail criterion's quantile of MSE over z's. The default is 95\%. 
#' @param skip_search_length	In the exhaustive search, how many designs are skipped? Default is 1 for 
#' 								full exhaustive search through all assignments provided for in \code{W_base_object}.
#' @param binary_search			If \code{TRUE}, a binary search is employed to find the optimal threshold instead of 
#' 								an exhaustive search. Default is \code{FALSE}.
#' @param excess_kurtosis_z		An estimate of the excess kurtosis in the measure on z. Default is 0.
#' @param frob_norm_sq_bias_correction_min_samples  This function requires computation of a Frobenius norm
#'                                          squared of a matrix. We employ a bias corrected estimate. 
#' However, this estimate suffers from high variance when there are not enough samples. Thus, we only implement
#' the correction beginning at this number of vectors. Default is 10.
#' @param dot_every_x_iters		Print out a dot every this many iterations. The default is 100. Set to
#' 								\code{NULL} for no printout.
#' @return 						A list containing the optimal design threshold, strategy, and
#' 								other information.
#' 
#' @author Adam Kapelner
#' @export
optimal_rerandomization_tail_approx = function(
		W_base_object,
		estimator = "linear",
		q = 0.95,
		skip_search_length = 1,
		binary_search = FALSE,
		excess_kurtosis_z = 0,
		frob_norm_sq_bias_correction_min_samples = 10,
		dot_every_x_iters = 100){
	
	optimal_rerandomization_argument_checks(W_base_object, estimator, q)
	
	c_val = qnorm(q)
	n = W_base_object$n
	X = W_base_object$X
	W_base_sort = W_base_object$W_base_sort
	max_designs = W_base_object$max_designs
	imbalance_by_w_sorted = W_base_object$imbalance_by_w_sorted
	
	if (estimator == "linear"){
	  Xt = t(X)
	  XtXinv = solve(Xt %*% X)
	  XtXinv_eigen <- eigen(XtXinv)
	  XtXinv_sqrt <- XtXinv_eigen$vectors %*% diag(sqrt(XtXinv_eigen$values)) %*% solve(XtXinv_eigen$vectors)
	  X_orth = X %*% XtXinv_sqrt
	  X_orth_T = t(X_orth)
	  P = X %*% XtXinv %*% Xt
	  I = diag(n)
	  I_min_P = I - P
	}


	Q_primes = array(NA, max_designs)
	#diagram all terms
	balances = array(NA, max_designs)
	frob_norm_sqs_debiased = array(NA, max_designs)
	tr_gds = array(NA, max_designs)
	tr_d_sqs = array(NA, max_designs)
	r_i_sqs = array(NA, max_designs)
	
	w_w_T_running_sum = matrix(0, n, n)
	if (estimator == "linear"){
	  w_w_T_P_w_w_T_running_sum = matrix(0, n, n)
	}
	for (s in seq(from = 1, to = max_designs, by = skip_search_length)){
		if (!is.null(dot_every_x_iters)){
			if (s %% dot_every_x_iters == 0){
				cat(".")
			}
		}

	  w_s = W_base_sort[s, , drop = FALSE]
	  w_s_w_s_T = t(w_s) %*% w_s
	  w_w_T_running_sum = w_w_T_running_sum + w_s_w_s_T
	  Sigma_W = 1 / (s / skip_search_length) * w_w_T_running_sum
	  if (estimator == "linear"){	    
	    w_w_T_P_w_w_T_running_sum = w_w_T_P_w_w_T_running_sum + w_s_w_s_T %*% P %*% w_s_w_s_T
	    D = 1 / (s / skip_search_length) * w_w_T_P_w_w_T_running_sum
	    G = I_min_P %*% Sigma_W %*% I_min_P
	    balances[s] = tr(X_orth_T %*% Sigma_W %*% X_orth)
	    frob_norm_sqs_debiased[s] = #see Appendix 5.9	 
	      frob_norm_sq_debiased(Sigma_W, s / skip_search_length, n, frob_norm_sq_bias_correction_min_samples)	-
	      frob_norm_sq_debiased_times_matrix(Sigma_W, P, s / skip_search_length, n, frob_norm_sq_bias_correction_min_samples)
	    tr_gds[s] = tr(G %*% D) / n
	    tr_d_sqs[s] = tr(D^2) / n^2
	    r_i_sqs[s] = sum((diag(G) + 2 * diag(D) / n)^2)
	    Q_primes[s] = balances[s] +
	      c_val * sqrt(
	        2 * frob_norm_sqs_debiased[s] + 
	          8 * tr_gds[s] +
	          8 * tr_d_sqs[s] + 
			  excess_kurtosis_z * r_i_sqs[s]
	      )
	  
	  }
	}
	cat("\n")
	
	
	#now find it by starting at the top and going down until it appears
	
	s_star = NULL
	Q_star = Inf
	for (s in seq(from = max_designs - skip_search_length, to = 1, by = -skip_search_length)){
	  if (Q_primes[s] > Q_primes[s + skip_search_length]){
	    Q_star = Q_primes[s + skip_search_length]
	    s_star = s + skip_search_length
	    break
	  }
	}


	all_data_from_run = data.frame(
	  imbalance_by_w_sorted = imbalance_by_w_sorted, 
	  Q_primes = Q_primes,
	  balances = balances,
	  frob_norm_sqs_debiased = frob_norm_sqs_debiased,
	  tr_gds = tr_gds,
	  tr_d_sqs = tr_d_sqs,
	  r_i_sqs = r_i_sqs
	)

	ll = list(
		type = "approx",
		estimator = estimator,
		q = q,
		W_base_object = W_base_object,
		W_star = W_base_sort[1 : s_star, ],
		W_star_size = s_star,
		a_star = imbalance_by_w_sorted[s_star],
		a_stars = imbalance_by_w_sorted[1 : s_star],
		all_data_from_run = all_data_from_run,
		Q_star = Q_star
	)
	class(ll) = "optimal_rerandomization_obj"
	ll
}

