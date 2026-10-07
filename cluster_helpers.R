as_sigma_vector <- function(sigma2_clust, K) {
  if (length(sigma2_clust) == 1) {
    return(rep(as.numeric(sigma2_clust), K))
  }
  
  if (length(sigma2_clust) != K) {
    stop("sigma_clust must be either a scalar or a vector of length K.")
  }
  
  as.numeric(sigma2_clust)
}

as_dirichlet_vector <- function(v_lambda, K) {
  if (length(v_lambda) == 1) {
    return(rep(as.numeric(v_lambda), K))
  }
  
  if (length(v_lambda) != K) {
    stop("v_lambda must be either a scalar or a vector of length K.")
  }
  
  as.numeric(v_lambda)
}


rdirichlet_simple <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}

rinvgamma <- function(n, shape, rate) {
  if (!is.finite(shape) || shape <= 0) stop("shape must be positive")
  if (!is.finite(rate) || rate <= 0) stop("rate must be positive")
  1 / rgamma(n, shape = shape, rate = rate)
}

cluster_tangent_ssq <- function(Z, mu_g, idx_g) {
  if (length(idx_g) == 0) return(0)
  
  mu_g <- as.numeric(mu_g)
  p <- length(mu_g) - 1
  mu0 <- c(1, rep(0, p))
  
  ssq <- 0
  
  for (i in idx_g) {
    u_i <- inv_exp_map(Z[i, ], mu_g)
    v_i <- parallel_transport(mu_g, mu0, u_i)
    v_hat <- v_i[-1]
    
    ssq <- ssq + sum(v_hat^2)
  }
  
  ssq
}


initial_labels_clust <- function(Z_init, mu_clust_init) {
  Z_init <- as.matrix(Z_init)
  mu_clust_init <- as.matrix(mu_clust_init)
  
  n <- nrow(Z_init)
  K <- nrow(mu_clust_init)
  
  labels_init <- integer(n)
  
  for (i in 1:n) {
    d_i <- rep(NA_real_, K)
    
    for (g in 1:K) {
      d_i[g] <- hyp_dist_lorentz(Z_init[i, ], mu_clust_init[g, ])
    }
    
    labels_init[i] <- which.min(d_i)
  }
  
  labels_init
}


normalise_lambda <- function(lambda) {
  lambda <- as.numeric(lambda)
  s <- sum(lambda)
  
  lambda / s
}


match_clusters_to_reference <- function(
    mu_curr,
    mu_ref
){
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Please install the 'clue' package for solve_LSAP().")
  }
  
  mu_curr <- as.matrix(mu_curr)
  mu_ref  <- as.matrix(mu_ref)
  
  K <- nrow(mu_curr)
  C <- matrix(0, nrow = K, ncol = K)
  
  for (g in 1:K) {
    for (h in 1:K) {
      C[g, h] <- hyp_dist_lorentz(mu_curr[g, ], mu_ref[h, ])^2
      
    }
  }
  
  assignment <- as.integer(clue::solve_LSAP(C))
  assignment
}

apply_label_permutation <- function(
    perm,
    mu_clust,
    clustering_labels,
    lambda = NULL,
    sigma2_clust = NULL
){
  ord <- order(perm)
  
  out <- list(
    mu_clust = mu_clust[ord, , drop = FALSE],
    clustering_labels = perm[clustering_labels]
  )
  
  if (!is.null(lambda)) {
    out$lambda <- lambda[ord]
  }
  
  if (!is.null(sigma2_clust)) {
    out$sigma2_clust <- sigma2_clust[ord]
  }
  
  out
}

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

