
# Cluster Sampling algorithms ---------------------------------------------

PHN_sampling <- function(mu, sigma2){
  dimension <- length(mu)
  mu0 <- c(1, rep(0, dimension - 1))
  
  samp <- MASS::mvrnorm(
    1,
    mu = rep(0, dimension - 1),
    Sigma = sigma2 * diag(dimension - 1)
  )
  
  v <- c(0, samp)
  v_mu <- parallel_transport(mu0, mu, v)
  z <- exp_map_hyp(mu, v_mu)
  
  return(z)
}

log_pdf_PHN <- function(z, mu, sigma2) {
  n <- length(mu) - 1L
  r <- hyp_dist_lorentz(z, mu)
  lg_pv <-
    -(n / 2) * log(2 * pi) -
    (n / 2) * log(sigma2) -
    (r^2 / (2 * sigma2))
  
  lg_jac <- if (r < 1e-10) {
    0
  } else {
    (n - 1) * log(sinh(r) / r)
  }
  
  lg_pv - lg_jac
}

scalar_update_step <- function(alpha, 
                               B0, 
                               B1, 
                               X,
                               Z, 
                               mu_al, 
                               mu_beta, 
                               sigma2_al, 
                               sigma2_beta, 
                               loglik_old, 
                               sigma2_prop_al, 
                               sigma2_prop_B,
                               stratified_network,
                               sampled_controls,
                               working_set,
                               D,
                               mode){
  
  alpha_old <- alpha
  alpha_new <- rnorm(1, mean = alpha_old, sd = sqrt(sigma2_prop_al))
  logprior_old <- dnorm(alpha_old, mean = mu_al, sd = sqrt(sigma2_al), log = TRUE)
  logprior_new <- dnorm(alpha_new, mean = mu_al, sd = sqrt(sigma2_al), log = TRUE)
  
  loglik_new <- fast_case_control_loglik(working_set = working_set,
                                         alpha = alpha_new,
                                         B0 = B0,
                                         B1 = B1
  )
  
  log_ratio <- (loglik_new + logprior_new) -
    (loglik_old + logprior_old)
  
  if (isTRUE(is.finite(log_ratio)) && log(runif(1)) < log_ratio) {
    alpha_hat <- alpha_new
    loglik_hat <- loglik_new
    accepted <- 1
  } else {
    alpha_hat <- alpha_old
    loglik_hat <- loglik_old
    accepted <- 0
  }
  
  return(list(
    scalar_hat = alpha_hat,
    loglik_hat = loglik_hat,
    accepted = accepted,
    log_ratio = log_ratio
  ))
}

MV_update_step <- function(
    alpha,
    B0,
    B1,
    X = NULL,
    Z,
    muB,
    sigma2,
    sigma2_prop,
    loglik_old,
    learn_B0 = TRUE,
    learn_B1 = TRUE,
    stratified_network,
    sampled_controls,
    working_set,
    D,
    mode
) {
  
  if (is.null(X)) {
    B_old <- c(B1)
    muB_use <- if (length(muB) == 1L) {
      muB
    } else {
      muB[length(muB)]
    }
    
    sigma2_use <- if (length(sigma2) == 1L) {
      sigma2
    } else {
      sigma2[length(sigma2)]
    }
    
    sigma2_prop_use <- if (length(sigma2_prop) == 1L) {
      sigma2_prop
    } else {
      sigma2_prop[length(sigma2_prop)]
    }
    
  } else {
    B_old <- c(B0, B1)
    muB_use <- muB
    sigma2_use <- sigma2
    sigma2_prop_use <- sigma2_prop
  }
  
  d <- length(B_old)
  
  Sigma_prior <- if (length(sigma2_use) == 1L) {
    diag(sigma2_use, d)
  } else {
    as.matrix(sigma2_use)
  }
  
  Sigma_prop <- if (length(sigma2_prop_use) == 1L) {
    diag(sigma2_prop_use, d)
  } else {
    as.matrix(sigma2_prop_use)
  }
  
  active <- if (is.null(X)) {
    isTRUE(learn_B1)
  } else {
    c(
      rep(isTRUE(learn_B0), length(B0)),
      isTRUE(learn_B1)
    )
  }

    B_new <- B_old
  
  if (any(active)) {
    B_new[active] <- as.numeric(
      MASS::mvrnorm(
        n = 1,
        mu = B_old[active],
        Sigma = Sigma_prop[active, active, drop = FALSE]
      )
    )
  }
  
  if (length(muB_use) == 1L) {
    muB_use <- rep(muB_use, d)
  }
  
  if (is.null(X)) {
    
    B0_new <- B0
    B1_new <- B_new[1L]
    
  } else {
    
    p <- length(B0)
    
    B0_new <- B_new[seq_len(p)]
    B1_new <- B_new[p + 1L]
  }
  
  logprior_old <- mvtnorm::dmvnorm(
    B_old,
    mean = muB_use,
    sigma = Sigma_prior,
    log = TRUE
  )
  
  logprior_new <- mvtnorm::dmvnorm(
    B_new,
    mean = muB_use,
    sigma = Sigma_prior,
    log = TRUE
  )
  
  loglik_new <- fast_case_control_loglik(
    working_set = working_set,
    alpha = alpha,
    B0 = B0_new,
    B1 = B1_new
  )
  
  log_ratio <- (
    loglik_new + logprior_new
  ) - (
    loglik_old + logprior_old
  )
  
  if (
    isTRUE(is.finite(log_ratio)) &&
    log(runif(1)) < log_ratio
  ) {
    
    B0_hat <- B0_new
    B1_hat <- B1_new
    loglik_hat <- loglik_new
    accepted <- 1L
    
  } else {
    
    B0_hat <- B0
    B1_hat <- B1
    loglik_hat <- loglik_old
    accepted <- 0L
  }
  
  list(
    B0_hat = B0_hat,
    B1_hat = B1_hat,
    loglik_hat = loglik_hat,
    accepted = accepted,
    log_ratio = log_ratio
  )
}



# Accept/reject schemes ---------------------------------------------------



sigma2_update_step_clust <- function(
    Z,
    mu_clust,
    clustering_labels,
    sigma2_clust,
    g,
    ig_shape,
    ig_rate,
    sigma2_min = 1e-8
){
  K <- nrow(mu_clust)
  mu_g <- as.numeric(mu_clust[g, ])
  idx_g <- which(clustering_labels == g)
  p <- length(mu_g) - 1
  
  ssq_g <- cluster_tangent_ssq(
    Z = Z,
    mu_g = mu_g,
    idx_g = idx_g
  )
  
  m_g <- length(idx_g)
  
  post_shape <- ig_shape + 0.5 * m_g * p
  post_rate  <- ig_rate  + 0.5 * ssq_g
  
  sigma2_draw <- rinvgamma(
    n = 1,
    shape = post_shape,
    rate = post_rate
  )
  
  sigma2_clust[g] <- max(sigma2_draw, sigma2_min)
  
  list(
    sigma2_clust = sigma2_clust,
    g = g,
    n_g = m_g,
    ssq_g = ssq_g,
    post_shape = post_shape,
    post_rate = post_rate
  )
}

Z_update_step_clust <- function(
    Z, 
    Y, 
    X, 
    alpha, 
    B0, 
    B1,
    mu_clust, 
    sigma2_clust, 
    clustering_labels,
    sigma_z_prop,
    i,
    stratified_network,
    sampled_controls,
    D,
    loglik_function,
    collect_delta = FALSE
) {
  z_i <- as.numeric(Z[i, ])
  k_i <- clustering_labels[i]
  mu_k <- as.numeric(mu_clust[k_i, ])
  sigma_k <- sigma2_clust[k_i]
  H <- stratified_network$strata_dimension
  
  prop_out <- 
    proposal_tangent_rw_h2(
      z_old = z_i,
      sigma_prop = sigma_z_prop
    )
  
  z_new <- prop_out$z_new
  
  if (!all(is.finite(z_new))) {
    return(list(
      Z = Z,
      accepted = 0L,
      i = i,
      log_r = -Inf,
      delta_ih = NULL,
      d_i = NULL
    ))
  }
  
  n <- nrow(Z)
  
  edge_idx <- stratified_network$index[[i]][[1L]]
  controls_i <- sampled_controls[[i]]
  
  if (is.list(controls_i)) {
    control_idx <- unlist(
      controls_i,
      use.names = FALSE
    )
  } else {
    control_idx <- controls_i
  }
  
  needed_idx <- unique(c(
    edge_idx,
    control_idx
  ))
  
  needed_idx <- needed_idx[needed_idx != i]
  d_prop <- rep(NA_real_, n)
  d_prop[i] <- 0
  
  if (length(needed_idx) > 0L) {
    d_prop[needed_idx] <- hyp_dist_to_indices(
      z_i = z_new,
      Z = Z,
      idx = needed_idx
    )
  }
  
  log_prior_old <- log_pdf_PHN(
    z_i,
    mu_k,
    sigma_k
  )
  
  log_prior_new <- log_pdf_PHN(
    z_new,
    mu_k,
    sigma_k
  )
  
  
  log_like_old <- loglik_function(
    i = i,
    X = X,
    alpha = alpha,
    B0 = B0,
    B1 = B1,
    stratified_network = stratified_network,
    sampled_controls = sampled_controls,
    D = D,
    d_i = NULL
  )
  
  log_like_new <- loglik_function(
    i = i,
    X = X,
    alpha = alpha,
    B0 = B0,
    B1 = B1,
    stratified_network = stratified_network,
    sampled_controls = sampled_controls,
    D = NULL,
    d_i = d_prop
  )
  
  
  log_likelihood_old <- log_like_old$loglik_i
  log_likelihood_new <- log_like_new$loglik_i
  
  delta_ih <- NULL
  
  if(collect_delta){
  H <- stratified_network$strata_dimension
  strata_old <- log_like_old$strata_contribution
  strata_new <- log_like_new$strata_contribution
  delta_ih <- numeric(H)
  delta_ih[2:H] <- strata_new[2:H] - strata_old[2:H]
  }
  
  log_q_forward <- tryCatch(
    log_q_tangent_rw_h2(
      z_to = z_new,
      z_from = z_i,
      sigma_prop = sigma_z_prop
    ),
    error = function(e) -Inf
  )
  
  log_q_reverse <- tryCatch(
    log_q_tangent_rw_h2(
      z_to = z_i,
      z_from = z_new,
      sigma_prop = sigma_z_prop
    ),
    error = function(e) -Inf
  )
  
  if (!is.finite(log_q_forward) || !is.finite(log_q_reverse)) {
    return(list(
      Z = Z,
      accepted = 0L,
      i = i,
      log_r = -Inf,
      delta_ih = delta_ih,
      d_i = NULL
    ))
  }
  
  log_r <-
    (log_likelihood_new +
       log_prior_new +
       log_q_reverse) -
    (log_likelihood_old +
       log_prior_old +
       log_q_forward)
  
  if (is.finite(log_r) && log(runif(1)) < min(0, log_r)) {
    remaining_idx <- setdiff(
      seq_len(n),
      c(i, needed_idx)
    )
    
    if (length(remaining_idx) > 0L) {
      
      d_prop[remaining_idx] <- hyp_dist_to_indices(
        z_i = z_new,
        Z = Z,
        idx = remaining_idx
      )
    }
    
    Z[i, ] <- z_new
    accepted <- 1L
    
  } else {
    accepted <- 0L
  }
  
  list(
    Z = Z,
    accepted = accepted,
    i = i,
    log_r = log_r,
    delta_ih = delta_ih,
    d_i = if (accepted == 1) d_prop else NULL
  )
}





mu_update_step_clust <- function(
    mu_clust, Z, clustering_labels, sigma2_clust,
    sigma_mu_prop,
    g,
    mu_prior_centre,
    sigma_mu_prior
){
  K <- nrow(mu_clust)
  mu_old <- as.numeric(mu_clust[g, ])
  sigma_g <- sigma2_clust[g]
  idx_g <- which(clustering_labels == g)
  
  prop_out <- tryCatch(
    proposal_tangent_rw_h2(z_old = mu_old, sigma_prop = sigma_mu_prop),
    error = function(e) NULL
  )
  
  if (is.null(prop_out)) {
    return(list(
      mu_clust = mu_clust,
      accepted = 0,
      g = g,
      log_r = -Inf
    ))
  }
  
  mu_new <- prop_out$z_new
  
  if (!all(is.finite(mu_new))) {
    return(list(
      mu_clust = mu_clust,
      accepted = 0,
      g = g,
      log_r = -Inf
    ))
  }
  
  log_prior_old <- log_pdf_PHN(mu_old, mu_prior_centre, sigma_mu_prior)
  log_prior_new <- log_pdf_PHN(mu_new, mu_prior_centre, sigma_mu_prior)
  
  if (length(idx_g) > 0) {
    log_like_old <- sum(vapply(
      idx_g,
      function(i) log_pdf_PHN(Z[i, ], mu_old, sigma_g),
      numeric(1)
    ))
    
    log_like_new <- sum(vapply(
      idx_g,
      function(i) log_pdf_PHN(Z[i, ], mu_new, sigma_g),
      numeric(1)
    ))
  } else {
    log_like_old <- 0
    log_like_new <- 0
  }
  
  log_q_forward <- tryCatch(
    log_q_tangent_rw_h2(z_to = mu_new, z_from = mu_old, sigma_prop = sigma_mu_prop),
    error = function(e) -Inf
  )
  
  log_q_reverse <- tryCatch(
    log_q_tangent_rw_h2(z_to = mu_old, z_from = mu_new, sigma_prop = sigma_mu_prop),
    error = function(e) -Inf
  )
  
  log_r <- (log_prior_new + log_like_new + log_q_reverse) -
    (log_prior_old + log_like_old + log_q_forward)
  
  accepted <- 0
  
  if (is.finite(log_r) && log(runif(1)) < log_r) {
    mu_clust[g, ] <- mu_new
    accepted <- 1
  }
  
  list(
    mu_clust = mu_clust,
    accepted = accepted,
    g = g,
    log_r = log_r
  )
}


label_update_step_clust <- function(
    Z, mu_clust, sigma2_clust, lambda,
    clustering_labels,
    i
){
  K <- nrow(mu_clust)
  z_i <- as.numeric(Z[i, ])
  
  log_probs <- rep(NA_real_, K)
  
  for (g in 1:K) {
    log_probs[g] <- log(lambda[g]) +
      log_pdf_PHN(z_i, mu_clust[g, ], sigma2_clust[g])
  }
  
  log_norm <- log_sum_exp(log_probs)
  probs <- exp(log_probs - log_norm)
  
  label_hat <- sample.int(K, size = 1, prob = probs)
  clustering_labels[i] <- label_hat
  
  list(
    clustering_labels = clustering_labels,
    label_hat = label_hat,
    probs = probs,
    log_probs = log_probs
  )
}



lambda_update_step_clust <- function(clustering_labels, v_lambda, K) {
  counts <- tabulate(clustering_labels, nbins = K)
  lambda_hat <- rdirichlet_simple(v_lambda + counts)
  
  list(
    lambda_hat = lambda_hat,
    counts = counts
  )
}

# Tangent Proposal Schemes ------------------------------------------------

log_jacobian_exp_h2 <- function(u, eps = 1e-10) {
  u <- as.numeric(u)
  r2 <- as.numeric(lorentz_inner(u, u))
  
  r <- sqrt(max(r2, 0))
  
  if (r < eps) {
    return(0)
  }
  
  log(sinh(r) / r)
}

proposal_tangent_rw_h2 <- function(z_old, sigma_prop) {
  z_old <- as.numeric(z_old)
  mu0 <- c(1, 0, 0)

  v0 <- c(0, rnorm(2, mean = 0, sd = sqrt(sigma_prop)))
  
  u_fwd <- parallel_transport(mu0, z_old, v0)
  
  z_new <- exp_map_hyp(z_old, u_fwd)
  
  list(
    z_new = as.numeric(z_new),
    u_fwd = as.numeric(u_fwd),
    v0_fwd = as.numeric(v0[-1])
  )
}

log_q_tangent_rw_h2 <- function(z_to, z_from, sigma_prop) {
  z_to   <- as.numeric(z_to)
  z_from <- as.numeric(z_from)
  mu0 <- c(1, 0, 0)
  
  u <- inv_exp_map(z_to, z_from)

  v0 <- parallel_transport(z_from, mu0, u)
  a  <- as.numeric(v0[-1])
  
  log_q_tan <- sum(dnorm(a, mean = 0, sd = sqrt(sigma_prop), log = TRUE))
  
  log_jac <- log_jacobian_exp_h2(u)
  
  log_q_tan - log_jac
}


