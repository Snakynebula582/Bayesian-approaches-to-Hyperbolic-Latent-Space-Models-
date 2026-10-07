mixture_bic <- function(model, K) {
  
  Z_fixed <- model$Z_bary_mean
  d <- ncol(Z_fixed) - 1
  n <- nrow(Z_fixed)
  
  min_var <- 0.01
  if (K == 1L) {
    mu_origin <- c(1, rep(0, d))
    sigma2_init <- mean(model$sigma2_chain[, 1L])
    fit <- stats::optim(
      par = log(sigma2_init),
      fn = negative_loglik_K1,
      Z = Z_fixed,
      method = "L-BFGS-B",
      lower = log(min_var),
      control = list(maxit = 2000)
    )
    
    if (fit$convergence != 0L) {
      warning("optim did not converge successfully for K = 1")
    }
    
    sigma2_hat <- exp(fit$par)
    n_parameters <- 1L
    
    bic <-
      2 * fit$value +
      n_parameters * log(n)
    
    return(list(
      K = 1L,
      bic = bic,
      loglik = -fit$value,
      convergence = fit$convergence,
      mu = matrix(mu_origin, nrow = 1L),
      sigma2 = sigma2_hat,
      weights = 1,
      n_parameters = n_parameters
    ))
  }
  init_pars <- initialise_hyperbolic_clusters(
    Z_fixed,
    K = K
  )
  
  mu_init <- init_pars$mu_clust
  sigma2_init <- init_pars$sigma2_clust
  lambda_init <- init_pars$lambda
  
  par_init <- parameter_helper_BIC(
    mu_init,
    sigma2_init,
    lambda_init,
    K = K
  )
  
  lower <- c(
    rep(-Inf, K * d),
    rep(log(min_var), K),
    rep(-Inf, K - 1L)
  )
  
  fit <- stats::optim(
    par = par_init,
    fn = negative_loglik_BIC,
    Z = Z_fixed,
    K = K,
    method = "L-BFGS-B",
    lower = lower,
    control = list(maxit = 2000)
  )
  
  if (fit$convergence != 0L) {
    warning("optim did not converge successfully")
  }
  
  n_parameters <- length(par_init)
  bic <-
    2 * fit$value +
    n_parameters * log(n)
  pars <- unpack_parameters_BIC(
    fit$par,
    K,
    d
  )
  
  return(list(
    K = K,
    bic = bic,
    loglik = -fit$value,
    convergence = fit$convergence,
    mu = pars$mu,
    sigma2 = pars$sigma2,
    weights = exp(pars$log_lambda),
    n_parameters = n_parameters
  ))
}

network_loglik_BIC <- function(
    par,
    model,
    Y,
    X = NULL,
    learn_alpha,
    learn_B0,
    learn_B1
) {
  Z <- model$Z_bary_mean
  q <- if (is.null(X)) {
    0L
  } else if (length(dim(X)) == 2L) {
    1L
  } else {
    dim(X)[3L]
  }
  
  alpha <- model$alpha_chain[1L]
  
  B0 <- if (q > 0L) {
    as.numeric(model$B0_chain[1L, ])
  } else {
    numeric(0L)
  }
  
  B1 <- model$B1_chain[1L]
  pos <- 1L
  
  if (isTRUE(learn_alpha)) {
    alpha <- par[pos]
    pos <- pos + 1L
  }
  
  if (isTRUE(learn_B0) && q > 0L) {
    B0 <- par[pos:(pos + q - 1L)]
    pos <- pos + q
  }
  
  if (isTRUE(learn_B1)) {
    B1 <- par[pos]
  }
  
  D <- D_t(Z)
  idx <- upper.tri(Y)
  eta <- alpha - B1 * D[idx]

  if (!is.null(X)) {
    if (length(dim(X)) == 2L) {
      eta <- eta + B0 * X[idx]
    } else {
      X_work <- sapply(
        seq_len(q),
        function(k) X[, , k][idx]
      )
      eta <- eta + drop(X_work %*% B0)
    }
  }
  y <- Y[idx]
  loglik <- sum(
    y * eta - softplus(eta)
  )
  -loglik
}

negative_loglik_K1 <- function(log_sigma2, Z) {
  sigma2 <- exp(log_sigma2)
  n <- nrow(Z)
  d <- ncol(Z) - 1
  mu_origin <- c(1,rep(0, d)
  )
  
  loglik <- sum(
    vapply(
      seq_len(n),
      function(i) {
        log_pdf_PHN(
          Z[i, ],
          mu_origin,
          sigma2
        )
      },
      numeric(1)
    )
  )
  
  -loglik
}

network_bic <- function(model, Y, X = NULL) {
  
  learn_alpha <- model$init$learn_alpha
  learn_B0 <- model$init$learn_B0
  learn_B1 <- model$init$learn_B1
  
  q <- if (is.null(X)) {
    0L
  } else if (length(dim(X)) == 2L) {
    1L
  } else {
    dim(X)[3L]
  }
  
  if (q == 0L) {
    learn_B0 <- FALSE
  }
  par_init <- numeric(0L)
  
  if (isTRUE(learn_alpha)) {
    
    par_init <- c(
      par_init,
      mean(model$alpha_chain)
    )
  }
  
  if (isTRUE(learn_B0) && q > 0L) {
    
    B0_init <- if (q == 1L) {
      mean(model$B0_chain)
    } else {
      colMeans(model$B0_chain)
    }
    
    par_init <- c(
      par_init,
      B0_init
    )
  }
  
  if (isTRUE(learn_B1)) {
    
    par_init <- c(
      par_init,
      mean(model$B1_chain)
    )
  }
  
  n_parameters <- length(par_init)

  if (n_parameters > 0L) {
    fit <- stats::optim(
      par = par_init,
      fn = network_loglik_BIC,
      model = model,
      Y = Y,
      X = X,
      learn_alpha = learn_alpha,
      learn_B0 = learn_B0,
      learn_B1 = learn_B1,
      method = "BFGS"
    )
    if (fit$convergence != 0L) {
      warning("Network likelihood optimisation did not converge")
    }
    neg_loglik <- fit$value
    convergence <- fit$convergence
    
  } else {
    neg_loglik <- network_loglik_BIC(
      par = numeric(0L),
      model = model,
      Y = Y,
      X = X,
      learn_alpha = FALSE,
      learn_B0 = FALSE,
      learn_B1 = FALSE
    )
    
    fit <- NULL
    convergence <- 0L
  }
  n_ties <- sum(Y) / 2
  bic <-
    2 * neg_loglik +
    n_parameters * log(n_ties)
  
  return(list(
    bic = bic,
    loglik = -neg_loglik,
    convergence = convergence,
    n_parameters = n_parameters,
    n_obs = n_ties,
    learn_alpha = learn_alpha,
    learn_B0 = learn_B0,
    learn_B1 = learn_B1,
    fit = fit
  ))
}

Approx_BIC <- function(model, Y, X = NULL, K = NULL, K_vector = NULL){
  if(!is.null(K) && is.null(K_vector)){
  network_BIC <- network_bic(model,Y, X = X)
  mixture_BIC <- mixture_bic(model, K)
  
  total_BIC <- network_BIC$bic + mixture_BIC$bic
  return(list(total_BIC = total_BIC,network_BIC = network_BIC,mixture_BIC = mixture_BIC))
  }
  else if(is.null(K) && !is.null(K_vector)){
    network_BIC <- network_bic(model,Y, X = X)
    n_K <- length(K_vector)
    total_BIC <- numeric(length = n_K)
    
    All_K_mixture_BIC <- vector("list", length = n_K)
    for (i in K_vector) {
      K_i <- K_vector[i]
      holder <- mixture_bic(model, K_i)
      All_K_mixture_BIC[[i]] <- holder
      total_BIC[i] <- network_BIC$bic + holder$bic
    }
    return(list(total_BIC=total_BIC,network_BIC = network_BIC,mixture_BIC = All_K_mixture_BIC))
  }
}


parameter_helper_BIC <- function(mu_init, sigma2_init, lambda_init, K) {
  mu_spatial <- mu_init[, -1, drop = FALSE]
  log_sigma2 <- log(sigma2_init)
  logits <- log(lambda_init[seq_len(K - 1)] / lambda_init[K])
  c(
    as.vector(mu_spatial),
    log_sigma2,
    logits
  )
}

unpack_parameters_BIC <- function(par, K, d) {
    mu_spatial <- matrix(
    par[seq_len(K * d)],
    nrow = K,
    ncol = d
  )
    mu <- cbind(
    sqrt(1 + rowSums(mu_spatial^2)),
    mu_spatial
  )
  sigma2 <- exp(par[K * d + seq_len(K)])
  logits <- par[K * d + K + seq_len(K - 1)]
  
  a <- c(logits, 0)
  a <- a - max(a)
  log_lambda <- a - log(sum(exp(a)))
  
  list(mu = mu, sigma2 = sigma2, log_lambda = log_lambda)
}

negative_loglik_BIC <- function(par, Z, K) {
  d <- ncol(Z) - 1
  pars <- unpack_parameters_BIC(par, K, d)
  loglik <- 0
  for (i in seq_len(nrow(Z))) {
    component_terms <- numeric(K)
    for (g in seq_len(K)) {
      component_terms[g] <-
      pars$log_lambda[g] +
      log_pdf_PHN(Z[i, ], pars$mu[g, ], pars$sigma2[g])
    }
    m <- max(component_terms)
    loglik <- loglik + m + log(sum(exp(component_terms - m)))
  }
  
  return(-loglik)
}



logistic <- function(alpha, B0, B1, X, z_i,z_j){
  if (is.null(alpha)) alpha <- 0; if (is.null(B0)) B0 <- 0; if (is.null(B1)) B1 <- 1;
  if (is.null(X)){
  eta <- alpha - B1*hyp_dist_lorentz(z_i,z_j)
  p_ij <- exp(eta)/(1 + exp(eta))
  return(p_ij)}
  else {
  eta <- alpha + B0*X - B1*hyp_dist_lorentz(z_i,z_j)
  p_ij <- exp(eta)/(1 + exp(eta))
  return(p_ij)
  }
}

mean_distance_trace <- function(model){
  Z <- model$Z_chain
  mean_distance <- numeric(dim(Z)[1])
  for (t in seq_len(dim(Z)[1])) {
    D <- D_t(Z[t, , ])
    mean_distance[t] <- mean(D[upper.tri(D)])
  }
  
  plot(
    mean_distance,
    type = "l",
    xlab = "Iteration",
    ylab = "Mean pairwise distance"
  )}


posterior_pred_check <- function(model, X = NULL) {
  
  n <- dim(model$Z_chain)[2]
  post_size <- dim(model$Z_chain)[1]
  
  t <- sample(seq_len(post_size), 1)
  
  al <- model$alpha_chain[t]
  B0 <- model$B0_chain[t, ]
  B1 <- model$B1_chain[t]
  Z  <- model$Z_chain[t, , ]
  
  Y_sim <- matrix(0, n, n)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      X_ij <- if (is.null(X)){NULL} else {X[i,j,]}
      p_ij <- logistic(
        al, B0, B1, X_ij,
        Z[i, ],
        Z[j, ]
      )
      
      Y_sim[i, j] <- rbinom(1, 1, p_ij)
      Y_sim[j, i] <- Y_sim[i, j]
    }
  }
  
  Y_sim
}

graph_statistics <- function(g = NULL, Y = NULL) {
  if(!is.null(g)){
  gen_stats <- c(
    density = igraph::edge_density(g),
    clustering = igraph::transitivity(g, type = "global"),
    mean_distance = igraph::mean_distance(g, unconnected = TRUE),
    diameter = igraph::diameter(g, unconnected = TRUE)
  )
  degree_dist <- igraph::degree(g)
  }
  else if(is.null(g) && !is.null(Y)){
    g_y <- igraph::graph_from_adjacency_matrix(Y, mode = "undirected",diag = FALSE)
    
    gen_stats <- c(
      density = igraph::edge_density(g_y),
      clustering = igraph::transitivity(g_y, type = "global"),
      mean_distance = igraph::mean_distance(g_y),
      diameter = igraph::diameter(g_y)
    )
    degree_dist <- igraph::degree(g_y)
  }
  return(list(gen_stats = gen_stats,degree_dist = degree_dist))
}


batch_posterior_pred_check <- function(model,Y,X = NULL, n_rep = 100){
  n <- dim(model$Z_chain)[2]
  Y_dim <- nrow(Y)
  simulated_adjacency_matrices <- array(NA_real_,dim = c(n_rep,Y_dim,Y_dim))
  simulated_graphs <- vector("list",length = n_rep)
  simulated_graphs_statistics <- matrix(NA_real_,nrow = n_rep, ncol = 4)
  degree_distribution <- matrix(NA_real_,nrow = n_rep, ncol = n)
  
  for (i in 1:n_rep) {
    simulated_adjacency_matrices[i, , ] <- posterior_pred_check(model, X = X)   
  }
  simulated_graphs <- lapply(seq_len(n_rep),function(i){
    igraph::graph_from_adjacency_matrix(simulated_adjacency_matrices[i,,],
                                mode = "undirected",
                                diag = FALSE)
  })
  for (j in 1:n_rep){
    hold <- graph_statistics(simulated_graphs[[j]])
    simulated_graphs_statistics[j,] <- hold$gen_stats
    degree_distribution[j,] <- hold$degree_dist
  }
  colnames(simulated_graphs_statistics) <- c(
                                             "Edge Density",
                                             "transitivity",
                                             "Mean Distance",
                                             "Diameter")
  return(list(simulated_graphs = simulated_graphs,
              simulated_graphs_statistics = simulated_graphs_statistics,
              degree_distributions = degree_distribution))
}



