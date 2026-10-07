mcmc_clust <- function(
    # Y the networks adjacency matrix, K the number of clusters to be fitted, 
    # X an (n,n,p) array of covariate information (optional).
    Y,
    K,
    X = NULL,
    
    # Optional manual initial configurations for parameters to be estimated.
    # If left "NULL", they are automatically initialised. 
    Z_init = NULL,
    mu_clust_init = NULL,
    sigma2_clust_init = NULL,
    lambda_init = NULL,
    clustering_labels_init = NULL,
    
    #Number of iterations and burn, total saved draws will be n_iter - burn.
    n_iter,
    burn = 0,
    
    # Alpha initialisation and hyperparameter selection.
    alpha_init = 0,
    learn_alpha = FALSE,
    mu_alpha = 0,
    sigma2_alpha = .5,
    sigma2_prop_alpha = .3,
    
    # Beta initialisation and hyperparameter selection, B0 corresponds to the covariate
    # information scale parameters and is only used when X is not NULL.
    # B1 is the hyperbolic distance scale parameter.
    learn_B0 = FALSE,
    learn_B1 = FALSE,
    B0_init = 0,
    B1_init = 1,
    mu_B = 0,
    sigma2_B = .5,
    sigma2_prop_B = .3,
    
    # Proposal variance for the latent positions Z.
    sigma_z_prop = .5,
    
    # Proposal variance for the cluster means mu.
    sigma_mu_prop = .3,
    
    # For K > 1; hyperparameter selection for the cluster variance.
    learn_sigma2 = TRUE,
    sigma2_ig_shape = NULL,
    sigma2_ig_rate = NULL,
    sigma2_min = 1e-8,
    v_lambda = 1,
    
    # For K > 1; hyperparameter selection for the cluster means prior.
    # For K = 1; mu_prior_centre and sigma_mu_prior is fixed.
    mu_prior_centre = c(1, 0, 0),
    sigma_mu_prior = 1,
    
    # Version control; 
    # "pilot" estimates the weights for the case control likelihood, leave weights NULL.
    # "production" requires weights not to be NULL and builds the case control likelihood
    # using the provided fixed weights from the pilot run.
    mode = c("pilot","production"),
    weights = NULL 
){
  start_time <- Sys.time()
  mode <- match.arg(mode)
  loglik_fun <- switch(
    mode,
    pilot      = pilot_loglik_i,
    production = production_loglik_i
  )
  
  if (is.null(Z_init)) {
      g <- igraph::graph_from_adjacency_matrix(
        Y,
        mode = "undirected",
        diag = FALSE
      )
      
      D <- as.matrix(igraph::distances(g))
    
    
    Z_init <- G_to_Z_MDS(
      D = D,
      dim = 2,
      curvature = 1,
      alpha = 1.1,
      equi.adj = 0.5,
      radial_scale = 1
    )
  }
  
  Z <- as.matrix(Z_init)
  n <- nrow(Z)
  d <- ncol(Z)
  
  cluster_init <- initialise_hyperbolic_clusters(
    Z = Z,
    K = K,
    mu_clust_init = mu_clust_init,
    clustering_labels_init = clustering_labels_init,
    sigma2_clust_init = sigma2_clust_init,
    lambda_init = lambda_init,
    mu_prior_centre = mu_prior_centre,
    sigma2_min = sigma2_min,
    nstart = 25L,
    lambda_pseudocount = 1
  )
  
  mu_clust <- cluster_init$mu_clust
  clustering_labels <- cluster_init$clustering_labels
  sigma2_clust <- cluster_init$sigma2_clust
  lambda <- cluster_init$lambda
  mu_ref_relab <- mu_clust
  
  init <- list(
    Z_mds = Z_init,
    mu_clust = mu_clust,
    clustering_labels = clustering_labels,
    sigma2_clust = sigma2_clust,
    lambda = lambda,
    learn_alpha = learn_alpha,
    learn_B0 = learn_B0,
    learn_B1 = learn_B1
  )
  
  alpha <- alpha_init
  B1 <- as.numeric(B1_init)
  
  if (length(B1) != 1L) {
    stop("B1_init must be a single scalar.")
  }
  
  if (is.null(X)) {
    
    # B0 is unused when there are no covariates
    B0 <- numeric(0L)
    n_covariates <- 0L
    covariate_names <- NULL
    
  } else {
    
    if (length(dim(X)) != 3L) {
      stop("X must be an n x n x q covariate array.")
    }
    
    if (dim(X)[1L] != n || dim(X)[2L] != n) {
      stop(
        "The first two dimensions of X must equal the number of actors."
      )
    }
    
    n_covariates <- dim(X)[3L]
    covariate_names <- dimnames(X)[[3L]]
    
    if (is.null(covariate_names)) {
      covariate_names <- paste0("covariate_", seq_len(n_covariates))
    }
    
    if (length(B0_init) == 1L) {
      B0 <- rep(as.numeric(B0_init), n_covariates)
    } else if (length(B0_init) == n_covariates) {
      B0 <- as.numeric(B0_init)
    } else {
      stop(
        paste0(
          "B0_init must have length 1 or length ",
          n_covariates,
          ", the number of covariates in X."
        )
      )
    }
    
    names(B0) <- covariate_names
  }
  

  Z_ref <- Z_init
  Z_ref <- as.matrix(Z_ref)
  
  n_save <- n_iter - burn
  
  alpha_chain <- rep(NA_real_, n_save)
  B0_chain <- matrix(
    NA_real_,
    nrow = n_save,
    ncol = n_covariates,
    dimnames = list(
      iteration = NULL,
      covariate = covariate_names
    )
  )
  B1_chain <- rep(NA_real_, n_save)
  
  Z_chain <- array(NA_real_, dim = c(n_save, n, d))
  Z_chain_raw <- array(NA_real_, dim = c(n_save, n, d))
  
  mu_chain <- array(NA_real_, dim = c(n_save, K, d))
  mu_chain_raw <- array(NA_real_, dim = c(n_save, K, d))
  
  lambda_chain <- matrix(NA_real_, nrow = n_save, ncol = K)
  sigma2_chain <- matrix(NA_real_, nrow = n_save, ncol = K)
  clustering_chain <- matrix(NA_integer_, nrow = n_save, ncol = n)
  relabel_perm_chain <- matrix(NA_integer_, nrow = n_save, ncol = K)
  
  alpha_accept <- rep(NA_real_, n_save)
  B_accept <- rep(NA_real_, n_save)
  z_accept <- rep(NA_real_, n_save)
  mu_accept <- matrix(NA_real_, nrow = n_save, ncol = K)
  
  last_iter <- 0
  save_iter <- 0
  
  stratified_network <- strata(Y)
  N_ih <- stratified_network$N_ih
  distance_mat <- stratified_network$distance_mat
  strata_dimension <- stratified_network$strata_dimension
  number_of_nodes <- stratified_network$number_nodes
  
  pilot_control_pool <- lapply(seq_len(number_of_nodes), function(i) {
    unlist(
      stratified_network$index[[i]][2:strata_dimension],
      use.names = FALSE
    )
  })
  
  control_budget <- n_i0(
    Y = Y,
    strata = stratified_network,
    r = 1
  )
  n_i0 <- control_budget
  
  if (mode == "pilot") {
    collect_delta <- TRUE
    weights <- matrix(0, nrow = n, ncol = strata_dimension)
    weight_sum <- matrix(0, nrow = n, ncol = strata_dimension)
    w_t_full <- matrix(0, nrow = n, ncol = strata_dimension)
    delta_mat <- matrix(0, nrow = n, ncol = strata_dimension)
    weight_count <- 0
  } else {
    collect_delta <- FALSE
    if (is.null(weights)) {
      stop("Supply pilot weights for the production run.")
    }
    n_ih <- allocate_controls(
      weights = weights,
      stratified_network = stratified_network,
      n_i0 = control_budget
    )
    
    weight_sum <- NULL
    w_t_full <- NULL
    delta_mat <- NULL
    weight_count <- 0
  }
  
  D_iter <- D_t(Z)
  
  for (iter in 1:n_iter) {
    z_acc_this_iter <- 0
    update_order <- sample.int(n)
    
    if (mode == "pilot") {
      sampled_non_neighbours_iter <- pilot_sampler(
        pilot_control_pool = pilot_control_pool,
        n_i0 = n_i0
      )
    } else {
      sampled_non_neighbours_iter <- stratified_sampler(
        strata = stratified_network,
        n_ih = n_ih
      )
    }
    
    for (s in seq_along(update_order)) {
      
      delta_ih <- numeric(strata_dimension)
      i <- update_order[s]
      
      z_step <- Z_update_step_clust(
        Z = Z,
        Y = Y,
        X = X,
        alpha = alpha,
        B0 = B0,
        B1 = B1,
        mu_clust = mu_clust,
        sigma2_clust = sigma2_clust,
        clustering_labels = clustering_labels,
        sigma_z_prop = sigma_z_prop,
        i = i,
        stratified_network = stratified_network,
        sampled_controls = sampled_non_neighbours_iter,
        D = D_iter,
        loglik_fun = loglik_fun,
        collect_delta = collect_delta
      )
      
      if(isTRUE(collect_delta)){
      delta_mat[i, ] <- z_step$delta_ih}
      Z <- z_step$Z
      if (z_step$accepted == 1L) {
        D_iter[i, ] <- z_step$d_i
        D_iter[, i] <- z_step$d_i
      }
      z_acc_this_iter <- z_acc_this_iter + z_step$accepted
    }
    
    if (mode == "pilot" && iter > burn){
    delta_non_neighbours <- delta_mat[, 2:strata_dimension, drop = FALSE]
    denom <- rowSums(delta_non_neighbours)
    w_t <- abs(
      sweep(delta_non_neighbours, 1, denom, "/")
    )
    
    w_t_full[, 2:strata_dimension] <- w_t
    weight_sum <- weight_sum + w_t_full
    weight_count <- weight_count + 1L
    }
    
    working_set <- build_case_control_set(
      X = X,
      stratified_network = stratified_network,
      sampled_controls = sampled_non_neighbours_iter,
      D = D_iter,
      mode = mode
    )
    
    loglik_current <- fast_case_control_loglik(
      working_set = working_set,
      alpha = alpha,
      B0 = B0,
      B1 = B1
    )

    if (isTRUE(learn_alpha)) {
      alpha_step <- scalar_update_step(
        alpha = alpha,
        B0 = B0,
        B1 = B1,
        mu_al = mu_alpha,
        mu_beta = mu_B,
        sigma2_al = sigma2_alpha,
        sigma2_beta = sigma2_B,
        Z = Z,
        X = X,
        loglik_old = loglik_current,
        sigma2_prop_al = sigma2_prop_alpha,
        sigma2_prop_B = sigma2_prop_B,
        stratified_network = stratified_network,
        sampled_controls = sampled_non_neighbours_iter,
        D = D_iter,
        working_set = working_set,
        mode = mode
        )
    } else {
      alpha_step <- list(
        scalar_hat = alpha,
        loglik_hat = loglik_current,
        accepted = NA_real_
      )
    }
    
    
    alpha <- alpha_step$scalar_hat
    loglik_current <- alpha_step$loglik_hat
    
   
    if ((isTRUE(learn_B0) && !is.null(X)) || isTRUE(learn_B1)) {
      B_step <- MV_update_step(
        alpha = alpha,
        B0 = B0,
        B1 = B1,
        X = X,
        Z = Z,
        muB = mu_B,
        sigma2 = sigma2_B,
        sigma2_prop = sigma2_prop_B,
        loglik_old = loglik_current,
        learn_B0 = learn_B0,
        learn_B1 = learn_B1,
        stratified_network = stratified_network,
        sampled_controls = sampled_non_neighbours_iter,
        working_set = working_set,
        D = D_iter,
        mode = mode
      )
    } else {
      B_step <- list(
        B0_hat = B0,
        B1_hat = B1,
        loglik_hat = loglik_current,
        accepted = NA_real_
      )}
    
    
    B0 <- B_step$B0_hat
    B1 <- B_step$B1_hat
    loglik_current <- B_step$loglik_hat
    
    
    
    mu_acc_this_iter <- rep(0, K)
    
    if (K > 1L){
    for (g in 1:K) {
      mu_step <- mu_update_step_clust(
        mu_clust = mu_clust,
        Z = Z,
        clustering_labels = clustering_labels,
        sigma2_clust = sigma2_clust,
        sigma_mu_prop = sigma_mu_prop,
        g = g,
        mu_prior_centre = mu_prior_centre,
        sigma_mu_prior = sigma_mu_prior
      )
      
      mu_clust <- mu_step$mu_clust
      mu_acc_this_iter[g] <- mu_step$accepted
    }
    }
    
    if(isTRUE(learn_sigma2)){
      for (g in 1:K) {
        sigma2_step <- sigma2_update_step_clust(
          Z = Z,
          mu_clust = mu_clust,
          clustering_labels = clustering_labels,
          sigma2_clust = sigma2_clust,
          g = g,
          ig_shape = sigma2_ig_shape,
          ig_rate = sigma2_ig_rate,
          sigma2_min = sigma2_min
        )
        
        sigma2_clust <- sigma2_step$sigma2_clust
      }
    }
    
    
    if(K > 1L){
    label_order <- sample.int(n)
    
    for (s in seq_along(label_order)) {
      i <- label_order[s]
      
      lab_step <- label_update_step_clust(
        Z = Z,
        mu_clust = mu_clust,
        sigma2_clust = sigma2_clust,
        lambda = lambda,
        clustering_labels = clustering_labels,
        i = i
      )
      
      clustering_labels <- lab_step$clustering_labels
    }
    
    lambda_step <- lambda_update_step_clust(
      clustering_labels = clustering_labels,
      v_lambda = v_lambda,
      K = K
    )
    
    
    lambda <- lambda_step$lambda_hat
    } else {
      clustering_labels <- rep(1L, n)
      lambda <- 1
    }
    
    last_iter <- iter
    
    if (iter > burn) {
      save_iter <- save_iter + 1
      
      Z_chain_raw[save_iter, , ] <- Z
      if (K == 1L) {
        mu_chain_raw[save_iter, 1L, ] <- as.numeric(mu_prior_centre)
      } else {
        mu_chain_raw[save_iter, , ] <- matrix(
          mu_clust,
          nrow = K,
          ncol = d
        )
      }
      
      if (K == 1L) {
        
        align_out <- align_fixed_origin(
          Z = Z,
          Z_ref = Z_ref,
          allow_reflection = TRUE
        )
        
      } else {
        
        align_out <- align_one_clustering_draw(
          Z = Z,
          mu_clust = mu_clust,
          Z_ref = Z_ref,
          allow_reflection = TRUE
        )
      }
      
      Z_store <- align_out$Z_store
      mu_store <- align_out$mu_store
      
      labels_store <- clustering_labels
      lambda_store <- lambda
      sigma2_store <- sigma2_clust
      perm_store <- seq_len(K)
      
          perm_store <- match_clusters_to_reference(
            mu_curr = mu_store,
            mu_ref = mu_ref_relab
          )
          
          relab_store <- apply_label_permutation(
            perm = perm_store,
            mu_clust = mu_store,
            clustering_labels = labels_store,
            lambda = lambda_store,
            sigma2_clust = sigma2_store
          )
          
          mu_store <- relab_store$mu_clust
          labels_store <- relab_store$clustering_labels
          lambda_store <- relab_store$lambda
          sigma2_store <- relab_store$sigma2_clust
        
      
      
      alpha_chain[save_iter] <- alpha
      if (n_covariates > 0L) {
        B0_chain[save_iter, ] <- B0
      }
      B1_chain[save_iter] <- B1
      B_accept[save_iter] <- B_step$accepted
      Z_chain[save_iter, , ] <- Z_store
      if (K == 1L) {
        mu_chain[save_iter, 1L, ] <- as.numeric(mu_prior_centre)
      } else {
        mu_chain[save_iter, , ] <- matrix(
          mu_store,
          nrow = K,
          ncol = d
        )
      }
      lambda_chain[save_iter, ] <- lambda_store
      sigma2_chain[save_iter, ] <- sigma2_store
      clustering_chain[save_iter, ] <- labels_store
      relabel_perm_chain[save_iter, ] <- perm_store
      
      alpha_accept[save_iter] <- alpha_step$accepted
      z_accept[save_iter] <- z_acc_this_iter / n
      mu_accept[save_iter, ] <- mu_acc_this_iter
    }
   
  }
  
  weights <- if(mode == "pilot" && weight_count >0){weight_sum/weight_count}else{NULL}
  
  if (save_iter > 0) {
    alpha_chain <- alpha_chain[1:save_iter]
    
    B0_chain <- B0_chain[
      seq_len(save_iter),
      ,
      drop = FALSE
    ]
    
    B1_chain <- B1_chain[1:save_iter]
    B_accept <- B_accept[1:save_iter]
    
    Z_chain <- Z_chain[1:save_iter, , , drop = FALSE]
    Z_chain_raw <- Z_chain_raw[1:save_iter, , , drop = FALSE]
    
    mu_chain <- mu_chain[1:save_iter, , , drop = FALSE]
    mu_chain_raw <- mu_chain_raw[1:save_iter, , , drop = FALSE]
    
    lambda_chain <- lambda_chain[1:save_iter, , drop = FALSE]
    sigma2_chain <- sigma2_chain[1:save_iter, , drop = FALSE]
    clustering_chain <- clustering_chain[1:save_iter, , drop = FALSE]
    relabel_perm_chain <- relabel_perm_chain[1:save_iter, , drop = FALSE]
    
    alpha_accept <- alpha_accept[1:save_iter]
    z_accept <- z_accept[1:save_iter]
    mu_accept <- mu_accept[1:save_iter, , drop = FALSE]
    
    alpha_mean <- mean(alpha_chain)
    alpha_median <- median(alpha_chain)
    lambda_mean <- colMeans(lambda_chain)
    sigma2_mean <- colMeans(sigma2_chain)
    
    label_posterior_probs <- matrix(0, nrow = n, ncol = K)
    for (i in 1:n) {
      label_posterior_probs[i, ] <-
        tabulate(clustering_chain[, i], nbins = K) / save_iter
    }
    
    clustering_map <- max.col(label_posterior_probs, ties.method = "first")
  } 
    
  Z_bary_mean <- NULL
  mu_bary_mean <- NULL
  
  if (save_iter > 0) {
    Z_bary_mean <- barycenters(
      Z_chain = Z_chain
    )
  }
  
  if (save_iter > 0) {
    mu_bary_mean <- matrix(NA_real_, nrow = K, ncol = d)
    
    for (g in 1:K) {
      mu_chain_g <- mu_chain[, g, , drop = FALSE]
      dim(mu_chain_g) <- c(save_iter, 1, d)
      
      mu_bary_mean[g, ] <- barycenters(
        Z_chain = mu_chain_g
      )[1, ]
    }
  }
  
  if (n_covariates > 0L && save_iter > 0L) {
    
    B0_mean <- colMeans(B0_chain)
    
    B0_median <- apply(
      B0_chain,
      MARGIN = 2L,
      FUN = median
    )
    
  } else {
    
    B0_mean <- numeric(0L)
    B0_median <- numeric(0L)
  }
  
  end_time <- Sys.time()
  
  runtime <- end_time - start_time
  
  return(list(
    alpha_chain = alpha_chain,
    B0_chain = B0_chain,
    B1_chain = B1_chain,
    B_accept = B_accept,
    B_accept_rate = mean(B_accept, na.rm = TRUE),
    B0_mean = B0_mean,
    B0_median = B0_median,
    B1_mean = mean(B1_chain),
    Z_chain = Z_chain,
    Z_chain_raw = Z_chain_raw,
    mu_chain = mu_chain,
    mu_chain_raw = mu_chain_raw,
    lambda_chain = lambda_chain,
    sigma2_chain = sigma2_chain,
    clustering_chain = clustering_chain,
    relabel_perm_chain = relabel_perm_chain,
    alpha_accept = alpha_accept,
    z_accept = z_accept,
    mu_accept = mu_accept,
    Z_ref = Z_ref,
    last_iter = last_iter,
    burn = burn,
    n_saved = save_iter,
    alpha_mean = alpha_mean,
    alpha_median = alpha_median,
    lambda_mean = lambda_mean,
    sigma2_mean = sigma2_mean,
    label_posterior_probs = label_posterior_probs,
    clustering_map = clustering_map,
    Z_bary_mean = Z_bary_mean,
    mu_bary_mean = mu_bary_mean,
    sigma2_clust = sigma2_clust,
    init = init,
    weights = weights,
    runtime = runtime
  ))
}