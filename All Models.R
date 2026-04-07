Radial_Angular_Model <- function(Y,
                                 d = 2,
                                 n_iter = 5000,
                                 burn = 1000,
                                 thin = 1,
                                 mu_alpha = 0,
                                 sigma_alpha = 0.3,
                                 step_r,
                                 step_u,
                                 step_alpha = 0.2,
                                 postproc = c("naive", "hpa"),
                                 allow_reflection = TRUE,
                                 r_init = NULL,
                                 U_init = NULL,
                                 alpha_init = 0,
                                 sigma_r_prior,
                                 learn_curvature = FALSE,
                                 K = 1,
                                 K_init = K,
                                 step_K = 0.15,
                                 mu_log_K = 0,
                                 sd_log_K = 0.15,
                                 compute_posterior_mean = TRUE,
                                 compute_barycenter = FALSE,
                                 bary_tol = 1e-6,
                                 bary_max_iter = 100,
                                 drop_bad = TRUE,
                                 bary_tol_hyp = 1e-5) {
  
  postproc <- match.arg(postproc)  
  Y <- as.matrix(Y)
  n <- nrow(Y)
  
  if (!is.null(r_init) && !is.null(U_init)) {
    r <- as.numeric(r_init)
    U <- as.matrix(U_init)
    alpha <- alpha_init
  } else {
    U <- matrix(rnorm(n * d), nrow = n, ncol = d)
    U <- U / sqrt(rowSums(U^2))
    r <- abs(rnorm(n, mean = 0, sd = sigma_r_prior))
    alpha <- alpha_init
  }
  
  K_curr <- if (learn_curvature) K_init else K
  row_norms <- sqrt(rowSums(U^2))
  U <- U / row_norms
  
  keep_idx <- seq(from = burn + 1, to = n_iter, by = thin)
  n_keep <- length(keep_idx)
  
  out_alpha <- numeric(n_keep)
  out_r_raw <- matrix(NA_real_, nrow = n_keep, ncol = n)
  out_U_raw <- array(NA_real_, dim = c(n_keep, n, d))
  out_K <- numeric(n_keep)
  
  acc_alpha <- 0
  acc_joint <- rep(0, n)
  acc_K <- 0
  
  kk <- 0
  
  for (t in 1:n_iter) {
    
    ll_curr <- loglik_full(Y, alpha, r, U, K = K_curr)
    lp_a_curr <- dnorm(alpha, mean = mu_alpha, sd = sigma_alpha, log = TRUE)
    
    alpha_prop <- propose_alpha(alpha, step_alpha)
    ll_prop <- loglik_full(Y, alpha_prop, r, U, K = K_curr)
    lp_a_prop <- dnorm(alpha_prop, mean = mu_alpha, sd = sigma_alpha, log = TRUE)
    
    log_acc_alpha <- (ll_prop + lp_a_prop) - (ll_curr + lp_a_curr)
    
    if (log(runif(1)) < log_acc_alpha) {
      alpha <- alpha_prop
      acc_alpha <- acc_alpha + 1
    }
    
    for (i in 1:n) {
      r_old <- r[i]
      u_old <- U[i, ]
      theta_old <- u_to_theta(u_old)
      
      ll_i_curr <- loglik_node(Y, alpha, r, U, i, K = K_curr)
      lp_r_curr <- log_halfnorm(r_old, sigma = sigma_r_prior)
      
      prop_i <- propose_r_theta_joint(
        r = r_old,
        theta = theta_old,
        step_r = step_r,
        step_theta = step_u
      )
      
      r[i] <- prop_i$r
      U[i, ] <- prop_i$u
      
      ll_i_prop <- loglik_node(Y, alpha, r, U, i, K = K_curr)
      lp_r_prop <- log_halfnorm(r[i], sigma = sigma_r_prior)
      
      log_acc_joint <- (ll_i_prop + lp_r_prop) - (ll_i_curr + lp_r_curr)
      
      if (log(runif(1)) < log_acc_joint) {
        acc_joint[i] <- acc_joint[i] + 1
      } else {
        r[i] <- r_old
        U[i, ] <- u_old
      }
    }
    
    if (learn_curvature) {
      log_K_old <- log(K_curr)
      log_K_prop <- propose_log_K(log_K_old, step_log_K = step_K)
      K_prop <- exp(log_K_prop)
      
      ll_K_curr <- loglik_full(Y, alpha, r, U, K = K_curr)
      ll_K_prop <- loglik_full(Y, alpha, r, U, K = K_prop)
      
      lp_K_curr <- log_prior_K(K_curr, mu_log_K, sd_log_K)
      lp_K_prop <- log_prior_K(K_prop, mu_log_K, sd_log_K)
      
      log_acc_K <- (ll_K_prop + lp_K_prop + log_K_prop) -
        (ll_K_curr + lp_K_curr + log_K_old)
      
      if (log(runif(1)) < log_acc_K) {
        K_curr <- K_prop
        acc_K <- acc_K + 1
      }
    }
    
    if (t %in% keep_idx) {
      kk <- kk + 1
      out_alpha[kk] <- alpha
      out_r_raw[kk, ] <- r
      out_U_raw[kk, , ] <- U
      out_K[kk] <- K_curr
    }
    
  }
  
  U_ref <- U_init
  r_ref <- r_init
  
  if (is.null(r_ref) || is.null(U_ref)) {
    r_ref <- out_r_raw[1, ]
    U_ref <- out_U_raw[1, , ]
    if (is.null(dim(U_ref))) {
      U_ref <- matrix(U_ref, nrow = n, ncol = d)
    }
  } else {
    U_ref <- as.matrix(U_ref)
  }
  
  fit <- list(
    alpha = out_alpha,
    r_raw = out_r_raw,
    U_raw = out_U_raw,
    K = out_K,
    acc = list(
      alpha = acc_alpha / n_iter,
      joint = acc_joint / n_iter,
      r = acc_joint / n_iter,
      u = acc_joint / n_iter,
      K = if (learn_curvature) acc_K / n_iter else NA_real_
    ),
    prior = list(
      r_prior = "half-normal",
      sigma_r_prior = sigma_r_prior,
      theta_prior = "uniform on circle",
      K_prior = if (learn_curvature) {
        list(name = "lognormal", mu_log_K = mu_log_K, sd_log_K = sd_log_K)
      } else {
        list(name = "fixed", value = K)
      }
    ),
    keep_idx = keep_idx,
    last = list(
      alpha = alpha,
      r = r,
      U = U,
      K = K_curr
    ),
    learn_curvature = learn_curvature,
    postproc = postproc
  )
  
  if (postproc == "naive") {
    out_r <- out_r_raw
    out_U <- array(NA_real_, dim = c(n_keep, n, d))
    
    for (m in 1:n_keep) {
      U_m <- out_U_raw[m, , ]
      if (is.null(dim(U_m))) {
        U_m <- matrix(U_m, nrow = n, ncol = d)
      }
      out_U[m, , ] <- align_directions(U_m, U_ref)
    }
    
    fit$r <- out_r
    fit$U <- out_U
    
    Z_chain <- array(NA_real_, dim = c(n_keep, n, d + 1))
    for (m in 1:n_keep) {
      U_m <- out_U[m, , ]
      if (is.null(dim(U_m))) {
        U_m <- matrix(U_m, nrow = n, ncol = d)
      }
      Z_chain[m, , ] <- polar_to_hyperboloid(
        r = out_r[m, ],
        u_mat = U_m,
        R = 1
      )
    }
    
    fit$Z <- Z_chain
    fit$Z_ref <- polar_to_hyperboloid(r = r_ref, u_mat = U_ref, R = 1)
    fit$Z_ref_bar <- NULL
    
    if (compute_posterior_mean) {
      U_mean <- apply(out_U, c(2, 3), mean)
      row_norms_mean <- sqrt(rowSums(U_mean^2))
      U_mean <- U_mean / row_norms_mean
      r_mean <- colMeans(out_r)
      Z_mean <- polar_to_hyperboloid(r = r_mean, u_mat = U_mean, R = 1)
      
      fit$posterior_mean <- list(
        alpha = mean(out_alpha),
        r = r_mean,
        U = U_mean,
        Z = Z_mean,
        K = mean(out_K)
      )
    }
    
  } else if (postproc == "hpa") {
    Z_ref <- polar_to_hyperboloid(r = r_ref, u_mat = U_ref, R = 1)
    
    Z_ref_bar <- posterior_barycenter_one_node(
      samples_node = Z_ref,
      mu_init = Z_ref[1, ],
      tol = bary_tol,
      max_iter = bary_max_iter,
      drop_bad = drop_bad,
      tol_hyp = bary_tol_hyp
    )
    Z_ref_bar <- project_to_hyperboloid(Z_ref_bar)
    
    Z_chain_raw <- array(NA_real_, dim = c(n_keep, n, d + 1))
    Z_chain <- array(NA_real_, dim = c(n_keep, n, d + 1))
    
    for (m in 1:n_keep) {
      U_m <- out_U_raw[m, , ]
      if (is.null(dim(U_m))) {
        U_m <- matrix(U_m, nrow = n, ncol = d)
      }
      
      Z_m <- polar_to_hyperboloid(
        r = out_r_raw[m, ],
        u_mat = U_m,
        R = 1
      )
      
      Z_chain_raw[m, , ] <- Z_m
      
      hpa_out <- HPA_align_draw(
        Z_ref = Z_ref,
        Z_new = Z_m,
        Z_ref_bar = Z_ref_bar,
        allow_reflection = allow_reflection,
        bary_tol = bary_tol,
        bary_max_iter = bary_max_iter,
        drop_bad = drop_bad,
        tol_hyp = bary_tol_hyp,
        tol_tangent = 1e-8,
        repair = TRUE
      )
      
      Z_chain[m, , ] <- hpa_out$Z_aligned
    }
    
    fit$Z_raw <- Z_chain_raw
    fit$Z <- Z_chain
    fit$Z_ref <- Z_ref
    fit$Z_ref_bar <- Z_ref_bar
    
    if (compute_posterior_mean) {
      Z_mean <- posterior_polar_means(Z_chain)
      fit$posterior_mean <- list(
        alpha = mean(out_alpha),
        Z = Z_mean,
        K = mean(out_K)
      )
    }
  }
  
  if (compute_barycenter) {
    fit$posterior_barycenter <- posterior_barycenters(
      Z_chain = fit$Z,
      tol = bary_tol,
      max_iter = bary_max_iter,
      drop_bad = drop_bad,
      tol_hyp = bary_tol_hyp
    )
  }
  
  fit
}



Centred_Wrapped_Normal <- function(
    Y,
    Z_init = NULL,
    alpha_init,
    n_iter,
    burn = 0,
    mu_alpha,
    sigma2_alpha,
    sigma2_prop_alpha,
    sigma2_z,
    sigma_prop_z,
    mu_z,
    mds_dim = 2,
    mds_curvature = 1,
    mds_alpha = 1.1,
    mds_equi_adj = 0,
    radial_scale = 1,
    alignment_method = c("none", "angular", "hpa"),
    Z_ref = NULL,
    allow_reflection = TRUE,
    compute_barycenter = TRUE,
    bary_tol = 1e-6,
    bary_max_iter = 100,
    bary_drop_bad = TRUE,
    bary_tol_hyp = 1e-5,
    compute_polar_mean = FALSE
){
  alignment_method <- match.arg(alignment_method)
  
  if (is.null(Z_init)) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("igraph is required to compute D from Y when Z_init isNULL.")
    }
    
    g <- igraph::graph_from_adjacency_matrix(
      Y,
      mode = "undirected",
      diag = FALSE
    )
    
    D <- as.matrix(igraph::distances(g))
    
    Z_init <- G_to_Z_MDS(
      D = D,
      dim = mds_dim,
      curvature = mds_curvature,
      alpha = mds_alpha,
      equi.adj = mds_equi_adj,
      radial_scale = radial_scale
    )
  }
  
  Z <- Z_init
  alpha <- alpha_init
  
  n <- nrow(Z)
  d <- ncol(Z)
  
  if (is.null(Z_ref)) {
    Z_ref <- Z_init
  }
  
  for (i in seq_len(nrow(Z_ref))) {
    Z_ref[i, ] <- project_to_hyperboloid(Z_ref[i, ])
  }
  
  Z_ref_bar <- NULL
  if (alignment_method == "hpa") {
    Z_ref_bar <- posterior_barycenter_one_node(
      samples_node = Z_ref,
      mu_init = Z_ref[1, ],
      tol = bary_tol,
      max_iter = bary_max_iter,
      drop_bad = bary_drop_bad,
      tol_hyp = bary_tol_hyp
    )
  }
  
  n_save <- n_iter - burn
  
  alpha_chain <- rep(NA_real_, n_save)
  Z_chain <- array(NA_real_, dim = c(n_save, n, d))
  Z_chain_raw <- array(NA_real_, dim = c(n_save, n, d))
  alpha_accept <- rep(NA_real_, n_save)
  z_accept <- rep(NA_real_, n_save)
  proc_angle <- rep(NA_real_, n_save)
  
  loglik_current <- Full_likelihood_initial(Z, Y, alpha)
  last_iter <- 0
  save_iter <- 0
  
  for (iter in 1:n_iter) {
    step_out <- tryCatch({
      
      alpha_step <- alpha_update_step(
        alpha = alpha,
        mu = mu_alpha,
        sigma2 = sigma2_alpha,
        Z = Z,
        Y = Y,
        loglik_old = loglik_current,
        sigma2_prop = sigma2_prop_alpha
      )
      
      alpha <- alpha_step$alpha_hat
      loglik_current <- alpha_step$loglik_hat
      
      z_acc_this_iter <- 0
      
      for (k in sample.int(n)) {
        z_step <- Z_update_step_noclust(
          Z = Z,
          Y = Y,
          alpha = alpha,
          sigma2_prior = sigma2_z,
          mu_z = mu_z,
          sigma_prop_z = sigma_prop_z,
          i = k,
          loglik_old = loglik_current
        )
        
        Z <- z_step$Z
        loglik_current <- z_step$loglik_hat
        z_acc_this_iter <- z_acc_this_iter + z_step$accepted
      }
      
      last_iter <- iter
      
      if (iter > burn) {
        save_iter <- save_iter + 1
        
        Z_chain_raw[save_iter, , ] <- Z
        
        if (alignment_method == "none") {
          Z_store <- Z
          proc_angle[save_iter] <- NA_real_
          
        } else if (alignment_method == "angular") {
          proc_out <- procrustes_transform(Z, Z_ref, renormalise = TRUE)
          Z_store <- proc_out$Z_aligned
          proc_angle[save_iter] <- proc_out$alpha
          
        } else if (alignment_method == "hpa") {
          hpa_out <- HPA_align_draw(
            Z_ref = Z_ref,
            Z_new = Z,
            Z_ref_bar = Z_ref_bar,
            allow_reflection = allow_reflection,
            bary_tol = bary_tol,
            bary_max_iter = bary_max_iter,
            drop_bad = bary_drop_bad,
            tol_hyp = bary_tol_hyp,
            repair = TRUE
          )
          
          Z_store <- hpa_out$Z_aligned
          proc_angle[save_iter] <- NA_real_
          
        } else {
          stop("Unknown alignment_method.")
        }
        
        alpha_chain[save_iter] <- alpha
        Z_chain[save_iter, , ] <- Z_store
        alpha_accept[save_iter] <- alpha_step$accepted
        z_accept[save_iter] <- z_acc_this_iter / n
      }
      
      TRUE
    }, error = function(e) {
      warning(paste("MCMC stopped early due to error, likely numerical"))
      FALSE
    })
    
    if (!step_out) break
  }
  
  if (save_iter > 0) {
    alpha_chain <- alpha_chain[1:save_iter]
    Z_chain <- Z_chain[1:save_iter, , , drop = FALSE]
    Z_chain_raw <- Z_chain_raw[1:save_iter, , , drop = FALSE]
    alpha_accept <- alpha_accept[1:save_iter]
    z_accept <- z_accept[1:save_iter]
    proc_angle <- proc_angle[1:save_iter]
    
    alpha_mean <- mean(alpha_chain)
    alpha_median <- median(alpha_chain)
  } else {
    alpha_chain <- numeric(0)
    Z_chain <- array(NA_real_, dim = c(0, n, d))
    Z_chain_raw <- array(NA_real_, dim = c(0, n, d))
    alpha_accept <- numeric(0)
    z_accept <- numeric(0)
    proc_angle <- numeric(0)
    alpha_mean <- NA_real_
    alpha_median <- NA_real_
  }
  
  if (compute_barycenter && save_iter > 0) {
    Z_barycenter <- posterior_barycenters(
      Z_chain = Z_chain,
      tol = bary_tol,
      max_iter = bary_max_iter,
      drop_bad = bary_drop_bad,
      tol_hyp = bary_tol_hyp
    )
  } else {
    Z_barycenter <- NULL
  }
  
  if (compute_polar_mean && save_iter > 0) {
    Z_polar_mean <- posterior_polar_means(Z_chain)
  } else {
    Z_polar_mean <- NULL
  }
  
  return(list(
    alpha_chain = alpha_chain,
    Z_chain = Z_chain,
    Z_chain_raw = Z_chain_raw,
    alpha_accept = alpha_accept,
    z_accept = z_accept,
    proc_angle = proc_angle,
    Z_ref = Z_ref,
    Z_ref_bar = Z_ref_bar,
    alignment_method = alignment_method,
    last_iter = last_iter,
    burn = burn,
    n_saved = save_iter,
    alpha_mean = alpha_mean,
    alpha_median = alpha_median,
    Z_barycenter = Z_barycenter,
    Z_polar_mean = Z_polar_mean
  ))
}



Wrapped_Normal_Clustering_model <- function(
    Y,
    K,
    Z_init = NULL,
    clustering_labels_init = NULL,
    mu_clust_init = NULL,
    lambda_init = NULL,
    alpha_init,
    n_iter,
    burn = 0,
    mu_alpha,
    sigma2_alpha,
    sigma2_prop_alpha,
    sigma_z_prop,
    sigma_mu_prop,
    sigma2_clust_init,
    sigma2_ig_shape,
    sigma2_ig_rate,
    sigma2_min = 1e-8,
    v_lambda = 1,
    mu_prior_centre = c(1, 0, 0),
    sigma_mu_prior = 1,
    mu_init_r = 1.5,
    align_method = c("none", "polar", "HPA"),
    Z_ref = NULL,
    relabel_by_angle = FALSE,
    relabel_reference = TRUE,   
    compute_barycenter = TRUE,
    bary_tol = 1e-6,
    bary_max_iter = 100,
    bary_drop_bad = TRUE,
    bary_tol_hyp = 1e-5,
    compute_polar_mean = TRUE
){
  align_method <- match.arg(align_method)
  sigma2_clust <- as_sigma_vector(sigma2_clust_init, K)
  v_lambda <- as_dirichlet_vector(v_lambda, K)
  
  if (is.null(Z_init)) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("igraph is required to compute D from Y when Z_init and D are both NULL.")
    }
    
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
      alpha = 1,
      equi.adj = 0,
      radial_scale = 1
    )
  }
  
  Z <- as.matrix(Z_init)
  
  n <- nrow(Z)
  d <- ncol(Z)
  
  if (is.null(mu_clust_init)) {
    mu_clust_init <- initial_cluster_means_clust(
      K = K,
      r = mu_init_r,
      R = 1
    )
  }
  
  mu_clust <- as.matrix(mu_clust_init)
  
  if (is.null(clustering_labels_init)) {
    clustering_labels_init <- initial_labels_clust(Z, mu_clust)
  }
  
  clustering_labels <- as.integer(clustering_labels_init)
  
  if (is.null(lambda_init)) {
    lambda_init <- tabulate(clustering_labels, nbins = K) / n
  } 
  
  lambda <- normalise_lambda(lambda_init)
  alpha <- alpha_init
  
  if (is.null(Z_ref)) {
    Z_ref <- Z_init
  }
  Z_ref <- as.matrix(Z_ref)
  
  n_save <- n_iter - burn
  
  alpha_chain <- rep(NA_real_, n_save)
  Z_chain <- array(NA_real_, dim = c(n_save, n, d))
  Z_chain_raw <- array(NA_real_, dim = c(n_save, n, d))
  
  mu_chain <- array(NA_real_, dim = c(n_save, K, d))
  mu_chain_raw <- array(NA_real_, dim = c(n_save, K, d))
  
  lambda_chain <- matrix(NA_real_, nrow = n_save, ncol = K)
  sigma2_chain <- matrix(NA_real_, nrow = n_save, ncol = K)
  clustering_chain <- matrix(NA_integer_, nrow = n_save, ncol = n)
  relabel_perm_chain <- matrix(NA_integer_, nrow = n_save, ncol = K)
  
  alpha_accept <- rep(NA_real_, n_save)
  z_accept <- rep(NA_real_, n_save)
  mu_accept <- matrix(NA_real_, nrow = n_save, ncol = K)
  proc_angle <- rep(NA_real_, n_save)
  
  loglik_current <- Full_likelihood_initial(Z, Y, alpha)
  last_iter <- 0
  save_iter <- 0
  
  mu_ref_relab <- NULL
  labels_ref_relab <- NULL
  lambda_ref_relab <- NULL
  
  for (iter in 1:n_iter) {
    
    z_acc_this_iter <- 0
    update_order <- sample.int(n)
    
    for (s in seq_along(update_order)) {
      i <- update_order[s]
      
      z_step <- Z_update_step_clust(
        Z = Z,
        Y = Y,
        alpha = alpha,
        mu_clust = mu_clust,
        sigma2_clust = sigma2_clust,
        clustering_labels = clustering_labels,
        sigma_z_prop = sigma_z_prop,
        i = i
      )
      
      Z <- z_step$Z
      z_acc_this_iter <- z_acc_this_iter + z_step$accepted
    }
    
    loglik_current <- Full_likelihood_initial(Z, Y, alpha)
    
    alpha_step <- alpha_update_step(
      alpha = alpha,
      mu = mu_alpha,
      sigma2 = sigma2_alpha,
      Z = Z,
      Y = Y,
      loglik_old = loglik_current,
      sigma2_prop = sigma2_prop_alpha
    )
    
    alpha <- alpha_step$alpha_hat
    loglik_current <- alpha_step$loglik_hat
    
    mu_acc_this_iter <- rep(0, K)
    
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
    
    if (relabel_by_angle && !relabel_reference) {
      relab <- relabel_by_mu_angle(
        mu_clust = mu_clust,
        clustering_labels = clustering_labels,
        lambda = lambda,
        sigma2_clust = sigma2_clust
      )
      
      mu_clust <- relab$mu_clust
      clustering_labels <- relab$clustering_labels
      lambda <- relab$lambda
      sigma2_clust <- relab$sigma2_clust
    }
    
    last_iter <- iter
    
    if (iter > burn) {
      save_iter <- save_iter + 1
      
      Z_chain_raw[save_iter, , ] <- Z
      mu_chain_raw[save_iter, , ] <- mu_clust
      
      align_out <- align_one_clustering_draw(
        Z = Z,
        mu_clust = mu_clust,
        Z_ref = Z_ref,
        align_method = align_method,
        bary_tol = bary_tol,
        bary_max_iter = bary_max_iter,
        bary_drop_bad = bary_drop_bad,
        bary_tol_hyp = bary_tol_hyp,
        allow_reflection = TRUE,
        repair = TRUE
      )
      
      Z_store <- align_out$Z_store
      mu_store <- align_out$mu_store
      proc_angle[save_iter] <- align_out$proc_angle
      
      labels_store <- clustering_labels
      lambda_store <- lambda
      sigma2_store <- sigma2_clust
      perm_store <- seq_len(K)
      
      if (relabel_reference) {
        if (is.null(mu_ref_relab)) {
          mu_ref_relab <- mu_store
          labels_ref_relab <- labels_store
          lambda_ref_relab <- lambda_store
        } else {
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
        }
      }
      
      alpha_chain[save_iter] <- alpha
      Z_chain[save_iter, , ] <- Z_store
      mu_chain[save_iter, , ] <- mu_store
      lambda_chain[save_iter, ] <- lambda_store
      sigma2_chain[save_iter, ] <- sigma2_store
      clustering_chain[save_iter, ] <- labels_store
      relabel_perm_chain[save_iter, ] <- perm_store
      
      alpha_accept[save_iter] <- alpha_step$accepted
      z_accept[save_iter] <- z_acc_this_iter / n
      mu_accept[save_iter, ] <- mu_acc_this_iter
    }
  }
  
  if (save_iter > 0) {
    alpha_chain <- alpha_chain[1:save_iter]
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
    proc_angle <- proc_angle[1:save_iter]
    
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
  } else {
    alpha_mean <- NA_real_
    alpha_median <- NA_real_
    lambda_mean <- rep(NA_real_, K)
    sigma2_mean <- rep(NA_real_, K)
    label_posterior_probs <- matrix(NA_real_, nrow = n, ncol = K)
    clustering_map <- rep(NA_integer_, n)
  }
  
  Z_bary_mean <- NULL
  Z_polar_mean <- NULL
  mu_bary_mean <- NULL
  
  if (compute_barycenter && save_iter > 0) {
    Z_bary_mean <- posterior_barycenters(
      Z_chain = Z_chain,
      tol = bary_tol,
      max_iter = bary_max_iter,
      drop_bad = bary_drop_bad,
      tol_hyp = bary_tol_hyp
    )
  }
  
  if (compute_polar_mean && save_iter > 0) {
    Z_polar_mean <- posterior_polar_means(Z_chain)
  }
  
  if (compute_barycenter && save_iter > 0) {
    mu_bary_mean <- matrix(NA_real_, nrow = K, ncol = d)
    
    for (g in 1:K) {
      mu_chain_g <- mu_chain[, g, , drop = FALSE]
      dim(mu_chain_g) <- c(save_iter, 1, d)
      
      mu_bary_mean[g, ] <- posterior_barycenters(
        Z_chain = mu_chain_g,
        tol = bary_tol,
        max_iter = bary_max_iter,
        drop_bad = bary_drop_bad,
        tol_hyp = bary_tol_hyp
      )[1, ]
    }
  }
  
  return(list(
    alpha_chain = alpha_chain,
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
    proc_angle = proc_angle,
    align_method = align_method,
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
    Z_polar_mean = Z_polar_mean,
    mu_bary_mean = mu_bary_mean,
    sigma2_clust = sigma2_clust
  ))
}