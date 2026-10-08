parallel_transport <- function(v, mu, x,
                               eps = 1e-10,
                               tol_hyp = 1e-6,
                               tol_tangent = 1e-8,
                               tol_same = 1e-8
                               ) {
  v  <- as.numeric(v)
  mu <- as.numeric(mu)
  x  <- as.numeric(x)
  tan_err <- lorentz_inner(v, x)

  if (abs(tan_err) > tol_tangent) {
    x <- project_to_tangent(v, x)
  }
  
  if (max(abs(v - mu)) < tol_same) {
    y <- x
    y <- project_to_tangent(mu, y)
    return(y)
  }
  
  alpha <- -lorentz_inner(v, mu)
  alpha <- max(alpha, 1)
  
  denom <- alpha + 1
  beta <- mu - alpha * v
  coeff_num <- lorentz_inner(beta, x)
  y <- x + (coeff_num / denom) * (v + mu)
  y <- project_to_tangent(mu, y)
  
  y
}


project_to_tangent <- function(mu, x) {
  x + lorentz_inner(mu, x) * mu
}

project_to_hyperboloid <- function(z, eps = 1e-12) {
  z <- as.numeric(z)
  spatial <- z[-1]
  t_new <- sqrt(max(1 + sum(spatial^2), eps))
  
  c(t_new, spatial)
}

exp_map_hyp <- function(mu, x,
                        eps = 1e-10,
                        tol_tangent = 1e-8,
                        tol_hyp = 1e-5
                        ) {
  mu <- as.numeric(mu)
  x  <- as.numeric(x)
  
  mu_norm <- lorentz_inner(mu, mu)
  if (abs(mu_norm + 1) > tol_hyp) {
    mu <- project_to_hyperboloid(mu)
    mu_norm <- lorentz_inner(mu, mu)
  } 
  else if (abs(mu_norm + 1) > 1e-12) {
    mu <- project_to_hyperboloid(mu)
  }
  
  tan_err <- lorentz_inner(mu, x)
  if (abs(tan_err) > tol_tangent) {
    x <- project_to_tangent(mu, x)
  }
  
  x_sqnorm <- lorentz_inner(x, x)
  x_sqnorm <- max(x_sqnorm, 0)
  x_norm <- sqrt(x_sqnorm)
  
  if (x_norm < 1e-8) {
    z <- mu + x
  } else {
    cosh_fac <- cosh(x_norm)
    sinh_over_x <- sinh(x_norm) / x_norm
    z <- cosh_fac * mu + sinh_over_x * x
  }
  z <- project_to_hyperboloid(z)
  
  z
}



inv_exp_map <- function(z, mu, eps = 1e-10, tol_close = 1e-8) {
  alpha <- -lorentz_inner(mu, z)
  
  if ((alpha - 1) < tol_close) {
    v <- z - mu
    v <- v + lorentz_inner(mu, v) * mu
    return(v)
  }
  
  if (alpha < 1) {
    alpha <- 1
  }
  
  denom <- sqrt(alpha^2 - 1)
  
  if (!is.finite(denom) || denom < eps) {
    v <- z - mu
    v <- v + lorentz_inner(mu, v) * mu
    return(v)
  }
  
  coeff <- acosh_safe(alpha) / denom
  v <- coeff * (z - alpha * mu)

  v + lorentz_inner(mu, v) * mu
}


lorentz_inner <- function(X, Y) {
  X <- as.numeric(X)
  Y <- as.numeric(Y)
  
  -X[1] * Y[1] + sum(X[-1] * Y[-1])
}

acosh_safe <- function(x, eps = 1e-10) {
  x <- pmax(x, 1)
  log(x + sqrt(pmax(x^2 - 1, 0)))
}


hyp_dist_lorentz <- function(X, Y, K = 1) {
  ip <- lorentz_inner(X, Y)
  acosh_safe(-ip) / sqrt(K)
}

hyp_dist_to_indices <- function(z_i, Z, idx) {
  if (length(idx) == 0L) {
    return(numeric(0L))
  }
  
  Z_sub <- Z[idx, , drop = FALSE]
  
  lorentz_products <-
    -Z_sub[, 1L] * z_i[1L] +
    rowSums(
      Z_sub[, -1L, drop = FALSE] *
        matrix(
          z_i[-1L],
          nrow = length(idx),
          ncol = length(z_i) - 1L,
          byrow = TRUE
        )
    )
  
  acosh(pmax(-lorentz_products, 1))
}

D_t <- function(Z){
  Z_time <- Z[,1]
  Z_spatial <- Z[,-1]
  
  IP <- -tcrossprod(Z_time) + tcrossprod(Z_spatial)
  D_t <- acosh_safe(-IP)
  return(D_t)
}

polar_to_hyperboloid <- function(r, u_mat, R = 1) {
  r <- as.numeric(r)
  u_mat <- as.matrix(u_mat)
  
  n <- nrow(u_mat)
  d <- ncol(u_mat)
  row_norms <- sqrt(rowSums(u_mat^2))
  
  u_mat <- u_mat / row_norms
  
  x0 <- R * cosh(r / R)
  spatial <- R * sinh(r / R) * u_mat
  
  Z <- cbind(x0, spatial)
  Z <- matrix(Z, nrow = n, ncol = d + 1)
  
  return(Z)
}



hyperboloid_to_polar <- function(x, R = 1, tol = 1e-10) {
  x <- as.numeric(x)
  x0 <- x[1]
  spatial <- x[-1]
  d_spatial <- length(spatial)
  r <- R * acosh(x0 / R)
  
  if (abs(r) < tol) {
    directional <- c(1, rep(0, d_spatial - 1))
  } else {
    denom <- R * sinh(r / R)
    directional <- spatial / denom
    
    dir_norm <- sqrt(sum(directional^2))
    directional <- directional / dir_norm
  }
  
  out <- list(
    r = r,
    directional = directional
  )
  
  if (d_spatial == 2) {
    out$theta <- atan2(directional[2], directional[1])
  }
  
  return(out)
}




hyperboloid_to_polar_matrix <- function(X, R = 1, tol = 1e-10) {
  X <- as.matrix(X)
  n <- nrow(X)
  d_spatial <- ncol(X) - 1
  
  r <- numeric(n)
  directional <- matrix(NA_real_, nrow = n, ncol = d_spatial)
  
  if (d_spatial == 2) {
    theta <- numeric(n)
  }
  
  for (i in 1:n) {
    out_i <- hyperboloid_to_polar(X[i, ], R = R, tol = tol)
    
    r[i] <- out_i$r
    directional[i, ] <- out_i$directional
    
    if (d_spatial == 2) {
      theta[i] <- out_i$theta
    }
  }
  
  out <- list(
    r = r,
    directional = directional
  )
  
  if (d_spatial == 2) {
    out$theta <- theta
  }
  
  return(out)
}

inv_exp_map_batch <- function(
    Z,
    mu,
    eps = 1e-10,
    tol_close = 1e-8
) {
  Z <- as.matrix(Z)
  mu <- as.numeric(mu)
  alpha <-
    Z[, 1L] * mu[1L] -
    drop(
      Z[, -1L, drop = FALSE] %*% mu[-1L]
    )
  
  alpha <- pmax(alpha, 1)
  denom <- sqrt(
    pmax(alpha^2 - 1, 0)
  )
  
  coeff <- numeric(length(alpha))
  
  close <-
    (alpha - 1 < tol_close) |
    !is.finite(denom) |
    denom < eps
    coeff[close] <- 1
  
  coeff[!close] <-
    acosh_safe(alpha[!close]) /
    denom[!close]
    V <- Z - alpha %o% mu
    V <- sweep(
    V,
    1L,
    coeff,
    "*"
  )
  
  V
}

barycenter_one_node <- function(
    samples_node,
    tol = 1e-6,
    max_iter = 100,
    tol_hyp = 1e-5,
    step_shrink = 0.5,
    min_step = 1 / 128
) {
  
  samples_node <- as.matrix(samples_node)
  spatial <- samples_node[, -1L, drop = FALSE]
  
  samples_node[, 1L] <-
    sqrt(1 + rowSums(spatial^2))
  
  mu <- samples_node[1L, ]
  mu <- project_to_hyperboloid(mu)
  
  for (iter in seq_len(max_iter)) {
    
    mu <- project_to_hyperboloid(mu)
    V <- inv_exp_map_batch(
      Z = samples_node,
      mu = mu
    )
    
    v_bar <- colMeans(V)
    
    v_bar <- project_to_tangent(
      mu,
      v_bar
    )
    
    step_sq <- lorentz_inner(
      v_bar,
      v_bar
    )
    
    step_norm <- sqrt(
      max(step_sq, 0)
    )
    
    if (step_norm < tol) {
      break
    }
    
    step_scale <- 1
    moved <- FALSE
    
    while (step_scale >= min_step) {
      
      mu_try <- try(
        exp_map_hyp(
          mu = mu,
          x = step_scale * v_bar,
          tol_hyp = tol_hyp
        ),
        silent = TRUE
      )
      
      if (
        !inherits(mu_try, "try-error") &&
        all(is.finite(mu_try))
      ) {
        
        mu <- project_to_hyperboloid(
          mu_try
        )
        
        moved <- TRUE
        break
      }
      
      step_scale <-
        step_scale * step_shrink
    }
    
    if (!moved) {
      stop(
        "Could not take stable barycenter step at iteration ",
        iter
      )
    }
  }
  
  mu
}

barycenters <- function(
    Z_chain,
    tol = 1e-6,
    max_iter = 100,
    tol_hyp = 1e-5
) {
  
  M <- dim(Z_chain)[1L]
  n <- dim(Z_chain)[2L]
  d <- dim(Z_chain)[3L]
  
  Z_mean <- matrix(
    NA_real_,
    nrow = n,
    ncol = d
  )
  
  for (i in seq_len(n)) {
    
    samples_node <- matrix(
      Z_chain[, i, ],
      nrow = M,
      ncol = d
    )
    
    Z_mean[i, ] <-
      barycenter_one_node(
        samples_node = samples_node,
        tol = tol,
        max_iter = max_iter,
        tol_hyp = tol_hyp
      )
  }
  
  Z_mean
}





G_to_Z_MDS <- function(
    D,
    dim = 2,
    curvature = 1,
    alpha = 1.1,
    equi.adj = 0,
    radial_scale = 1
) {
  initial_Z_MDS_data <- hydra::hydraPlus(
    D,
    dim = dim,
    curvature = curvature,
    alpha = alpha,
    equi.adj = equi.adj
  )
  
  r_init <- as.numeric(initial_Z_MDS_data$r)
  directional_init <- as.matrix(initial_Z_MDS_data$directional)
  if (nrow(directional_init) != length(r_init) &&
      ncol(directional_init) == length(r_init)) {
    directional_init <- t(directional_init)
  }
  r_scaled <- radial_scale * r_init
  
  Z_init <- polar_to_hyperboloid(r_scaled, directional_init, R = 1)
  Z_init
}

initialise_hyperbolic_clusters <- function(
    Z,
    K,
    mu_clust_init = NULL,
    clustering_labels_init = NULL,
    sigma2_clust_init = NULL,
    lambda_init = NULL,
    mu_prior_centre = c(1, 0, 0),
    nstart = 25L,
    sigma2_min = 1e-4,
    lambda_pseudocount = 1
) {
  
  Z <- as.matrix(Z)
  
  n <- nrow(Z)
  d <- ncol(Z)
  p <- d - 1L
  
  if (K < 1L || K > n) {
    stop("K must be between 1 and nrow(Z).")
  }
  
  if (K == 1L) {
    
    mu_clust <- matrix(
      mu_prior_centre,
      nrow = 1L,
      byrow = TRUE
    )
    
    clustering_labels <- rep(1L, n)
    lambda <- 1
    
    if (is.null(sigma2_clust_init)) {
      
      ssq <- cluster_tangent_ssq(
        Z = Z,
        mu_g = mu_clust[1L, ],
        idx_g = seq_len(n)
      )
      
      sigma2_clust <- max(
        ssq / (n * p),
        sigma2_min
      )
      
    } else {
      
      sigma2_clust <- as.numeric(sigma2_clust_init)[1L]
    }
    
    return(
      list(
        mu_clust = mu_clust,
        clustering_labels = clustering_labels,
        sigma2_clust = sigma2_clust,
        lambda = lambda,
        cluster_sizes = n,
        method = "fixed_origin"
      )
    )
  }
  
  kmeans_out <- NULL
  
  if (is.null(mu_clust_init)) {
    
    origin <- c(1, rep(0, p))
    
   tangent_ambient <- t(
      vapply(
        seq_len(n),
        function(i) {
          inv_exp_map(
            z = Z[i, ],
            mu = origin
          )
        },
        numeric(d)
      )
    )
    
    tangent_coordinates <-
      tangent_ambient[, -1L, drop = FALSE]
    
    kmeans_out <- stats::kmeans(
      x = tangent_coordinates,
      centers = K,
      nstart = nstart
    )
    
    tangent_centres <- cbind(
      0,
      kmeans_out$centers
    )
    
    mu_clust <- t(
      vapply(
        seq_len(K),
        function(g) {
          exp_map_hyp(
            mu = origin,
            x = tangent_centres[g, ]
          )
        },
        numeric(d)
      )
    )
    
  } else {
    
    mu_clust <- as.matrix(mu_clust_init)
    
    if (!all(dim(mu_clust) == c(K, d))) {
      stop("mu_clust_init must be a K x ncol(Z) matrix.")
    }
    
    for (g in seq_len(K)) {
      mu_clust[g, ] <-
        project_to_hyperboloid(mu_clust[g, ])
    }
  }

  if (!is.null(clustering_labels_init)) {
    
    clustering_labels <-
      as.integer(clustering_labels_init)
    
  } else if (!is.null(kmeans_out)) {
    
    clustering_labels <-
      as.integer(kmeans_out$cluster)
    
  } else {

    clustering_labels <- initial_labels_clust(
      Z_init = Z,
      mu_clust_init = mu_clust
    )
  }
  
  
  if (length(clustering_labels) != n) {
    stop("clustering_labels_init must have length nrow(Z).")
  }
  
  if (any(!clustering_labels %in% seq_len(K))) {
    stop("Cluster labels must be integers in 1:K.")
  }
  
  
  cluster_sizes <- tabulate(
    clustering_labels,
    nbins = K
  )
  
  if (!is.null(sigma2_clust_init)) {
    
    sigma2_clust <- as_sigma_vector(
      sigma2_clust_init,
      K
    )
    
  } else {
    sigma2_clust <- numeric(K)
    for (g in seq_len(K)) {
      
      idx_g <- which(
        clustering_labels == g
      )
      
      n_g <- length(idx_g)
      
      if (n_g > 0L) {
        
        ssq_g <- cluster_tangent_ssq(
          Z = Z,
          mu_g = mu_clust[g, ],
          idx_g = idx_g
        )
        
        sigma2_clust[g] <-
          ssq_g / (p * n_g)
        
      } else {
        
        sigma2_clust[g] <- NA_real_
      }
    }
    
    valid <-
      is.finite(sigma2_clust) &
      sigma2_clust > sigma2_min
    
    fallback <- if (any(valid)) {
      median(sigma2_clust[valid])
    } else {
      1
    }
    
    sigma2_clust[!valid] <- fallback
    
    sigma2_clust <- pmax(
      sigma2_clust,
      sigma2_min
    )
  }
  
  if (!is.null(lambda_init)) {
    
    lambda <- normalise_lambda(
      lambda_init
    )
    
  } else {
    
    lambda <-
      cluster_sizes +
      lambda_pseudocount
    
    lambda <- lambda / sum(lambda)
  }
  
  if (
    length(sigma2_clust) != K ||
    any(!is.finite(sigma2_clust)) ||
    any(sigma2_clust <= 0)
  ) {
    stop("Initial cluster variances must be finite and positive.")
  }
  
  if (
    length(lambda) != K ||
    any(!is.finite(lambda)) ||
    any(lambda < 0)
  ) {
    stop("Initial lambda is invalid.")
  }
  
  
  list(
    mu_clust = mu_clust,
    clustering_labels = clustering_labels,
    sigma2_clust = sigma2_clust,
    lambda = lambda,
    cluster_sizes = cluster_sizes,
    kmeans_object = kmeans_out,
    method = if (is.null(mu_clust_init)) {
      "tangent_kmeans"
    } else {
      "manual_centres"
    }
  )
}


poincare_from_hyperboloid <- function(Z, eps = 1e-10) {
  Z <- as.matrix(Z)
  
  if (ncol(Z) != 3) {
    stop("Z must be an n x 3 matrix.")
  }
  
  denom <- 1 + Z[, 1]
  
  if (any(!is.finite(denom)) || any(abs(denom) < eps)) {
    stop("Invalid hyperboloid points for Poincare projection.")
  }
  
  cbind(
    Z[, 2] / denom,
    Z[, 3] / denom
  )
}

