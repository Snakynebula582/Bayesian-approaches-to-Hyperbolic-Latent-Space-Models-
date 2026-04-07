library(MASS)
library(hydra)
library(latentnet)
library(intergraph)
library(igraph)
library(rgl)
library(plotly)

parallel_transport <- function(v, mu, x,
                               eps = 1e-10,
                               tol_hyp = 1e-6,
                               tol_tangent = 1e-8,
                               tol_same = 1e-8,
                               repair = TRUE) {
  v  <- as.numeric(v)
  mu <- as.numeric(mu)
  x  <- as.numeric(x)
  
  
  tan_err <- lorentz_inner(v, x)
  if (abs(tan_err) > tol_tangent) {
    if (!repair) {
      stop(paste("parallel_transport x is not tangent at v"))
    }
    x <- project_to_tangent(v, x)
  }
  
  if (max(abs(v - mu)) < tol_same) {
    y <- x
    y <- project_to_tangent(mu, y)
    return(y)
  }
  
  alpha <- -lorentz_inner(v, mu)
  if (alpha < 1 - tol_hyp) {
    stop(paste("parallel_transport invalid alpha < 1"))
  }
  alpha <- max(alpha, 1)
  denom <- alpha + 1
  beta <- mu - alpha * v    
  coeff_num <- lorentz_inner(beta, x)
  y <- x + (coeff_num / denom) * (v + mu)
  y <- project_to_tangent(mu, y)
  
  out_tan_err <- lorentz_inner(mu, y)
  if (abs(out_tan_err) > 1e-6 && !repair) {
    stop(paste("parallel_transport transported vector not tangent at mu"))
  }
  
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
                        tol_hyp = 1e-5,
                        repair = TRUE) {
  mu <- as.numeric(mu)
  x  <- as.numeric(x)
  
  mu_norm <- lorentz_inner(mu, mu)
  
  if (abs(mu_norm + 1) > tol_hyp) {
    if (!repair) {
      stop(paste("exp_map_hyp mu is not on hyperboloid"))
    }
    mu <- project_to_hyperboloid(mu)
    mu_norm <- lorentz_inner(mu, mu)
  } else if (repair && abs(mu_norm + 1) > 1e-12) {
    mu <- project_to_hyperboloid(mu)
  }
  
  tan_err <- lorentz_inner(mu, x)
  
  if (abs(tan_err) > tol_tangent) {
    if (!repair) {
      stop(paste("exp_map_hyp x is not tangent"))
    }
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
  
  z_norm <- lorentz_inner(z, z)
  if (repair) {
    z <- project_to_hyperboloid(z)
  } else {
    if (abs(z_norm + 1) > 1e-4) {
      stop(paste("exp_map_hyp: output drift too large"))
    }
  }
  
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

G_to_Z_MDS <- function(
    D,
    dim = 2,
    curvature = 1,
    alpha = 1.1,
    equi.adj = 0,
    radial_scale = 1
) {
  initial_Z_MDS_data <- hydraPlus(
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


posterior_barycenter_one_node <- function(samples_node, mu_init = NULL,
                                          tol = 1e-6,
                                          max_iter = 100,
                                          drop_bad = TRUE,
                                          tol_hyp = 1e-5,
                                          step_shrink = 0.5,
                                          min_step = 1/128) {
  samples_node <- as.matrix(samples_node)
  
  if (nrow(samples_node) == 0) {
    stop("posterior_barycenter_one_node: empty samples_node")
  }
  
  if (drop_bad) {
    keep <- apply(samples_node, 1, function(x) {
      all(is.finite(x)) && is.finite(lorentz_inner(x, x))
    })
    samples_node <- samples_node[keep, , drop = FALSE]
  }
  
  M <- nrow(samples_node)
  
  for (m in 1:M) {
    nn <- lorentz_inner(samples_node[m, ], samples_node[m, ])
    samples_node[m, ] <- project_to_hyperboloid(samples_node[m, ])
  }
  
  if (is.null(mu_init)) {
    mu <- samples_node[1, ]
  } else {
    mu <- as.numeric(mu_init)
    mu <- project_to_hyperboloid(mu)
  }
  
  for (iter in 1:max_iter) {
    mu <- project_to_hyperboloid(mu)
    V <- matrix(NA_real_, nrow = M, ncol = ncol(samples_node))
    
    for (m in 1:M) {
      v_m <- inv_exp_map(z = samples_node[m, ], mu = mu)
      v_m <- project_to_tangent(mu, v_m)
      V[m, ] <- v_m
    }
    
    v_bar <- colMeans(V)
    v_bar <- project_to_tangent(mu, v_bar)
    
    li <- lorentz_inner(v_bar, v_bar)
    
    if (li < -tol_hyp) {
      stop(paste("barycenter step has strongly negative tangent norm"))
    }
    
    li <- max(li, 0)
    step_norm <- sqrt(li)
    
    if (step_norm < tol) {
      break
    }
    
    step_scale <- 1
    moved <- FALSE
    
    while (step_scale >= min_step) {
      mu_try <- try(
        exp_map_hyp(mu = mu,
                    x = step_scale * v_bar,
                    tol_hyp = tol_hyp,
                    repair = TRUE),
        silent = TRUE
      )
      
      if (!inherits(mu_try, "try-error") && all(is.finite(mu_try))) {
        mu <- project_to_hyperboloid(mu_try)
        moved <- TRUE
        break
      }
      
      step_scale <- step_scale * step_shrink
    }
    
    if (!moved) {
      stop(paste("posterior_barycenter_one_node: could not take stable step at iteration", iter))
    }
  }
  
  mu
}


posterior_barycenters <- function(Z_chain, tol = 1e-6, max_iter = 100,
                                  drop_bad = TRUE, tol_hyp = 1e-5) {
  M <- dim(Z_chain)[1]
  n <- dim(Z_chain)[2]
  d <- dim(Z_chain)[3]
  
  Z_mean <- matrix(NA_real_, nrow = n, ncol = d)
  
  for (i in 1:n) {
    samples_node <- Z_chain[, i, , drop = FALSE]
    samples_node <- matrix(samples_node, nrow = M, ncol = d)
    
    Z_mean[i, ] <- posterior_barycenter_one_node(
      samples_node = samples_node,
      mu_init = samples_node[1, ],
      tol = tol,
      max_iter = max_iter,
      drop_bad = drop_bad,
      tol_hyp = tol_hyp
    )
  }
  
  Z_mean
}


circular_mean <- function(theta) {
  s <- mean(sin(theta))
  c <- mean(cos(theta))
  atan2(s, c)
}


posterior_polar_mean_one_node <- function(samples_node) {
  samples_node <- as.matrix(samples_node)
  polar <- hyperboloid_to_polar_matrix(samples_node)
  
  r <- polar$r
  
  if (!is.null(polar$theta)) {
    theta <- polar$theta
  } else if (!is.null(polar$directional)) {
    directional <- polar$directional
    
    if (is.vector(directional)) {
      theta <- directional
    } else if (is.matrix(directional) && ncol(directional) == 2) {
      theta <- atan2(directional[, 2], directional[, 1])
    } 
  } 
  
  r_mean <- mean(r)
  theta_mean <- circular_mean(theta)
  
  dir_mean <- matrix(c(cos(theta_mean), sin(theta_mean)), nrow = 1)
  z_mean <- polar_to_hyperboloid(r_mean, dir_mean)
  
  as.numeric(z_mean)
}


posterior_polar_means <- function(Z_chain) {
  dims <- dim(Z_chain)  
  M <- dims[1]
  n <- dims[2]
  d <- dims[3]  
  Z_mean <- matrix(NA_real_, nrow = n, ncol = d)
  
  for (i in 1:n) {
    samples_node <- Z_chain[, i, , drop = FALSE]
    samples_node <- matrix(samples_node, nrow = M, ncol = d)
    Z_mean[i, ] <- posterior_polar_mean_one_node(samples_node)
  }
  
  Z_mean
}

fit_ru_to_hyperboloid <- function(fit, draw = 1, use_mean = FALSE, R = 1) {
  r_array <- fit$r
  U_array <- fit$U  
  n_keep <- nrow(r_array)
  n <- ncol(r_array)
  
  if (use_mean) {
    r <- colMeans(r_array)
    U <- apply(U_array, c(2, 3), mean)
    
    row_norms <- sqrt(rowSums(U^2))
    U <- U / row_norms
  } else {    
    r <- as.numeric(r_array[draw, ])
    U <- U_array[draw, , ]
    
    if (is.null(dim(U))) {
      U <- matrix(U, ncol = dim(U_array)[3])
    }
  }
  
  x0 <- R * cosh(r / R)
  spatial <- R * sinh(r / R) * U
  
  Z <- cbind(x0, spatial)
  colnames(Z) <- c("x0", paste0("x", 1:ncol(U)))
  
  return(Z)
}

PHN_sampling <- function(mu, sigma2){
  dimension <- length(mu)
  mu0 <- c(1, rep(0, dimension - 1))
  
  samp <- mvrnorm(
    1,
    mu = rep(0, dimension - 1),
    Sigma = sigma2 * diag(dimension - 1)
  )
  
  
  v <- c(0, samp)
  v_mu <- parallel_transport(mu0, mu, v)
  tangent_check <- lorentz_inner(mu, v_mu)
  if (abs(tangent_check) > 1e-5) {
    stop(paste("transported vector is not tangent"))
  }
  if (lorentz_inner(v_mu, v_mu) < -1e-8) stop("transported vector has negative Lorentz norm")
  
  z <- exp_map_hyp(mu, v_mu)
  
  return(z)
}

log_pdf_PHN <- function(z, mu, sigma2){
  dimension <- length(mu)
  n <- dimension - 1
  
  Sigma <- diag(rep(sigma2, n))
  mu0 <- c(1, rep(0, n))
  
  u <- inv_exp_map(z, mu)
  v <- parallel_transport(mu, mu0, u)
  v_hat <- v[-1]
  
  r <- sqrt(max(as.numeric(lorentz_inner(u, u)), 0))
  
  lg_pv <- -(n / 2) * log(2 * pi) -
    0.5 * log(det(Sigma)) -
    0.5 * as.numeric(t(v_hat) %*% solve(Sigma, v_hat))
  
  lg_jac <- ifelse(abs(r) < 1e-10, 0, (n - 1) * log(sinh(r) / r))
  
  lg_pz <- lg_pv - lg_jac
  return(lg_pz)
}


Full_likelihood_initial <- function(Z, Y, alpha) {
  n <- nrow(Z)
  total <- 0
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      d_ij <- hyp_dist_lorentz(Z[i, ], Z[j, ])
      eta_ij <- alpha - d_ij
      total <- total + Y[i, j] * eta_ij - log1p(exp(eta_ij))
    }
  }
  
  total
}


loglik_i <- function(i, z_i, Z, Y, alpha) {
  n <- nrow(Z)
  total <- 0
  
  for (j in seq_len(n)) {
    if (j == i) next
    y_ij <- if (i < j) Y[i, j] else Y[j, i]
    
    d_ij <- hyp_dist_lorentz(z_i, Z[j, ])
    eta_ij <- alpha - d_ij
    
    total <- total + y_ij * eta_ij - log1p(exp(eta_ij))
  }
  
  total
}

procrustes_transform <- function(Z_new, Z_ref, renormalise = TRUE) {
  t_new <- Z_new[, 1]
  X_new <- Z_new[, 2:3, drop = FALSE]
  X_ref <- Z_ref[, 2:3, drop = FALSE]
  
  theta_new <- atan2(X_new[, 2], X_new[, 1])
  theta_ref <- atan2(X_ref[, 2], X_ref[, 1])
  
  alpha <- Arg(sum(exp(1i * (theta_ref - theta_new))))
  
  R <- matrix(c(cos(alpha), -sin(alpha),
                sin(alpha),  cos(alpha)),
              nrow = 2, byrow = TRUE)
  
  X_aligned <- X_new %*% R
  
  if (renormalise) {
    t_aligned <- sqrt(1 + rowSums(X_aligned^2))
  } else {
    t_aligned <- t_new
  }
  
  Z_aligned <- cbind(t_aligned, X_aligned)
  
  list(
    Z_aligned = Z_aligned,
    alpha = alpha,
    R = R
  )
}


all_permutations_int <- function(x) {
  x <- as.integer(x)
  
  if (length(x) == 1L) {
    return(matrix(x, nrow = 1L))
  }
  
  out <- vector("list", length(x))
  
  for (i in seq_along(x)) {
    rest <- x[-i]
    rest_perm <- all_permutations_int(rest)
    out[[i]] <- cbind(x[i], rest_perm)
  }
  
  do.call(rbind, out)
}


match_clusters_to_reference <- function(mu_curr, mu_ref, w_mu = 1) {
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Package 'clue' is required for cluster matching")
  }
  
  mu_curr <- as.matrix(mu_curr)
  mu_ref  <- as.matrix(mu_ref)
  
  K <- nrow(mu_curr)  
  C <- matrix(0, nrow = K, ncol = K)
  
  for (g in 1:K) {
    for (h in 1:K) {
      C[g, h] <- w_mu * hyp_dist_lorentz(mu_curr[g, ], mu_ref[h, ])^2
    }
  }
  
  C <- C - min(C)
  
  as.integer(clue::solve_LSAP(C))
}


apply_reference_relabel <- function(
    perm,
    mu_clust,
    clustering_labels,
    lambda = NULL
) {
  
  ord <- order(perm)
  
  out <- list(
    mu_clust = mu_clust[ord, , drop = FALSE],
    clustering_labels = perm[clustering_labels],
    ord = ord
  )
  
  if (!is.null(lambda)) {
    out$lambda <- lambda[ord]
  }
  
  out
}

HPA_translation <- function(Z_ref_bar, Z_new,
                            mu_init_new = NULL,
                            bary_tol = 1e-6,
                            bary_max_iter = 100,
                            drop_bad = TRUE,
                            tol_hyp = 1e-5,
                            repair = TRUE) {
  Z_new <- as.matrix(Z_new)
  Z_ref_bar <- as.numeric(Z_ref_bar)
  
  if (repair) {
    Z_ref_bar <- project_to_hyperboloid(Z_ref_bar)
    for (i in seq_len(nrow(Z_new))) {
      Z_new[i, ] <- project_to_hyperboloid(Z_new[i, ])
    }
  }
  
  if (is.null(mu_init_new)) {
    mu_init_new <- Z_new[1, ]
  }
  
  Z_new_bar <- posterior_barycenter_one_node(
    samples_node = Z_new,
    mu_init = mu_init_new,
    tol = bary_tol,
    max_iter = bary_max_iter,
    drop_bad = drop_bad,
    tol_hyp = tol_hyp
  )
  
  if (repair) {
    Z_new_bar <- project_to_hyperboloid(Z_new_bar)
  }
  
  alpha_val <- -lorentz_inner(Z_ref_bar, Z_new_bar)
  alpha_val <- max(alpha_val, 1)
  denom <- alpha_val + 1
  
  a_plus_b_over <- (Z_ref_bar + Z_new_bar) / denom
  lam_vec <- (Z_ref_bar - (2 * alpha_val + 1) * Z_new_bar) / denom
  
  npoints <- nrow(Z_new)
  Z_new_aligned <- matrix(NA_real_, nrow = npoints, ncol = ncol(Z_new))
  
  for (i in seq_len(npoints)) {
    z <- Z_new[i, ]
    
    beta_i <- -lorentz_inner(a_plus_b_over, z)
    lambda_i <- lorentz_inner(lam_vec, z)
    
    z_aligned <- z - beta_i * Z_new_bar + lambda_i * Z_ref_bar
    
    if (repair) {
      z_aligned <- project_to_hyperboloid(z_aligned)
    }
    
    Z_new_aligned[i, ] <- z_aligned
  }
  
  Z_new_aligned
}


HPA_rotation <- function(Z_ref, Z_new, mu,
                         allow_reflection = TRUE,
                         tol_hyp = 1e-5,
                         tol_tangent = 1e-8,
                         repair = TRUE) {
  
  Z_ref <- as.matrix(Z_ref)
  Z_new <- as.matrix(Z_new)
  mu <- as.numeric(mu)
  mu <- project_to_hyperboloid(mu)
  
  for (i in seq_len(nrow(Z_ref))) {
    Z_ref[i, ] <- project_to_hyperboloid(Z_ref[i, ])
    Z_new[i, ] <- project_to_hyperboloid(Z_new[i, ])
  }
  
  P_mu <- function(v) {
    v <- as.numeric(v)
    as.numeric(v[-1])
  }
  
  P_mu_inv <- function(s, mu) {
    s <- as.numeric(s)
    v0 <- sum(mu[-1] * s) / mu[1]
    v <- c(v0, s)
    project_to_tangent(mu, v)
  }
  
  n <- nrow(Z_ref)
  V_ref <- matrix(NA_real_, nrow = n, ncol = 3)
  V_new <- matrix(NA_real_, nrow = n, ncol = 3)
  
  for (i in seq_len(n)) {
    v_ref_i <- inv_exp_map(z = Z_ref[i, ], mu = mu)
    v_new_i <- inv_exp_map(z = Z_new[i, ], mu = mu)
    
    v_ref_i <- project_to_tangent(mu, v_ref_i)
    v_new_i <- project_to_tangent(mu, v_new_i)
    V_ref[i, ] <- v_ref_i
    V_new[i, ] <- v_new_i
  }
  
  S_ref_raw <- t(apply(V_ref, 1, P_mu))
  S_new_raw <- t(apply(V_new, 1, P_mu))
  
  sbar_ref <- colMeans(S_ref_raw)
  sbar_new <- colMeans(S_new_raw)
  
  S_ref <- sweep(S_ref_raw, 2, sbar_ref, "-")
  S_new <- sweep(S_new_raw, 2, sbar_new, "-")
  
  C <- t(S_new) %*% S_ref
  sv <- svd(C)
  R <- sv$u %*% t(sv$v)
  
  if (!allow_reflection && det(R) < 0) {
    D_fix <- diag(c(1, -1))
    R <- sv$u %*% D_fix %*% t(sv$v)
  }
  
  Z_rotated <- matrix(NA_real_, nrow = n, ncol = 3)
  
  for (i in seq_len(n)) {
    s_i <- S_new_raw[i, ]
    
    s_i_rot <- as.numeric(R %*% (s_i - sbar_new)) + sbar_ref
    
    v_i_rot <- P_mu_inv(s_i_rot, mu)
    v_i_rot <- project_to_tangent(mu, v_i_rot)
    
    z_i_rot <- exp_map_hyp(mu = mu, x = v_i_rot,
                           tol_tangent = tol_tangent,
                           tol_hyp = tol_hyp,
                           repair = repair)
    
    Z_rotated[i, ] <- project_to_hyperboloid(z_i_rot)
  }
  
  list(
    Z_rotated = Z_rotated,
    R = R,
    S_ref = S_ref,
    S_new = S_new,
    sbar_ref = sbar_ref,
    sbar_new = sbar_new
  )
}

HPA_align_draw <- function(Z_ref, Z_new, Z_ref_bar = NULL,
                           allow_reflection = TRUE,
                           bary_tol = 1e-6,
                           bary_max_iter = 100,
                           drop_bad = TRUE,
                           tol_hyp = 1e-5,
                           tol_tangent = 1e-8,
                           repair = TRUE) {
  
  Z_ref <- as.matrix(Z_ref)
  Z_new <- as.matrix(Z_new)
  n <- nrow(Z_ref)
  
  if (repair) {
    for (i in seq_len(n)) {
      Z_ref[i, ] <- project_to_hyperboloid(Z_ref[i, ])
      Z_new[i, ] <- project_to_hyperboloid(Z_new[i, ])
    }
  }
  
  if (is.null(Z_ref_bar)) {
    Z_ref_bar <- posterior_barycenter_one_node(
      samples_node = Z_ref,
      mu_init = Z_ref[1, ],
      tol = bary_tol,
      max_iter = bary_max_iter,
      drop_bad = drop_bad,
      tol_hyp = tol_hyp
    )
  } else {
    Z_ref_bar <- as.numeric(Z_ref_bar)
    if (repair) {
      Z_ref_bar <- project_to_hyperboloid(Z_ref_bar)
    }
  }
  
  Z_trans <- HPA_translation(
    Z_ref_bar = Z_ref_bar,
    Z_new = Z_new,
    mu_init_new = Z_new[1, ],
    bary_tol = bary_tol,
    bary_max_iter = bary_max_iter,
    drop_bad = drop_bad,
    tol_hyp = tol_hyp,
    repair = repair
  )
  
  rot_out <- HPA_rotation(
    Z_ref = Z_ref,
    Z_new = Z_trans,
    mu = Z_ref_bar,
    allow_reflection = allow_reflection,
    tol_hyp = tol_hyp,
    tol_tangent = tol_tangent,
    repair = repair
  )
  
  Z_aligned <- rot_out$Z_rotated
  
  if (repair) {
    for (i in seq_len(n)) {
      Z_aligned[i, ] <- project_to_hyperboloid(Z_aligned[i, ])
    }
  }
  
  list(
    Z_aligned = Z_aligned,
    Z_trans = Z_trans,
    Z_ref_bar = Z_ref_bar,
    R = rot_out$R,
    S_ref = rot_out$S_ref,
    S_new = rot_out$S_new,
    sbar_ref = rot_out$sbar_ref,
    sbar_new = rot_out$sbar_new
  )
}



HPA_translation_apply <- function(X, Z_ref_bar, Z_new_bar,
                                  tol_hyp = 1e-5,
                                  repair = TRUE) {
  X <- as.matrix(X)
  Z_ref_bar <- as.numeric(Z_ref_bar)
  Z_new_bar <- as.numeric(Z_new_bar)
  
  if (repair) {
    Z_ref_bar <- project_to_hyperboloid(Z_ref_bar)
    Z_new_bar <- project_to_hyperboloid(Z_new_bar)
    for (i in seq_len(nrow(X))) {
      X[i, ] <- project_to_hyperboloid(X[i, ])
    }
  }
  
  alpha_val <- -lorentz_inner(Z_ref_bar, Z_new_bar)
  alpha_val <- max(alpha_val, 1)
  denom <- alpha_val + 1
  
  a_plus_b_over <- (Z_ref_bar + Z_new_bar) / denom
  lam_vec <- (Z_ref_bar - (2 * alpha_val + 1) * Z_new_bar) / denom
  
  X_out <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(X))
  
  for (i in seq_len(nrow(X))) {
    x <- X[i, ]
    
    beta_i <- -lorentz_inner(a_plus_b_over, x)
    lambda_i <- lorentz_inner(lam_vec, x)
    
    x_aligned <- x - beta_i * Z_new_bar + lambda_i * Z_ref_bar
    
    if (repair) {
      x_aligned <- project_to_hyperboloid(x_aligned)
    }
    
    X_out[i, ] <- x_aligned
  }
  
  X_out
}


HPA_rotation_apply <- function(X, mu, R, sbar_from, sbar_to,
                               tol_hyp = 1e-5,
                               tol_tangent = 1e-8,
                               repair = TRUE) {
  X <- as.matrix(X)
  mu <- as.numeric(mu)
  R <- as.matrix(R)
  mu <- project_to_hyperboloid(mu)
  
  P_mu <- function(v) {
    as.numeric(v[-1])
  }
  
  P_mu_inv <- function(s, mu) {
    v0 <- sum(mu[-1] * s) / mu[1]
    v <- c(v0, s)
    project_to_tangent(mu, v)
  }
  
  X_out <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(X))
  
  for (i in seq_len(nrow(X))) {
    x <- if (repair) project_to_hyperboloid(X[i, ]) else X[i, ]
    
    v_i <- inv_exp_map(z = x, mu = mu)
    v_i <- project_to_tangent(mu, v_i)
    
    s_i <- P_mu(v_i)
    s_i_rot <- as.numeric(R %*% (s_i - sbar_from)) + sbar_to
    
    v_i_rot <- P_mu_inv(s_i_rot, mu)
    v_i_rot <- project_to_tangent(mu, v_i_rot)
    
    x_rot <- exp_map_hyp(
      mu = mu,
      x = v_i_rot,
      tol_tangent = tol_tangent,
      tol_hyp = tol_hyp,
      repair = repair
    )
    
    X_out[i, ] <- project_to_hyperboloid(x_rot)
  }
  
  X_out
}


align_one_clustering_draw <- function(
    Z,
    mu_clust,
    Z_ref,
    align_method = c("none", "polar", "HPA"),
    bary_tol = 1e-6,
    bary_max_iter = 100,
    bary_drop_bad = TRUE,
    bary_tol_hyp = 1e-5,
    allow_reflection = TRUE,
    repair = TRUE
) {
  align_method <- match.arg(align_method)
  
  if (align_method == "none") {
    return(list(
      Z_store = Z,
      mu_store = mu_clust,
      proc_angle = NA_real_,
      proc_R = NULL,
      Z_ref_bar = NULL
    ))
  }
  
  if (align_method == "polar") {
    proc_out <- procrustes_transform(Z, Z_ref, renormalise = TRUE)
    
    return(list(
      Z_store = proc_out$Z_aligned,
      mu_store = apply_spatial_rotation_h2(mu_clust, proc_out$R, renormalise = TRUE),
      proc_angle = proc_out$alpha,
      proc_R = proc_out$R,
      Z_ref_bar = NULL
    ))
  }
  
  Z_ref_bar <- posterior_barycenter_one_node(
    samples_node = Z_ref,
    mu_init = Z_ref[1, ],
    tol = bary_tol,
    max_iter = bary_max_iter,
    drop_bad = bary_drop_bad,
    tol_hyp = bary_tol_hyp
  )
  Z_ref_bar <- project_to_hyperboloid(Z_ref_bar)
  
  hpa_out <- HPA_align_draw(
    Z_ref = Z_ref,
    Z_new = Z,
    Z_ref_bar = Z_ref_bar,
    allow_reflection = allow_reflection,
    bary_tol = bary_tol,
    bary_max_iter = bary_max_iter,
    drop_bad = bary_drop_bad,
    tol_hyp = bary_tol_hyp,
    tol_tangent = 1e-8,
    repair = repair
  )
  
  Z_bar <- posterior_barycenter_one_node(
    samples_node = Z,
    mu_init = Z[1, ],
    tol = bary_tol,
    max_iter = bary_max_iter,
    drop_bad = bary_drop_bad,
    tol_hyp = bary_tol_hyp
  )
  Z_bar <- project_to_hyperboloid(Z_bar)
  
  mu_trans <- HPA_translation_apply(
    X = mu_clust,
    Z_ref_bar = Z_ref_bar,
    Z_new_bar = Z_bar,
    tol_hyp = bary_tol_hyp,
    repair = repair
  )
  
  mu_store <- HPA_rotation_apply(
    X = mu_trans,
    mu = Z_ref_bar,
    R = hpa_out$R,
    sbar_from = hpa_out$sbar_new,
    sbar_to = hpa_out$sbar_ref,
    tol_hyp = bary_tol_hyp,
    tol_tangent = 1e-8,
    repair = repair
  )
  
  list(
    Z_store = hpa_out$Z_aligned,
    mu_store = mu_store,
    proc_angle = NA_real_,
    proc_R = hpa_out$R,
    Z_ref_bar = Z_ref_bar
  )
}


alpha_sample_noclust <- function(mu =1 , sigma2 =1){
  al_new <- rnorm(1, mean = mu, sd = sqrt(sigma2))
  prob_al_new <- dnorm(al_new, mean = mu, sd = sqrt(sigma2))
  return(list(al_new = al_new, prob_al_new = prob_al_new))
}

PHN_sampling_noclust <- function(z, sigma2,Y){
  z_new <- PHN_sampling(z, sigma2)
  prob_z_new <- log_pdf_PHN(z_new,z,sigma2)
  return(list(z_new = z_new, prob_z_new = prob_z_new))
}

alpha_update_step <- function(alpha, mu, sigma2, Z, Y, loglik_old, sigma2_prop){
  
  alpha_old <- alpha
  
  # random-walk proposal
  alpha_new <- rnorm(1, mean = alpha_old, sd = sqrt(sigma2_prop))
  
  # log prior densities
  logprior_old <- dnorm(alpha_old, mean = mu, sd = sqrt(sigma2), log = TRUE)
  logprior_new <- dnorm(alpha_new, mean = mu, sd = sqrt(sigma2), log = TRUE)
  
  # new log-likelihood
  loglik_new <- Full_likelihood_initial(Z, Y, alpha_new)
  
  # symmetric proposal, so proposal terms cancel
  log_ratio <- (loglik_new + logprior_new) -
    (loglik_old + logprior_old)
  
  # accept / reject
  if (is.finite(log_ratio) && log(runif(1)) < log_ratio) {
    alpha_hat <- alpha_new
    loglik_hat <- loglik_new
    accepted <- 1
  } else {
    alpha_hat <- alpha_old
    loglik_hat <- loglik_old
    accepted <- 0
  }
  
  return(list(
    alpha_hat = alpha_hat,
    loglik_hat = loglik_hat,
    accepted = accepted,
    log_ratio = log_ratio
  ))
}


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



Z_update_step_noclust <- function(
    Z, Y, alpha,
    sigma2_prior, mu_z,
    sigma_prop_z,
    i,
    loglik_old
){
  n <- nrow(Z)
  
  z_i <- as.numeric(Z[i, ])
  
  prop_out <- tryCatch(
    proposal_tangent_rw_h2(z_old = z_i, sigma_prop = sigma_prop_z),
    error = function(e) NULL
  )
  
  if (is.null(prop_out)) {
    return(list(
      Z = Z,
      accepted = 0,
      log_r = -Inf,
      loglik_hat = loglik_old
    ))
  }
  
  z_new <- prop_out$z_new
  
  if (!all(is.finite(z_new))) {
    return(list(
      Z = Z,
      accepted = 0,
      log_r = -Inf,
      loglik_hat = loglik_old
    ))
  }
  
  log_prior_old <- log_pdf_PHN(z_i,   mu_z, sigma2_prior)
  log_prior_new <- log_pdf_PHN(z_new, mu_z, sigma2_prior)
  
  log_like_old_i <- loglik_i(i, z_i,   Z, Y, alpha)
  log_like_new_i <- loglik_i(i, z_new, Z, Y, alpha)
  
  log_q_forward <- tryCatch(
    log_q_tangent_rw_h2(z_to = z_new, z_from = z_i, sigma_prop = sigma_prop_z),
    error = function(e) -Inf
  )
  
  log_q_reverse <- tryCatch(
    log_q_tangent_rw_h2(z_to = z_i, z_from = z_new, sigma_prop = sigma_prop_z),
    error = function(e) -Inf
  )
  
  log_r <- (log_like_new_i + log_prior_new + log_q_reverse) -
    (log_like_old_i + log_prior_old + log_q_forward)
  
  accepted <- 0
  loglik_hat <- loglik_old
  
  if (is.finite(log_r) && log(runif(1)) < log_r) {
    Z[i, ] <- z_new
    accepted <- 1
    loglik_hat <- loglik_old + (log_like_new_i - log_like_old_i)
  }
  
  list(
    Z = Z,
    accepted = accepted,
    log_r = log_r,
    loglik_hat = loglik_hat
  )
}


postprocess_hpa_hyp_lsm <- function(fit,
                                    r_ref,
                                    U_ref,
                                    allow_reflection = TRUE,
                                    bary_tol = 1e-6,
                                    bary_max_iter = 100,
                                    drop_bad = TRUE,
                                    tol_hyp = 1e-5,
                                    tol_tangent = 1e-8,
                                    R = 1) {
  r_draws <- fit$r
  U_draws <- fit$U
  n_keep <- nrow(r_draws)
  n <- ncol(r_draws)
  d <- dim(U_draws)[3]
  
  Z_ref <- polar_to_hyperboloid(r = r_ref, u_mat = U_ref, R = R)
  
  Z_ref_bar <- posterior_barycenter_one_node(
    samples_node = Z_ref,
    mu_init = Z_ref[1, ],
    tol = bary_tol,
    max_iter = bary_max_iter,
    drop_bad = drop_bad,
    tol_hyp = tol_hyp
  )
  Z_ref_bar <- project_to_hyperboloid(Z_ref_bar)
  
  Z_hpa <- array(NA_real_, dim = c(n_keep, n, d + 1))
  
  r_hpa <- matrix(NA_real_, nrow = n_keep, ncol = n)
  U_hpa <- array(NA_real_, dim = c(n_keep, n, d))
  
  for (m in 1:n_keep) {
    U_m <- U_draws[m, , ]
    if (is.null(dim(U_m))) {
      U_m <- matrix(U_m, nrow = n, ncol = d)
    }
    
    Z_m <- polar_to_hyperboloid(
      r = r_draws[m, ],
      u_mat = U_m,
      R = R
    )
    
    hpa_out <- HPA_align_draw(
      Z_ref = Z_ref,
      Z_new = Z_m,
      Z_ref_bar = Z_ref_bar,
      allow_reflection = allow_reflection,
      bary_tol = bary_tol,
      bary_max_iter = bary_max_iter,
      drop_bad = drop_bad,
      tol_hyp = tol_hyp,
      tol_tangent = tol_tangent,
      repair = TRUE
    )
    
    Z_hpa[m, , ] <- hpa_out$Z_aligned
    
    polar_m <- hyperboloid_to_polar_matrix(hpa_out$Z_aligned, R = R)
    r_hpa[m, ] <- polar_m$r
    U_hpa[m, , ] <- polar_m$directional
  }
  
  fit$Z_raw <- fit$Z
  fit$Z <- Z_hpa
  
  fit$r_raw <- fit$r
  fit$U_raw <- fit$U
  
  fit$r <- r_hpa
  fit$U <- U_hpa
  
  fit$alignment <- list(
    method = "hpa_postprocess",
    reference = "initialisation",
    Z_ref = Z_ref,
    Z_ref_bar = Z_ref_bar
  )
  
  fit
}


init_hyperbolic_mds_polar <- function(D, d = 2) {
  emb <- hydraPlus(
    D,
    dim = d,
    curvature = 1,
    alpha = 1.1
  )
  
  r <- emb$r
  U <- emb$directional
  
  list(r = r, U = U)
}

log_r_prior <- function(r) {
  dgamma(r, shape = 4, rate = 2, log = TRUE)
}

inv_logit <- function(x) 1 / (1 + exp(-x))


log_halfnorm <- function(r, sigma) {
  if (r < 0) return(-Inf)
  log(sqrt(2) / (sigma * sqrt(pi))) - (r^2) / (2 * sigma^2)
}

log_prior_K <- function(K, mu_log_K = 0, sd_log_K = 0.15) {
  if (K <= 0 || !is.finite(K)) return(-Inf)
  dnorm(log(K), mean = mu_log_K, sd = sd_log_K, log = TRUE) - log(K)
}

cosh_dist_ij <- function(ri, ui, rj, uj) {
  cosh(ri) * cosh(rj) - sinh(ri) * sinh(rj) * sum(ui * uj)
}

dist_ij <- function(ri, ui, rj, uj, K = 1) {
  if (!is.finite(K) || K <= 0) {
    stop("K must be a positive finite number.")
  }
  
  acosh_safe(cosh_dist_ij(ri, ui, rj, uj)) / sqrt(K)
}

loglik_pair <- function(yij, alpha, dij) {
  eta <- alpha - dij
  yij * eta - log1p(exp(eta))
}

loglik_full <- function(Y, alpha, r, U, K = 1) {
  n <- nrow(Y)
  ll <- 0
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      dij <- dist_ij(r[i], U[i, ], r[j], U[j, ], K = K)
      ll <- ll + loglik_pair(Y[i, j], alpha, dij)
    }
  }
  
  ll
}

loglik_node <- function(Y, alpha, r, U, i, K = 1) {
  n <- nrow(Y)
  ll <- 0
  
  for (j in 1:n) {
    if (j == i) next
    
    yij <- if (i < j) Y[i, j] else Y[j, i]
    dij <- dist_ij(r[i], U[i, ], r[j], U[j, ], K = K)
    ll <- ll + loglik_pair(yij, alpha, dij)
  }
  
  ll
}

align_directions <- function(U, U_ref) {
  C <- t(U_ref) %*% U
  sv <- svd(C)
  R <- sv$v %*% t(sv$u)
  U %*% R
}


propose_u <- function(u, step_u) {
  v <- u + rnorm(length(u), 0, step_u)
  v / sqrt(sum(v^2))
}

propose_r_reflect <- function(r, step_r) {
  rp <- r + rnorm(1, 0, step_r)
  if (rp < 0) rp <- -rp
  rp
}

propose_alpha <- function(alpha, step_a) {
  alpha + rnorm(1, 0, step_a)
}


propose_log_K <- function(log_K, step_log_K) {
  log_K + rnorm(1, 0, step_log_K)
}

wrap_angle <- function(theta) {
  atan2(sin(theta), cos(theta))
}

u_to_theta <- function(u) {
  u <- as.numeric(u)
  nu <- sqrt(sum(u^2))
  u <- u / nu
  atan2(u[2], u[1])
}

theta_to_u <- function(theta) {
  c(cos(theta), sin(theta))
}

propose_r_theta_joint <- function(r, theta, step_r, step_theta) {
  r_prop <- r + rnorm(1, mean = 0, sd = step_r)
  if (r_prop < 0) {
    r_prop <- -r_prop
  }
  
  theta_prop <- wrap_angle(theta + rnorm(1, mean = 0, sd = step_theta))
  
  list(
    r = r_prop,
    theta = theta_prop,
    u = theta_to_u(theta_prop)
  )
}


default_mu_prior_centre <- function(d = 2) {
  c(1, rep(0, d))
}

as_sigma_vector <- function(sigma2_clust, K) {
  if (length(sigma2_clust) == 1) {
    return(rep(as.numeric(sigma2_clust), K))
  }
  
  as.numeric(sigma2_clust)
}

as_dirichlet_vector <- function(v_lambda, K) {
  if (length(v_lambda) == 1) {
    return(rep(as.numeric(v_lambda), K))
  }  
  as.numeric(v_lambda)
}

rdirichlet_simple <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

rinvgamma <- function(n, shape, rate) {
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
  
  ssq_g <- cluster_tangent_ssq(Z = Z, mu_g = mu_g, idx_g = idx_g)
  m_g <- length(idx_g)
  
  post_shape <- ig_shape + 0.5 * m_g * p
  post_rate  <- ig_rate  + 0.5 * ssq_g
  
  sigma2_draw <- rinvgamma(1, shape = post_shape, rate = post_rate)
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



apply_spatial_rotation_h2 <- function(Z, R, renormalise = TRUE) {
  Z <- as.matrix(Z)
  
  X <- Z[, 2:3, drop = FALSE] %*% R
  
  if (renormalise) {
    t_coord <- sqrt(1 + rowSums(X^2))
  } else {
    t_coord <- Z[, 1]
  }
  
  cbind(t_coord, X)
}


log_prior_mu_radial <- function(mu, r_max = 5) {
  mu <- as.numeric(mu)
  
  polar_mu <- hyperboloid_to_polar(mu)
  r_mu <- as.numeric(polar_mu$r)
  
  if (!is.finite(r_mu)) {
    return(-Inf)
  }
  
  if (r_mu < 0 || r_mu > r_max) {
    return(-Inf)
  }
  
  0
}

initial_cluster_means_clust <- function(K, r = 1.5, R = 1) {
  mu_clust_init <- matrix(NA_real_, nrow = K, ncol = 3)
  
  for (k in 1:K) {
    theta <- 2 * pi * (k - 1) / K
    
    mu_clust_init[k, ] <- c(
      R * cosh(r / R),
      R * sinh(r / R) * cos(theta),
      R * sinh(r / R) * sin(theta)
    )
  }
  
  mu_clust_init
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


relabel_by_mu_angle <- function(mu_clust, clustering_labels, lambda, sigma2_clust) {
  mu_clust <- as.matrix(mu_clust)
  
  angles <- atan2(mu_clust[, 3], mu_clust[, 2])
  ord <- order(angles)
  
  inv_ord <- integer(length(ord))
  inv_ord[ord] <- seq_along(ord)
  
  list(
    mu_clust = mu_clust[ord, , drop = FALSE],
    clustering_labels = inv_ord[clustering_labels],
    lambda = lambda[ord],
    sigma2_clust = sigma2_clust[ord],
    ord = ord
  )
}

match_clusters_to_reference <- function(mu_curr, mu_ref) {
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Please install the 'clue' package")
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
  
  as.integer(clue::solve_LSAP(C))
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



Z_update_step_clust <- function(
    Z, Y, alpha,
    mu_clust, sigma2_clust, clustering_labels,
    sigma_z_prop,
    i
){
  n <- nrow(Z)
  z_i <- as.numeric(Z[i, ])
  k_i <- clustering_labels[i]
  mu_k <- as.numeric(mu_clust[k_i, ])
  sigma_k <- sigma2_clust[k_i]
  
  prop_out <- tryCatch(
    proposal_tangent_rw_h2(z_old = z_i, sigma_prop = sigma_z_prop),
    error = function(e) NULL
  )
  
  if (is.null(prop_out)) {
    return(list(
      Z = Z,
      accepted = 0,
      i = i,
      log_r = -Inf
    ))
  }
  
  z_new <- prop_out$z_new
  
  if (!all(is.finite(z_new))) {
    return(list(
      Z = Z,
      accepted = 0,
      i = i,
      log_r = -Inf
    ))
  }
  
  log_prior_old <- log_pdf_PHN(z_i,   mu_k, sigma_k)
  log_prior_new <- log_pdf_PHN(z_new, mu_k, sigma_k)
  
  log_like_old <- loglik_i(
    i = i,
    z_i = z_i,
    Z = Z,
    Y = Y,
    alpha = alpha
  )
  
  log_like_new <- loglik_i(
    i = i,
    z_i = z_new,
    Z = Z,
    Y = Y,
    alpha = alpha
  )
  
  log_q_forward <- tryCatch(
    log_q_tangent_rw_h2(z_to = z_new, z_from = z_i, sigma_prop = sigma_z_prop),
    error = function(e) -Inf
  )
  
  log_q_reverse <- tryCatch(
    log_q_tangent_rw_h2(z_to = z_i, z_from = z_new, sigma_prop = sigma_z_prop),
    error = function(e) -Inf
  )
  
  log_r <- (log_like_new + log_prior_new + log_q_reverse) -
    (log_like_old + log_prior_old + log_q_forward)
  
  accepted <- 0
  
  if (is.finite(log_r) && log(runif(1)) < log_r) {
    Z[i, ] <- z_new
    accepted <- 1
  }
  
  list(
    Z = Z,
    accepted = accepted,
    i = i,
    log_r = log_r
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
  
  if (is.null(mu_new) || !all(is.finite(mu_new))) {
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

poincare_from_hyperboloid <- function(Z) {
  Z <- as.matrix(Z)
  cbind(
    Z[, 2] / (1 + Z[, 1]),
    Z[, 3] / (1 + Z[, 1])
  )
}

plot_poincare_embedding <- function(
    Z,
    Y = NULL,
    clustering_labels = NULL,
    point_cex = 1.2,
    point_pch = 19,
    edge_col = "grey80",
    node_col = "black",
    main = "Poincare disk embedding"
) {
  Z <- as.matrix(Z)
  U <- poincare_from_hyperboloid(Z)
  n <- nrow(U)
  
  if (!is.null(clustering_labels)) {
    cols <- rainbow(max(clustering_labels))[as.integer(clustering_labels)]
  } else {
    cols <- rep(node_col, n)
  }
  
  plot(
    U[, 1], U[, 2],
    type = "n",
    asp = 1,
    xlim = c(-1.05, 1.05),
    ylim = c(-1.05, 1.05),
    xlab = "",
    ylab = "",
    main = main
  )
  
  theta <- seq(0, 2 * pi, length.out = 400)
  lines(cos(theta), sin(theta), lwd = 1.2)
  
  if (!is.null(Y)) {
    Y <- as.matrix(Y)
    idx <- which(upper.tri(Y) & Y != 0, arr.ind = TRUE)
    
    for (k in 1:nrow(idx)) {
      i <- idx[k, 1]
      j <- idx[k, 2]
      segments(U[i, 1], U[i, 2], U[j, 1], U[j, 2], col = edge_col)
    }
  }
  
  points(U[, 1], U[, 2], pch = point_pch, cex = point_cex, col = cols)
  
  invisible(U)
}
