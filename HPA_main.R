prepare_hyperboloid_matrix <- function(X, tol_hyp = 1e-5) {
  
  X <- as.matrix(X)
  
  for (i in seq_len(nrow(X))) {
    norm_i <- lorentz_inner(X[i, ], X[i, ])
    if (
      !is.finite(norm_i) ||
      abs(norm_i + 1) > tol_hyp ||
      X[i, 1L] <= 0
    ) {
      X[i, ] <- project_to_hyperboloid(X[i, ])
    }
  }
  
  X
}

HPA_centroid_parameter <- function(X, weights = NULL) {
  X <- as.matrix(X)
  n <- nrow(X)
  
  if (is.null(weights)) {
    
    weights <- rep(1 / n, n)
    
  } else {
    
    weights <- as.numeric(weights)
    
    if (length(weights) != n) {
      stop("weights must have length nrow(X).")
    }
    
    if (
      any(!is.finite(weights)) ||
      any(weights < 0) ||
      sum(weights) <= 0
    ) {
      stop("weights must be finite, non-negative and have positive sum.")
    }
    
    weights <- weights / sum(weights)
  }
  
  x_bar <- colSums(
    sweep(X, 1L, weights, "*")
  )
  
  x_bar_norm <- lorentz_inner(x_bar, x_bar)
  
  if (!is.finite(x_bar_norm) || x_bar_norm >= 0) {
    stop(
      paste0(
        "The ambient mean is not timelike: ",
        "<x_bar, x_bar>_L = ",
        signif(x_bar_norm, 6),
        "."
      )
    )
  }
  
  scale <- sqrt(-x_bar_norm)
  
  m <- x_bar[-1L] / scale
  
  list(
    m = as.numeric(m),
    x_bar = as.numeric(x_bar),
    scale = scale,
    weights = weights
  )
}

HPA_translation_matrix <- function(b) {
  
  b <- as.numeric(b)
  
  if (length(b) < 1L || any(!is.finite(b))) {
    stop("b must be a finite vector.")
  }
  
  d <- length(b)
  
  gamma <- sqrt(1 + sum(b^2))
  
  spatial_block <-
    diag(d) +
    tcrossprod(b) / (1 + gamma)
  
  R_b <- matrix(
    0,
    nrow = d + 1L,
    ncol = d + 1L
  )
  
  R_b[1L, 1L] <- gamma
  R_b[1L, -1L] <- b
  R_b[-1L, 1L] <- b
  R_b[-1L, -1L] <- spatial_block
  
  R_b
}

HPA_rotation_matrix <- function(U) {
  
  U <- as.matrix(U)
  
  if (nrow(U) != ncol(U)) {
    stop("U must be square.")
  }
  
  d <- nrow(U)
  
  orth_error <- max(
    abs(crossprod(U) - diag(d))
  )
  
  if (!is.finite(orth_error) || orth_error > 1e-6) {
    stop("U is not numerically orthogonal.")
  }
  
  R_U <- matrix(
    0,
    nrow = d + 1L,
    ncol = d + 1L
  )
  
  R_U[1L, 1L] <- 1
  R_U[-1L, -1L] <- U
  
  R_U
}

HPA_apply_matrix <- function(X, A, repair = FALSE) {
  
  X <- as.matrix(X)
  A <- as.matrix(A)
  
  if (
    nrow(A) != ncol(A) ||
    ncol(X) != nrow(A)
  ) {
    stop("A must be square and compatible with X.")
  }
  
  X_out <- X %*% t(A)
  
  if (repair) {
    for (i in seq_len(nrow(X_out))) {
      X_out[i, ] <- project_to_hyperboloid(X_out[i, ])
    }
  }
  
  X_out
}

HPA_align_draw <- function(
    Z_ref,
    Z_new,
    allow_reflection = TRUE,
    weights = NULL,
    tol_hyp = 1e-5,
    repair_output = FALSE
) {
  
  Z_ref <- prepare_hyperboloid_matrix(
    Z_ref,
    tol_hyp = tol_hyp
  )
  
  Z_new <- prepare_hyperboloid_matrix(
    Z_new,
    tol_hyp = tol_hyp
  )
  
  if (!all(dim(Z_ref) == dim(Z_new))) {
    stop("Z_ref and Z_new must have matching dimensions.")
  }
  
  n <- nrow(Z_ref)
  d <- ncol(Z_ref) - 1L
  
  if (n < d) {
    warning(
      "There are fewer points than spatial dimensions; ",
      "the Procrustes rotation may not be identifiable."
    )
  }
  
  ref_centroid <- HPA_centroid_parameter(
    X = Z_ref,
    weights = weights
  )
  
  new_centroid <- HPA_centroid_parameter(
    X = Z_new,
    weights = weights
  )
  
  m_ref <- ref_centroid$m
  m_new <- new_centroid$m
  
  R_minus_m_ref <- HPA_translation_matrix(-m_ref)
  R_minus_m_new <- HPA_translation_matrix(-m_new)
  
  Z_ref_centred <- HPA_apply_matrix(
    X = Z_ref,
    A = R_minus_m_ref
  )
  
  Z_new_centred <- HPA_apply_matrix(
    X = Z_new,
    A = R_minus_m_new
  )
  
  P_ref <- Z_ref_centred[, -1L, drop = FALSE]
  P_new <- Z_new_centred[, -1L, drop = FALSE]
  
  if (is.null(weights)) {
    
    C <- crossprod(P_ref, P_new)
    
  } else {
    
    weights_use <- as.numeric(weights)
    weights_use <- weights_use / sum(weights_use)
    
    C <- crossprod(
      P_ref,
      sweep(P_new, 1L, weights_use, "*")
    )
  }
  
  sv <- svd(C)
  
  U_hat <- sv$u %*% t(sv$v)
  
  if (!allow_reflection && det(U_hat) < 0) {
    
    D_fix <- diag(d)
    D_fix[d, d] <- -1
    
    U_hat <- sv$u %*% D_fix %*% t(sv$v)
  }
  
  R_U <- HPA_rotation_matrix(U_hat)
  R_m_ref <- HPA_translation_matrix(m_ref)
  
  A_full <-
    R_m_ref %*%
    R_U %*%
    R_minus_m_new
  
  Z_aligned <- HPA_apply_matrix(
    X = Z_new,
    A = A_full,
    repair = repair_output
  )
  
  list(
    Z_aligned = Z_aligned,
    A = A_full,
    U = U_hat,
    m_ref = m_ref,
    m_new = m_new,
    Z_ref_centred = Z_ref_centred,
    Z_new_centred = Z_new_centred,
    P_ref = P_ref,
    P_new = P_new,
    singular_values = sv$d,
    reflection_used = det(U_hat) < 0
  )
}


align_one_clustering_draw <- function(
    Z,
    mu_clust,
    Z_ref,
    allow_reflection = TRUE,
    tol_hyp = 1e-5
) {
  
  Z <- as.matrix(Z)
  mu_clust <- as.matrix(mu_clust)
  Z_ref <- as.matrix(Z_ref)
  
  if (ncol(mu_clust) != ncol(Z)) {
    stop("mu_clust and Z must have the same number of columns.")
  }
  
  hpa_out <- HPA_align_draw(
    Z_ref = Z_ref,
    Z_new = Z,
    allow_reflection = allow_reflection,
    tol_hyp = tol_hyp
  )
  
  mu_store <- HPA_apply_matrix(
    X = mu_clust,
    A = hpa_out$A
  )
  
  list(
    Z_store = hpa_out$Z_aligned,
    mu_store = mu_store,
    proc_R = hpa_out$U,
    proc_A = hpa_out$A,
    m_ref = hpa_out$m_ref,
    m_new = hpa_out$m_new,
    
    singular_values = hpa_out$singular_values,
    reflection_used = hpa_out$reflection_used
  )
}


align_fixed_origin <- function(
    Z,
    Z_ref,
    allow_reflection = TRUE,
    tol_hyp = 1e-5
) {
  
  Z <- prepare_hyperboloid_matrix(
    Z,
    tol_hyp = tol_hyp
  )
  
  Z_ref <- prepare_hyperboloid_matrix(
    Z_ref,
    tol_hyp = tol_hyp
  )
  
  if (!all(dim(Z) == dim(Z_ref))) {
    stop("Z and Z_ref must have matching dimensions.")
  }
  
  d <- ncol(Z) - 1L
  
  P_new <- Z[, -1L, drop = FALSE]
  P_ref <- Z_ref[, -1L, drop = FALSE]
  C <- crossprod(P_ref, P_new)
  
  sv <- svd(C)
  
  U_hat <- sv$u %*% t(sv$v)
  
  if (!allow_reflection && det(U_hat) < 0) {
    
    D_fix <- diag(d)
    D_fix[d, d] <- -1
    
    U_hat <- sv$u %*% D_fix %*% t(sv$v)
  }
  
  A_full <- HPA_rotation_matrix(U_hat)
  
  Z_store <- HPA_apply_matrix(
    X = Z,
    A = A_full
  )
  
  list(
    Z_store = Z_store,
    mu_store = matrix(
      c(1, rep(0, d)),
      nrow = 1L
    ),
    rotation = U_hat,
    proc_R = U_hat,
    proc_A = A_full,
    reflection_used = det(U_hat) < 0
  )
}

