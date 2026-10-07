neighbours <- function(Y) {
  dimen <- nrow(Y)
  
  N <- integer(dimen)
  non_neighbours_ind <- vector("list", dimen)
  neighbours_ind <- vector("list", dimen)
  
  for (i in seq_len(dimen)) {
    possible <- seq_len(dimen) < i
    
    non_neighbours_ind[[i]] <- which(Y[i, ] == 0 & possible)
    neighbours_ind[[i]] <- which(Y[i, ] == 1 & possible)
    
    N[i] <- length(non_neighbours_ind[[i]])
  }
  
  
  list(
    N = N,
    non_neighbours = non_neighbours_ind,
    neighbours = neighbours_ind
  )
}

strata <- function(Y){
 n <- nrow(Y)
 graph <- igraph::graph_from_adjacency_matrix(Y)
 graph <- igraph::as_undirected(graph = graph)
 distance_mat <- igraph::distances(graph = graph, v = igraph::V(graph), to = igraph::V(graph))
 
 max_dist <- max(distance_mat)

 N_ih <- matrix(0,nrow = n, ncol = max_dist)
 index <- vector("list",n)
 
 for (i in 1:n){
     for(k in 1:max_dist){
       N_ih[i,k] <- length(which(distance_mat[i,] == k))
       index[[i]][[k]] <- which(distance_mat[i,] == k)
     }
 }
 
 return(list(
   N_ih = N_ih,
   index = index,
   distance_mat = distance_mat,
   strata_dimension = max_dist,
   number_nodes = n
 ))
}

stratified_sampler <- function(strata, n_ih){
  n <- strata$number_nodes
  no_bins <- strata$strata_dimension
  N_ih <- strata$N_ih
  
  sampled_non_neighbours_ind <- vector("list",n)
  
  for (i in 1:n){
    for (k in 2:no_bins){
    sample_size <- min(N_ih[i,k],n_ih[i,k])
    sampled_non_neighbours_ind[[i]][[k]] <- sample(strata$index[[i]][[k]], sample_size, replace = FALSE)
    }
  }
  return(sampled_non_neighbours_ind)
}

n_i0 <- function(Y, strata,r = 1){
  H <- strata$strata_dimension
  n <- strata$number_nodes

  d <- sum(Y)/nrow(Y)
  
  n_i0 <- rep(ceiling(r * d), n)
  return(n_i0)
}

softplus <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

eta_indices <- function(
    i,
    inds,
    X = NULL,
    alpha,
    B0,
    B1,
    D = NULL,
    d_i = NULL
) {
  
  if (length(inds) == 0L) {
    return(numeric(0L))
  }
  inds <- as.numeric(unlist(inds))
  bins <- length(inds)
  
  d <- if (is.null(d_i)) D[i, inds] else d_i[inds]
  
  if (is.null(X)) {
    return(alpha - B1 * d)
  }
  
  X_i <- X[i, inds, , drop = FALSE]
  
  X_i <- matrix(
    X_i,
    nrow = length(inds),
    ncol = length(B0)
  )
  
  if (anyNA(X_i)) {
    stop("Missing covariate information for node ", i)
  }
  
  alpha + drop(X_i %*% B0) - B1 * d
}

pilot_sampler <- function(
    pilot_control_pool,
    n_i0
) {
  n <- length(pilot_control_pool)
  if (length(n_i0) == 1L) {
    n_i0 <- rep(n_i0, n)
  }
  sampled <- vector("list", n)
  for (i in seq_len(n)) {
    controls_i <- pilot_control_pool[[i]]
    
    sample_size <- min(
      length(controls_i),
      n_i0[i]
    )
    
    sampled[[i]] <- sample(
      controls_i,
      size = sample_size,
      replace = FALSE
    )
  }
  
  sampled
}


pilot_loglik_i <- function(
    i,
    X = NULL,
    alpha,
    B0,
    B1,
    stratified_network,
    sampled_controls,
    D = NULL,
    d_i = NULL
) {
  H <- stratified_network$strata_dimension
  N_ih <- stratified_network$N_ih
  edge_inds <- stratified_network$index[[i]][[1]]
  eta_edge <- eta_indices(
    i = i,
    inds = edge_inds,
    X = X,
    alpha = alpha,
    B0 = B0,
    B1 = B1,
    D = D,
    d_i = d_i
  )
  
  lik_edges <- sum(
    eta_edge - softplus(eta_edge)
  )
  
  control_inds <- sampled_controls[[i]]
  N_i0 <- sum(N_ih[i, 2:H])
  n_i0 <- length(control_inds)
  strata_contribution <- numeric(H)
  
  if (N_i0 > 0L) {
    if (n_i0 == 0L) {
      stop("Node ", i, " has controls but none were sampled.")
    }
    eta_control <- eta_indices(
      i = i,
      inds = control_inds,
      X = X,
      alpha = alpha,
      B0 = B0,
      B1 = B1,
      D = D,
      d_i = d_i
    )
    
    control_terms <- -softplus(eta_control)
    sampled_h <- as.integer(
    stratified_network$distance_mat[i, control_inds]
    )
    
    scale <- N_i0 / n_i0
      tmp <- tapply(
      control_terms,
      sampled_h,
      sum
    )
    
    if (length(tmp) > 0L) {
      strata_contribution[
        as.integer(names(tmp))
      ] <- scale * as.numeric(tmp)
    }
  }
  
  list(
    loglik_i = lik_edges + sum(strata_contribution),
    strata_contribution = strata_contribution,
    lik_edges = lik_edges
  )
}

production_loglik_i <- function(
    i,
    X = NULL,
    alpha,
    B0,
    B1,
    stratified_network,
    sampled_controls,
    D = NULL,
    d_i = NULL
) {
  H <- stratified_network$strata_dimension
  N_ih <- stratified_network$N_ih
  edge_inds <- stratified_network$index[[i]][[1]]
  eta_edge <- eta_indices(
    i = i,
    inds = edge_inds,
    X = X,
    alpha = alpha,
    B0 = B0,
    B1 = B1,
    D = D,
    d_i = d_i
  )
  
  lik_edges <- sum(
    eta_edge - softplus(eta_edge)
  )
  lik_controls <- numeric(H)
  
  for (h in 2:H) {
    
    N_h <- N_ih[i, h]
    
    if (N_h == 0L) {
      next
    }
    
    inds <- sampled_controls[[i]][[h]]
    n_h <- length(inds)
    
    if (n_h == 0L) {
      stop(
        "Stratum ", h,
        " for node ", i,
        " contains controls but none were sampled."
      )
    }
    
    eta_h <- eta_indices(
      i = i,
      inds = inds,
      X = X,
      alpha = alpha,
      B0 = B0,
      B1 = B1,
      D = D,
      d_i = d_i
    )
    
    lik_controls[h] <-
      (N_h / n_h) *
      sum(-softplus(eta_h))
  }
  
  list(
    loglik_i =
      lik_edges + sum(lik_controls),
    
    strata_contribution =
      lik_controls,
    
    lik_edges =
      lik_edges
  )
}


full_case_control_loglik <- function(
    Z,
    X = NULL,
    alpha,
    B0,
    B1,
    stratified_network,
    sampled_controls,
    D,
    mode = c("pilot", "production")
) {
  mode <- match.arg(mode)
  node_loglik_fun <- switch(
    mode,
    pilot      = pilot_loglik_i,
    production = production_loglik_i
  )
  n <- nrow(Z)
  full_loglik <- 0
  
  for (i in seq_len(n)) {
    out <- node_loglik_fun(
      i = i,
      X = X,
      alpha = alpha,
      B0 = B0,
      B1 = B1,
      stratified_network = stratified_network,
      sampled_controls = sampled_controls,
      D = D
    )
    full_loglik <- full_loglik + out$loglik_i
  }
  return(0.5*full_loglik)
}


allocate_controls <- function(weights,stratified_network,n_i0){
  N_ih <- stratified_network$N_ih
  stopifnot(
    is.matrix(weights),
    is.matrix(N_ih),
    identical(dim(weights), dim(N_ih)),
    all(is.finite(N_ih)),
    all(N_ih >= 0),
    all(N_ih == floor(N_ih))
  )
  n <- nrow(N_ih)
  H <- ncol(N_ih)
  
  if (length(n_i0) == 1L) {
    n_i0 <- rep(n_i0, n)
  }
  
  stopifnot(
    length(n_i0) == n,
    all(is.finite(n_i0)),
    all(n_i0 >= 0),
    all(n_i0 == floor(n_i0))
  )
  
  n_ih <- matrix(
    0L, nrow = n, ncol = H,
    dimnames = dimnames(N_ih)
  )
  
  for (i in seq_len(n)) {
    bins <- which(seq_len(H) > 1L & N_ih[i, ] > 0)
    if (!length(bins)) next
    
    capacity <- N_ih[i, bins]
    w <- weights[i, bins]
    
    if (any(!is.finite(w)) || any(w < 0)) {
      stop("Invalid weights for node ", i)
    }
    
    target <- min(n_i0[i], sum(capacity))
    
    if (target < length(bins)) {
      stop(
        "n_i0 is too small for node ", i,
        ": need at least ", length(bins),
        " controls to cover every nonempty stratum."
      )
    }
    counts <- rep(1L, length(bins))
    remaining <- target - sum(counts)
    
    while (remaining > 0) {
      eligible <- which(counts < capacity)
      p <- w[eligible]
      if (all(p == 0)) {
        p[] <- 1
      }
      
      p <- p / max(p)
      p <- p / sum(p)
      
      quota <- remaining * p
      extra <- pmin(
        floor(quota),
        capacity[eligible] - counts[eligible]
      )
      if (sum(extra) == 0) {
        extra[which.max(quota)] <- 1L
      }
      
      counts[eligible] <- counts[eligible] + extra
      remaining <- remaining - sum(extra)
    }
    
    n_ih[i, bins] <- as.integer(counts)
  }
  
  n_ih
}


build_case_control_set <- function(
    X = NULL,
    stratified_network,
    sampled_controls,
    D,
    mode = c("pilot", "production")
) {
  mode <- match.arg(mode)
  n <- stratified_network$number_nodes
  H <- stratified_network$strata_dimension
  N_ih <- stratified_network$N_ih
  
  i_list <- list()
  j_list <- list()
  y_list <- list()
  weight_list <- list()
  
  counter <- 0L
  for (i in seq_len(n)) {
    edge_inds <- stratified_network$index[[i]][[1L]]
    if (length(edge_inds) > 0L) {
      counter <- counter + 1L
      
      i_list[[counter]] <- rep.int(i, length(edge_inds))
      j_list[[counter]] <- edge_inds
      
      y_list[[counter]] <- rep.int(1L, length(edge_inds))
      weight_list[[counter]] <- rep(1, length(edge_inds))
    }
    
    if (mode == "pilot") {
      inds <- sampled_controls[[i]]
      if (length(inds) > 0L) {
        N_i0 <- if (H >= 2L) {
          sum(N_ih[i, 2:H])
        } else {
          0
        }
      
        n_i0 <- length(inds)
        counter <- counter + 1L
        i_list[[counter]] <- rep.int(i, n_i0)
        j_list[[counter]] <- inds
        y_list[[counter]] <- rep.int(0L, n_i0)
        weight_list[[counter]] <-
          rep(N_i0 / n_i0, n_i0)
      }
    } else {
      if (H >= 2L) {
        for (h in 2:H) {
          N_h <- N_ih[i, h]
          if (N_h == 0L) {
            next
          }
          inds <- sampled_controls[[i]][[h]]
          n_h <- length(inds)
          if (n_h == 0L) {
            stop(
              "Stratum ", h,
              " for node ", i,
              " contains controls but none were sampled."
            )
          }
          
          counter <- counter + 1L
          i_list[[counter]] <- rep.int(i, n_h)
          j_list[[counter]] <- inds
          y_list[[counter]] <- rep.int(0L, n_h)
          weight_list[[counter]] <-
            rep(N_h / n_h, n_h)
        }
      }
    }
  }
  
  i_idx <- unlist(i_list, use.names = FALSE)
  j_idx <- unlist(j_list, use.names = FALSE)
  y <- unlist(y_list, use.names = FALSE)
  weights <- unlist(weight_list, use.names = FALSE)
  d <- D[cbind(i_idx, j_idx)]
  if (is.null(X)) {
    X_work <- NULL
  } else {
    q <- dim(X)[3L]
    X_work <- matrix(
      NA_real_,
      nrow = length(i_idx),
      ncol = q
    )
    
    for (k in seq_len(q)) {
      X_work[, k] <-
        X[cbind(
          i_idx,
          j_idx,
          rep.int(k, length(i_idx))
        )]
    }
    
    if (anyNA(X_work)) {
      stop("Missing values found in case-control covariates.")
    }
    cov_names <- dimnames(X)[[3L]]
    
    if (!is.null(cov_names)) {
      colnames(X_work) <- cov_names
    }
  }
  list(
    y = y,
    weights = weights,
    d = d,
    X = X_work,
    i = i_idx,
    j = j_idx
  )
}

fast_case_control_loglik <- function(
    working_set,
    alpha,
    B0,
    B1
) {
  eta <-
    alpha -
    B1 * working_set$d
  if (!is.null(working_set$X)) {
    eta <-
      eta +
      drop(working_set$X %*% B0)
  }
  
  loglik <-
    sum(
      working_set$weights *
        (
          working_set$y * eta -
            softplus(eta)
        )
    )
  
  0.5 * loglik
}

