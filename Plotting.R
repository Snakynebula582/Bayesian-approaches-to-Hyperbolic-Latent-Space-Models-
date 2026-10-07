plot_poincare_embedding <- function(
    Z_summary,                   
    Y = NULL,                     
    node_labels = NULL,          
    clustering_labels = NULL,     
    mu_summary = NULL,            
    posterior_nodes = NULL, 
    posterior_clusters = NULL,    
    
    show_edges = TRUE,
    show_labels = TRUE,
    show_cluster_means = TRUE,
    show_posterior_clouds = FALSE,
    posterior_for = c("nodes", "clusters", "both"),
    
    point_cex = 1.2,
    point_pch = 19,
    cluster_mean_cex = 2,
    cluster_mean_pch = 8,
    posterior_cex = 0.35,
    posterior_pch = 16,
    posterior_alpha = 0.08,
    
    edge_col = "grey80",
    node_col = "lightskyblue",
    cluster_mean_col = "red",
    main = "Poincare disk embedding"
) {
  
  posterior_for <- match.arg(posterior_for)
  
  Z_summary <- as.matrix(Z_summary)
  n <- nrow(Z_summary)
  
  if (ncol(Z_summary) != 3) {
    stop("Z_summary must be an n x 3 matrix.")
  }
  
  U <- poincare_from_hyperboloid(Z_summary)

  # colours
  if (!is.null(clustering_labels)) {
    clustering_labels <- as.integer(clustering_labels)
    
    if (length(clustering_labels) != n) {
      stop("clustering_labels must have length nrow(Z_summary).")
    }
    
    K <- max(clustering_labels)
    cluster_cols <- rainbow(K)
    node_cols <- cluster_cols[clustering_labels]
  } else {
    node_cols <- rep(node_col, n)
    cluster_cols <- NULL
  }
  
  # blank plot
  plot(
    U[,1], U[,2],
    type = "n",
    asp = 1,
    xlim = c(-1.05, 1.05),
    ylim = c(-1.05, 1.05),
    xlab = "",
    ylab = "",
    main = main
  )
  
  theta <- seq(0, 2*pi, length.out = 400)
  lines(cos(theta), sin(theta), lwd = 1.2)
  
  # edges
  if (!is.null(Y) && show_edges) {
    Y <- as.matrix(Y)
    
    if (!all(dim(Y) == c(n, n))) {
      stop("Y must be an n x n matrix matching Z_summary.")
    }
    
    idx <- which(upper.tri(Y) & Y != 0, arr.ind = TRUE)
    
    if (nrow(idx) > 0) {
      for (k in 1:nrow(idx)) {
        i <- idx[k, 1]
        j <- idx[k, 2]
        segments(U[i,1], U[i,2], U[j,1], U[j,2], col = edge_col)
      }
    }
  }
  
  # node posterior clouds
  if (show_posterior_clouds && posterior_for %in% c("nodes", "both") && !is.null(posterior_nodes)) {
    
    # case 1: posterior_nodes is array(M, n, 3)
    if (is.array(posterior_nodes) && length(dim(posterior_nodes)) == 3) {
      M <- dim(posterior_nodes)[1]
      
      if (dim(posterior_nodes)[2] != n || dim(posterior_nodes)[3] != 3) {
        stop("posterior_nodes array must have dimension (M, n, 3).")
      }
      
      for (i in 1:n) {
        Z_post_i <- matrix(posterior_nodes[, i, ], nrow = M, ncol = 3)
        U_post_i <- poincare_from_hyperboloid(Z_post_i)
        
        points(
          U_post_i[,1], U_post_i[,2],
          pch = posterior_pch,
          cex = posterior_cex,
          col = grDevices::adjustcolor(node_cols[i], alpha.f = posterior_alpha)
        )
      }
    }
    
    # case 2: posterior_nodes is list of length n, each entry M_i x 3
    else if (is.list(posterior_nodes)) {
      if (length(posterior_nodes) != n) {
        stop("posterior_nodes list must have length n.")
      }
      
      for (i in 1:n) {
        Z_post_i <- as.matrix(posterior_nodes[[i]])
        
        if (ncol(Z_post_i) != 3) {
          stop("Each posterior_nodes[[i]] must be an M_i x 3 matrix.")
        }
        
        U_post_i <- poincare_from_hyperboloid(Z_post_i)
        
        points(
          U_post_i[,1], U_post_i[,2],
          pch = posterior_pch,
          cex = posterior_cex,
          col = grDevices::adjustcolor(node_cols[i], alpha.f = posterior_alpha)
        )
      }
    }
    
    else {
      stop("posterior_nodes must be either an array(M, n, 3) or a list of n matrices.")
    }
  }
  
  # cluster mean posterior clouds
  if (show_posterior_clouds && posterior_for %in% c("clusters", "both") && !is.null(posterior_clusters)) {
    
    if (is.null(mu_summary)) {
      stop("If posterior_clusters is supplied, mu_summary must also be supplied.")
    }
    
    K_mu <- nrow(mu_summary)
    
    # case 1: array(M, K, 3)
    if (is.array(posterior_clusters) && length(dim(posterior_clusters)) == 3) {
      M <- dim(posterior_clusters)[1]
      
      if (dim(posterior_clusters)[2] != K_mu || dim(posterior_clusters)[3] != 3) {
        stop("posterior_clusters array must have dimension (M, K, 3).")
      }
      
      for (k in 1:K_mu) {
        Z_post_k <- matrix(posterior_clusters[, k, ], nrow = M, ncol = 3)
        U_post_k <- poincare_from_hyperboloid(Z_post_k)
        
        col_k <- if (!is.null(cluster_cols)) cluster_cols[k] else cluster_mean_col
        
        points(
          U_post_k[,1], U_post_k[,2],
          pch = posterior_pch,
          cex = posterior_cex,
          col = grDevices::adjustcolor(col_k, alpha.f = posterior_alpha)
        )
      }
    }
    
    # case 2: list of K matrices
    else if (is.list(posterior_clusters)) {
      if (length(posterior_clusters) != K_mu) {
        stop("posterior_clusters list must have length nrow(mu_summary).")
      }
      
      for (k in 1:K_mu) {
        Z_post_k <- as.matrix(posterior_clusters[[k]])
        
        if (ncol(Z_post_k) != 3) {
          stop("Each posterior_clusters[[k]] must be an M_k x 3 matrix.")
        }
        
        U_post_k <- poincare_from_hyperboloid(Z_post_k)
        
        col_k <- if (!is.null(cluster_cols)) cluster_cols[k] else cluster_mean_col
        
        points(
          U_post_k[,1], U_post_k[,2],
          pch = posterior_pch,
          cex = posterior_cex,
          col = grDevices::adjustcolor(col_k, alpha.f = posterior_alpha)
        )
      }
    }
    
    else {
      stop("posterior_clusters must be either an array(M, K, 3) or a list of K matrices.")
    }
  }
  
  # nodes
  points(U[,1], U[,2], pch = point_pch, cex = point_cex, col = node_cols)
  
  # labels
  if (show_labels) {
    if (is.null(node_labels)) {
      node_labels <- 1:n
    }
    
    text(U[,1], U[,2], labels = node_labels, pos = 3, cex = 0.75, col = node_cols)
  }
  
  # cluster means
  if (show_cluster_means && !is.null(mu_summary)) {
    U_mu <- poincare_from_hyperboloid(mu_summary)
    
    if (!is.null(cluster_cols) && nrow(mu_summary) <= length(cluster_cols)) {
      mu_cols <- cluster_cols[1:nrow(mu_summary)]
    } else {
      mu_cols <- rep(cluster_mean_col, nrow(mu_summary))
    }
    
    points(
      U_mu[,1], U_mu[,2],
      pch = cluster_mean_pch,
      cex = cluster_mean_cex,
      lwd = 2,
      col = mu_cols
    )
  }
  
  invisible(list(U = U))
}


