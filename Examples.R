# This file is meant to serve an example on how to call these models, the main input is 
# the networks adjacency matrix, the functions will compute the MDS initialisations by default
# Here is a toy tree to run this model on, note that these models require a lot of iterations to 
# mix properly so the resulting outputs are not very representative of the model fits and are meant to 
# just demonstrate that these models run


# Generate a random graph
set.seed(67)
n_nodes <- 20
tree_graph <- sample_pa(
  n = n_nodes,
  power = 1,
  m = 1,
  directed = FALSE
)

toy_adj <- as.matrix(as_adjacency_matrix(tree_graph, sparse = FALSE))
plot(tree_graph)



# Radial Angular model
fit_radial <- Radial_Angular_Model(
  Y = toy_adj, # netowrk adjacency matrix
  d = 2, # dimension, (currently only works for 2)
  n_iter = 2000,  # number of iterations, burn and thin
  burn = 1000,
  thin = 1,
  mu_alpha = 0, # alpha is normal a priori, this is its mean and variance
  sigma_alpha = 0.3, 
  step_r = 0.15, # proposal parameters for r, theta and alpha
  step_u = 0.20,
  step_alpha = 0.20, 
  postproc = "hpa", # type of post processing allignment, choose "naive" (polar) or "hpa"
  allow_reflection = TRUE, # Allows the post processing alignment to deal with reflections
  alpha_init = 0, 
  sigma_r_prior = 3, # r is halfnormal a priori, this is its variance parameter
  learn_curvature = FALSE, # learn curvature switch
  K = 1, # initial/ fixed curvature 
  compute_posterior_mean = TRUE, # Compute naive polar means (fast)
  compute_barycenter = TRUE, #Compute posterior barycentre (slow)
)

plot_poincare_embedding(fit_radial$posterior_barycenter, Y = toy_adj)




# The Centred Wrapped Normal example call
fit_centred <- Centred_Wrapped_Normal(
  Y = toy_adj, 
  Z_init = NULL,
  alpha_init = 0,
  n_iter = 2000,
  burn = 1000,
  mu_alpha = 0, 
  sigma2_alpha = 1,
  sigma2_prop_alpha = 0.5,
  sigma2_z = 1, 
  sigma_prop_z = 0.15, #isotrpoic variance for the tangent random walk proposal   
  mu_z = c(1, 0, 0), #z are assumed to draw from a hyperbolid wrapped normal distribution with this centre
  alignment_method = "hpa", # post processing method, choose "none" or "angular" or "hpa"
  Z_ref = NULL, # Optional initial reference input, uses MDS if null
  allow_reflection = TRUE, 
  compute_barycenter = TRUE,
)

plot_poincare_embedding(fit_centred$Z_barycenter, Y = toy_adj)


# The wrapped normal clustering model example call
fit_clust <- Wrapped_Normal_Clustering_model(
  Y = toy_adj,
  K = 2,  # number of clusters
  Z_init = NULL,  # optional initial/reference configuration, uses MDS if null
  clustering_labels_init = NULL, #optional, assigns labels based on proximty to cluster centres if null
  mu_clust_init = NULL, #optional, assume evenly spaced clusters on hyperboloid if null
  lambda_init = NULL, #initial cluster probabilities, assigns even probabilities if null
  alpha_init = 0, 
  n_iter = 3000,
  burn = 1500,
  mu_alpha = 0,
  sigma2_alpha = 1,
  sigma2_prop_alpha = 0.5,
  sigma_z_prop = 0.15, # random walk proposal variance for z and mu
  sigma_mu_prop = 0.15,
  sigma2_clust_init = c(0.5, 0.5), # initial cluster variances
  sigma2_ig_shape = 3, # clusters assumed inv gamma a priori, shape and rate, make sure both are positive!
  sigma2_ig_rate = 1,
  sigma2_min = 1e-8, # lower bound on cluster varinces so your computer doesnt explode numerically
  v_lambda = c(1, 1), # lambda are assumed to be drawn from a dirichlet distribution, this is its prior 
  mu_prior_centre = c(1, 0, 0), # cluster means are assumed wrapped normal around the hyperboloid with this centre and variance
  sigma_mu_prior = 1,
  mu_init_r = 1.5, # initial radial distance of the cluster means from the centre of the hyperboloid
  align_method = "HPA", # post processing method, choose "none" or "polar" or "HPA"
  Z_ref = NULL, 
  relabel_by_angle = FALSE, # Initial attempt at relabbelling clusters by angle, performs poorly 
  relabel_reference = TRUE, # relabel by minimal proximity cost 
  compute_barycenter = TRUE, 
)

plot_poincare_embedding(fit_clust$Z_bary_mean, Y = toy_adj, clustering_labels = fit_clust$clustering_map)


