### Bayesian Additive Semi-Structured Regression Trees ###
### Bayesian additive model class                      ###

library(R6)
library(igraph)
library(deldir)

### Model class ----

# function to get log likelihood ratio
EvalLogLikeRatio <- function(move, e_m, ids_1, ids_2, sigmasq_y, sigmasq_mu) {
  sigma_ratio = sigmasq_y / sigmasq_mu
  csize_1 = length(ids_1); csize_2 = length(ids_2)
  sum_e_1 = sum(e_m[ids_1]); sum_e_2 = sum(e_m[ids_2])
  
  if(move == 'split') {
    logdetdiff = -0.5 * (log(csize_1+sigma_ratio) + log(csize_2+sigma_ratio) - 
                           log(csize_1+csize_2+sigma_ratio) - log(sigma_ratio))
    quaddiff = 0.5 * (sum_e_1^2/(csize_1+sigma_ratio) + sum_e_2^2/(csize_2+sigma_ratio) - 
                        (sum_e_1+sum_e_2)^2/(csize_1+csize_2+sigma_ratio)) / sigmasq_y
  }
  
  if(move == 'merge') {
    logdetdiff = -0.5 * (log(csize_1+csize_2+sigma_ratio) - log(csize_1+sigma_ratio) - 
                           log(csize_2+sigma_ratio) + log(sigma_ratio))
    quaddiff = 0.5 * ((sum_e_1+sum_e_2)^2/(csize_1+csize_2+sigma_ratio) - sum_e_1^2/(csize_1+sigma_ratio) - 
                        sum_e_2^2/(csize_2+sigma_ratio)) / sigmasq_y
  }
  return(logdetdiff + quaddiff)
}

# function to copy decision trees
CopyTree <- function(decision_root) {
  decision_root_new = decision_root$clone(deep = TRUE)
  # DFS
  stack = collections::stack()
  stack$push(list(decision_root_new, decision_root))
  while (stack$size() > 0) {
    node_new = stack$peek()[[1]]
    node = stack$peek()[[2]]
    stack$pop()
    if (!node$is_terminal) {
      # copy children
      node_new$children[[1]] = node$children[[1]]$clone(deep = TRUE)
      node_new$children[[2]] = node$children[[2]]$clone(deep = TRUE)
      
      stack$push(list(node_new$children[[2]], node$children[[2]]))
      stack$push(list(node_new$children[[1]], node$children[[1]]))
    }
  }
  return(decision_root_new)
}

# function to deep copy partitions
CopyPartitions <- function(partitions, keep_spanning_tree = FALSE) {
  partitions_copy = list()
  for(i in 1:length(partitions)) {
    partitions_copy = c(partitions_copy, partitions[[i]]$clone(deep = TRUE))
    partitions_copy[[i]]$decision_root = CopyTree(partitions[[i]]$decision_root)
    if(!keep_spanning_tree)
      partitions_copy[[i]]$spanning_tree = NULL
  }
  return(partitions_copy)
}

Model <- R6Class('Model', list(
  Y = NULL,
  X = NULL,
  coords = NULL,
  graphs = NULL,
  std_par = NULL,       # parameters to unstandardize Y
  X_cutpoints = NULL,   # cutting points of X (p * numcut matrix)
  projections = NULL,   # mapping from observed locations to knots
  n_knots = NULL,       # number of knots
  
  hyperpar = NULL,  # hyperparameters
  init_val = NULL,  # initial values
  
  ## for out-of-sample prediction
  X_new = NULL,
  coords_new = NULL,
  projections_new = NULL,  # array with dimension (nrow(X_new), nn, M)
  
  ## MCMC samples
  partition_out = NULL,
  sigmasq_y_out = NULL,
  g_out = NULL,
  # log_post_out = NULL,
  importance_out = NULL,
  # log_like_out = NULL,
  Y_new_out = NULL,
  
  ## MCMC settings
  MCMC = NULL,
  BURNIN = NULL,
  THIN = NULL,
  
  initialize = function(Y, X, graphs, projections, hyperpar,
                        X_new = NULL, projections_new = NULL) {
    self$Standardize(Y)
    self$hyperpar = hyperpar
    self$graphs = graphs
    self$projections = projections
    self$n_knots = sapply(graphs, vcount)
    
    self$X_new = X_new
    self$projections_new = projections_new
    
    for (m in 1:length(graphs)) {
      if('name' %in% names(vertex_attr(self$graphs[[m]]))) {
        self$graphs[[m]] = delete_vertex_attr(self$graphs[[m]], 'name')
      }
    }
    
    # find lambda_s
    nu = self$hyperpar['nu']; q = self$hyperpar['q']
    quant = qchisq(1-q, nu)
    self$hyperpar['lambda_s'] = quant * var(self$Y) / nu
    
    # get cutting points of X
    numcut = self$hyperpar['numcut']
    res = BART::bartModelMatrix(X, numcut, rm.const = T)
    self$X = res$X
    self$X_cutpoints = res$xinfo
    
    # remove constant covariate in test data
    if (!is.null(X_new))
      X_new = X_new[, res$rm.const]
  },
  
  # function to fit model by MCMC
  Fit = function(init_val, MCMC, BURNIN, THIN, seed = 1234, save_partitions = F) {
    set.seed(seed)
    
    self$init_val = init_val
    self$MCMC = MCMC
    self$BURNIN = BURNIN
    self$THIN = THIN
    
    # check if we need to do prediction
    prediction_flag = !(is.null(self$X_new) | is.null(self$projections_new))
    
    # set initial spanning tree as minimal spanning tree / forest
    if(is.null(self$init_val$spanning_tree)) {
      self$init_val$spanning_tree = list()
      for (m in 1:length(self$graphs))
        self$init_val$spanning_tree[[m]] = mst(self$graphs[[m]])
    }
    for (m in 1:length(self$graphs)) {
      self$graphs[[m]] = delete_edge_attr(self$graphs[[m]], 'weight')
      self$init_val$spanning_tree[[m]] = delete_edge_attr(self$init_val$spanning_tree[[m]], 'weight')
    }
    
    n = nrow(self$X)
    p = ncol(self$X)
    # hyper-parameters
    M = self$hyperpar['M']  # number of trees
    sigmasq_mu = self$hyperpar['sigmasq_mu']
    lambda_s = self$hyperpar['lambda_s']
    nu = self$hyperpar['nu']
    alpha = self$hyperpar['alpha']
    beta = self$hyperpar['beta']
    prob_split_by_x_0 = self$hyperpar['prob_split_by_x']
    
    # initialize
    sigmasq_y = init_val[['sigmasq_y']]
    partitions = list()
    numcut = apply(self$X_cutpoints, 1, function(x) max(which(!is.na(x))))  # effective number of cutpoints
    for(m in 1:M) {
      partitions[[m]] = Partition$new(self$init_val$spanning_tree[[m]],
                                      ncol(self$X), numcut,
                                      self$projections[, m], self$n_knots[m])
    }
    g = matrix(0, nrow = n, ncol = M) # n*M matrix of fitted mu's
    Y_hat = rep(0, n)
    
    cnt_split = rep(0, p + 1)  # number of splits by each covariate and spatial coords
    
    # setup result containers
    if (save_partitions)
      self$partition_out = list()
    self$sigmasq_y_out = numeric((MCMC-BURNIN)/THIN)
    self$g_out = array(0, dim = c((MCMC-BURNIN)/THIN, n, M))
    # self$log_like_out = numeric((MCMC-BURNIN)/THIN)
    self$importance_out = matrix(0, nrow = (MCMC-BURNIN)/THIN, ncol = p + 1)
    if (prediction_flag) {
      n_new = nrow(self$X_new)
      self$Y_new_out = matrix(0, nrow = (MCMC-BURNIN)/THIN, ncol = n_new)
    }
    
    ## Run MCMC 
    for(iter in 1:MCMC) {
      for(m in 1:M) {
        e_m = self$Y - (Y_hat - g[, m])
        
        # get splittable terminal nodes
        splittables = partitions[[m]]$GetSplittableNodes(type = 'all')
        
        # propose a move
        rs = 0; rhy = 0  # depreciated
        if(partitions[[m]]$n_terminals == 1) {
          rb = 1 - (rs + rhy)
          rd = 0
        } else if( length(splittables$x) == 0 & length(splittables$loc) == 0 ) {
          rb = 0
          rd = 1 - (rs + rhy)
        } else {
          rb = (1 - (rs + rhy)) / 2
          rd = rb
        }
        move = sample(4, 1, prob = c(rb, rd, rs, rhy))
        
        if(move == 1) { ## birth move
          # split a terminal node
          if(length(splittables$x) > 0 & length(splittables$loc) > 0) {
            prob_split_by_x = prob_split_by_x_0
          } else if(length(splittables$x) == 0) {
            prob_split_by_x = 0
          } else {
            prob_split_by_x = 1
          }
          
          split_proposal = partitions[[m]]$ProposeSplit(splittables, prob_split_by_x, self$X, self$X_cutpoints, 
                                                        alpha, beta, rb, rs, rhy, e_m, sigmasq_y, sigmasq_mu, 
                                                        self$projections[, m])
          acc_prob = split_proposal$acc_prob
          if(runif(1) < acc_prob){
            # accept
            partitions[[m]]$Split(split_proposal, self$projections[, m])
            
            # update covariate split counts
            split_rule = split_proposal$split_rule
            if (class(split_rule)[1] == 'SplitRuleX') {
              cnt_split[split_rule$xid] = cnt_split[split_rule$xid] + 1
            } else {
              cnt_split[p + 1] = cnt_split[p + 1] + 1
            }
          }
        }
        
        if(move == 2) { ## death move
          # merge a node
          merge_proposal = partitions[[m]]$ProposeMerge(alpha, beta, rd, rs, rhy, e_m, sigmasq_y, sigmasq_mu,
                                                        length(splittables$x), length(splittables$loc))
          acc_prob = merge_proposal$acc_prob
          if(runif(1) < acc_prob){
            # update covariate split counts
            split_rule = merge_proposal$node_merge$split_rule
            if (class(split_rule)[1] == 'SplitRuleX') {
              cnt_split[split_rule$xid] = cnt_split[split_rule$xid] - 1
            } else {
              cnt_split[p + 1] = cnt_split[p + 1] - 1
            }
            
            # accept
            partitions[[m]]$Merge(merge_proposal)
          }
        }
        
        # update mu_m and partial residual
        Y_hat = Y_hat - g[, m]
        g[, m] = partitions[[m]]$UpdateMu(e_m, sigmasq_y, sigmasq_mu)
        Y_hat = Y_hat + g[, m]
      }
      
      # update sigmasq_y
      rate = 0.5*(nu*lambda_s + sum((self$Y - Y_hat)^2))
      sigmasq_y = 1/rgamma(1, shape = (n+nu)/2, rate = rate)
      
      # save result
      if(iter > BURNIN & (iter - BURNIN) %% THIN == 0) {
        if (save_partitions)
          self$partition_out[[(iter-BURNIN)/THIN]] = CopyPartitions(partitions)
        self$sigmasq_y_out[(iter-BURNIN)/THIN] = sigmasq_y
        self$g_out[(iter-BURNIN)/THIN, , ] = g
        self$importance_out[(iter-BURNIN)/THIN, ] = cnt_split
        
        # out-os-sample prediction
        if (prediction_flag) {
          g_new = matrix(0, nrow = n_new, ncol = M)
          for(m in 1:M)
            g_new[, m] = partitions[[m]]$Predict2(self$X_new, self$projections_new[, , m])
          self$Y_new_out[(iter-BURNIN)/THIN, ] = rowSums(g_new)
        }
      }
      
      if(iter %% 100 == 0) 
        cat('Iteration', iter, 'done\n')
    }
    
    # post-processing
    colnames(self$importance_out) = c(colnames(self$X), 'spatial')
    self$sigmasq_y_out = self$Unstandardize(self$sigmasq_y_out, nomean = T, s2 = T)
    if (prediction_flag)
      self$Y_new_out = self$Unstandardize(self$Y_new_out)
    
  },
  
  # function to standardize Y
  Standardize = function(Y) {
    ymean = mean(Y)
    Y = Y - ymean
    yscale = 2 * max(abs(Y))
    self$Y = Y / yscale
    self$std_par = c('mean' = ymean, 'scale' = yscale)
  },
  
  # function to unstandardize
  Unstandardize = function(x, nomean = F, s2 = F) {
    if(s2) {
      x = x * self$std_par['scale'] ^ 2
    } else {
      x = x * self$std_par['scale']
    }
    if(!nomean) x = x + self$std_par['mean']
    return(x)
  },
  
  # function to compute WAIC
  ComputeWAIC = function(penalty = 2) {
    Y_post = apply(self$g_out, c(1, 2), sum)
    Y_post = self$Unstandardize(Y_post)
    Y = self$Unstandardize(self$Y)
    n = ncol(Y_post); n_post = nrow(Y_post)
    
    # compute log likelihood for each observation at each posterior sample
    residual_sq = t(apply(Y_post, 1, function(Y_hat) (Y - Y_hat) ^ 2))  # n_post * n
    log_like = -log(self$sigmasq_y_out)/2 - residual_sq / (2*self$sigmasq_y_out) - log(2*pi)/2
    
    # compute WAIC
    log_mean_like = apply(log_like, 2, logSumExp) - log(n_post)
    mean_log_like = colMeans(log_like)
    sum_log_mean_like = sum(log_mean_like)
    
    p_WAIC1 = 2 * sum(log_mean_like - mean_log_like)
    WAIC1 = -2 * sum_log_mean_like + penalty * p_WAIC1
    
    p_WAIC2 = sum( apply(log_like, 2, var) )
    WAIC2 = -2 * sum_log_mean_like + penalty * p_WAIC2
    
    return(c('WAIC1' = WAIC1, 'WAIC2' = WAIC2))
  },
  
  
  # function to predict
  # when arguments are NULL, return in-sample prediction
  # using knots for out-of-sample prediction
  Predict2 = function(X_new = NULL, projections_new = NULL, return_all = F, seed = 12345) {
    if(is.null(X_new) | is.null(projections_new)) {
      ## in-sample prediction using posterior means
      Y_post = apply(self$g_out, c(1, 2), sum)
      Y_hat = self$Unstandardize(Y_post)
      if (!return_all)
        Y_hat = colMeans(Y_hat)
      return(Y_hat)
    } else {
      ## out-of-sample prediction using posterior means
      if (length(self$partitions_out) == 0)
        stop("Partitions were not saved")
      set.seed(seed)
      
      M = self$hyperpar['M']
      n_post = length(self$sigmasq_y_out)
      g_new = array(0, dim = c(n_post, nrow(X_new), M))
      
      for(i in 1:n_post) {
        partitions = self$partition_out[[i]]
        for(m in 1:M)
          g_new[i, , m] = partitions[[m]]$Predict2(X_new, projections_new[, , m])
      }
      Y_new_post = apply(g_new, c(1, 2), sum)
      Y_new_hat = self$Unstandardize(Y_new_post)
      if (!return_all)
        Y_new_hat = colMeans(Y_new_hat)
      return(Y_new_hat)
    }
  }
))

### Utils -----

# function to compute log-sum-exp
logSumExp <- function(x) {
  x_max = max(x)
  return(x_max + log( sum(exp(x - x_max)) ))
}

