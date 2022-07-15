##### Demo code for BAMDT ######

rm(list=ls())
setwd('D:\\Documents\\_TAMU\\Projects\\Github\\BAMDT')
library(ggplot2)
source('ComplexDomainFun.R')
source('Trees.R')
source('Model.R')

### Get model inputs -----

load("input_U.RData")

# get initial values
init_val = list()
init_val[['sigmasq_y']] = 1  # initial value for noise variance

nu = 3; q = 0.9
M = 50
hyperpar = c()
hyperpar['sigmasq_mu'] = (0.5/(2*sqrt(M)))^2          # variance of prior for mu
hyperpar['q'] = q                                     # quantile in prior for noise variance
hyperpar['nu'] = nu                                   # degree of freedom in prior for noise variance
hyperpar['M'] = M                                     # number of trees
hyperpar['alpha'] = 0.95                              # hyperparameter in tree generating process
hyperpar['beta'] = 2                                  # hyperparameter in tree generating process
hyperpar['numcut'] = 100                              # number of candidate cutting points for x                       
hyperpar['prob_split_by_x'] = min(p / (p + 2), 0.85)  # probability for a split on x

# make graph on spatial knots
n_knot = rep(100, M)  # numver of knots for each tree
coords_knots = list() # knot coordinates for each tree
graphs = list()       # spatial graph for each tree
knot_idx = list()     # indices of knots in training locations
set.seed(1234)
for (m in 1:M) {
  # subsample training locations as knots
  knot_idx[[m]] = sample.int(n, n_knot[m])
  coords_knots[[m]] = coords[knot_idx[[m]], ]
  dist_mat_knot = dist_mat_all[knot_idx[[m]], knot_idx[[m]]]
  
  # get CDT graph on knots
  mesh = gen2dMesh(coords_knots[[m]], ubnd)
  graph0 = constrainedDentri(n_knot[m], mesh, 
                             gaurantee_connected = T, dist_mat = dist_mat_knot)
  E(graph0)$eid = as.integer(1:ecount(graph0))  # edge id
  V(graph0)$vid = as.integer(1:vcount(graph0))  # vertex id
  graphs[[m]] = graph0
}

# ensuring spatial graphs are connected
for (m in 1:M) {
  if (components(graphs[[m]])$no != 1)
    stop(paste("Disconnected graph:", m))
}

# assign observations to their nearest knots
# projections[i, j] is the index of the nearest knot of obs. i in weak learner j
# similirly, projections_ho is for hold-out locations
projections = array(0, dim = c(n, M))
projections_ho = array(0, dim = c(n_ho, M))
for (m in 1:M) {
  # get distance between training locations and knots
  cdist_mat_knots = dist_mat_all[ 1:n, knot_idx[[m]] ]
  projections[, m] = apply(cdist_mat_knots, 1, which.min)
  
  # get distance between hold-out locations and knots
  cdist_mat_ho_knots = dist_mat_all[ (n + 1):(n + n_ho), knot_idx[[m]] ]
  projections_ho[, m] = apply(cdist_mat_ho_knots, 1, which.min)
}


### Fit BAMDT model ------

Y = Y_rep[1, ]; Y_ho = Y_ho_rep[1, ]

# MCMC iterations
MCMC = 1000
BURNIN = 500
THIN = 5

# train BAMDT and predict for hold-out data
# NOTE: this may take a while
model = Model$new(Y, X, graphs, projections, hyperpar, X_ho, projections_ho)
model$Fit(init_val, MCMC, BURNIN, THIN, seed = 12345)

# get posterior mean prediction
Y_ho_hat = colMeans(model$Y_new_out)

# plot prediction
ggplot() + 
  geom_boundary(ubnd) +
  geom_point(aes(x = lon, y = lat, col = Y_ho_hat),
             data = as.data.frame(coords_ho)) +
  scale_color_gradientn(colors = rainbow(5)) +
  ggtitle('Prediction of BAMDT')

# check feature importance metric
colMeans(model$importance_out)
