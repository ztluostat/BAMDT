### Simulation data on U-shape domain ###

rm(list=ls())
setwd('D:\\Documents\\_TAMU\\Projects\\Github\\BAMDT')
library(ggplot2)
source('ComplexDomainFun.R')

### generate data ---- 
set.seed(1234)

# generate U-shaped domain
rot_angle = 45
ubnd = genURegion(angle = rot_angle, n_line = 14)

# generate n uniform training locations
n = 500
coords = genLocations(n, ubnd)

# generate n_ho uniform hold-out locations
n_ho = 200
coords_ho = genLocations(n_ho, ubnd)

# estimate geodesic distance between observed/testing locations
# NOTE: this may take a while
coords_all = rbind(coords, coords_ho)
dist_res = gdist(coords_all, ubnd)
dist_mat_all = dist_res[1:(n + n_ho), 1:(n + n_ho)]
rm(dist_res)

# generate X from GP
set.seed(1234)
p = 10
X_all = matrix(0, nrow = n + n_ho, ncol = p)
X_all = simGpU(p, coords_all, angle = rot_angle)
colnames(X_all) = paste("X", 1:p, sep = "")

X = X_all[1:n, , drop = F]
X_ho = X_all[-(1:n), , drop = F]

# check what X[, 1] looks like
ggplot() + 
  geom_boundary(ubnd) +
  geom_point(aes(x = lon, y = lat, color = X[, 1]), data = as.data.frame(coords)) +
  scale_color_gradientn(colours = rainbow(5))  +
  labs(x = 's_1', y = 's_2', title = 'X[, 1]')

# generate true function
f_true = evalFunU(coords, X, angle = rot_angle)
f_ho_true = evalFunU(coords_ho, X_ho, angle = rot_angle)

# true partition
cluster_true = ifelse(coords[, 1] ^ 2 + coords[, 2] ^ 2 <= 0.9^2, 1, 0)
idx = cluster_true == 0
cluster_true[idx] = ifelse(coords[idx, 1] > coords[idx, 2], 2, 3)

cluster_ho_true = ifelse(coords_ho[, 1] ^ 2 + coords_ho[, 2] ^ 2 <= 0.9^2, 1, 0)
idx = cluster_ho_true == 0
cluster_ho_true[idx] = ifelse(coords_ho[idx, 1] > coords_ho[idx, 2], 2, 3)

# modify true function to create discontinuity
f_true[cluster_true == 3] = f_true[cluster_true == 3] - 4
f_true[cluster_true == 2] = f_true[cluster_true == 2] + 4
f_true[cluster_true == 1] = -0.5 * f_true[cluster_true == 1]

f_ho_true[cluster_ho_true == 3] = f_ho_true[cluster_ho_true == 3] - 4
f_ho_true[cluster_ho_true == 2] = f_ho_true[cluster_ho_true == 2] + 4
f_ho_true[cluster_ho_true == 1] = -0.5 * f_ho_true[cluster_ho_true == 1]

# check what the true function looks like
ggplot() + 
  geom_boundary(ubnd) +
  geom_ellipse(0.9, 0.9, 0, 0) +
  geom_point(aes(x = lon, y = lat, col = f_true), data = as.data.frame(coords)) +
  scale_color_gradientn(colours = rainbow(5), name = 'True f') +
  labs(x = 's_1', y = 's_2', title = 'True function')


# generate replicate data sets
n_rep = 50
Y_rep = matrix(0, nrow = n_rep, ncol = n)
Y_ho_rep = matrix(0, nrow = n_rep, ncol = n_ho)

for (i in 1:n_rep) {
  Y = f_true + rnorm(n, 0, 0.1)
  Y_ho = f_ho_true + rnorm(n_ho, 0, 0.1)
  
  Y_rep[i, ] = Y
  Y_ho_rep[i, ] = Y_ho
}

# save data
save(Y_rep, Y_ho_rep, X, coords, n, p, ubnd,
     X_ho, coords_ho, n_ho, dist_mat_all,
     file = "input_U.RData", compress = "bzip2")
