### Functions related to complex domains ###

library(fdaPDE)
library(igraph)
library(sf)

### Functions for 2-d constrained domain -----

# function to rotate locations counterclockwise (in degrees)
rotate <- function(x, y, angle = 0) {
  angle = angle / 180 * pi
  x_rot = x * cos(angle) - y * sin(angle)
  y_rot = x * sin(angle) + y * cos(angle)
  return(list(x = x_rot, y = y_rot))
}

# function to rotate locations counterclockwise (in degrees)
rotate2 <- function(coords, angle = 0) {
  res = rotate(coords[, 1], coords[, 2], angle)
  coords[, 1] = res$x; coords[, 2] = res$y
  return(coords)
}

# function to check if points are inside a polygon
insidePolygon <- function(bnd, x, y) {
  bnd$x = round(bnd$x, 9)
  bnd$y = round(bnd$y, 9)
  pg = st_polygon(list( cbind(bnd$x, bnd$y) ))
  
  coords = as.data.frame(cbind(x, y))
  pt = st_as_sf(coords, coords = c(1, 2))
  return( st_intersects(pt, pg, sparse = F)[, 1] )
}

# function to generate U-shape region (adapted from R package "mgcv")
# allowing for rotation counterclockwise (in degrees)
genURegion <- function (r0 = 0.1, r = 0.5, l = 3, n_theta = 20, n_line = 30, angle = 0) {
  rr = r + (r - r0)
  theta <- seq(pi, pi/2, length.out = n_theta)
  x = rr * cos(theta)
  y = rr * sin(theta)
  
  x = c(x, seq(l/n_line, l - l/n_line, length.out = n_line))
  y = c(y, rep(rr, n_line))
  
  theta = seq(pi/2, -pi/2, length = n_theta)
  x = c(x, (r - r0) * cos(theta) + l)
  y = c(y, (r - r0) * sin(theta) + r)
  
  x = c(x, seq(l - l/n_line, l/n_line, length.out = n_line))
  y = c(y, rep(r0, n_line))
  
  theta = seq(pi/2, pi, length.out = round(n_theta * r0 / (2*r - r0)))
  x = c(x, r0 * cos(theta))
  y = c(y, r0 * sin(theta))
  n = length(x)
  x = c(x, x[n:1])
  y = c(y, -y[n:1])
  
  # rotate
  rot = rotate(x, y, angle)
  return(rot)
}

# function to generate uniform locations in a region
genLocations <- function(n, bnd = NULL) {
  if(is.null(bnd)) {
    # generate locations in a unit square
    coords = cbind(runif(n), runif(n))
  } else {
    x_min = min(bnd$x); x_max = max(bnd$x)
    y_min = min(bnd$y); y_max = max(bnd$y)
    
    x = runif(3*n, x_min, x_max); y = runif(3*n, y_min, y_max)
    idx_in = insidePolygon(bnd, x, y)
    coords = cbind(x[idx_in], y[idx_in])
    
    if(nrow(coords) >= n) {
      coords = coords[1:n, ]
    } else {
      for(i in 1:(n - nrow(coords))) {
        x = runif(1, x_min, x_max); y = runif(1, y_min, y_max)
        while(!insidePolygon(bnd, x, y)) 
          x = runif(1, x_min, x_max); y = runif(1, y_min, y_max)
        coords = rbind(coords, c(x, y))
      }
    }
  }
  colnames(coords) = c('lon', 'lat')
  return(coords)
}

# function to generate grids in a region
genGrids <- function(nx, ny, bnd = NULL) {
  if(is.null(bnd)) {
    # generate locations in a unit square
    coords_grid = expand.grid(seq(0, 1, length.out = nx), seq(0, 1, length.out = ny))
    coords_grid = as.matrix(coords_grid)
  } else {
    x_min = min(bnd$x); x_max = max(bnd$x)
    y_min = min(bnd$y); y_max = max(bnd$y)
    
    # to avoid points on boundaries
    x_buffer = (x_max - x_min) * 1e-5
    y_buffer = (y_max - y_min) * 1e-5
    
    grid_x = seq(x_min + x_buffer, x_max - x_buffer, length.out = nx)
    grid_y = seq(y_min + y_buffer, y_max - y_buffer, length.out = ny)
    coords_grid = as.matrix(expand.grid(grid_x, grid_y))
    
    x = coords_grid[, 1]; y = coords_grid[, 2]
    idx_in = insidePolygon(bnd, x, y)
    coords_grid = coords_grid[idx_in, ]
  }
  colnames(coords_grid) = c('lon', 'lat')
  return(coords_grid)
}

# function to create a 2d mesh
gen2dMesh <- function(coords, bnd, ...) {
  # note the first and last boundary nodes are the same
  n = nrow(coords); n_bnd = length(bnd$x) - 1
  coords_all = rbind(coords, cbind(bnd$x, bnd$y)[1:n_bnd, ])
  
  # get boundary segments
  segments = cbind( (n+1):(n+n_bnd), c((n+2):(n+n_bnd), n+1) )
  
  mesh = create.mesh.2D(coords_all, segments = segments, ...)
  mesh$n_int = n  # number of interior nodes
  # edges that connect boundary nodes
  mesh$bnd_edges = apply(mesh$edges, 1, FUN = function(x) any(x > n))
  return(mesh)
}

# function to obtain a constrained Delaunay triangulation graph from a mesh
constrainedDentri <- function(n, mesh, threshold = 1e6, crs = NULL,
                              gaurantee_connected = F, dist_mat = NULL) {
  coords = mesh$nodes[1:n, ]
  colnames(coords) = c("X", "Y")
  
  # drop edges that connect boundary nodes
  rid_drop = mesh$bnd_edges
  edge_list = mesh$edges[!rid_drop, ]
  
  # compute edge length
  if (!is.null(dist_mat)) {
    distance = dist_mat[edge_list]
  } else {
    if (is.null(crs)) {
      distance = sqrt( rowSums((coords[edge_list[, 1], ] - coords[edge_list[, 2], ]) ^ 2) )
    } else {
      coords_sf = st_as_sf(as.data.frame(coords), coords = c("X", "Y"), crs = crs)
      coords_sf = st_geometry(coords_sf)
      distance = st_distance(coords_sf[edge_list[, 1]], coords_sf[edge_list[, 2]], 
                             by_element = T)
      distance = as.numeric(distance)
    }
  }
  
  rid_drop = distance > threshold
  edge_list = edge_list[!rid_drop, ]
  distance = distance[!rid_drop]
  
  adj_mat = matrix(0, nrow = n, ncol = n)
  adj_mat[edge_list] = distance
  adj_mat = adj_mat + t(adj_mat)
  graph0 = graph_from_adjacency_matrix(adj_mat, "undirected", weighted = "weight")
  
  if (gaurantee_connected & !is.null(dist_mat)) {
    graph0 = connectGraph(graph0, dist_mat)
    E(graph0)$weight = dist_mat[as_edgelist(graph0)]
  }
    
  return(graph0)
}

# function to connect each component to its nearest neighbor such that
# we can obtain a connected graph from a disconnected one
connectGraph <- function(graph, dist_mat) {
  n = vcount(graph)
  connected_comp = components(graph)
  memberships = connected_comp$membership
  n_comp = connected_comp$no
  if (n_comp == 1)
    return(graph)
  
  edge_list = as_edgelist(graph)
  for (i in 1:(n_comp - 1)) {
    # merge one component to its nearest component
    idx_1 = memberships == 1
    idx_2 = memberships != 1
    cdist_mat = dist_mat[idx_1, idx_2, drop = F]
    idx_min = which(cdist_mat == min(cdist_mat), arr.ind = TRUE)[1, ]
    
    # find the pair of vertices that has minimum distance
    vid_1 = c(1:n)[idx_1][ idx_min[1] ]
    vid_2 = c(1:n)[idx_2][ idx_min[2] ]
    
    # connect vid_1 with vid_2
    edge_new = sort(c(vid_1, vid_2))
    edge_list = rbind(edge_list, edge_new)
    
    # merge two components where vid_1 and vid_2 lie in
    comp_id_merged = memberships[vid_2]
    memberships[memberships == comp_id_merged] = 1
  }
  
  return(graph_from_edgelist(edge_list, directed = F))
}

# function to generate equally spaced points along a given line
refineLine <- function(start, end, n) {
  grids_x = seq(start[1], end[1], length.out = n)
  grids_y = seq(start[2], end[2], length.out = n)
  grids = cbind(grids_x, grids_y)
  return(grids)
}


# function to estimate geodesic distance on a 2-d constrained domain
# by using shortest path length between two vertices on a dense graph
gdist <- function(coords, bnd, nx = 50, ny = 50, k_nn = 8, 
                  coords2 = NULL, crs = NULL, return_both = F) {
  require(igraph)
  require(fields)
  require(FNN)
  
  # generate grids
  coords_grids = genGrids(nx, ny, bnd)
  n_grids = nrow(coords_grids)
  
  # get nearest neighbor graph
  coords_all = rbind(coords_grids, coords)
  if(!is.null(coords2))
    coords_all = rbind(coords_all, coords2)
  if (is.null(crs)) {
    nn_res_grid = get.knn(coords_grids, k = k_nn)
    nn_res_nongrid = get.knnx(coords_grids, coords_all[-(1:n_grids), ], k = k_nn)
  } else {
    coords_all_sf = st_as_sf(as.data.frame(coords_all), coords = c('lon', 'lat'), crs = crs)
    coords_all_sf = st_geometry(coords_all_sf)
    nn_res_grid = st_knn(coords_all_sf[1:n_grids], k = k_nn)
    nn_res_nongrid = st_knn(coords_all_sf[1:n_grids], coords_all_sf[-(1:n_grids)], k = k_nn)
  }
  knngraph = constrainedKNN(nn_res_grid, nn_res_nongrid, coords_all, bnd)
  
  # check if the graph is connected
  if(igraph::components(knngraph)$no != 1)
    stop("Disconnected kNN graph; try larger 'k'")
  
  # compute geodesic distance matrix
  from = (n_grids + 1):(n_grids + nrow(coords))
  n1 = ifelse(is.null(coords2), 0, nrow(coords))
  to = (n_grids + n1 + 1):nrow(coords_all)
  cdist_mat = igraph::distances(knngraph, v = from, to = to)
  
  if(return_both & !is.null(coords2)) {
    # return distance matrix of coords as well
    from = (n_grids + 1):(n_grids + nrow(coords))
    dist_mat = igraph::distances(knngraph, v = from, to = from)
    return(list('dist' = dist_mat, 'cdist' = cdist_mat))
  } else {
    # return cross distance matrix only
    return(cdist_mat)
  }
}

# function to obtain intrinsic coordinates on a rotated U-shape domain
# adapted from R package "mgcv"
intrinsicCoordsU <- function(coords, r0 = 0.1, r = 0.5, l = 3, angle = 0) {
  # rotate to horizontal
  coords = rotate2(coords, -angle)
  
  n = nrow(coords)
  q = pi * r/2 # 1/2 length of semi-circle part of centre curve
  a = rep(0, n); d = rep(0, n) # along and distance to arrays
  x = coords[, 1]; y = coords[, 2]
  
  ## convert x,y to along curve and distance to curve (a,d) 
  ## co-ordinates. 0 distance along is at (x = -r, y = 0)  
  
  ind = (x >= 0) & (y > 0)
  a[ind] = q + x[ind]
  d[ind] = y[ind]-r
  
  ind = (x >= 0) & (y <= 0)
  a[ind] = -q - x[ind]
  d[ind] = -r - y[ind]
  
  ind = (x < 0) 
  a[ind] = -atan(y[ind] / x[ind]) * r
  d[ind] = sqrt(x[ind]^2 + y[ind]^2) - r
  
  return(cbind(a, d))
}

# function to simulate a GP on a U-shape domain using intrinsic coordinates
simGpU <- function(n, coords, r0 = 0.1, r = 0.5, l = 3, angle = 0) {
  # get intrinsic coordinate (a, d)
  coords_in = intrinsicCoordsU(coords, r0, r, l, angle)
  
  # length scale
  phi_a = 1; phi_d = 1
  coords_in[, 1] = coords_in[, 1] / phi_a
  coords_in[, 2] = coords_in[, 2] / phi_d
  # get covariance matrix based on (a, d)
  dist_mat = as.matrix(dist(coords_in))
  corr = exp(-dist_mat)
  varcov = corr
  chol_varcov = t(chol(varcov))
  
  # simulate GP
  X = matrix(0, nrow = nrow(coords_in), ncol = n)
  for (j in 1:n)
    X[, j] = chol_varcov %*% rnorm(nrow(coords_in))
  return(X)
}

# evaluate test function on a horizontal U-shape domain
evalFunU <- function(coords, X, r0 = 0.1, r = 0.5, l = 3, angle = 0) {
  coords_in = intrinsicCoordsU(coords, r0, r, l, angle)
  a = coords_in[, 1]; d = coords_in[, 2]
  f = a * X[, 1] + d ^ 2
  return(f)
}

### Constrained KNN/RNN on 2-d constrained domain -----

# helper function to check if 3 points are in counterclockwise order
ccw <- function(A, B, C) {
  return( (C[2] - A[2]) * (B[1] - A[1]) > (B[2] - A[2]) * (C[1] - A[1]) )
}

# function to ckeck if two segments AB and CD intersect
segmentIntersect <- function(A, B, C, D) {
  return( ccw(A, C, D) != ccw(B, C, D) & ccw(A, B, C) != ccw(A, B, D) )
}

# helper function to get RNN list for one location
# return idx in d that is less than r
# if no idx is found, return which.min(d)
getRNNList <- function(d, r) {
  idx = which(d <= r)
  if(length(idx) == 0)
    idx = which.min(d)
  return(idx)
}

# function to obtain constrained RNN given geodesic distances
constrainedRNN <- function(dist_mat, r = 1, bnd = NULL, return_graph = F) {
  n1 = nrow(dist_mat); n2 = ncol(dist_mat)
  nn_list = apply(dist_mat, 1, getRNNList, r = r)
  
  if (return_graph & n1 == n2) {
    require(igraph)
    require(Matrix)
    i_all = c(); j_all = c(); x_all = c()
    
    for(i in 1:n1) {
      for(idx in 1:length(nn_list[[i]])) {
        j = nn_list[[i]][idx]
        if (i == j) next
        i_all = c(i_all, i)
        j_all = c(j_all, j)
        x_all = c(x_all, dist_mat[i, j])
      }
    }
    adj_mat = sparseMatrix(i = i_all, j = j_all, x = x_all, dims = c(n1, n2))
    # get rnn graph
    rnngraph = graph_from_adjacency_matrix(adj_mat, mode = 'directed', weighted = T)
    rnngraph = as.undirected(rnngraph, mode = 'collapse', edge.attr.comb = 'first')
    
    # check if the graph is connected
    if(components(rnngraph)$no > 1)
      warning("Disconnected RNN graph; try larger 'r'")
    
    return(rnngraph)
  } else {
    return(nn_list)
  }
}

# function to get a KNN graph given a distance matrix
KNNGraph <- function(dist_mat, k_nn = 5, cross_dist = F, return_graph = T) {
  
  n1 = nrow(dist_mat); n2 = ncol(dist_mat)
  if(cross_dist) {
    # dist_mat is cross distance matrix
    adj_list = apply(dist_mat, 1, function(x) order(x)[1:k_nn])
  } else {
    adj_list = apply(dist_mat, 1, function(x) order(x)[2:(k_nn+1)])
  }
  adj_list = t(adj_list)
  
  if(return_graph & n1 == n2) {
    require(igraph)
    require(Matrix)
    
    i_all = c(); j_all = c(); x_all = c()
    for(i in 1:n1) {
      for(cidx in 1:k_nn) {
        i_all = c(i_all, i)
        j_all = c(j_all, adj_list[i, cidx])
        x_all = c(x_all, dist_mat[i, adj_list[i, cidx] ])
      }
    }
    adj_mat = sparseMatrix(i = i_all, j = j_all, x = x_all, dims = c(n1, n2))
    # get knn graph
    knngraph = graph_from_adjacency_matrix(adj_mat, mode = 'directed', weighted = T)
    knngraph = as.undirected(knngraph, mode = 'collapse', edge.attr.comb = 'first')
    return(knngraph)
  } else {
    return(lapply(1:n1, function(i) adj_list[i, ]))
  }
}

# function to build constrained kNN graph on a 2-d constrained domain given NNs
# remove NN pairs if the corresponding edge across the boundary
# nn_res: list returned by FNN::get.knn()
constrainedKNN <- function(nn_res_grid, nn_res_nongrid, coords, bnd) {
  require(igraph)
  require(Matrix)
  
  # get boundary segments
  n_bnd = length(bnd$x)
  bnd_segments = cbind(bnd$x[-n_bnd], bnd$y[-n_bnd], bnd$x[-1], bnd$y[-1])
  
  # combine nn_res
  nn_res = list(nn.index = rbind(nn_res_grid$nn.index, nn_res_nongrid$nn.index),
                nn.dist = rbind(nn_res_grid$nn.dist, nn_res_nongrid$nn.dist))
  
  # get sparse adjacency matrix
  # and remove edges across the boundary
  i_all = c(); j_all = c(); x_all = c()
  n = nrow(nn_res$nn.index)
  k_nn = ncol(nn_res$nn.index)  # number of nn
  for(i in 1:n) {
    coords_i = coords[i, ]
    cnt_nn = 0  # count how many valid nn (i.e., not crossing boundary)
    for(cidx in 1:k_nn) {
      j = nn_res$nn.index[i, cidx]
      coords_j = coords[j, ]
      
      # check if segment (i, j) crosses domain boundary
      intersect_bnd = apply(bnd_segments, 1, function(seg) segmentIntersect(
        coords_i, coords_j, seg[1:2], seg[3:4]
      ))
      
      if(!any(intersect_bnd)) {
        # add to adj matrix
        i_all = c(i_all, i); j_all = c(j_all, j)
        x_all = c(x_all, nn_res$nn.dist[i, cidx])
        cnt_nn = cnt_nn + 1
      }
    }
    
    # check if the current node has any neighbor
    if(cnt_nn == 0)
      stop("Disconnected kNN graph; try larger 'k'")
  }
  
  adj_mat = sparseMatrix(i = i_all, j = j_all, x = x_all, dims = c(n, n))
  
  # get knn graph
  knngraph = graph_from_adjacency_matrix(adj_mat, mode = 'directed', weighted = T)
  knngraph = as.undirected(knngraph, mode = 'collapse', edge.attr.comb = 'first')
  return(knngraph)
}


### Plotting functions -----

# function to plot complex domain
geom_boundary <- function(bnd, ...) {
  n = length(bnd$x)
  segments = cbind(bnd$x[-n], bnd$y[-n], bnd$x[-1], bnd$y[-1])
  segments = data.frame(segments)
  names(segments) = c('x1', 'y1', 'x2', 'y2')
  return(geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = segments, ...))
}

# function to plot ellipse
geom_ellipse <- function(rx, ry, xc, yc, color = "red", size = 0.3, ...) {
  x = xc + rx * cos(seq(0, pi, length.out=200))
  ymax = yc + ry * sin(seq(0, pi, length.out=200))
  ymin = yc + ry * sin(seq(0, -pi, length.out=200))
  annotate("ribbon", x = x, ymin = ymin, ymax = ymax, color = color, size = size, fill = NA, ...)
}

# plotting functions from SCC
plotGraph <- function(coords, graph, title = NULL){
  require(ggplot2)
  edgelist = get.edgelist(graph) 
  edgedata = data.frame(coords[edgelist[,1 ], ], coords[edgelist[, 2], ])
  colnames(edgedata) = c("X1","Y1","X2","Y2")
  
  ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edgedata, size = 0.5, colour="grey")+
    labs(title = title, x = "lon", y = "lat")+
    theme(plot.title = element_text(hjust = 0.5))
}

plotGraphData <- function(coords, graph, Data, title = NULL, col_lim = NULL){
  require(ggplot2)
  
  edgelist = get.edgelist(graph) 
  edgedata = data.frame(coords[edgelist[, 1], ], coords[edgelist[, 2], ])
  colnames(edgedata) = c("X1","Y1","X2","Y2")
  if(missing(col_lim)) col_lim = NA
  mtsub = data.frame(coords[index, ]);colnames(mtsub) = c('X1', 'Y1');
  ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2), data = edgedata, size = 0.5, colour="grey") + 
    geom_point(data=data.frame(coords), aes(lon, lat, color = Data)) +
    scale_color_gradientn(colours = rainbow(5)) +
    ggtitle(title)+
    theme(legend.title =element_blank(),plot.title = element_text(hjust = 0.5))+
    geom_text(data = mtsub, aes(x = X1,y=Y1, label =
                                  rownames(mtsub)), hjust = 0, size = 4)
    
}

### Util functions -----

# function to get knn of sfc objects
st_knn <- function(data, query = NULL, k = 1) {
  require(sf)
  if (is.null(query))
    dist_mat = st_distance(data)
  else
    dist_mat = st_distance(query, data)
  dist_mat = as.matrix(dist_mat)
  nn_index = t(apply(dist_mat, 1, WhichNMin, n = k))
  nn_dist = t(sapply(1:nrow(nn_index), function(i) dist_mat[i, nn_index[i, ]]))
  return(list('nn.index' = nn_index, 'nn.dist' = nn_dist))
}


# function to obtain indices of n smallest values in a vector
WhichNMin <- function(x, n = 1) {
  n = min(n, length(x))
  threshold = sort(x, partial = n)[n]
  which_n_min = which(x <= threshold)
  # deal with ties
  if(length(which_n_min) > n) {
    x = x[which_n_min]
    which_n_min = which_n_min[order(x)[1:n]]
  }
  return(which_n_min)
}
