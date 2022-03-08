# This file is part of metanetwork

# metanetwork is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# metanetwork is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with metanetwork.  If not, see <http://www.gnu.org/licenses/>


#' Default configuration for the diffusion kernel based t-sne
#'
#' A list with parameters customizing configuration for the diffusion kernel based t-sne (see 'tsne' R package documentation)
#'
#' @examples
#' # display all default settings
#' TL_tsne.default
#'
#' # create a new settings object with n_neighbors set to 5
#' TL_tsne.custom = TL_tsne.default
#' TL_tsne.custom$max_iter = 5
#' TL_tsne.custom
#' 
#' @export
TL_tsne.default = list(
  max_iter = 300,
  min_cost = 0,
  epoch_callback = NULL,
  epoch = 100,
  momentum = 0.01,
  final_momentum = 0.3,
  mom_switch_iter = 250,
  epsilon = 2,
  min_gain = 0.01,
  initial_P_gain = 4,
  eps = 2^(-52)
)

#to do: print method for the class metanetwork.config
class(TL_tsne.default) = 'metanetwork_config'

#computing diffusion kernel
compute_diffusion_kernel = function(g,beta){
  n = length(igraph::V(g))
  if(is.null(igraph::E(g)$weight)){
    A = igraph::get.adjacency(g)
    u  = igraph::degree(g)
  }else{
    A = igraph::get.adjacency(g,attr = "weight")
    u  = igraph::strength(g)
  }
  H = A + Matrix::t(A) - diag(u) 
  eig_H = eigen(H)
  eig_H_values = eig_H$values
  eig_H_vectors = eig_H$vectors
  exp_values_vec = exp(beta*eig_H_values)
  K = eig_H_vectors %*% diag(exp_values_vec) %*% t(eig_H_vectors)
  return(K)
}

#compute TL-tsne layout
get_nodes_position_TL_tsne <- function(g,TL,TL_tsne.config,beta){
  max_iter = TL_tsne.config$max_iter
  min_cost = TL_tsne.config$min_cost
  epoch_callback = TL_tsne.config$epoch_callback
  epoch = TL_tsne.config$epoch
  momentum = TL_tsne.config$momentum
  final_momentum = TL_tsne.config$final_momentum
  mom_switch_iter = TL_tsne.config$mom_switch_iter
  epsilon = TL_tsne.config$epsilon
  min_gain = TL_tsne.config$min_gain
  initial_P_gain = TL_tsne.config$initial_P_gain
  eps = TL_tsne.config$eps
  message(paste0('beta = ',beta))
  
  if(is.null(igraph::E(g)$weight)){
    X = igraph::get.adjacency(g)
  }else{
    X = igraph::get.adjacency(g,attr = "weight")
  }

  #adapted from tsne function (R package 'tsne')
  # initial_config = NULL
  k = 2
  # initial_dims = 10#30
  # whiten = TRUE
  n = nrow(X)
  ydata = matrix(rnorm(k * n), n)
  ydata[,1] = TL
  P = compute_diffusion_kernel(g,beta)
  P = abs(P)
  grads = matrix(0, nrow(ydata), ncol(ydata))
  incs = matrix(0, nrow(ydata), ncol(ydata))
  gains = matrix(1, nrow(ydata), ncol(ydata))
  for (iter in 1:max_iter) {
    if (iter%%epoch == 0) {
      cost = sum(apply(P * log((P + eps)/(Q + eps)), 1, 
                       sum))
      message("Epoch: Iteration #", iter, " error is: ", 
              cost)
      if (cost < min_cost) 
        break
      if (!is.null(epoch_callback)) 
        epoch_callback(ydata)
    }
    sum_ydata = apply(ydata^2, 1, sum)
    num = 1/(1 + sum_ydata + sweep(-2 * ydata %*% t(ydata), 
                                   2, -t(sum_ydata)))
    diag(num) = 0
    Q = num/sum(num)
    if (any(is.nan(num))) 
      message("NaN in grad. descent")
    Q[Q < eps] = eps
    stiffnesses = 4 * (P - Q) * num
    for (i in 1:n) {
      grads[i, ] = apply(sweep(-ydata, 2, -ydata[i, ]) * 
                           stiffnesses[, i], 2, sum)
    }
    gains = ((gains + 0.2) * abs(sign(grads) != sign(incs)) + 
               gains * 0.8 * abs(sign(grads) == sign(incs)))
    gains[gains < min_gain] = min_gain
    incs = momentum * incs - epsilon * (gains * grads)
    ydata[,2] = ydata[,2] + incs[,2]
    ydata = sweep(ydata, 2, apply(ydata, 2, mean))
    if (iter == mom_switch_iter) 
      momentum = final_momentum
    if (iter == 100)
      P = P/4
  }
  return(ydata)
}

# get TL-kpco layout
get_nodes_position_TL_kpco <- function(g,TL,beta){
  if(igraph::is.connected(g,mode = "weak")){
    if(igraph::vcount(g)>2){
      K = compute_diffusion_kernel(g,beta)
      nf_loc = dim(K)[1] - 1
      pco1 = ade4::dudi.pco(dist(K),scannf = FALSE, nf = nf_loc)
      pco2 = ade4::pcaivortho(pco1,TL, scan = F,nf = nf_loc)
      coords = cbind(TL,rowMeans(pco2$li))
      return(coords)
    } else{
      coords = cbind(TL,rep(0,igraph::vcount(g)))
      return(coords)
    }
  } else{
    coords_list = list() #stock coords for each connected component
    membership_loc = igraph::components(g)$membership
    for(comp in unique(membership_loc)){
      #nodes belonging to comp
      nodes_comp = igraph::V(g)[names(which(membership_loc == comp))]
      g_comp_loc = igraph::induced_subgraph(g,nodes_comp)
      coords_comp = get_nodes_position_TL_kpco(g = g_comp_loc,
                                               TL = igraph::V(g_comp_loc)$TL ,beta = beta)
      rownames(coords_comp) = igraph::V(g_comp_loc)$name
      #avoiding overlay of different components
      if(comp>1){
        coords_comp[,2] = coords_comp[,2] - min(coords_comp[,2]) + 
          max(coords_list[[comp-1]][,2]) #+ sd(coords_list[[comp-1]][,2])
      }
      coords_list = c(coords_list,list(coords_comp))
    }
    coords = do.call(rbind,coords_list)
    return(coords[igraph::V(g)$name,])
  }
}

#save TL-tsne layout as node attribute
attach_layout_g <- function(g,metanetwork,mode = 'TL-tsne',
                       beta,TL_tsne.config = TL_tsne.default){
  if(mode == 'TL-tsne'){
    TL = igraph::V(g)$TL
    coords = get_nodes_position_TL_tsne(g = g,TL = TL,
                                        TL_tsne.config,beta)
    #check if the layout has been computed already
    if(is.null(igraph::get.vertex.attribute(g,paste0("layout_beta",beta)))){
      g = igraph::set_vertex_attr(graph = g, name = paste0("layout_beta",beta),
                          value = as.vector(coords[,2]))
    } else{
      #add the different layout (for beta value)
      attr_names = igraph::vertex_attr_names(g)
      nrep = length(grep(paste0("layout_beta",beta),attr_names))
      g = igraph::set_vertex_attr(graph = g, name = paste0("layout_beta",beta,"_",nrep),
                          value = as.vector(coords[,2]))
    }

  }
  return(g)
}
  

#' compute and attach 'TL-tsne' layout
#'
#' Method to compute 'TL-tsne' layout and save it as node attributes of the focal network. Each node of the focal network has an attribute \code{layout_beta_VALUE}.
#' If this function is run several times for a given beta value, repetitions of the layout algorithm will be stored as node attributes
#'
#' @param metanetwork object of class 'metanetwork'
#' @param g character indicating the name of the network for which the 'TL-tsne' layout is computed,
#'  default is 'metaweb'
#' @param beta the diffusion parameter of the diffusion kernel, a positive scalar controlling the 
#' squeezing of the network
#' @param mode only 'TL-tsne' is supported for this function
#' @param TL_tsne.config configuration list for mode 'TL-tsne', default is TL_tsne.default
#' @return 'metanetwork' object with layout added as node attribute of the considered network
#'
#' @examples
#' library(metanetwork)
#' library(igraph)
#' # on angola dataset (metaweb)
#' meta_angola = attach_layout(meta_angola,beta = 0.05)
#' V(meta_angola$metaweb)$layout_beta0.05
#' # on a local network
#' meta_angola = attach_layout(meta_angola,g = meta_angola$X1986,beta = 0.05)
#' 
#' # getting repetitions
#' meta_angola = attach_layout(meta_angola,beta = 0.05)
#' V(meta_angola$metaweb)$layout_beta0.05
#' V(meta_angola$metaweb)$layout_beta0.05_1
#'
#' @export
attach_layout <- function(metanetwork,g = NULL,beta = 0.1,
                          mode = 'TL-tsne',TL_tsne.config = TL_tsne.default){
UseMethod("attach_layout",metanetwork)
}


#' @return \code{NULL}
#'
#' @rdname attach_layout
#' @exportS3Method attach_layout metanetwork
attach_layout.metanetwork <- function(metanetwork,g = NULL,beta = 0.1,
				mode = 'TL-tsne',TL_tsne.config = TL_tsne.default){
  if(is.null(g)){
    g = metanetwork$metaweb
  }
  #simplify the network (remove self loops)
  if(!(igraph::is.simple(g))){
    g = igraph::simplify(g)
  }
  g_lay = attach_layout_g(g,metanetwork,mode,beta,TL_tsne.config)
  #identification of g
  if(is.null(g$res)){
    metanetwork[[which(names(metanetwork) == g$name)]] = g_lay
  } else{
    #check if the network is at original resolution
    if(g$res == colnames(metanetwork$trophicTable)[1]){
      metanetwork[[which(names(metanetwork) == g$name)]] = g_lay
    } else{ #if not, combination of name and res is the identifiant
      metanetwork[[which(names(metanetwork) == paste0(g$name,"_",g$res))]] = g_lay
    }
  }
  return(metanetwork)
}

    
    
