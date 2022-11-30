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

#' Default configuration for group-TL-tsne layout
#'
#' A list with parameters customizing group-TL-tsne layout
#'
#' @examples
#' # display all default settings
#' group_layout.default
#'
#' # create a new settings object with n_neighbors set to 5
#' group_layout.custom = group_layout.default
#' group_layout.custom$group_height = 10
#' group_layout.custom
#' 
#' @importFrom stats rnorm sd quantile dist
#' @export
group_layout.default = list(
  nbreaks_group = 1,
  group_height = 5,
  group_width = 5
)

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
get_coord_TL_tsne <- function(g,TL,TL_tsne.config,beta){
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


#get group-TL-tsne layout
get_coord_group_TL_tsne <- function(g,metanetwork,res,beta,group_layout.config){
  #need for a res
  if(is.null(res)){
    stop("no resolution provided, group layout is impossible !")
  }
  #need for a trophic Table
  if(is.null(metanetwork$trophicTable)){
    stop("no trophicTable provided, group layout is impossible !")
  }
  networks = metanetwork[lapply(metanetwork,class) == "igraph"]
  if(!(res %in% sapply(networks,function(g) g$res))){
    stop("to use group-TL-tsne layout, you need to compute aggregated networks first, see append_agg_nets")
  }
  #group-TL-tsne layout is available only at the original resolution
  if(!(g$res == colnames(metanetwork$trophicTable)[1])){
    stop("group layout is only available at the orginal resolution")
  }
  
  #get the network at the desired res
  name_res_tab = sapply(networks,
                        function(g) list(res = g$res, name = g$name))
  g_agg = networks[[which(colSums(name_res_tab == c(res,g$name)) == 2)]]
  if(is.null(igraph::V(g_agg)$TL)){
    stop("to use group-TL-tsne layout, you need to compute trophic levels first, see compute_TL")
  }
  #check if TL-tsne layout is already computed for g_agg
  if(is.null(igraph::get.vertex.attribute(g_agg,paste0("layout_beta",beta)))){
    message(paste0("computing TL-tsne layout for ",g$name,"at resolution: ",res,". See attach_layout to store it."))
    g_agg = attach_layout_g(g_agg,metanetwork,mode = 'TL-tsne',beta,TL_tsne.config = TL_tsne.default)
  }
  nbreaks_group = group_layout.config$nbreaks_group
  group_height = group_layout.config$group_height
  group_width = group_layout.config$group_width
  
  groups = igraph::V(g_agg)$name
  if(nbreaks_group == 1){
    coords_list = c()
    for(k in 1:length(groups)){
      sp_loc = metanetwork$trophicTable[
        which(metanetwork$trophicTable[,res] == groups[k]),
        1]
      g_loc = igraph::induced_subgraph(metanetwork$metaweb,sp_loc)
      #setting coords to 0 if single node
      if(igraph::vcount(g_loc) == 1){
        layout_loc = matrix(0,ncol = 2,nrow = 1)
      }else{
        layout_loc = igraph::layout_with_graphopt(g_loc)
        #centering and scaling
        layout_loc[,1] = group_width*(layout_loc[,1] - mean(layout_loc[,1]))/(sd(layout_loc[,1]))
        layout_loc[,2] = group_height*(layout_loc[,2] - mean(layout_loc[,2]))/(sd(layout_loc[,2]))
      }
      rownames(layout_loc) = igraph::V(g_loc)$name
      #adding the coordinate of the current group in the aggregated layout
      layout_loc[,1] = layout_loc[,1] + 100*igraph::V(g_agg)$TL[k]/max(abs(igraph::V(g_agg)$TL))
      layout_loc[,2] = layout_loc[,2] + 100*igraph::get.vertex.attribute(
        g_agg,paste0("layout_beta",beta))[k]/
        max(abs(igraph::get.vertex.attribute(
          g_agg,paste0("layout_beta",beta))))
      coords_list = c(coords_list,list(layout_loc))
    }
    names(coords_list) = groups
    #rbinding coordinates
    coords_mat = do.call(rbind,coords_list)
    coords_mat = coords_mat[igraph::V(g)$name,]
  } else{
    #check if height and length are of appropriate length.
    if(!((length(group_height) == nbreaks_group) && (length(group_width) == nbreaks_group))){
      stop("group_height and group_length must be of length equal to nbreaks_group, see group_layout.config argument")
    } else{
      #get group sizes
      group_sizes = table(metanetwork$trophicTable[,res])[igraph::V(g_agg)$name]
      group_sizes_cut = cut(group_sizes, 
                             quantile(group_sizes,probs = seq(0,1,length.out = nbreaks_group+1) ) , 
                             include.lowest=TRUE) %>% as.numeric()
      coords_list = c()
      for(k in 1:length(groups)){
        ind_loc = group_sizes_cut[k]
        sp_loc = metanetwork$trophicTable[
            which(metanetwork$trophicTable[,res] == groups[k]),
            1]
        g_loc = igraph::induced_subgraph(metanetwork$metaweb,sp_loc)
        #setting coords to 0 if single node
        if(igraph::vcount(g_loc) == 1){
            layout_loc = matrix(0,ncol = 2,nrow = 1)
        }else{
            layout_loc = igraph::layout_with_graphopt(g_loc)
            #centering and scaling
            layout_loc[,1] = group_width[ind_loc]*(layout_loc[,1] - mean(layout_loc[,1]))/(sd(layout_loc[,1]))
            layout_loc[,2] = group_height[ind_loc]*(layout_loc[,2] - mean(layout_loc[,2]))/(sd(layout_loc[,2]))
        }
        rownames(layout_loc) = igraph::V(g_loc)$name
        #adding the coordinate of the current group in the aggregated layout
        layout_loc[,1] = layout_loc[,1] + 100*igraph::V(g_agg)$TL[k]/max(abs(igraph::V(g_agg)$TL))
        layout_loc[,2] = layout_loc[,2] + 100*igraph::get.vertex.attribute(
          g_agg,paste0("layout_beta",beta))[k]/
          max(abs(igraph::get.vertex.attribute(
          g_agg,paste0("layout_beta",beta))))
        coords_list = c(coords_list,list(layout_loc))
        }
        names(coords_list) = groups
        #rbinding coordinates
        coords_mat = do.call(rbind,coords_list)
        coords_mat = coords_mat[igraph::V(g)$name,]
      }
  }
  return(coords_mat)
}

# # get TL-kpco layout
# get_nodes_position_TL_kpco <- function(g,TL,beta){
#   if(igraph::is.connected(g,mode = "weak")){
#     if(igraph::vcount(g)>2){
#       K = compute_diffusion_kernel(g,beta)
#       nf_loc = dim(K)[1] - 1
#       pco1 = ade4::dudi.pco(dist(K),scannf = FALSE, nf = nf_loc)
#       pco2 = ade4::pcaivortho(pco1,TL, scan = F,nf = nf_loc)
#       coords = cbind(TL,rowMeans(pco2$li))
#       return(coords)
#     } else{
#       coords = cbind(TL,rep(0,igraph::vcount(g)))
#       return(coords)
#     }
#   } else{
#     coords_list = list() #stock coords for each connected component
#     membership_loc = igraph::components(g)$membership
#     for(comp in unique(membership_loc)){
#       #nodes belonging to comp
#       nodes_comp = igraph::V(g)[names(which(membership_loc == comp))]
#       g_comp_loc = igraph::induced_subgraph(g,nodes_comp)
#       coords_comp = get_nodes_position_TL_kpco(g = g_comp_loc,
#                                                TL = igraph::V(g_comp_loc)$TL ,beta = beta)
#       rownames(coords_comp) = igraph::V(g_comp_loc)$name
#       #avoiding overlay of different components
#       if(comp>1){
#         coords_comp[,2] = coords_comp[,2] - min(coords_comp[,2]) + 
#           max(coords_list[[comp-1]][,2]) #+ sd(coords_list[[comp-1]][,2])
#       }
#       coords_list = c(coords_list,list(coords_comp))
#     }
#     coords = do.call(rbind,coords_list)
#     return(coords[igraph::V(g)$name,])
#   }
# }

#save TL-tsne layout as node attribute
attach_layout_g <- function(g,metanetwork,mode = 'TL-tsne',
                       beta,TL_tsne.config = TL_tsne.default,res = NULL,
                       group_layout.config){
  if(mode == 'TL-tsne'){
    TL = igraph::V(g)$TL
    coords = get_coord_TL_tsne(g = g,TL = TL,
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

  } else if (mode == "group-TL-tsne"){
    coords = get_coord_group_TL_tsne(g = g,metanetwork = metanetwork,
                                     group_layout.config = group_layout.config,res = res,beta = beta)
    g = igraph::set_vertex_attr(graph = g, name = paste0("group_layout_x_beta",beta),
                                value = as.vector(coords[,1]))
    g = igraph::set_vertex_attr(graph = g, name = paste0("group_layout_y_beta",beta),
                                value = as.vector(coords[,2]))
  }
  return(g)
}
  

#' compute and attach metanetwork layouts
#'
#' Method to compute `'TL-tsne'` and `'group-TL-tsne'` layouts and save it as node attributes of the focal network.
#' 
#' The `'TL-tsne'` layout is a diffusion based layout algorithm specifically designed for trophic networks.
#' In metanetwork, first axis is the trophic level (see `compute_TL` method) whereas the second axis is computed using a diffusion graph kernel (Kondor & Lafferty 2002) 
#' and tsne dimension reduction algorithm to (see van der Maaten & Hinton (2008) and 'tsne' R package). 
#' Let \eqn{A} be the adjacency matrix of the considered network and \eqn{D} its degree diagonal matrix. 
#' The Laplacian matrix of the symmetrised network is defined by: 
#' 
#' \deqn{L = D - A - t(A)}
#' 
#' The diffusion graph kernel is:
#' 
#' \deqn{K = exp(-beta*L)}
#' 
#' It is a similarity matrix between nodes according to a diffusion process. `beta` is the diffusion constant,it must be provided by the user.
#' `beta` parameter influences the layout by grouping together similar paths (see `pyramid` vignette).
#' Each node of the focal network has an attribute \code{layout_beta_VALUE}.
#' If this function is run several times for a given beta value, repetitions of the layout algorithm will be stored as node attributes.
#' 
#' The `'group-TL-tsne'` layout is a variation of `'TL-tsne` layout. For a focal network, it mixes `'TL-tsne'` layout at the desired aggregated level
#' with the `layout_with_graphopt` function from `igraph`. It clusters nodes belonging to the same group. 
#' `'group-TL-tsne'` layout is recommended for large networks since you only need to compute `'TL-tsne'` at the aggregated network
#'  that is much smaller than the focal network. `group_layout.config` allows controlling the overall size of the groups.
#'
#' @param metanetwork object of class 'metanetwork'
#' @param g character indicating the name of the network for which the 'TL-tsne' layout is computed,
#'  default is 'metaweb'
#' @param beta the diffusion parameter of the diffusion kernel, a positive scalar controlling the 
#' squeezing of the network, default is 0.1
#' @param mode 'TL-tsne' or 'group-TL-tsne', default is 'TL-tsne'.
#' @param res resolution for the 'group-TL-tsne' layout
#' @param TL_tsne.config configuration list for mode 'TL-tsne', default is TL_tsne.default
#' @param group_layout.config configuration list for mode 'group-TL-tsne', default is group_layout.default
#' @return an object of class 'metanetwork', with the computed layout added as node attribute of the considered network
#' 
#' @seealso [ggmetanet()], [vismetaNetwork()],[group_layout.default]
#'
#' @references Kondor, R. I., & Lafferty, J. (2002, July). Diffusion kernels on graphs and other discrete structures.
#'  In Proceedings of the 19th international conference on machine learning (Vol. 2002, pp. 315-322). 
#'  Van der Maaten, L., & Hinton, G. (2008). Visualizing data using t-SNE. Journal of machine learning research, 9(11).
#'
#' @examples
#' library(metanetwork)
#' library(igraph)
#' # on angola dataset (metaweb)
#' data("meta_angola")
#' meta_angola = attach_layout(meta_angola,beta = 0.05)
#' V(meta_angola$metaweb)$layout_beta0.05
#' @export
attach_layout <- function(metanetwork,g = NULL,beta = 0.1,
                          mode = 'TL-tsne',TL_tsne.config = TL_tsne.default,
                          res = NULL,group_layout.config = group_layout.default){
UseMethod("attach_layout",metanetwork)
}

#' @return \code{NULL}
#'
#' @rdname attach_layout
#' @exportS3Method attach_layout metanetwork
attach_layout.metanetwork <- function(metanetwork,g = NULL,beta = 0.1,
				                              mode = 'TL-tsne',TL_tsne.config = TL_tsne.default,
				                              res = NULL,group_layout.config = group_layout.default){
  if(is.null(g)){
    g = metanetwork$metaweb
  }
  if(is.null(igraph::V(g)$TL)){
    stop("you must compute trophic levels first, see compute_TL")  
  }
  if(is.null(res)){
    message(paste0("attaching ",mode," layout for ",g$name,"_",g$res,"\n
                  beta = ",beta))
  }else{
    if(!(res %in% colnames(metanetwork$trophicTable))){
      stop(paste0("res must be one the the available resolutions:",
                  paste(colnames(metanetwork$trophicTable),collapse = ' ')))
    }else{
      message(paste0("attaching ",mode," layout for ",g$name, " at resolution: ",res,"\n
                  beta = ",beta))
    }
  }

  #simplify the network (remove self loops)
  if(!(igraph::is.simple(g))){
    g = igraph::simplify(g)
  }
  g_lay = attach_layout_g(g = g,metanetwork = metanetwork,mode = mode,beta = beta,
                          TL_tsne.config = TL_tsne.config,res = res,
                          group_layout.config = group_layout.config)
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

    
    
