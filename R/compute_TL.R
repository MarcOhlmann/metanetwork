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

compute_TL_laplacian <- function(G,metanetwork){
  #remove self loops
  if(!(igraph::is.simple(G))){
    G = igraph::simplify(G)
  }
  #recursive function on connected components of G
  if(is.null(E(G)$weight)){
    A = as.matrix(igraph::get.adjacency(G))
    names_loc = rownames(A)
    u  = igraph::degree(G)
    v =  igraph::degree(G,mode ='in') - igraph::degree(G,mode = 'out')
  } else{
    A = as.matrix(igraph::get.adjacency(G,attr = "weight"))
    names_loc = rownames(A)
    u  = igraph::strength(G)
    v =  igraph::strength(G,mode ='in') - igraph::strength(G,mode = 'out')
  }
  A = A[-1,-1]
  u = u[-1]
  v = v[-1]
  L = diag(u) - A - t(A) #matrix is made symmetric!

  
  if(igraph::is.connected(G,mode = 'weak')){
    if(igraph::vcount(G) == 1){
      #assigning metaweb trophic value if single node    
      TL_vec = igraph::V(metanetwork$metaweb)[igraph::V(G)[1]$name]$TL
    } else{
      TL_vec = Matrix::solve(L,v)
      TL_vec = c(0,TL_vec)
      TL_vec = TL_vec - min(TL_vec)
      names(TL_vec) = igraph::V(G)$name
      if(G$name != "metaweb"){
        #metaweb at the current resolution
        networks = metanetwork[sapply(metanetwork,class) == 'igraph']
        metaweb_loc = networks[sapply(networks,function(g) g$name == 'metaweb')]
        #if several resolutions, get the current one
        if(length(metaweb_loc)>1){
          metaweb_loc = networks[[which(sapply(networks,function(g)
            g$res == G$res && g$name == 'metaweb'))]]
        } else {metaweb_loc = metaweb_loc$metaweb}
        #assigning trophic level of metaweb for the
        #reference species in the local network
        TL_vec = TL_vec + mean(igraph::V(metaweb_loc)[names(TL_vec[TL_vec == 0])]$TL)
      }
      names(TL_vec) = names_loc
    }
    return(TL_vec)
  } else{
    TL_list = list() #stock the results for all connected components
    membership_loc = igraph::components(G)$membership
    for(comp in unique(membership_loc)){
      #nodes belonging to comp
      nodes_comp = igraph::V(G)[names(which(membership_loc == comp))]
      G_comp_loc = igraph::induced_subgraph(G,nodes_comp)
      TL_comp = compute_TL_laplacian(G = G_comp_loc,metanetwork = metanetwork)
      names(TL_comp) = igraph::V(G_comp_loc)$name
      TL_list = c(TL_list,list(TL_comp))
    }
    TL_vec = unlist(TL_list)
    TL_vec = TL_vec[igraph::V(G)$name]
    return(TL_vec)
  }
}

#' compute trophic levels
#'
#' Compute trophic levels using graph laplacian
#'
#' @param metanetwork object of class 'metanetwork'
#' @return object of class 'metanetwork' with a node attribude \code{TL}
#'
#' @examples
#' library(metanetwork)
#' library(igraph)
#' 
#' #on angola dataset
#' meta_angola  = compute_TL(meta_angola)
#' V(meta_angola$metaweb)$TL
#'
#' @export
compute_TL <- function(metanetwork){
  UseMethod("compute_TL", metanetwork)
}

#' @return \code{NULL}
#'
#' @rdname compute_TL
#' @exportS3Method compute_TL metanetwork
compute_TL.metanetwork <- function(metanetwork){ 
    TL_loc = compute_TL_laplacian(G = metanetwork$metaweb,
                                  metanetwork = metanetwork)
    metanetwork$metaweb = igraph::set_vertex_attr(
    metanetwork$metaweb, name = "TL", value = TL_loc)
    #TL for metaweb at different aggregation levels
    metaweb_names = names(metanetwork)[grep('metaweb_',x = names(metanetwork))]
    for(metaweb_name in metaweb_names){
      eval(parse(text=paste0('TL_loc = compute_TL_laplacian(metanetwork$',metaweb_name,')')))
      eval(parse(text=paste0('metanetwork$',metaweb_name, '= igraph::set_vertex_attr(metanetwork$',metaweb_name,', name = "TL", value = TL_loc)' )))
    }
    #TL for local networks
    local_networks = metanetwork[sapply(metanetwork,class) == 'igraph']
    local_networks = local_networks[-grep("metaweb",sapply(local_networks,function(g) g$name))]
    local_networks_ids = which(sapply(metanetwork,class) == 'igraph')
    local_networks_ids = local_networks_ids[-grep("metaweb",names(local_networks_ids))]
    
    TL_list = lapply(local_networks,function(g)
      compute_TL_laplacian(G = g,metanetwork = metanetwork))
    # for(g in local_networks){
    #   compute_TL_laplacian(G = g,metanetwork = metanetwork)
    # }
    loc_net_TL_list = list(locnet = local_networks,TL = TL_list)
    
    #setting trophic level as node attribute for local networks
    if(length(local_networks)){
    for(k in 1:length(local_networks)){
      g = local_networks[[k]]
      local_networks[[k]] = igraph::set_vertex_attr(g, name = "TL",
                                                value = TL_list[[k]])
    }
    metanetwork[local_networks_ids] = local_networks
    }
  return(metanetwork)
} 

#compute trophic levels of the difference network
compute_TL_diff <- function(metanetwork_diff,metanetwork){
  G = metanetwork_diff$metaweb
  #recursive function on connected components of G
  A = as.matrix(igraph::get.adjacency(G))
  names_loc = rownames(A)
  u  = igraph::degree(G)
  v =  igraph::degree(G,mode ='in') - igraph::degree(G,mode = 'out')
  
  A = A[-1,-1]
  u = u[-1]
  v = v[-1]
  L = diag(u) - A - t(A) #matrix is made symmetric!
  
  if(igraph::is.connected(G,mode = 'weak')){
    if(igraph::vcount(G) == 1){
      #assigning metaweb trophic value if single node
      TL_vec = igraph::V(metanetwork$metaweb)[igraph::V(G)[1]$name]$TL
    } else{
      TL_vec = Matrix::solve(L,v)
      TL_vec = c(0,TL_vec)
      TL_vec = TL_vec - min(TL_vec)
      names(TL_vec) = igraph::V(G)$name
      if(G$name != "metaweb"){
        #metaweb at the current resolution
        networks = metanetwork[sapply(metanetwork,class) == 'igraph']
        metaweb_loc = networks[sapply(networks,function(g) g$name == 'metaweb')]
        #if several resolutions, get the current one
        if(length(metaweb_loc)>1){
          metaweb_loc = networks[[which(sapply(networks,function(g)
            g$res == G$res && g$name == 'metaweb'))]]
        } else {metaweb_loc = metaweb_loc$metaweb}
        #assigning trophic level of metaweb for the
        #reference species in the local network
        TL_vec = TL_vec + mean(igraph::V(metaweb_loc)[names(TL_vec[TL_vec == 0])]$TL)
      }
      names(TL_vec) = names_loc
    }
    return(TL_vec)
  } else{
    TL_list = list() #stock the results for all connected components
    membership_loc = igraph::components(G)$membership
    for(comp in unique(membership_loc)){
      #nodes belonging to comp
      nodes_comp = igraph::V(G)[names(which(membership_loc == comp))]
      G_comp_loc = igraph::induced_subgraph(G,nodes_comp)
      TL_comp = compute_TL_laplacian(G = G_comp_loc,metanetwork = metanetwork)
      names(TL_comp) = igraph::V(G_comp_loc)$name
      TL_list = c(TL_list,list(TL_comp))
    }
    TL_vec = unlist(TL_list)
    TL_vec = TL_vec[igraph::V(G)$name]
    return(TL_vec)
  }
}
