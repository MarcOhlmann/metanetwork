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

#' append mean networks
#'
#' Append mean networks per (discrete) value of covariable
#'
#' @param metanetwork object of class 'metanetwork'
#' @return object of class 'metanetwork' with mean networks appended per value of the covariable at all resolutions
#'
#' @examples
#' library(igraph)
#' library(metanetwork)
#' 
#'
#' @export
append_mean_nets <- function(metanetwork){
  UseMethod("append_mean_nets",metanetwork)
}

#' @return \code{NULL}
#'
#' @rdname append_mean_nets
#' @exportS3Method append_mean_nets metanetwork
append_mean_nets.metanetwork <- function(metanetwork){
if(is.null(metanetwork$covariable)){
  stop('covariable is null!')
} else{
  covariable = metanetwork$covariable
  covariable.df = as.data.frame(covariable)
  networks = metanetwork[sapply(metanetwork,class) == 'igraph']
  #remove metawebs
  networks = networks[sapply(networks,function(g) g$name) != "metaweb"]
  cov_values = unique(covariable)
  res_vec  = unique(sapply(networks, function(g) g$res))
  origin_res = colnames(metanetwork$trophicTable)[1]
  for(cov_value in cov_values){
    for(res in res_vec){
      if(res == origin_res){
        #indices of networks at origin_res that have cov_value as covariable values
        ind_cov_res = which(sapply(networks,function(g) g$res) %in% origin_res)
        networks_loc = networks[ind_cov_res]
        ind_cov_value = which(sapply(networks_loc,function(g) g$name) %in% 
          rownames(covariable.df)[which(covariable.df$covariable == cov_value)])
        gList = networks_loc[ind_cov_value]
        g_union = get_mean_network(gList,metanetwork,cov_value)
        g_union$res = res
        g_union$name = paste0('mean_net_',cov_value)
        eval(parse(text = paste0("metanetwork$",'mean_net_',cov_value,"= g_union")))
      } else{
        #indices of corresponding networks at current res
        ind_cov_value_res = intersect(
          which(sapply(networks,function(g) g$name) %in% sapply(gList,function(g) g$name)),
          which(sapply(networks,function(g) g$res) %in% res)
        )
        gList = networks[ind_cov_value_res]
        g_union = get_mean_network(gList,metanetwork,cov_value)
        g_union$res = res
        g_union$name = paste0('mean_net_',cov_value)
        eval(parse(text = paste0("metanetwork$",'mean_net_',cov_value,'_',res,"= g_union")))
      }
    }
  }
  return(metanetwork)
}
} 

get_mean_network <- function(gList,metanetwork,cov_value){
   #check if all networks have the same resolution
  res_loc = unique(sapply(gList,function(g) g$res))
  if(length(res_loc)>1){
  stop("networks must have the same resolutions")  
  }
  n = length(gList)
  #easy case: finest resolutions, edges weights is 0 or 1
  if(res_loc == colnames(metanetwork$trophicTable)[1]){
    # g_union = do.call(what = union,args =  gList) #NOT WORKING IN THE PACKAGE TO FIX
    g_union = gList[[1]]
    for(k in 2:length(gList)){
      g_union = igraph::union(g_union,gList[[k]])
    }
    #set nodes in alphabetic order (to avoid confusion)
    g_union = igraph::permute(g_union,order(order(igraph::V(g_union)$name)))
    #get the abundances (NA when group is absent), ab_1, ab_2,... in the union graph
    g_union_abs = igraph::get.vertex.attribute(g_union)[
                   grep('^ab_*',names(igraph::get.vertex.attribute(g_union)))
                  ]
    #set 0 instead of NA
    g_union_abs = lapply(g_union_abs,function(x) ifelse(is.na(x),0,x))
    ab_mean_loc = Reduce('+',g_union_abs)/n
    #deleting useless vertiex and edge attributes
    igraph::graph.attributes(g_union) = igraph::graph.attributes(g_union)[1]
    names(igraph::graph_attr(g_union)) = "res"
    igraph::vertex_attr(g_union) = igraph::vertex_attr(g_union)["name"]
    #setting new vertices and edges attr
    g_union = igraph::set_graph_attr(g_union,name = "name",value = cov_value)
    g_union = igraph::set_vertex_attr(g_union, name = "ab", value = ab_mean_loc)
    igraph::edge_attr(g_union) = igraph::edge_attr(g_union)[1]
    names(igraph::edge_attr(g_union)) = "weight"
    igraph::E(g_union)$weight = 1
    return(g_union)
  } else{
    #at another resolution: get the mean network at the finest resolution and aggregate it
    networks = metanetwork[sapply(metanetwork,class) %in% 'igraph']
    #remove metawebs
    networks = networks[sapply(networks,function(g) g$name) != "metaweb"]
    #index of elements of metanetwork at the original resolution corresponding to networks at current resolution
    ind_res_origin = intersect(
      which(sapply(networks,function(g) g$name) %in% sapply(gList,function(g) g$name)),
      which(sapply(networks,function(g) g$res) %in% colnames(metanetwork$trophicTable)[1])
                               )
    #corresponding networks at original resolution
    gList_origin_res = networks[ind_res_origin]
    #recursive call (g_union_origin_res is at the original resolution)
    g_union_origin_res = get_mean_network(gList_origin_res,metanetwork,cov_value)
    #aggregate g_union
    g_union = get_aggregated_network(g = g_union_origin_res,scale = res_loc,metanetwork)
    return(g_union)
  }
}

