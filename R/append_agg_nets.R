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

#' append aggregated networks
#'
#' Method to append aggregated metawebs and local networks using 
#' the hierarchy described in `trophicTable`
#' 
#' It uses the network aggregation method developed in Ohlmann et al. 2019.
#' It computes group abundances and edge probabilities of the aggregated networks.
#'
#'
#' @param metanetwork object of class 'metanetwork'
#' @return an object of class 'metanetwork', with aggregated networks appended to the network list.
#' 
#' @seealso [plot_trophicTable()]
#' 
#' @references Ohlmann, M., Miele, V., Dray, S., Chalmandrier, L., O connor, L., & Thuiller, W. 2019.
#'  Diversity indices for ecological networks: a unifying framework using Hill numbers. Ecology letters, 22 4 , 737-747.
#'
#' @examples
#' library(metanetwork)
#' data(meta_angola)
#' meta_angola = append_agg_nets(meta_angola)
#' names(meta_angola)
#' @export
append_agg_nets <- function(metanetwork){
  UseMethod("append_agg_nets",metanetwork)
}

#' @return \code{NULL}
#'
#' @rdname append_agg_nets
#' @exportS3Method append_agg_nets metanetwork
append_agg_nets.metanetwork <- function(metanetwork){
  if(is.null(metanetwork$trophicTable)||(ncol(metanetwork$trophicTable)==1)){
    stop('single resolution, no aggregation possible')
    return(metanetwork)
  } else {
    N_resolution = ncol(metanetwork$trophicTable)
    #check if aggregated networks are already computed
    networks = extract_networks(metanetwork)
    if(unique(sapply(networks,function(g) g$res)) %>% length() == N_resolution){
      return(metanetwork)
    } else{
      for(scale in colnames(metanetwork$trophicTable)[2:N_resolution]){
        g = metanetwork$metaweb
        groups_loc =  metanetwork$trophicTable[,scale]
        names(groups_loc) = metanetwork$trophicTable[,1]
        array_ag = sbmParams(g = g,groups = groups_loc[igraph::V(g)$name])
        g_agg = igraph::graph_from_adjacency_matrix(t(array_ag$pi),weighted = T)
        g_agg = igraph::set_vertex_attr(g_agg, name = "ab", value = array_ag$alpha)
        g_agg = igraph::set_graph_attr(g_agg, name = "res", value = scale)
        g_agg = igraph::set_graph_attr(g_agg, name = "name", value = "metaweb")
        eval(parse(text = paste0("metanetwork$metaweb_",scale,"= g_agg")))
        
        local_networks_names = rownames(metanetwork$abTable)
        for(local_networks_name in local_networks_names){
          eval(parse(text = paste0('g = metanetwork$',local_networks_name) ))
          groups_loc =  metanetwork$trophicTable[,scale]
          names(groups_loc) = metanetwork$trophicTable[,1]
          array_ag = sbmParams(g,groups_loc[igraph::V(g)$name])
          g_agg = igraph::graph_from_adjacency_matrix(t(array_ag$pi),weighted = T)
          g_agg = igraph::set_graph_attr(g_agg, name = "res", value = scale)
          g_agg = igraph::set_graph_attr(g_agg, name = "name",
                                         value = local_networks_name)
          g_agg = igraph::set_vertex_attr(g_agg, name = "ab", value = array_ag$alpha)
          eval(parse(text = paste0("metanetwork$",local_networks_name,"_",scale,"= g_agg")))
        }
      }
      return(metanetwork) 
    }
  }
}

# get_aggregated_network <- function(g,scale,metanetwork){
#   groups_loc <- metanetwork$trophicTable[,scale]
#   names(groups_loc) = rownames(metanetwork$trophicTable)
#   array_ag = sbmParams(g,groups_loc[igraph::V(g)$name])
#   g_agg = igraph::graph_from_adjacency_matrix(t(array_ag$pi),weighted = T)
#   g_agg = igraph::set_vertex_attr(g_agg, name = "ab", value = array_ag$alpha)
#   return(g_agg)
# }
# 
# get_aggregated_networks <- function(gList,scale,metanetwork){
#   get_aggregated_network_loc <- function(g) {return(get_aggregated_network(g,scale,metanetwork))}
#   return(lapply(gList,FUN = get_aggregated_network_loc))
# }
