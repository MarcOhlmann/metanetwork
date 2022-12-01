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

#'Plot trophic groups hierarchy
#'
#'Function to represent trophic groups hierarchy provided by trophicTable
#'
#' @param metanetwork object of class 'metanetwork'
#' @param res resolutions included in the hierarchy representation. Default is "all" (all resolutions are then included)
#'  but can be also a vector of given resolutions
#' @param ggnet.config configuration list for ggnet representation, default is ggnet.default
#'
#' @return object of class 'ggnet', representation of group hierarchy
#'
#' @examples
#' library(metanetwork)
#' 
#' #on Angola data_set
#' data("meta_angola")
#' plot_trophicTable(meta_angola)
#'
#' @export
plot_trophicTable = function(metanetwork,res = "all",ggnet.config = ggnet.default){
  if(length(res) == 1){
    if(res == "all"){
      res = colnames(metanetwork$trophicTable)
    } else{
      stop("if not equal to 'all', resolutions should be at least of length 2")
    }
  } else {
    if(!(prod(res %in% colnames(metanetwork$trophicTable)))){
      stop("resolutions should be a vector (character) 
        indicating the resolutions considered (see available resolutions using print(metanetwork))")
      }
  }
  
  #make res in the same order as trophicTable
  res = res[order(match(res,colnames(metanetwork$trophicTable)))]
  N_res = length(res)
  
  trophicTable_loc = metanetwork$trophicTable 
  trophicTable_loc = trophicTable_loc[,res]
  
  graph_data_frame = data.frame(from = "origin", to = unique(trophicTable_loc[,N_res]))
  for(k in 2:N_res){
      graph_data_frame = rbind(graph_data_frame,
                               data.frame(from = trophicTable_loc[,N_res+1-(k-1)],
                                          to = trophicTable_loc[,N_res+1-k]))
    }
  group_tree = igraph::graph_from_data_frame(graph_data_frame)
  color_loc = c("origin",
                as.vector(unlist(sapply(res[seq(N_res,1,length.out = N_res)],
                                        function(x) rep(x,length(unique(metanetwork$trophicTable[,x]))))))
                )
  lay_loc = igraph::layout_as_tree(group_tree,circular = T)
  group_tree_Network = intergraph::asNetwork(igraph::simplify(group_tree))
  network::network.vertex.names(group_tree_Network) = igraph::V(group_tree)$name
  
  return(GGally::ggnet2(group_tree_Network,mode = lay_loc, color = color_loc,
         label = ggnet.config$label, arrow.size = ggnet.config$arrow.size,
         alpha = ggnet.config$alpha, size = ggnet.config$max_size,
         label.size = ggnet.config$label.size,palette = ggnet.config$palette,
         edge.alpha = ggnet.config$edge.alpha,edge.size = ggnet.config$edge.size,
         arrow.gap = ggnet.config$arrow.gap))
}
