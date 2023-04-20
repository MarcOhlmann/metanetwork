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


#' compute network metrics
#'
#' Function to compute (some) network metrics on the metaweb and local networks
#' 
#' XXXXXXXXXXXXXXXXXXXXXx
#' 
#' @references Ohlmann, M., Miele, V., Dray, S., Chalmandrier, L., O'connor, L., & Thuiller, W. (2019). Diversity indices for ecological networks: a unifying framework using Hill numbers. Ecology letters, 22(4), 737-747.
#'
#' @param metanetwork object of class 'metanetwork'
#' @param res a vector containing the resolutions at which the metrics are computed
#' @return a \code{data.frame}
#'
#' @examples
#' library(metanetwork)
#' library(igraph)
#' 
#' #on angola dataset
#' data("meta_angola")
#' compute_metrics(meta_angola)
#' 
#' #computing metrics only at Phylum level
#' compute_metrics(meta_angola,res = "Phylum")
#'
#' @export
compute_diversities <- function(metanetwork,res = NULL){
  # get the local networks
  networks = metanetwork[lapply(metanetwork,class) == "igraph"]
  metaweb_names = names(metanetwork)[grep('metaweb',x = names(metanetwork))]
  networks_loc = networks[!(names(networks) %in% metaweb_names)]
  metaweb_loc = networks[names(networks) %in% metaweb_names]

  if(is.null(res)){
    res_loc = unique(sapply(networks_loc,function(g) g$res))
  }else{
    if(!(res %in% colnames(metanetwork$trophicTable))){
      stop(paste0("res must be a vector with available resolutions: ",  paste(colnames(metanetwork$trophicTable),collapse = " ")))
    }else{
      res_loc = res
    }
  }
  
  if(is.null(igraph::V(metanetwork$metaweb)$TL)){
    stop("to use compute_metrics, you need to compute trophic levels first. See compute_TL")
  }
  
  for(res in res_loc){
    metrics_df = matrix(NA,nrow = 1 + nrow(metanetwork$abTable),ncol = 5)
    colnames(metrics_df) = c("connectance","mean_TL","max_TL","shortest_path_length","modularity")
    rownames(metrics_df) = c("metaweb",rownames(metanetwork$abTable))
    metrics_df = as.data.frame(metrics_df)
    
    metaweb_loc_loc = metaweb_loc[[which(sapply(metaweb_loc, function(g) g$res) == res)]]
    #metrics for the metaweb
    metrics_df$connectance[1] = get_connectance(metaweb_loc_loc)
    metrics_df$mean_TL = metaweb_loc_loc
    
    
    
  }
  return(list(nodes = diversites_df_P,links = diversites_df_L))
}

get_connectance <- function(g){
  adj = igraph::get.adjacency(g,attr = "weight") %>% as.matrix()
  C = t(igraph::V(g)$ab) %*% adj %*% igraph::V(g)$ab
  return(C)
}

