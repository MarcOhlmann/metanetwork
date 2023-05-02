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
#' This function computes network metrics for the metaweb and local networks at the available resolutions. This function computes the following metrics:
#'  - weighted connectance: connectance weighted by abundances and link probability (see Ohlmann et al. 2019)
#'  Denoting \eqn{p_q} abundance of node \eqn{q} and \eqn{\pi_{ql}} the interaction probability between nodes \eqn{q} and \eqn{l}, then the weighted connectance \eqn{C} is:
#'  \deqn{C = \sum_{q,l} \pi_{ql}p_q p_l}
#'  - mean_TL and max_TL: mean and max trophic level of each network, using values computed from compute_TL() method.
#'  - mean_shortest_path_length: mean shortest path lenght of each network computed using igraph mean_distance() function.
#'
#' 
#' @references 
#' - Csardi, G., & Nepusz, T. (2006). The igraph software package for complex network research. InterJournal, complex systems, 1695(5), 1-9.
#'
#' - Ohlmann, M., Miele, V., Dray, S., Chalmandrier, L., O'connor, L., & Thuiller, W. (2019). Diversity indices for ecological networks: a unifying framework using Hill numbers. Ecology letters, 22(4), 737-747.

#'
#' @param metanetwork object of class 'metanetwork'
#' @param res a vector containing the resolutions at which the metrics are computed
#' @return a \code{data.frame}
#'
#' @seealso [compute_TL()],[compute_div()]
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
compute_metrics <- function(metanetwork,res = NULL){
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
  
  #if several resolution
  if(!(is.null(metanetwork$trophicTable))){
    metrics_df_list = list()
    for(res in res_loc){
      if(nrow(metanetwork$abTable) == 1){ #singe network case
        metrics_df = get_metrics(metaweb_loc[[1]])
        names(metrics_df) = c("connectance","mean_TL","max_TL","mean_shortest_path_length")
      }else{
      metrics_df = matrix(NA,nrow = 1 + nrow(metanetwork$abTable),ncol = 4)
      colnames(metrics_df) = c("connectance","mean_TL","max_TL","mean_shortest_path_length")
      rownames(metrics_df) = c("metaweb",rownames(metanetwork$abTable))
      metrics_df = as.data.frame(metrics_df)
      
      metaweb_loc_loc = metaweb_loc[[which(sapply(metaweb_loc, function(g) g$res) == res)]]
      #metrics for the metaweb
      metrics_df[1,] = get_metrics(metaweb_loc_loc)
      #metrics for local networks
      networks_loc_loc = networks_loc[which(sapply(networks_loc, function(g) g$res) == res)]
      metrics_df[2:nrow(metrics_df),] = do.call(rbind, lapply(networks_loc_loc,get_metrics))
      metrics_df_list = c(metrics_df_list,list(metrics_df))
      }
    }
    names(metrics_df_list) = res_loc
    return(metrics_df_list)
  }else{ #single resolution case
    if(nrow(metanetwork$abTable) == 1){ #singe network case
      metrics_vec = get_metrics(metaweb_loc[[1]])
      names(metrics_vec) = c("connectance","mean_TL","max_TL","mean_shortest_path_length")
      return(metrics_vec)
    }else{
      metrics_df = matrix(NA,nrow = 1 + nrow(metanetwork$abTable),ncol = 4)
      colnames(metrics_df) = c("connectance","mean_TL","max_TL","mean_shortest_path_length")
      rownames(metrics_df) = c("metaweb",rownames(metanetwork$abTable))
      metrics_df = as.data.frame(metrics_df)
      #metrics for the metaweb
      metrics_df[1,] = get_metrics(metaweb_loc[[1]])
      #metrics for local networks
      metrics_df[2:nrow(metrics_df),] = do.call(rbind, lapply(networks_loc,get_metrics))
      return(metrics_df)
    }

  }
}

#get weighted connectance
get_connectance <- function(g){
  adj = igraph::get.adjacency(g,attr = "weight") %>% as.matrix()
  C = t(igraph::V(g)$ab) %*% adj %*% igraph::V(g)$ab
  return(C)
}


get_metrics <- function(g){
  C = get_connectance(g)
  mean_TL = mean(igraph::V(g)$TL)
  max_TL = max(igraph::V(g)$TL)
  short_path = igraph::mean_distance(g,directed = TRUE)
  return(c(C,mean_TL,max_TL,short_path))
}


