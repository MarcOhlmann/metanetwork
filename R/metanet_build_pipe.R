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

#' Build and execute 'metanetwork' pipeline
#' 
#' Method executing the whole metanetwork pipeline, including building 'metanetwork' object (\code{build_metanet},\code{append_agg_nets}, \code{compute_TL}, 
#' \code{attach_layout}) 
#'
#' @param metaweb metaweb of the metanetwork, object of class 'graph', 'matrix', 'data.frame' or 'dgCMatrix'. 
#' Metaweb needs to be directed and connected. This parameter must be non-null.
#' @param abTable abundances of nodes in local networks, matrix of class 'matrix',
#'  columns must have names corresponding to node labels of the metaweb,
#'  rows are node abundances in local networks.
#'  Default is null, in that case, uniform abundances are assigned
#' @param trophicTable a 'matrix' or 'data.frame' indicating hierarchy of the nodes.
#'  Names of the columns correspond to the different resolutions.
#'  It indicates the membership of each node of the metaweb. Default is null.
#' @param compute_local_nets a boolean, indicates whether local networks must be computed or not.
#' Default is \code{TRUE}
#' @param beta the diffusion parameter of the diffusion kernel, a positive scalar controlling the 
#' squeezing of the network
#' @param verbose a boolean indicating whether message along the pipeline should be printed
#'
#' @return object of class 'metanetwork', with computed layout stored as node attribute
#'
#' @examples
#' library(metanetwork)
#' library(igraph)
#' 
#' g = make_lattice(dimvector = c(4,4),2,3,directed = TRUE)
#' meta0 = metanet_build_pipe(g)
#' ggmetanet(meta0)
#' 
#'
#' @export
metanet_build_pipe <- function(metaweb,abTable = NULL, trophicTable = NULL,
                               compute_local_nets = TRUE,verbose = TRUE,beta = 0.1){
 
  if(verbose){
    message("building metanetwork")
    if(!(is.null(trophicTable))){
      message("appending aggregated networks")    
    }
    message("computing trophic levels")
    message(paste0("attaching layout for beta= ",beta))
  }
  #piping the different metanetwork operations
  meta_loc = build_metanet(metaweb,abTable = abTable, trophicTable = trophicTable,
                           compute_local_nets = compute_local_nets)
  if(is.null(meta_loc$trophicTable)){
    meta_loc = meta_loc %>% compute_TL() %>%
      attach_layout(beta = beta)
  }else{
    meta_loc = meta_loc %>% append_agg_nets() %>%  compute_TL() %>%
      attach_layout(beta = beta)
  }
  return(meta_loc)
}



