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
#' @param metanetwork object of class 'metanetwork'
#' @param beta the diffusion parameter of the diffusion kernel, a positive scalar controlling the 
#' squeezing of the network
#' @param verbose a boolean indicating wheter message along the pipeline should be printed
#'
#' @return object of class 'metanetwork'
#'
#' @examples
#' library(metanetwork)
#' library(igraph)
#' 
#' g = make_lattice(dimvector = c(4,4),2,3,directed = TRUE)
#' meta0 = build_metanet(g)
#' meta0 = metanet_pipe(meta0)
#' ggmetanet(meta0)
#' 
#' #on angola data set
#' meta_angola = metanet_pipe(meta_angola,beta = 0.05)
#'
#' @export
metanet_build_pipe <- function(metanetwork,beta = 0.1,verbose = T){
  UseMethod("metanet_build_pipe",metanetwork)
}

#' @return \code{NULL}
#'
#' @rdname metanet_build_pipe
#' @exportS3Method metanet_build_pipe metanetwork
metanet_build_pipe.metanetwork <- function(metanetwork,beta = 0.1,verbose = T){
  if(verbose){
    message("building metanetwork")
    if(!(is.null(metanetwork$trophicTable))){
      message("append aggregated networks")    
    }
    message("computing trophic levels")
    message(paste0("attaching layout for beta= ",beta))
  }
  #piping the different metanetwork operations
  meta_loc = build_metanetwork(metaweb,abTable = abTable, trophicTable = trophicTable,
                               compute_local_networks = compute_local_networks,
                               covariable = covariable)
  if(is.null(meta_loc$trophicTable)){
    meta_loc = meta_loc %>% compute_TL() %>%
      attach_layout(beta = beta)
  }else{
    meta_loc = meta_loc %>% append_agg_nets(metanetwork) %>%  compute_TL() %>%
      attach_layout(beta = beta)
  }
  return(meta_loc)
}



