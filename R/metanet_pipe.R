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

#' Execute 'metanetwork' pipeline
#' 
#' Method executing the whole metanetwork pipeline for the initial metanetwork object (\code{append_agg_nets}, \code{compute_TL}, 
#' \code{attach_layout}) 
#'
#' @param metanetwork object of class 'metanetwork'
#' @param beta the diffusion parameter of the diffusion kernel, a positive scalar controlling the 
#' squeezing of the network
#' @param verbose a boolean indicating whether message along the pipeline should be printed
#'
#' @return object of class 'metanetwork', with computed trophic levels and layout stored as node attribute
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
#' @export
metanet_pipe <- function(metanetwork,beta = 0.1,verbose = TRUE){
  UseMethod("metanet_pipe",metanetwork)
}

#' @return \code{NULL}
#'
#' @rdname metanet_pipe
#' @exportS3Method metanet_pipe metanetwork
metanet_pipe.metanetwork <- function(metanetwork,beta = 0.1,verbose = TRUE){
  if(verbose){
    if(!(is.null(metanetwork$trophicTable))){
      message("append aggregated networks")    
    }
    message("computing trophic levels")
    message(paste0("attaching layout for beta= ",beta))
  }
  #piping the different metanetwork operations
  if(is.null(metanetwork$trophicTable)){
    meta_loc = metanetwork %>% compute_TL() %>%
      attach_layout(beta = beta)
  }else{
    meta_loc = metanetwork %>% append_agg_nets() %>%  compute_TL() %>%
      attach_layout(beta = beta)
  }
  return(meta_loc)
}



