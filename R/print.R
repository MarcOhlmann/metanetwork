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

#' print metanetwork
#'
#' Print method for class \code{metanetwork}
#'
#' @param metanetwork object of class 'metanetwork'
#' @return character indicating number of nodes and edges of the metaweb, available resolutions 
#' and number of local networks
#'
#' @examples
#' library(metanetwork)
#' library(igraph)
#' 
#' g = make_ring(5,directed = TRUE)
#' meta = build_metanet(g)
#' print(meta)
#' 
#' #on Angola dataset
#' data("meta_angola")
#' print(meta_angola)
#'
#' #on Norway dataset
#' data("meta_norway")
#' print(meta_norway)
#' @export
print <- function(metanetwork){
  UseMethod("print",metanetwork)
}

#' @return \code{NULL}
#'
#' @rdname print
#' @exportS3Method  print metanetwork
print.metanetwork <- function(metanetwork){
  message("object of class metanetwork")
  cat("metaweb has",length(igraph::V(metanetwork$metaweb)),
      "nodes and",length(igraph::E(metanetwork$metaweb)),"edges","\n")
  if(nrow(metanetwork$abTable) == 1){
    cat("single network","\n")
  }else{
      cat(nrow(metanetwork$abTable),"local networks","\n")
  }
  if(is.null(metanetwork$trophicTable)){
    cat("single resolution available","\n")
  }else{
    networks = metanetwork[which(sapply(metanetwork,class) == 'igraph')]
    if(setequal(unique(sapply(networks,function(g) g$res)),colnames(metanetwork$trophicTable))){
      cat("available resolutions are:",colnames(metanetwork$trophicTable),"\n")
    }else{
      cat("available resolutions (not computed) are:",colnames(metanetwork$trophicTable),"\n")
    }
  }
}



