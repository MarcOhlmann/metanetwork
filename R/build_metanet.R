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

#' Build metanetwork object
#'
# Build an object of S3 class \code{metanetwork}
#'
#' @param metaweb metaweb of the metanetwork, object of class 'graph', 'matrix', 'data.frame' or 'dgCMatrix'. 
#' Metaweb needs to be directed and connected. This argument must be non-null.
#' @param abTable node abundances in local networks, matrix of class 'matrix',
#'  columns must have names corresponding to node labels of the metaweb,
#'  rows are node abundances in local networks.
#'  Default is null, in that case, uniform abundances are assigned
#' @param trophicTable a 'matrix' or 'data.frame' indicating hierarchy of the nodes.
#'  Names of the columns correspond to the different resolutions.
#'  It indicates the membership of each node of the metaweb. Default is null.
#' @param compute_local_nets a boolean, indicates whether local networks must be computed or not.
#' Default is \code{TRUE}
#'
#' @return an object of S3 class 'metanetwork'
#'
#' @examples
#' library(metanetwork)
#' library(igraph)
#' #with a single metaweb
#' g = igraph::make_ring(5,directed = TRUE)
#' meta = build_metanet(g)
#' 
#' #on Angola dataset (re-building the dataset)
#' data("meta_angola")
#' metaweb = meta_angola$metaweb
#' abTable = meta_angola$abTable
#' trophicTable = meta_angola$trophicTable 
#' meta_angola = build_metanet(metaweb,abTable,trophicTable)
#' print(meta_angola)
#' 
#' @importFrom igraph V get.adjacency vcount
#' @export
build_metanet <- function(metaweb,abTable = NULL,trophicTable = NULL,
                              compute_local_nets = TRUE){
  
  if(inherits(metaweb,'igraph')){
    if(is.null(igraph::V(metaweb)$name)){
      warning('nodes of metaweb do not have names. Assigning integers as names')
      igraph::V(metaweb)$name = as.character(1:length(igraph::V(metaweb)))
    }
  }else if(inherits(metaweb, c("matrix","data.frame","dgCMatrix"))){
    if(is.null(colnames(metaweb))){
      warning('colnames of metaweb do not have names. Assigning integers as names')
      colnames(metaweb) = as.character(1:ncol(metaweb))
    }
    metaweb = igraph::graph_from_adjacency_matrix(as.matrix(metaweb),mode = "directed")
    metaweb = igraph::permute(metaweb,order(order(igraph::V(metaweb)$name)))
  }else{
   stop("metaweb should be a network of class igraph or a matrix/dataframe") 
  }
  if(!(igraph::is.connected(metaweb,mode = "weak"))){
    stop("metaweb must be connected, consider extracting the largest connected component")
  }
  if(!(igraph::is.directed(metaweb))){
    stop("metaweb must be directed")
  }
  if(is.null(abTable)){
    if(is.null(igraph::V(metaweb)$ab)){
      abTable = rep(1,length(igraph::V(metaweb)))
    }else{
      abTable = igraph::V(metaweb)$ab
    }
    abTable = t(as.matrix(abTable))
    colnames(abTable) = igraph::V(metaweb)$name
  } else{
    if(is.null(dim(abTable))){
      stop("abTable must be a matrix, not a vector")
    } else if(is.null(colnames(abTable))){
      stop("abTable must have colnames")
    } else if(!(setequal(colnames(abTable),igraph::V(metaweb)$name))) {
      stop("colnames of abTalbe must correspond to metaweb node names")
    } else{
      abTable = abTable[,igraph::V(metaweb)$name]
    }
  }
  if(!(is.null(trophicTable))){
    if(is.null(dim(trophicTable))){
      stop("trophicTable must be a matrix, not a vector")
    } else if(is.null(colnames(trophicTable))){
      stop("trophicTable must have colnames")
    } else if(!setequal(trophicTable[,1],igraph::V(metaweb)$name)) {
      stop("first column of trophicTable must correspond to metaweb node names") 
    } else{
      rownames(trophicTable) = trophicTable[,1]
      trophicTable = trophicTable[order(rownames(trophicTable)),]
      }
    # trophicTable = as.data.frame(trophicTable)
    # if(!(is.dag(graph_from_edgelist(rbind(trophicTable[,1:2],trophicTable[,2:3]))))){
    #   stop('trophic table must represent a hierarchy, thus group dependencies must be acyclic')
    # }
  }

  #attributing relative mean abundance at the nodes of the metaweb
  if(is.null(igraph::V(metaweb)$ab)){
    igraph::V(metaweb)$ab = colMeans(abTable)/sum(colMeans(abTable))
  }
  #attributing weights of 1 if null weights
  if(is.null(igraph::E(metaweb)$weight)){
    igraph::E(metaweb)$weight = rep(1,igraph::ecount(metaweb))
  }
  
  if(!(is.null(trophicTable))){
    metaweb = igraph::set_graph_attr(graph = metaweb,name = 'res',value = colnames(trophicTable)[1])
  }
  metaweb$name = 'metaweb'
  
  metanetwork = list(metaweb = metaweb, abTable = abTable,
                      trophicTable = trophicTable)
  
  if(nrow(abTable) > 1 && compute_local_nets){
   metanetwork = c(metanetwork, 
                   get_local_networks(metanetwork))  
  }
  class(metanetwork) = "metanetwork"
  return(metanetwork)
}

#' Test of belonging to class metanetwork
#'
#' Return a boolean indicating whether the object belongs to class \code{metanetwork}
#'
#' @param metanetwork the object to test
#' @return a boolean indicating whether the object belongs to class \code{metanetwork}
#'
#' @examples
#' library(metanetwork)
#' library(igraph)
#'  
#' g = make_ring(5,directed = TRUE)
#' meta = build_metanet(g)
#' is.metanetwork(meta)
#' #on Angola dataset
#' data("meta_angola")
#' is.metanetwork(meta_angola)
#' @export
is.metanetwork <- function(metanetwork){
  UseMethod("is.metanetwork",metanetwork)
}

#' @return \code{NULL}
#'
#' @rdname is.metanetwork
#' @exportS3Method is.metanetwork metanetwork
is.metanetwork.metanetwork <- function(metanetwork) inherits(metanetwork, "metanetwork")

get_local_networks <- function(metanetwork){
  gList = list()
  abTable = sweep(metanetwork$abTable, 1, rowSums(metanetwork$abTable), "/")
  for(k in 1:nrow(abTable)){
    gLoc = igraph::induced_subgraph(metanetwork$metaweb,colnames(abTable)[which(metanetwork$abTable[k,]>0)])
    gLoc = igraph::permute(gLoc,order(order(igraph::V(gLoc)$name)))
    gLoc = igraph::set_vertex_attr(gLoc, name = "ab",
                                   value = as.numeric(abTable[k,igraph::V(gLoc)$name]))
    gLoc = igraph::set_graph_attr(gLoc,
           name = "res", value = colnames(metanetwork$trophicTable)[1])
    gLoc$name = rownames(abTable)[k]
    gList = c(gList,list(gLoc))
  }
  names(gList) = rownames(abTable)
  return(gList)
}

#' extract networks from a metanetwork object
#'
#' Function to extract metawebs and local networks from a metanetwork object
#'
#' Return a list of 'igraph' objects
#'
#' @param metanetwork the object whose networks need to be extracted
#' @return a list of \code{igraph} objects with attributes computed by \code{metanetwork}
#'
#' @examples
#' library(metanetwork)
#' data("meta_angola")
#' nets = extract_networks(meta_angola)  
#' sapply(nets,class)
#' @export
extract_networks <- function(metanetwork){
  return(metanetwork[which(lapply(metanetwork,class) == 'igraph')])
}



#' @import igraph
