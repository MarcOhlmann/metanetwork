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


#' compute pairwise network dissimilarity indices
#'
#' Function to compute pairwise network dissimilarity indices based on Hill numbers following the method described in Ohlmann et al. 2019
#' 
#' This function compute pairwise dissimilarity indices using Hill numbers on node and link abundances. Importantly, a viewpoint parameters \eqn{q} allows giving more weigth to abundant nodes/links. Given a network, we note \eqn{p_q} the abundance of node \eqn{q}
#'  (stored as node attribute \code{ab}) and \eqn{\pi_{ql}} interaction probability between nodes \eqn{q} and \eqn{l} (stored as edge attribute \code{weight}). The link abundance \eqn{L_{ql}} between nodes 
#' \eqn{q} and \eqn{l} is then:
#' \deqn{L_{ql} = \pi_{ql}p_q p_l}
#' Node diversity (for \eqn{q = 1}) is then computed as:
#' \deqn{D(p) = \exp (\sum_q - p_q \log p_q)}
#' Link diversity is computed as:
#' \deqn{D(L) = \exp (\sum_{ql} - \frac{L_{ql}}{C} \log L_{ql}{C})}
#' where \eqn{C} is the weighted connectance
#' \deqn{C = \sum_{ql} \pi_{ql}p_q p_l}
#' 
#' For \eqn{q = 1}, the pairwise node dissimilarity indices are computed from pairwise diversities as:
#' \deqn{\delta_P=\frac{\log(G_P)-log(A_P)}{\log 2}}
#' where \eqn{G_P} is node \eqn{\gamma}-diversity and \eqn{A_P} the node \eqn{\alpha}-diversity
#' 
#' Pairwise link diversity is:
#' \deqn{\delta_L=\frac{\log(G_L)-log(A_L)}{\log 2}}
#' where \eqn{G_P} is the link \eqn{\gamma}-diversity and \eqn{A_P} the link \eqn{\alpha}-diversity
#' 
#' For more details on \eqn{\alpha}-,\eqn{\beta}- and \eqn{\gamma}-diversity, see Ohlmann et al. 2019.
#' 
#' @references Ohlmann, M., Miele, V., Dray, S., Chalmandrier, L., O'connor, L., & Thuiller, W. (2019). Diversity indices for ecological networks: a unifying framework using Hill numbers. Ecology letters, 22(4), 737-747.
#'
#' @param metanetwork object of class 'metanetwork'
#' @param q viewpoint parameter controlling the weight given to abundant species/groups and links, default is 1
#' @param res a vector containing the resolutions at which the diversities are computed
#' @param ncores number of cores used for the computation, default is 4
#' @return a list of \code{data.frame} containing node and link pairwise dissimilarities 
#' 
#' @seealso [compute_div()]
#'
#' @examples
#' library(metanetwork)
#' library(igraph)
#' 
#' #on angola dataset
#' data(meta_angola)
#' compute_dis(meta_angola,q = 1,ncores = 1)
#' 
#' #computing dissimilarities only at Phylum level
#' compute_dis(meta_angola,q = 1,res = "Phylum",ncores = 1)
#' 
#' @export
compute_dis <- function(metanetwork,q = 1,res = NULL,ncores = 4){
  # get the local networks
  networks = metanetwork[lapply(metanetwork,class) == "igraph"]
  metaweb_names = names(metanetwork)[grep('metaweb',x = names(metanetwork))]
  networks_loc = networks[!(names(networks) %in% metaweb_names)]
  metaweb_loc = networks[names(networks) %in% metaweb_names]
  
  if(length(networks_loc) == 0){
    stop("To compute network diversity indices, you must have local networks")
  }
  
  if(is.null(res)){
    res_loc = unique(sapply(networks_loc,function(g) g$res))
  }else{
    if(!(res %in% colnames(metanetwork$trophicTable))){
      stop(paste0("res must be a vector with available resolutions: ",  paste(colnames(metanetwork$trophicTable),collapse = " ")))
    }else{
      res_loc = res
    }
  }
  
  
  if(!(is.null(metanetwork$trophicTable))){ #if several resolutions are available
    #get node and link abundances
    P_L_list = metawebParams(metaweb_loc,networks_loc,res_loc)
    N = nrow(metanetwork$abTable)
    #list to store the results
    dis_df_list_list = list()
    
    #loop over resolutions
    for(res in res_loc){
      dis_df_P = matrix(0,nrow = N, ncol = N)
      dis_df_L = matrix(0,nrow = N, ncol = N)
      colnames(dis_df_P) = rownames(metanetwork$abTable)
      rownames(dis_df_P) = rownames(metanetwork$abTable)
      dis_df_P = as.data.frame(dis_df_P)
      colnames(dis_df_L) = rownames(metanetwork$abTable)
      rownames(dis_df_L) = rownames(metanetwork$abTable)
      dis_df_L = as.data.frame(dis_df_L)
      #inde table for parallelisation
      ind_table = t(combn(1:N, 2))
      message(paste0("computing node pairwise dissimilarities at resolution: ",res))
      res_P_list = parallel::mclapply(1:nrow(ind_table),function(k)  
        compute_dis_loc(ind_table[k,],P_L_list = P_L_list,type = "P",q = q,res = res),mc.cores = ncores)
      
      message(paste0("computing link pairwise dissimilarities at resolution: ",res))
      res_L_list = parallel::mclapply(1:nrow(ind_table),function(k)  
        compute_dis_loc(ind_table[k,],P_L_list = P_L_list,type = "L",q = q,res = res),mc.cores = ncores)
      
      #filling the dis matrix
      dis_df_P[lower.tri(dis_df_P)] = unlist(res_P_list)
      dis_df_L[lower.tri(dis_df_L)] = unlist(res_L_list)
      
      dis_df_P = dis_df_P + t(dis_df_P)
      dis_df_L = dis_df_L + t(dis_df_L)
      
      dis_df_list = list(nodes = dis_df_P, links = dis_df_L)
      dis_df_list_list = c(dis_df_list_list,list(dis_df_list))
    }
    names(dis_df_list_list) = res_loc
    return(dis_df_list_list)
  }else{#single resolution case
    #get node and link abundances
    P_L_list = metawebParams(metaweb_loc,networks_loc)
    N = nrow(metanetwork$abTable)
    dis_df_P = matrix(0,nrow = N, ncol = N)
    dis_df_L = matrix(0,nrow = N, ncol = N)
    colnames(dis_df_P) = rownames(metanetwork$abTable)
    rownames(dis_df_P) = rownames(metanetwork$abTable)
    dis_df_P = as.data.frame(dis_df_P)
    colnames(dis_df_L) = rownames(metanetwork$abTable)
    rownames(dis_df_L) = rownames(metanetwork$abTable)
    dis_df_L = as.data.frame(dis_df_L)
    #inde table for parallelisation
    ind_table = t(combn(1:N, 2))
    message(paste0("computing node pairwise dissimilarities"))
    res_P_list = parallel::mclapply(1:nrow(ind_table),function(k)  
        compute_dis_loc(ind_table[k,],P_L_list = P_L_list,type = "P",q = q),mc.cores = ncores)
      
    message(paste0("computing link pairwise dissimilarities"))
    res_L_list = parallel::mclapply(1:nrow(ind_table),function(k)  
        compute_dis_loc(ind_table[k,],P_L_list = P_L_list,type = "L",q = q),mc.cores = ncores)
      
    #filling the dis matrix
    dis_df_P[lower.tri(dis_df_P)] = unlist(res_P_list)
    dis_df_L[lower.tri(dis_df_L)] = unlist(res_L_list)
    
    dis_df_P = dis_df_P + t(dis_df_P)
    dis_df_L = dis_df_L + t(dis_df_L)
      
    dis_df_list = list(nodes = dis_df_P, links = dis_df_L)
    return(dis_df_list)
    }
}

#function to compute pairwise dissimilarity (to execute in parallel)
compute_dis_loc <- function(index_vec,P_L_list,type = c("P","L"),q,res = NULL){
  ind_i = index_vec[1]
  ind_j = index_vec[2]
  if(type == "P"){
    if(is.null(res)){
      spxp.dummy= P_L_list$P[,c(ind_i,ind_j)]
    }else{
      spxp.dummy= P_L_list$P[[res]][,c(ind_i,ind_j)]
    }
  }
  if(type == "L"){
    if(is.null(res)){
      spxp.dummy = P_L_list$L[,,c(ind_i,ind_j)]
      spxp.dummy = aperm(spxp.dummy,c(2,1,3))  
      dim(spxp.dummy) = c(nrow(P_L_list$P)*nrow(P_L_list$P),2) 
      colnames(spxp.dummy) = colnames(P_L_list$P[,c(ind_i,ind_j)]) 
    }else{ 
      spxp.dummy = P_L_list$L[[res]][,,c(ind_i,ind_j)]
      spxp.dummy = aperm(spxp.dummy,c(2,1,3))  
      dim(spxp.dummy) = c(nrow(P_L_list$P[[res]])*nrow(P_L_list$P[[res]]),2) 
      colnames(spxp.dummy) = colnames(P_L_list$P[[res]][,c(ind_i,ind_j)]) 
    }
    #removing empty rows
    if(sum(rowSums(spxp.dummy)>0)<nrow(spxp.dummy)){
      spxp.dummy=spxp.dummy[-which(rowSums(spxp.dummy)==0),]
    }
  }
  div = abgDecompQ(t(spxp.dummy),q = q)
  if(q != 1){
    res_loc = 1-((1/div$Beta)^(q-1)-(1/2)^(q-1))/(1-(1/2)^(q-1))
  }
  if(q == 1){
    res_loc = (log(div$Gamma)-log(div$mAlpha))/(log(2))
  }
  return(res_loc)
}



# 
# # Functions
# ## divLeinster calculates the diversity of each site of a site by species matrix according to the q parameter according to Chao
# ## abgDecompQ performs a alpha, beta, gamma multiplicative decomposition using Chao diversity indices. 
# ##Allows a parametrization of the dominance effect
# 
# ## chaoObjects is a data preparation function. It returns adequate arguments for abgDecompQ, BetaDisQ and divLeinster to perform a diversity analysis using Chao's diversity index.
# 
# # Arguments
# ## spxp : sites (row) by species (cols) matrix with or without rownames and colnames.
# ## check : arguments specifying if the arguments should be checked.
# 
# divLeinster <- function(spxp, Z = NULL,q = 1, check = TRUE){
#   #Calcul the diversity of each site of sites by species matrix. 
#   #spxp columns and Z rows and columns are assumed to be in the same order.
#   Z <- diag(ncol(spxp))
#   spxp <- sweep(spxp, 1, rowSums(spxp), "/")
#   Zp <- Z %*% t(spxp)
#   if (q != 1 & q != Inf){
#     mat <- t(spxp) * (Zp)^(q-1)
#     mat[is.na(mat)] <- 0
#     D <- colSums(mat) ^ (1/(1-q))
#   }
#   if (q==Inf)  {
#     D <- 1/ apply(Zp, 2, max)
#   }  
#   if (q == 1){
#     D <- apply(Zp^t(spxp), 2, function(x) 1/prod(x))
#   }
#   return(D)
# }
# 
# abgDecompQ <- function(spxp, Z=NULL, q=2, check=TRUE) {
#   #Calcul the diversity of each site of sites by species matrix. 
#   #spxp columns and Z rows/cols are assumed to be in the same order.
#   Z <- diag(ncol(spxp))
#   
#   site.weight <- rep(1/nrow(spxp), nrow(spxp))
#   spxp <- sweep(spxp, 1, rowSums(spxp), "/")
#   
#   gamma.ab <- colSums(sweep(spxp, 1, site.weight, "*"))
#   
#   Gamma <- divLeinster(t(as.matrix(gamma.ab)), Z=Z , q=q, check = FALSE)
#   Alphas <- divLeinster(spxp, Z=Z , q=q, check = FALSE) 
#   
#   if (q != 1 & q != Inf) {
#     mAlpha <- (sum(site.weight * (Alphas ^ (1 - q))))^(1 / (1 - q))
#   }
#   if (q==1){
#     mAlpha <- exp(sum(site.weight * log(Alphas)))
#   }
#   if (q==Inf){
#     mAlpha <- min(Alphas)
#   }
#   Beta <- Gamma / mAlpha
#   
#   names(Alphas) <- row.names(spxp)
#   res <- list(Gamma=Gamma, Beta=Beta, mAlpha=mAlpha, Alphas=Alphas)
#   return(res)
# }
# #function to extract group and link abundances
# metawebParams <- function(metaweb_loc,networks_loc,res_loc){
#   P_mat_list = list()
#   L_array_list = list()
#   ## get the L array and P mat for a list of graph
#   for(res in res_loc){
#     #get the metaweb at the current res
#     metaweb_loc_res = metaweb_loc[sapply(metaweb_loc,function(g) g$res) == res][[1]]
#     #get the local networks at the current resolution
#     networks_loc_res = networks_loc[sapply(networks_loc,function(g) g$res) == res]
#     n = igraph::vcount(do.call(igraph::union,networks_loc_res))
#     #build abundance matrix
#     P_mat = matrix(0, nrow = n, ncol = length(networks_loc_res))
#     rownames(P_mat) = igraph::V(metaweb_loc_res)$name
#     colnames(P_mat) = names(networks_loc_res)
#     for(net_loc_name in names(networks_loc_res)){
#       P_mat[igraph::V(networks_loc_res[[net_loc_name]])$name,net_loc_name] =
#         igraph::V(networks_loc_res[[net_loc_name]])$ab
#     }
#     #build link abundance matrix
#     L_array = array(0, dim=c(n,n,length(networks_loc_res))) #stacked adjacency matrix at a group level
#     dimnames(L_array)[[1]] = if(n>1)  igraph::V(metaweb_loc_res)$name else list(igraph::V(metaweb_loc_res)$name)
#     dimnames(L_array)[[2]] = if(n>1)  igraph::V(metaweb_loc_res)$name else list(igraph::V(metaweb_loc_res)$name)
#     dimnames(L_array)[[3]] = names(networks_loc_res)
#     for(net_loc_name in dimnames(L_array)[[3]]){
#       L_array[,,net_loc_name] = (P_mat[,net_loc_name] %*% t(P_mat[,net_loc_name])) * 
#         as.matrix(igraph::get.adjacency(metaweb_loc_res))
#     }
#     P_mat_list = c(P_mat_list,list(P_mat))
#     L_array_list = c(L_array_list,list(L_array))
#   }
#   names(P_mat_list) = res_loc
#   names(L_array_list) = res_loc
#   return(list(P = P_mat_list, L = L_array_list))
# }
