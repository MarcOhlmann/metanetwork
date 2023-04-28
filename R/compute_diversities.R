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


#' compute network diversity indices
#'
#' Function to compute network diversity indices based on Hill numbers following the method described in Ohlmann et al. 2019
#' 
#' This function compute diversity indices using Hill numbers on node and link abundances. Importantly, a viewpoint parameters \eqn{q} allows giving more weigth to abundant nodes/links. Given a network, we note \eqn{p_q} the abundance of node \eqn{q}
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
#' The overall \eqn{\alpha}-diversity in node abundances is  (for \eqn{q = 1}):
#' \deqn{A_{P} =  \exp(\sum_{q=1}^{Q} \sum_{k} - \frac{P_{qk}}{P_{++}} \log(\frac{P_{qk}}{P_{++}}) - \log(K))}
#' The overall \eqn{\alpha}-diversity in node abundances is  (for \eqn{q = 1}):
#' \deqn{A_{P} =  exp(\sum_{q=1}^{Q} \sum_{k=1}^{K} - \frac{P_{qk}}{P_{++}} log (  \frac{P_{qk}}{P_{++}}) - log(K) )}
#' where \eqn{P_{++}=\sum_{k} \sum_{q} P_{qk}}
#' The overall \eqn{\alpha}-diversities in link abundances and link probabilities are equal to:
#' \deqn{A_{L} =  exp (\sum_{q,l=1}^{Q} \sum_{k=1}^{K} - \frac{L_{qlk}}{L_{+++}} \log (\frac{L_{qlk}}{L_{+++}}) - log(K))}
#' where \eqn{L_{+++}=\sum_{k} \sum_{q,l} L_{qlk}}.
#' The \eqn{\gamma}-diversity in node  and link abundances are  (for \eqn{q = 1}):
#' \deqn{G_{P}=exp(\sum_{q=1}^{Q} -\frac{P_{q+}}{P_{++}} log  (\frac{P_{q+}}{P_{++}}))}
#' \deqn{G_{L}=exp(\sum_{q,l=1}^{Q} -\frac{L_{ql+}}{L_{+++}} log(\frac{L_{ql+}}{L_{+++}}))}
#' The \eqn{\beta}-diversity is then defined in a multiplicative way:
#' \deqn{B_{P}=\frac{G_P}{A_P}}
#' \deqn{B_{L}=\frac{G_P}{A_L}}
#' 
#' 
#' For more details on \eqn{\alpha}-,\eqn{\beta}- and \eqn{\gamma}-diversity, see Ohlmann et al. 2019.
#'   
#' @references Ohlmann, M., Miele, V., Dray, S., Chalmandrier, L., O'connor, L., & Thuiller, W. (2019). Diversity indices for ecological networks: a unifying framework using Hill numbers. Ecology letters, 22(4), 737-747.
#'
#' @param metanetwork object of class 'metanetwork'
#' @param q viewpoint parameter controlling the weight given to abundant species/groups and links, default is 1
#' @param res a vector containing the resolutions at which the diversities are computed
#' @return a \code{data.frame}
#' 
#' @seealso [compute_dis()]
#'
#' @examples
#' library(metanetwork)
#' library(igraph)
#' 
#' #on angola dataset
#' data("meta_angola")
#' compute_div(meta_angola,q = 1)
#' 
#' #computing diversities only at Phylum level
#' compute_div(meta_angola,q = 1,res = "Phylum")
#'
#' @export
compute_div <- function(metanetwork,q = 1,res = NULL){
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
  
  if(!(is.null(metanetwork$trophicTable))){ #if several resolutions
  diversites_df_P = matrix(0,nrow = 3 + nrow(metanetwork$abTable),
                         ncol = length(res_loc))
  diversites_df_L = matrix(0,nrow = 3 + nrow(metanetwork$abTable),
                           ncol = length(res_loc))
  colnames(diversites_df_P) = res_loc
  colnames(diversites_df_L) = res_loc
  rownames(diversites_df_P) = c("Gamma_P","mean_Alpha_P","Beta_P",paste0("Alpha_",rownames(metanetwork$abTable),"_P"))
  rownames(diversites_df_L) = c("Gamma_L","mean_Alpha_L","Beta_L",paste0("Alpha_",rownames(metanetwork$abTable),"_L"))
  diversites_df_P = as.data.frame(diversites_df_P)
  diversites_df_L = as.data.frame(diversites_df_L)
  #get node and link abundances
  P_L_list = metawebParams(metaweb_loc,networks_loc,res_loc)
  #node diversity at the different resolutions
  for(res in res_loc){
    #node diversity
    abg_P_loc = abgDecompQ(spxp = t(P_L_list$P[[res]]),q = q)
    diversites_df_P["Gamma_P",res] = abg_P_loc$Gamma
    diversites_df_P["mean_Alpha_P",res] = abg_P_loc$mAlpha
    diversites_df_P["Beta_P",res] = abg_P_loc$Beta
    diversites_df_P[paste0("Alpha_",rownames(metanetwork$abTable),"_P"),res] = abg_P_loc$Alphas
    #link diversity
    L_array = P_L_list$L[[res]]
    dim(L_array) = c(nrow(P_L_list$P[[res]])*nrow(P_L_list$P[[res]]),
                     ncol(P_L_list$P[[res]])) 
    colnames(L_array) = colnames(P_L_list$P[[res]])
    abg_L_loc = abgDecompQ(spxp = t(L_array),q = q)
    diversites_df_L["Gamma_L",res] = abg_L_loc$Gamma
    diversites_df_L["mean_Alpha_L",res] = abg_L_loc$mAlpha
    diversites_df_L["Beta_L",res] = abg_L_loc$Beta
    diversites_df_L[paste0("Alpha_",rownames(metanetwork$abTable),"_L"),res] = abg_L_loc$Alphas
  }
  return(list(nodes = diversites_df_P,links = diversites_df_L))
  }else{#single resolution case
    #get node and link abundances
    P_L_list = metawebParams(metaweb_loc,networks_loc)
    #node diversity
    abg_P_loc = abgDecompQ(spxp = t(P_L_list$P),q = q)
    #link diversity
    L_array = P_L_list$L
    dim(L_array) = c(nrow(P_L_list$P)*nrow(P_L_list$P),
                     ncol(P_L_list$P)) 
    colnames(L_array) = colnames(P_L_list$P)
    abg_L_loc = abgDecompQ(spxp = t(L_array),q = q)
    return(list(nodes = list(Gamma = abg_P_loc$Gamma, mean_Alpha = abg_P_loc$mAlpha,
                             Beta = abg_P_loc$Beta,Alphas = abg_P_loc$Alphas),
                links = list(Gamma = abg_L_loc$Gamma, mean_Alpha = abg_L_loc$mAlpha,
                             Beta = abg_L_loc$Beta,Alphas = abg_L_loc$Alphas)))
  }
}



# Functions
## divLeinster calculates the diversity of each site of a site by species matrix according to the q parameter according to Chao
## abgDecompQ performs a alpha, beta, gamma multiplicative decomposition using Chao diversity indices. 
##Allows a parametrization of the dominance effect

## chaoObjects is a data preparation function. It returns adequate arguments for abgDecompQ, BetaDisQ and divLeinster to perform a diversity analysis using Chao's diversity index.

# Arguments
## spxp : sites (row) by species (cols) matrix with or without rownames and colnames.
## check : arguments specifying if the arguments should be checked.

divLeinster <- function(spxp, Z = NULL,q = 1, check = TRUE){
  #Calcul the diversity of each site of sites by species matrix. 
  #spxp columns and Z rows and columns are assumed to be in the same order.
  Z <- diag(ncol(spxp))
  spxp <- sweep(spxp, 1, rowSums(spxp), "/")
  Zp <- Z %*% t(spxp)
  if (q != 1 & q != Inf){
    mat <- t(spxp) * (Zp)^(q-1)
    mat[is.na(mat)] <- 0
    D <- colSums(mat) ^ (1/(1-q))
  }
  if (q==Inf)  {
    D <- 1/ apply(Zp, 2, max)
  }  
  if (q == 1){
    D <- apply(Zp^t(spxp), 2, function(x) 1/prod(x))
  }
  return(D)
}

abgDecompQ <- function(spxp, Z=NULL, q=1, check=TRUE) {
  #Calcul the diversity of each site of sites by species matrix. 
  #spxp columns and Z rows/cols are assumed to be in the same order.
  Z <- diag(ncol(spxp))
  
  site.weight <- rep(1/nrow(spxp), nrow(spxp))
  spxp <- sweep(spxp, 1, rowSums(spxp), "/")
  
  gamma.ab <- colSums(sweep(spxp, 1, site.weight, "*"))
  
  Gamma <- divLeinster(t(as.matrix(gamma.ab)), Z=Z , q=q, check = FALSE)
  Alphas <- divLeinster(spxp, Z=Z , q=q, check = FALSE) 
  
  if (q != 1 & q != Inf) {
    mAlpha <- (sum(site.weight * (Alphas ^ (1 - q))))^(1 / (1 - q))
  }
  if (q==1){
    mAlpha <- exp(sum(site.weight * log(Alphas)))
  }
  if (q==Inf){
    mAlpha <- min(Alphas)
  }
  Beta <- Gamma / mAlpha
  
  names(Alphas) <- row.names(spxp)
  res <- list(Gamma=Gamma, Beta=Beta, mAlpha=mAlpha, Alphas=Alphas)
  return(res)
}
#function to extract group and link abundances
metawebParams <- function(metaweb_loc,networks_loc,res_loc = NULL){
  if(!(is.null(res_loc))){
    P_mat_list = list()
    L_array_list = list()
    ## get the L array and P mat for a list of graph
    for(res in res_loc){
      #get the metaweb at the current res
      metaweb_loc_res = metaweb_loc[sapply(metaweb_loc,function(g) g$res) == res][[1]]
      #get the local networks at the current resolution
      networks_loc_res = networks_loc[sapply(networks_loc,function(g) g$res) == res]
      n = igraph::vcount(do.call(igraph::union,networks_loc_res))
      #build abundance matrix
      P_mat = matrix(0, nrow = n, ncol = length(networks_loc_res))
      rownames(P_mat) = igraph::V(metaweb_loc_res)$name
      colnames(P_mat) = names(networks_loc_res)
      for(net_loc_name in names(networks_loc_res)){
        P_mat[igraph::V(networks_loc_res[[net_loc_name]])$name,net_loc_name] =
          igraph::V(networks_loc_res[[net_loc_name]])$ab
      }
      #build link abundance matrix
      L_array = array(0, dim=c(n,n,length(networks_loc_res))) #stacked adjacency matrix at a group level
      dimnames(L_array)[[1]] = if(n>1)  igraph::V(metaweb_loc_res)$name else list(igraph::V(metaweb_loc_res)$name)
      dimnames(L_array)[[2]] = if(n>1)  igraph::V(metaweb_loc_res)$name else list(igraph::V(metaweb_loc_res)$name)
      dimnames(L_array)[[3]] = names(networks_loc_res)
      for(net_loc_name in dimnames(L_array)[[3]]){
        L_array[,,net_loc_name] = (P_mat[,net_loc_name] %*% t(P_mat[,net_loc_name])) * 
          as.matrix(igraph::get.adjacency(metaweb_loc_res))
      }
      P_mat_list = c(P_mat_list,list(P_mat))
      L_array_list = c(L_array_list,list(L_array))
    }
    names(P_mat_list) = res_loc
    names(L_array_list) = res_loc
    return(list(P = P_mat_list, L = L_array_list))
  }else{
    n = igraph::vcount(metaweb_loc[[1]])
    #build abundance matrix
    P_mat = matrix(0, nrow = n, ncol = length(networks_loc))
    rownames(P_mat) = igraph::V(metaweb_loc[[1]])$name
    colnames(P_mat) = names(networks_loc)
    for(net_loc_name in names(networks_loc)){
        P_mat[igraph::V(networks_loc[[net_loc_name]])$name,net_loc_name] =
        igraph::V(networks_loc[[net_loc_name]])$ab
    }
    #build link abundance matrix
    L_array = array(0, dim=c(n,n,length(networks_loc))) #stacked adjacency matrix at a group level
    dimnames(L_array)[[1]] = if(n>1)  igraph::V(metaweb_loc[[1]])$name else list(igraph::V(metaweb_loc[[1]])$name)
    dimnames(L_array)[[2]] = if(n>1)  igraph::V(metaweb_loc[[1]])$name else list(igraph::V(metaweb_loc[[1]])$name)
    dimnames(L_array)[[3]] = names(networks_loc)
    for(net_loc_name in dimnames(L_array)[[3]]){
      L_array[,,net_loc_name] = (P_mat[,net_loc_name] %*% t(P_mat[,net_loc_name])) * 
        as.matrix(igraph::get.adjacency(metaweb_loc[[1]]))
    }
    return(list(P = P_mat, L = L_array))
  }
}



