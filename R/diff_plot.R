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

#' plot difference network
#'
#' Function to represent difference between two networks belonging to a metanetwork 
#' with specific layout ('TL-tsne' or group 'TL-tsne') using either 'ggnet' or 'visNetwork' visualisation.
#'  This function represent the difference between g1 and g2 (g1-g2). 
#'
#' @param metanetwork object of class 'metanetwork'
#' @param g1 network (of class 'igraph') of metanetwork
#' @param g2 network (of class 'igraph') of metanetwork
#' @param beta the diffusion parameter of the diffusion kernel, a positive scalar controlling the 
#' squeezing of the network
#' @param mode mode used for layout, either 'TL-tsne' or 'group-TL-tsne' (see \code{attach_layout()}). Default is 'TL-tsne'
#' @param vis_tool a character indicating the visualisation tool, either 'ggnet' or visNetwork
#' @param edge_thrs if non-null, a numeric (between 0 and 1) indicating an edge threshold for the representation
#' @param layout_metaweb a boolean indicating whether the layout of the metaweb should be used to represent the difference network.
#' to use metaweb layout = T, you need first to compute 'TL-tsne' layout for the metaweb for this beta value using \code{attach_layout()}
#' @param nrep_ly If several layouts for this beta value are attached to the metaweb 
#' (if \code{layout_metaweb = T}), index of the layout to use, see \code{attach_layout()}
#' @param flip_coords a boolean indicating whether coordinates should be flipped. 
#' In that case, y-axis is the trophic level and x-axis is the layout axis
#' @param alpha_per_group controlling alpha per group (only for 'ggnet' vis), a list of format 
#' \code{list(resolutions = "XX",groups = XX,alpha_focal = XX,alpha_hidden = XX)}, see example
#' @param alpha_per_node controlling alpha per node (only for 'ggnet' vis), a list of format 
#' \code{list(nodes = XX,alpha_focal = XX,alpha_hidden = XX)}, see example
#' @param TL_tsne.config configuration list for mode 'TL-tsne', default is TL_tsne.default
#' @param ggnet.config configuration list for ggnet representation, default is ggnet.default
#' @param visNetwork.config configuration list for visNetwork representation, default is visNetwork.default
#' 
#' @seealso [attach_layout()]
#' 
#' @return an object of class \code{ggplot} or \code{visNetwork}, representation of the difference network
#'
#' @examples
#' #on Angola dataset
#' library(igraph)
#' library(metanetwork)
#' 
#' data(meta_angola)
#'
#' diff_plot(g1 = meta_angola$X2003,g2 = meta_angola$X1986,metanetwork = meta_angola,
#' beta = 0.05)
#'
#' 
#' 
#' @export
#' 
#' 

diff_plot <- function(metanetwork,g1,g2,beta = 0.1,mode ='TL-tsne',
                      vis_tool = "ggnet",
                      edge_thrs = NULL,layout_metaweb = FALSE,flip_coords = FALSE,
                      alpha_per_group = NULL,alpha_per_node = NULL,
                      TL_tsne.config = TL_tsne.default,nrep_ly = 1,
                      ggnet.config = ggnet.default,visNetwork.config = visNetwork.default){
  
  if(!is.metanetwork(metanetwork)){
    stop("metanetwork is an object of class metanetwork, see build_metanet")
  }
  message(paste0("mode is ",mode))
  if(!(mode %in% c('TL-tsne','group-TL-tsne','fr','kk','circle'))){
    stop("mode must be one of: \n
         'TL-tsne','fr','kk','circle")
  }
  
  networks = extract_networks(metanetwork)
  if(!(prod(c(g1$name,g2$name) %in% sapply(networks,function(g) g$name)))){
    stop('g1 and g2 must belongs to metanetwork')
  } else{
      #special case : uniform abundance (binary data)
      if(length(unique(igraph::V(g1)$ab)) == 1 && length(unique(igraph::V(g2)$ab)) == 1){
        igraph::V(g1)$ab = 1
        igraph::V(g2)$ab = 1
      }
      g_union = igraph::union(g1,g2)
      #test if g_union is conneted
      if(!(igraph::is.connected(g_union)) && !(layout_metaweb)){
        stop("the union network is not connected, you must use layout_metaweb=T to represent it (you then also have to compute layout for 
             the metaweb using attach_layout)")
      }
      g_union = igraph::permute(g_union,order(order(igraph::V(g_union)$name)))
      #ab of gUnion: V(g1)$ab - V(g2)$ab if the node is present in both networks
      # -10 if node is absent
      igraph::V(g_union)$ab_1[is.na(igraph::V(g_union)$ab_1)] = -10
      igraph::V(g_union)$ab_2[is.na(igraph::V(g_union)$ab_2)] = -10
      
      igraph::V(g_union)$ab = igraph::V(g_union)$ab_1 - igraph::V(g_union)$ab_2

      
      #weights of gUnion: E(g1)$weight - E(g2)$weight if the edge is present in both networks
      # -10 if edge is absent
      igraph::E(g_union)$weight_1[is.na(igraph::E(g_union)$weight_1)] = -10
      igraph::E(g_union)$weight_2[is.na(igraph::E(g_union)$weight_2)] = -10
      igraph::E(g_union)$weight_col = igraph::E(g_union)$weight_1 - igraph::E(g_union)$weight_2
      igraph::E(g_union)$weight_col_bis = igraph::E(g_union)$weight_col
      #when weights are equal, E(g_union)$weight_col == 0, 
      #setting this value to weight in g1 (or g2) for representaiton purpose
      igraph::E(g_union)$weight_col_bis = ifelse(igraph::E(g_union)$weight_col_bis == 0,
                                                 igraph::E(g_union)$weight_1,
                                                 igraph::E(g_union)$weight_col)
      #when edge is absent from g2 or g1, setting weight_col_bis to value in g1 or g2
      igraph::E(g_union)$weight_col_bis = ifelse(igraph::E(g_union)$weight_col_bis > 5,
                                         igraph::E(g_union)$weight_1,
                                         igraph::E(g_union)$weight_col_bis)
      igraph::E(g_union)$weight_col_bis = ifelse(igraph::E(g_union)$weight_col_bis < -5,
                                                 igraph::E(g_union)$weight_2,
                                                 igraph::E(g_union)$weight_col_bis)
      if(is.null(metanetwork$trophicTable)){
        metanetwork_diff = build_metanet(g_union)  
      }else{
        if(g_union$res_1 == colnames(metanetwork$trophicTable)[1]){
          #append trophicTable to g_union, keep only focal nodes
          metanetwork_diff = build_metanet(metaweb = g_union,
                                               trophicTable = metanetwork$trophicTable[
                                                 igraph::V(g_union)$name,])
        }else{
          metanetwork_diff = build_metanet(g_union)
        }
      }
      #associate layout of the metaweb to metanetwork if layout_metaweb = T
      if(layout_metaweb){
        #get the attributes of the metaweb at the resolution of the networks
        if(is.null(metanetwork$metaweb$res)){
          ind_metaweb = which(sapply(networks,function(g) g$name) == "metaweb")
          current_metaweb = networks[[ind_metaweb]]
        }else{
          ind_metaweb = intersect(which(sapply(networks,function(g) g$res) == g1$res),
                                  which(sapply(networks,function(g) g$name) == "metaweb"))
          current_metaweb = networks[[ind_metaweb]]
        }
        attr_names = igraph::vertex_attr_names(current_metaweb)
        if(length(grep(paste0("beta",beta),attr_names)) == 0){
          stop("to use 'layout_metaweb = T', you need to attach a layout to the metaweb
        for the desired beta value and resolution, see attach_layout function")
        } else{
        #set trophic levels as metaweb trophic levels
        metaweb_TL = igraph::V(current_metaweb)$TL
        names(metaweb_TL) = igraph::V(current_metaweb)$name
        igraph::V(metanetwork_diff$metaweb)$TL =
          metaweb_TL[igraph::V(metanetwork_diff$metaweb)$name]
        if(mode == "TL-tsne"){
          layout_loc = igraph::get.vertex.attribute(graph = current_metaweb,
                                                    name = paste0("layout_beta",beta))
          names(layout_loc) = V(current_metaweb)$name
          metanetwork_diff$metaweb =  
            igraph::set_vertex_attr(graph = metanetwork_diff$metaweb,
                                    name = paste0("layout_beta",beta),
                                    value = layout_loc[V(metanetwork_diff$metaweb)$name])
        } else if(mode == "group-TL-tsne"){
          layout_loc = sapply(c("x","y"),
                              function(k) igraph::get.vertex.attribute(graph = current_metaweb,
                                                                       name = paste0("group_layout_",k,"_beta",beta)))
          rownames(layout_loc) = V(current_metaweb)$name
          metanetwork_diff$metaweb =  
            igraph::set_vertex_attr(graph = metanetwork_diff$metaweb,
                                    name = paste0("group_layout_x_beta",beta),
                                    value = layout_loc[V(metanetwork_diff$metaweb)$name,"x"]) 
          metanetwork_diff$metaweb =  
            igraph::set_vertex_attr(graph = metanetwork_diff$metaweb,
                                    name = paste0("group_layout_y_beta",beta),
                                    value = layout_loc[V(metanetwork_diff$metaweb)$name,"y"])
        }

        }
      } else{
        #compute trophic level of the difference network if layout_metaweb = F
        if(mode == "TL-tsne"){
          igraph::V(metanetwork_diff$metaweb)$TL = compute_TL_diff(metanetwork_diff,metanetwork)
        } else if(mode == "group-TL-tsne"){
          stop("to use 'group-TL-tsne' layout, you need to set metaweb_layout = T")
        }
      } 
      if(is.null(g1$res)){
        message(paste0('plotting: ',g1$name,' - ',g2$name))
      }else{message(paste0('plotting: ',g1$name,'_',g1$res,' - ',g2$name,'_',g2$res))}
      
      #choice of visualisation tool
      if(vis_tool == "ggnet"){
        return(ggmetanet(metanetwork = metanetwork_diff,mode = mode,diff_plot_bool = T,
                         ggnet.config = ggnet.config,TL_tsne.config = TL_tsne.config,
                         beta = beta,edge_thrs = edge_thrs,layout_metaweb,flip_coords = flip_coords,
                         alpha_per_node = alpha_per_node,alpha_per_group = alpha_per_group))
      }else if(vis_tool == "visNetwork"){
        return(vismetaNetwork(metanetwork = metanetwork_diff,mode = mode,diff_plot_bool = T,
                         visNetwork.config = visNetwork.config,TL_tsne.config = TL_tsne.config,
                         beta = beta,edge_thrs = edge_thrs,layout_metaweb = layout_metaweb,flip_coords = flip_coords,x_y_range = c(100,100)))
      }
  }
}

