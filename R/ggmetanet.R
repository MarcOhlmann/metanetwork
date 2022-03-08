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


#' Default configuration for ggnet
#'
#' A list with parameters customizing ggmetanet representation (see ggnet documentations)
#'
#' @examples
#' # display all default settings
#' ggnet.default
#'
#' # create a new settings
#' ggnet.custom = ggnet.default
#' ggnet.custom$edge.size = 2
#' ggnet.custom
#' 
#' @export
ggnet.default = list(
  label = TRUE,
  label.size = 3,
  max_size = 5,
  edge.size = 0.5,
  arrow.size = 6, 
  arrow.gap = 0.015,
  alpha = 0.8,
  edge.alpha = 0.5,
  alpha_diff = 0.8,
  edge.alpha_diff = 0.8,
  size.cut = 5,
  palette = "Set2",
  default.color = "grey75",
  legend.position = "bottom")

class(ggnet.default) = 'metanetwork_config'



#' ggmetanet
#'
#' Function that provides network static representation  (using 'ggnet') from a 
#' 'metanetwork' object with a layout based on a diffusion kernel
#'
#' @param metanetwork object of class metanetwork
#' @param g network (igraph object) to represent, default is metaweb
#' @param beta the diffusion parameter of the diffusion kernel, a positive scalar controlling the 
#' vertical squeezing of the network
#' @param legend resolution for the legend, legend resolution must be a coarser resolution than the resolution of g, default is NULL
#' @param mode mode used for layout, 'TL-tsne' for trophic level t-sne and 'TL-kpco' for trophic level kernel based pco. Default is 'TL-tsne'
#' @param edge_thrs if non-null, a numeric (between 0 and 1) indicating an edge threshold for the representation
#' @param layout_metaweb a boolean indicating whether the layout of the metaweb should be used to represent the network
#' to use metaweb layout = T, you need first to compute metaweb layout for this beta value using \code{attach_layout()}
#' @param nrep_ly If several layouts for this beta value are attached to the metaweb (if \code{layout_metaweb = T}), index of the layout to use, see \code{attach_layout()}
#' @param flip_coords a boolean indicating wheter coordinates should be flipped. 
#' @param alpha_per_group controlling alpha per group (only for 'ggnet' vis), a list of format 
#' \code{list(resolutions = "XX",groups = XX,alpha_focal = XX,alpha_hidden = XX)}, see example
#' @param alpha_per_node controlling alpha per node (only for 'ggnet' vis), a list of format 
#' \code{list(nodes = XX,alpha_focal = XX,alpha_hidden = XX)}, see example
#' In that case, y-axis is the trophic level and x-axis is the layout axis
#' @param alpha_interactive a boolean indicating whether alpha (that is node transparency) 
#' should be asked in interactive mode to the user
#' @param ggnet.config configuration list for ggnet representation, default is ggnet.default
#' @param TL_tsne.config configuration list for mode 'TL-tsne', default is TL_tsne.default
#' @param diff_plot_bool boolean, do not edit by hand
#' 
#' @return object of class 'ggplot' 
#'
#' @examples
#' library(metanetwork)
#' library(igraph)
#' 
#' g = make_ring(5,directed = TRUE)
#' meta0 = build_metanet(g)
#' meta0 = compute_TL(meta0)
#' ggmetanet(meta0)
#'
#'# angola dataset
#'meta_angola = compute_TL(meta_angola)
#'ggmetanet(meta_angola,legend = 'Phylum',beta = 0.05)
#'
#' @export
ggmetanet <- function(metanetwork,g = NULL,beta = 0.1,
                           legend = NULL,mode = 'TL-tsne',
                           edge_thrs = NULL,
                           layout_metaweb = F,nrep_ly = 1,
                           flip_coords = F,
                           diff_plot_bool = F,
                           alpha_per_group = NULL,alpha_per_node = NULL,
                           alpha_interactive = F,
                           ggnet.config = ggnet.default,
                           TL_tsne.config = TL_tsne.default){

  if(!is.metanetwork(metanetwork)){
    stop("metanetwork is an object of class metanetwork, see build_metanetwork")
  }
  message(paste0("mode is ",mode))
  if(!(mode %in% c('TL-tsne','TL-kpco','fr','kk','circle'))){
    stop("mode must be one of: \n
         'TL-tsne','TL-kpco','fr','kk','circle")
  }
  
  if(!(diff_plot_bool)){
    if(is.null(g)){
      g = metanetwork$metaweb
    } 
    #simplify the network (remove self loops)
    if(!(igraph::is.simple(g))){
      g = igraph::simplify(g)
    }
    #edge weight threshold
    if(!(is.null(edge_thrs))){
      g = igraph::delete.edges(graph = g,edges =  which(igraph::E(g)$weight < edge_thrs))
    }
    if(!(is.null(legend))){
      # get the resolution of g
      res_local = g$res
      #legend must be one of the available resolution
      if(!(legend %in% colnames(metanetwork$trophicTable))){
        if(is.null(metanetwork$trophicTable)){
          stop('single resolution available, legend is not possible')
        }else{
          stop(paste0('legend must be one of the available resolutions : ',
                      Reduce(paste,colnames(metanetwork$trophicTable))))
        }
       }
      #legend must be a coarser resolution than resolution of the current network
      if(!(which(colnames(metanetwork$trophicTable) == res_local) <
           which(colnames(metanetwork$trophicTable) == legend))){
        stop("legend must be a coarser resolution than resolution of the current network")
      } else{
        # if(!(is.null(groups))){
          #get the largest connected components only
          # g = igraph::induced_subgraph(g,
          #                              which(metanetwork$trophicTable[igraph::V(g)$name,legend] %in% groups))
          # components_loc = igraph::components(g,mode='weak')
          # g = igraph::induced_subgraph(g,names(which(components_loc$membership == 1)))
          # trophic_table_loc = unique(metanetwork$trophicTable[,c(res_local,legend)])
          # rownames(trophic_table_loc) = trophic_table_loc[,1]
          # color_loc = trophic_table_loc[igraph::V(g)$name,legend]
        # } else {
          trophic_table_loc = unique(metanetwork$trophicTable[,c(res_local,legend)])
          rownames(trophic_table_loc) = trophic_table_loc[,1]
          color_loc = trophic_table_loc[igraph::V(g)$name,legend]
        # }
      }
    }
  }else{
    #metaweb of the two networks
    g = metanetwork$metaweb
    #simplify the network (remove self loops) 
    if(!(igraph::is.simple(g))){
      g = igraph::simplify(g,edge.attr.comb = "concat")
    }
    #edge weight threshold
    if(!(is.null(edge_thrs))){
      g = igraph::delete.edges(graph = g,
                               edges =  which(abs(igraph::E(g)$weight_col) < edge_thrs))
    }
    color_loc =  ifelse(igraph::V(g)$ab>0,'more abundant in g1','more abundant in g2')
    color_loc[which(igraph::V(g)$ab > 2)] = 'only present in g1'
    color_loc[which(igraph::V(g)$ab < -2)] = 'only present in g2'
    color_loc[which(igraph::V(g)$ab ==0)] = 'shared'
    edge_color_loc = ifelse(igraph::E(g)$weight_col>0,'#a1d99b','#fc9272')
    edge_color_loc[which(igraph::E(g)$weight_col == 0)] = 'black'
    edge_color_loc[which(igraph::E(g)$weight_col > 2)] = '#31a354'
    edge_color_loc[which(igraph::E(g)$weight_col < -2)] = '#de2d26'
    if(length(igraph::E(g)$weight_col[which(igraph::E(g)$weight_col == 0)])>0){
      igraph::E(g)$weight_col[which(igraph::E(g)$weight_col == 0)] =  
        igraph::E(g)$weight_col_bis[which(igraph::E(g)$weight_col == 0)]
    }
    if(length(igraph::E(g)$weight_col[which(igraph::E(g)$weight_col < 2)])>0){
      igraph::E(g)$weight_col[which(igraph::E(g)$weight_col < 2)] =
        igraph::E(g)$weight_col_bis[which(igraph::E(g)$weight_col < 2)]
    }
    if(length(igraph::E(g)$weight_col[which(igraph::E(g)$weight_col > 2)])>0){
      igraph::E(g)$weight_col[which(igraph::E(g)$weight_col > 2)] =
        igraph::E(g)$weight_col_bis[which(igraph::E(g)$weight_col > 2)]
    }
    igraph::E(g)$weight = abs(igraph::E(g)$weight_col)
    igraph::V(g)$ab = abs(igraph::V(g)$ab)
  }
  
  # if(!(is.null(vertex_names))){
  #   if(length(vertex_names) == length(V(g))){
  #     network.vertex.names(g_Network) = vertex_names
  #   }
  #   if(!(vertex_names)){network.vertex.names(g_Network) = rep("",length(V(g)))}
  # }

  
  if(mode %in% c("fr","kk","circle")){
    if(mode == "fr"){
      mode_loc = "fruchtermanreingold"
    } else if(mode == "kk"){
      mode_loc = "kamadakawai"
    } else{
      mode_loc = mode
    }
  }
  if(mode == "TL-kpco"){
    if(is.null(igraph::V(g)$TL)){
      if(!diff_plot_bool){
        stop("you must compute trophic levels first, see compute_TL")  
      }
      igraph::V(g)$TL = compute_TL_laplacian(g,metanetwork)
    }
    mode_loc =  get_nodes_position_TL_kpco(g = g,TL = igraph::V(g)$TL,beta = beta)
    rownames(mode_loc) = igraph::V(g)$name
  }
  if(mode == 'TL-tsne'){
    if(is.null(igraph::V(g)$TL)){
      stop("you must compute trophic levels first, see compute_TL")
    }
    #layout_metaweb option: get metaweb layout at the current resolution
    if(layout_metaweb){
      if(!(diff_plot_bool)){
        networks = extract_networks(metanetwork)
        if(is.null(g$res)){
          metaweb = networks[[which(sapply(networks,function(g) g$name) == "metaweb")]]
        } else{
          metaweb = networks[[intersect(which(sapply(networks,function(g) g$res) %in% g$res),
                                        which(sapply(networks,function(g) g$name) == "metaweb"))]]
        }
        attr_names = igraph::vertex_attr_names(metaweb)
        if(length(grep(paste0("beta",beta),attr_names)) == 0){
          stop("to use 'metaweb_layout = T', you need to attach a layout to the metaweb
        for the desired beta value and resolution, see attach_layout function")
        }else{
          mode_loc_metaweb = cbind(igraph::V(metaweb)$TL,
                                   igraph::get.vertex.attribute(
                                     metaweb,attr_names[grep(paste0("beta",beta),attr_names)[nrep_ly]]))
          rownames(mode_loc_metaweb) = igraph::V(metaweb)$name
          mode_loc = mode_loc_metaweb[V(g)$name,]
        }
      } else{ #diff_plot case
        mode_loc = cbind(igraph::V(metanetwork$metaweb)$TL,
                         igraph::get.vertex.attribute(metanetwork$metaweb,paste0("layout_beta",beta)))
        rownames(mode_loc) = igraph::V(metanetwork$metaweb)$name
      }
    }else{
      attr_names = igraph::vertex_attr_names(g)
      #check if a layout is attached to the focal network
      if(length(grep(paste0("beta",beta),attr_names))>0){
        mode_loc = cbind(igraph::V(g)$TL,
                         igraph::get.vertex.attribute(g,attr_names[grep(paste0("beta",beta),attr_names)[nrep_ly]]))
        rownames(mode_loc) = igraph::V(g)$name
      }else{
        mode_loc = get_nodes_position_TL_tsne(g = g,TL = igraph::V(g)$TL,beta = beta,
                                              TL_tsne.config = TL_tsne.config)
        rownames(mode_loc) = igraph::V(g)$name
      }
    }
  }
  
  igraph::E(g)$weight = ggnet.config$edge.size*igraph::E(g)$weight
  g_Network <- intergraph::asNetwork(g)
  #possible to flip coordinates
  if(flip_coords){
    mode_loc = cbind(mode_loc[,2],mode_loc[,1])
  }
  
  if(diff_plot_bool){
    #for legend
    network::set.vertex.attribute(x = g_Network,attrname = "color",value = color_loc)
    alpha = ifelse(color_loc == 'shared',ggnet.config$alpha,ggnet.config$alpha_diff)
    #assigning max_size if single abundance value
    size_loc = as.numeric(igraph::V(g)$ab/max(igraph::V(g)$ab))
    if(length(unique(size_loc)) == 1){
      size_loc = ggnet.config$max_size
    }
    #enhance groups or nodes
    if(!(is.null(alpha_per_group) && is.null(alpha_per_node))){
      if(!is.null(alpha_per_group)){
        if(is.null(alpha_per_group$resolution)){
          if(ncol(metanetwork$trophicTable) == 2){
            legend = colnames(metanetwork$trophicTable)[2]
          } else {stop("you must provide resolution in alpha_per_group")}
        } else{
          legend = alpha_per_group$resolution
        }
        groups_focal = alpha_per_group$groups
        alpha_focal = alpha_per_group$alpha_focal
        alpha_hidden = alpha_per_group$alpha_hidden
        alpha_bis = ifelse(metanetwork$trophicTable[rownames(mode_loc),legend] %in% groups_focal,
                       as.numeric(alpha_focal),as.numeric(alpha_hidden))
        alpha = alpha*alpha_bis
        edge_alpha = (alpha %*% t(alpha)) * igraph::get.adjacency(g) * ggnet.config$edge.alpha
        g_Network = network::set.edge.value(g_Network,"edge_alpha",edge_alpha)
        edge_alpha_vec = network::get.edge.attribute(g_Network,'edge_alpha')
        
        return(GGally::ggnet2(g_Network,mode = mode_loc,color = "color",edge.color = edge_color_loc,
                              size = size_loc,edge.size = "weight",
                              label = ggnet.config$label, label.size = ggnet.config$label.size,
                              max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                              arrow.size = ggnet.config$arrow.size,
                              arrow.gap = ggnet.config$arrow.gap,
                              alpha = alpha,
                              edge.alpha = edge_alpha_vec*ifelse(edge_color_loc == 'black',ggnet.config$edge.alpha,ggnet.config$edge.alpha_diff),
                              legend.position = ggnet.config$legend.position,
                              palette = c("only present in g1" = "#31a354", "more abundant in g1" = "#a1d99b","more abundant in g2" = "#fc9272",
                                          "only present in g2" = "#de2d26","shared" = "grey75")))
      } else{
        nodes_focal = alpha_per_node$nodes
        alpha_focal = alpha_per_node$alpha_focal
        alpha_hidden = alpha_per_node$alpha_hidden
        alpha_bis = ifelse(V(g)$name %in% nodes_focal,
                       as.numeric(alpha_focal),as.numeric(alpha_hidden))
        alpha = alpha_bis*alpha
        edge_alpha = (alpha %*% t(alpha)) * igraph::get.adjacency(g) * ggnet.config$edge.alpha
        g_Network = network::set.edge.value(g_Network,"edge_alpha",edge_alpha)
        edge_alpha_vec = network::get.edge.attribute(g_Network,'edge_alpha')
        return(GGally::ggnet2(g_Network,mode = mode_loc,color = "color",edge.color = edge_color_loc,
                              size = as.numeric(igraph::V(g)$ab/max(igraph::V(g)$ab)),
                              label = ggnet.config$label, 
                              label.size = ggnet.config$label.size*alpha,
                              max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                              edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                              arrow.gap = ggnet.config$arrow.gap,
                              alpha = alpha,
                              edge.alpha = edge_alpha_vec,
                              legend.position = ggnet.config$legend.position,
                              palette = c("only present in g1" = "#31a354", "more abundant in g1" = "#a1d99b","more abundant in g2" = "#fc9272",
                                          "only present in g2" = "#de2d26","shared" = "grey75")))
      }
    } else{
      return(GGally::ggnet2(g_Network,mode = mode_loc,color = "color",edge.color = edge_color_loc,
                     size = size_loc,edge.size = "weight",
                     label = ggnet.config$label, label.size = ggnet.config$label.size,
                     max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                     arrow.size = ggnet.config$arrow.size,
                     arrow.gap = ggnet.config$arrow.gap,
                     alpha = alpha,
                     edge.alpha = ifelse(edge_color_loc == 'black',ggnet.config$edge.alpha,ggnet.config$edge.alpha_diff),
                     legend.position = ggnet.config$legend.position,
                     palette = c("only present in g1" = "#31a354", "more abundant in g1" = "#a1d99b","more abundant in g2" = "#fc9272",
                                 "only present in g2" = "#de2d26","shared" = "grey75")))
    }

  }else{
    if(!(is.null(alpha_per_group) && is.null(alpha_per_node))){
      if(!is.null(alpha_per_group)){
        groups_focal = alpha_per_group$groups
        alpha_focal = alpha_per_group$alpha_focal
        alpha_hidden = alpha_per_group$alpha_hidden
        alpha = ifelse(metanetwork$trophicTable[rownames(mode_loc),legend] %in% groups_focal,
                       as.numeric(alpha_focal),as.numeric(alpha_hidden))
        
        edge_alpha = (alpha %*% t(alpha)) * igraph::get.adjacency(g) * ggnet.config$edge.alpha
        g_Network = network::set.edge.value(g_Network,"edge_alpha",edge_alpha)
        edge_alpha_vec = network::get.edge.attribute(g_Network,'edge_alpha')
        
        nb_cols = length(unique(color_loc))
        #use colors only if reasonible number of colors
        if(nb_cols < 15){
          mycolors = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, ggnet.config$palette))(nb_cols)
          names(mycolors) = unique(color_loc)
          
          return(GGally::ggnet2(g_Network,mode = mode_loc,color = color_loc,
                                size = as.numeric(igraph::V(g)$ab/max(igraph::V(g)$ab)),
                                label = ggnet.config$label, 
                                label.size = ggnet.config$label.size,
                                max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                                edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                                arrow.gap = ggnet.config$arrow.gap,
                                alpha = alpha,
                                edge.alpha = edge_alpha_vec,
                                palette = mycolors,
                                legend.position = ggnet.config$legend.position)
                )
        } else{ #use colors and shapes
          #assigning shapes (four types 15,16,17,18) and colors
          
          groups_loc = unique(color_loc)
          n_groups_loc = length(groups_loc)
          
          #shpes
          k = floor(n_groups_loc/4)
          shapes_loc = c(rep(15,k),rep(16,k),rep(17,k),
                         rep(18,n_groups_loc-3*k)) 
          names(shapes_loc) = groups_loc
          shapes_loc_nodes = sapply(metanetwork$trophicTable[V(g)$name,legend],
                                    function(x) shapes_loc[x])
          names(shapes_loc_nodes) = V(g)$name
          #colors
          mycolors = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, ggnet.config$palette))(floor(nb_cols/4)+3)
          colors_loc = c(mycolors[1:k],mycolors[1:k],mycolors[1:k],mycolors[1:(n_groups_loc-3*k)])
          names(colors_loc) = groups_loc
          colors_loc_nodes = sapply(metanetwork$trophicTable[V(g)$name,legend],
                                    function(x) colors_loc[x])
          names(colors_loc_nodes) = V(g)$name
          
          # shape_color_table = data.frame(shape = shapes_loc,color = shapes_loc)
          # rownames(shape_color_table) = groups_loc
          return(GGally::ggnet2(g_Network,mode = mode_loc,
                                color = colors_loc_nodes,
                                size = as.numeric(igraph::V(g)$ab/max(igraph::V(g)$ab)),
                                label = ggnet.config$label, 
                                label.size = ggnet.config$label.size,
                                max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                                edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                                arrow.gap = ggnet.config$arrow.gap,
                                alpha = alpha,
                                edge.alpha = edge_alpha_vec,
                                legend.position = ggnet.config$legend.position,
                                shape = shapes_loc_nodes)
                 )
        }
      }else{
        nodes_focal = alpha_per_node$nodes
        alpha_focal = alpha_per_node$alpha_focal
        alpha_hidden = alpha_per_node$alpha_hidden
        alpha = ifelse(V(g)$name %in% nodes_focal,
                       as.numeric(alpha_focal),as.numeric(alpha_hidden))
        
        edge_alpha = (alpha %*% t(alpha)) * igraph::get.adjacency(g) * ggnet.custom$edge.alpha
        g_Network = network::set.edge.value(g_Network,"edge_alpha",edge_alpha)
        edge_alpha_vec = network::get.edge.attribute(g_Network,'edge_alpha')
        
        nb_cols = length(unique(color_loc))
        #use colors only if reasonible number of colors
        if(nb_cols < 15){
          mycolors = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, ggnet.config$palette))(nb_cols)
          names(mycolors) = unique(color_loc)
          
          return(GGally::ggnet2(g_Network,mode = mode_loc,color = color_loc,
                                size = size_loc,
                                label = ggnet.config$label, 
                                label.size = ggnet.config$label.size,
                                max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                                edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                                arrow.gap = ggnet.config$arrow.gap,
                                alpha = alpha,
                                edge.alpha = edge_alpha_vec,
                                palette = mycolors,
                                legend.position = ggnet.config$legend.position)
          )
        } else{ #use colors and shapes
          #assigning shapes (four types 15,16,17,18) and colors
          
          groups_loc = unique(color_loc)
          n_groups_loc = length(groups_loc)
          
          #shpes
          k = floor(n_groups_loc/4)
          shapes_loc = c(rep(15,k),rep(16,k),rep(17,k),
                         rep(18,n_groups_loc-3*k)) 
          names(shapes_loc) = groups_loc
          shapes_loc_nodes = sapply(metanetwork$trophicTable[V(g)$name,legend],
                                    function(x) shapes_loc[x])
          names(shapes_loc_nodes) = V(g)$name
          #colors
          mycolors = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, ggnet.config$palette))(floor(nb_cols/4)+3)
          colors_loc = c(mycolors[1:k],mycolors[1:k],mycolors[1:k],mycolors[1:(n_groups_loc-3*k)])
          names(colors_loc) = groups_loc
          colors_loc_nodes = sapply(metanetwork$trophicTable[V(g)$name,legend],
                                    function(x) colors_loc[x])
          names(colors_loc_nodes) = V(g)$name
          
          # shape_color_table = data.frame(shape = shapes_loc,color = shapes_loc)
          # rownames(shape_color_table) = groups_loc
          return(GGally::ggnet2(g_Network,mode = mode_loc,
                                color = colors_loc_nodes,
                                size = size_loc,
                                label = ggnet.config$label, 
                                label.size = ggnet.config$label.size,
                                max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                                edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                                arrow.gap = ggnet.config$arrow.gap,
                                alpha = alpha,
                                edge.alpha = edge_alpha_vec,
                                legend.position = ggnet.config$legend.position,
                                shape = shapes_loc_nodes)
          )
        }
      }
    }else{
      if(!(alpha_interactive)){
        #assigning max_size if single abundance value
        size_loc = as.numeric(igraph::V(g)$ab/max(igraph::V(g)$ab))
        if(length(unique(size_loc)) == 1){
          size_loc = ggnet.config$max_size
        }
        
        if(!(is.null(legend))){
          nb_cols = length(unique(color_loc))
          
          #use colors only if reasonible number of colors
          if(nb_cols < 15){
            mycolors = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, ggnet.config$palette))(nb_cols)
            names(mycolors) = unique(color_loc)
            
            return(GGally::ggnet2(g_Network,mode = mode_loc,
                                  color = color_loc,
                                  size = size_loc,
                                  label = ggnet.config$label, 
                                  label.size = ggnet.config$label.size,
                                  max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                                  edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                                  arrow.gap = ggnet.config$arrow.gap,
                                  alpha = ggnet.config$alpha,
                                  edge.alpha = ggnet.config$edge.alpha,
                                  legend.position = ggnet.config$legend.position,
                                  palette = mycolors)
            )
          } else{ #use colors and shapes
            #assigning shapes (four types 15,16,17,18) and colors
            
            groups_loc = unique(color_loc)
            n_groups_loc = length(groups_loc)
            
            #shpes
            k = floor(n_groups_loc/4)
            shapes_loc = c(rep(15,k),rep(16,k),rep(17,k),
                           rep(18,n_groups_loc-3*k)) 
            names(shapes_loc) = groups_loc
            shapes_loc_nodes = sapply(metanetwork$trophicTable[V(g)$name,legend],
                                      function(x) shapes_loc[x])
            names(shapes_loc_nodes) = V(g)$name
            #colors
            mycolors = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, ggnet.config$palette))(floor(nb_cols/4)+3)
            colors_loc = c(mycolors[1:k],mycolors[1:k],mycolors[1:k],mycolors[1:(n_groups_loc-3*k)])
            names(colors_loc) = groups_loc
            colors_loc_nodes = sapply(metanetwork$trophicTable[V(g)$name,legend],
                                      function(x) colors_loc[x])
            names(colors_loc_nodes) = V(g)$name
            
            # shape_color_table = data.frame(shape = shapes_loc,color = shapes_loc)
            # rownames(shape_color_table) = groups_loc
            return(GGally::ggnet2(g_Network,mode = mode_loc,
                                  color = colors_loc_nodes,
                                  size = size_loc,
                                  label = ggnet.config$label, 
                                  label.size = ggnet.config$label.size,
                                  max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                                  edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                                  arrow.gap = ggnet.config$arrow.gap,
                                  alpha = ggnet.config$alpha,
                                  edge.alpha = ggnet.config$edge.alpha,
                                  legend.position = ggnet.config$legend.position,
                                  shape = shapes_loc_nodes)
            )
          }
    } else {
          return(GGally::ggnet2(g_Network,mode = mode_loc,node.color = ggnet.config$default.color,
                        size = size_loc,
                        label = ggnet.config$label, label.size = ggnet.config$label.size,
                        max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                        edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                        arrow.gap = ggnet.config$arrow.gap, alpha = ggnet.config$alpha,
                        edge.alpha = ggnet.config$edge.alpha, 
                        legend.position = ggnet.config$legend.position))
          
        }
      } else {
        if(is.null(legend)){
          message("Alpha interactive mode, enter nodes to enhance")
          message(paste0("Choose nodes among: "),Reduce(paste,V(g)$name))
          nodes_loc = readline(prompt = "nodes (format: name1,name2,...): ")       
          nodes_focal = strsplit(nodes_loc,',')[[1]]
          alpha_focal = readline(prompt = "node alpha of focal nodes: ")
          alpha_hidden = readline(prompt = "node alpha of hidden nodes: ")
          alpha = ifelse(V(g)$name %in% nodes_focal,
                         as.numeric(alpha_focal),as.numeric(alpha_hidden))
          edge_alpha = (alpha %*% t(alpha)) * igraph::get.adjacency(g) * ggnet.config$edge.alpha
          g_Network = network::set.edge.value(g_Network,"edge_alpha",edge_alpha)
          edge_alpha_vec = network::get.edge.attribute(g_Network,'edge_alpha')
          
          #formatting nodes_focal (string) for message
          names_loc = sapply(nodes_focal,function(x) paste0("'",x,"'"))
          names_to_print = gsub(" ",",",Reduce(x = as.character(names_loc),f = paste))
          message(paste0('you can reproduce this result in non-interactive mode by setting: \n',
                         'alpha_per_node = ','list(nodes = c(',names_to_print,'),alpha_focal =',alpha_focal,
                         ',alpha_hidden =',alpha_hidden,')'))
          
          return(GGally::ggnet2(g_Network,mode = mode_loc,
                                size = as.numeric(igraph::V(g)$ab/max(igraph::V(g)$ab)),
                                label = ggnet.config$label, 
                                label.size = ggnet.config$label.size*alpha,
                                max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                                edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                                arrow.gap = ggnet.config$arrow.gap,
                                alpha = alpha,
                                edge.alpha = edge_alpha_vec,
                                legend.position = ggnet.config$legend.position)
          )
          
          
          
        } else{
          message("Alpha interactive mode, enter groups to enhance")
          message(paste0("Choose groups among: "),Reduce(paste,unique(metanetwork$trophicTable[,legend])))
          groups_loc = readline(prompt="groups (format: group1,group2,...): ")
          groups_focal = strsplit(groups_loc,',')[[1]]
          alpha_focal = readline(prompt="node alpha of focal groups: ")
          alpha_hidden = readline(prompt="node alpha of hidden groups: ")
          # edge_alpha_focal = readline(prompt="edge alpha of focal groups: ")
          # edge_alpha_hidden = readline(prompt="edge alpha of hidden groups: ")
          alpha = ifelse(metanetwork$trophicTable[rownames(mode_loc),legend] %in% groups_focal,
                         as.numeric(alpha_focal),as.numeric(alpha_hidden))
          
          edge_alpha = (alpha %*% t(alpha)) * igraph::get.adjacency(g) * ggnet.config$edge.alpha
          g_Network = network::set.edge.value(g_Network,"edge_alpha",edge_alpha)
          edge_alpha_vec = network::get.edge.attribute(g_Network,'edge_alpha')
          
          #formatting groups_focal (string) for message
          names_loc = sapply(groups_focal,function(x) paste0("'",x,"'"))
          names_to_print = gsub(" ",",",Reduce(x = as.character(names_loc),f = paste))
          message(paste0('you can reproduce this result in non-interactive mode by setting: \n',
                         'alpha_per_group = ','list(groups = c(',names_to_print,'),alpha_focal =',alpha_focal,
                         ',alpha_hidden =',alpha_hidden,')'))
          
          nb_cols = length(unique(color_loc))
          mycolors = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, ggnet.config$palette))(nb_cols)
          names(mycolors) = unique(color_loc)
          return(GGally::ggnet2(g_Network,mode = mode_loc,color = color_loc,
                                size = as.numeric(igraph::V(g)$ab/max(igraph::V(g)$ab)),
                                label = ggnet.config$label, 
                                label.size = ggnet.config$label.size,
                                max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                                edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                                arrow.gap = ggnet.config$arrow.gap,
                                alpha = alpha,
                                edge.alpha = edge_alpha_vec,
                                palette = mycolors,
                                legend.position = ggnet.config$legend.position)
          )
        }
      }
    }
  }
}    
     
    
    
    
