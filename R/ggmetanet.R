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
  legend.position = "bottom",
  legend.big.nets = T,
  img_PATH = NULL)

class(ggnet.default) = 'metanetwork_config'



#' ggmetanet
#'
#' Function that provides network static representation  (using 'ggnet') from a 
#' 'metanetwork' object using 'TL-tsne' or 'group-TL-tsne' layout.
#' 
#' At each call of the function with 'TL-tsne' layout, it computes a layout for the current beta value.
#' If a layout is already attached to the current network, it uses directly this layout (without computing). 
#' This function provides many static visualisation tools:
#'  - customising ggnet parameters wrapped in `ggnet.config`
#'  - legending using the trophicTable
#'  - playing on group transparency (alpha)
#'  - using the metaweb layout
#'  - building a legend for large networks.
#'
#' @param metanetwork object of class metanetwork
#' @param g network (igraph object) to represent, default is metaweb
#' @param beta the diffusion parameter of the diffusion kernel, a positive scalar controlling the 
#' vertical squeezing of the network
#' @param legend resolution for the legend, legend resolution must be a coarser resolution than the resolution of g, default is NULL
#' @param mode mode used for layout, 'TL-tsne' or 'group-TL-tsne' Default is 'TL-tsne'.
#' This argument can also be a two-column matrix for custom layout.
#' @param edge_thrs if non-null, a numeric (between 0 and 1) indicating an edge threshold for the representation
#' @param layout_metaweb a boolean indicating whether the layout of the metaweb should be used to represent the network
#' to use metaweb layout = TRUE, you need first to compute metaweb layout for this beta value using \code{attach_layout()}
#' @param nrep_ly If several layouts for this beta value are attached to the metaweb (if \code{layout_metaweb = T}), index of the layout to use, see \code{attach_layout()}
#' @param flip_coords a boolean indicating whether coordinates should be flipped. 
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
#' @return an object of class \code{ggplot}, the current network representation
#' 
#' @seealso [attach_layout()],[ggnet.default]
#' 
#' @examples
#' library(metanetwork)
#' library(igraph)
#' 
#' #lattice example
#' g = make_lattice(dim = 2,length = 4,directed = TRUE)
#' #building metanetwork and computing trophic levels
#' meta0 = build_metanet(g) 
#' meta0 = compute_TL(meta0)
#' ggmetanet(meta0)
#' #storing layout
#' meta0 = attach_layout(meta0)
#' ggmetanet(meta0)
#' 
#'#custom ggnet parameters
#'ggnet.custom = ggnet.default
#'ggnet.custom$label = TRUE
#'ggnet.custom$edge.alpha = 0.5
#'ggnet.custom$alpha = 0.7
#'ggnet.custom$arrow.size = 1
#'ggnet.custom$max_size = 12
#'
#'# using pre-computed layout and custom ggnet parametersfor vertebrates metaweb
#'data("meta_vrtb")
#'#custom ggnet parameters
#'ggnet.custom = ggnet.default
#'ggnet.custom$label = TRUE
#'ggnet.custom$edge.alpha = 0.5
#'ggnet.custom$alpha = 0.7
#'ggnet.custom$arrow.size = 1
#'ggnet.custom$max_size = 12
#' #at SBM group level
#'beta = 0.005
#'ggmetanet(meta_vrtb,g = meta_vrtb$metaweb_group,flip_coords = TRUE,
#'          beta = beta,legend = "group",
#'          ggnet.config = ggnet.custom,edge_thrs = 0.1)
#'
#' @importFrom ggplot2 aes element_blank element_line element_rect xlim ylim
#' @importFrom intergraph asNetwork
#' @importFrom graphics image
#' @importFrom rlang .data
#' @import sna
#' 
#' @export
#
#
#
ggmetanet <- function(metanetwork,g = NULL,beta = 0.1,
                           legend = NULL,mode = 'TL-tsne',
                           edge_thrs = NULL,
                           layout_metaweb = FALSE,nrep_ly = 1,
                           flip_coords = FALSE,
                           diff_plot_bool = FALSE,
                           alpha_per_group = NULL,alpha_per_node = NULL,
                           alpha_interactive = FALSE,
                           ggnet.config = ggnet.default,
                           TL_tsne.config = TL_tsne.default){

  if(!is.metanetwork(metanetwork)){
    stop("metanetwork is an object of class metanetwork, see build_metanetwork")
  }
  if(sum(class(mode) == "character")>0){
    message(paste0("mode is ",mode))
    if(!(mode %in% c('TL-tsne','fr','kk','circle','group-TL-tsne'))){
      stop("mode must be one of: \n
         'TL-tsne','fr','kk','circle', 'group-TL-tsne' \n
           or, alternatively, a 2-column matrix with coordinates" )
    }
  } else if (sum(class(mode) %in% c("matrix","data.frame"))>0){
    message(paste0("mode is custom"))
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
      if(!(which(colnames(metanetwork$trophicTable) == res_local) <=
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
          order_loc = order(color_loc)
          #permuting in the order of colors for legend
          g = igraph::permute(g,order(order_loc))
          #ordering colors for legend
          color_loc = color_loc[order(color_loc)]
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
    color_loc =  ifelse(igraph::V(g)$ab>0,'more ab in g1','more ab in g2')
    color_loc[which(igraph::V(g)$ab > 2)] = 'only pres in g1'
    color_loc[which(igraph::V(g)$ab < -2)] = 'only pres in g2'
    color_loc[which(igraph::V(g)$ab ==0)] = 'shared'
    edge_color_loc = ifelse(igraph::E(g)$weight_col>0,'#d55e00','#56b4e9')
    edge_color_loc[which(igraph::E(g)$weight_col == 0)] = 'black'
    edge_color_loc[which(igraph::E(g)$weight_col > 2)] = '#8c2a00'
    edge_color_loc[which(igraph::E(g)$weight_col < -2)] = '#1a2d4d'
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

  if(sum(class(mode) == "character")>0){
  if(mode %in% c("fr","kk","circle")){
    if(mode == "fr"){
      mode_loc = "fruchtermanreingold"
    } else if(mode == "kk"){
      mode_loc = "kamadakawai"
    } else{
      mode_loc = mode
    }
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
        #check if nrep_ly indicates computed layout
        if(nrep_ly>length(grep(paste0("beta",beta),attr_names))){
          stop(paste0("nrep_ly argument is:",nrep_ly,
                      "whereas number of computed layout for this beta value are:",
                      length(grep(paste0("beta",beta),attr_names))," .You must decrease nrep_ly or
                      compute more layouts for this beta value"))
        }
        mode_loc = cbind(igraph::V(g)$TL,
                         igraph::get.vertex.attribute(g,attr_names[grep(paste0("beta",beta),attr_names)[nrep_ly]]))
        rownames(mode_loc) = igraph::V(g)$name
      }else{
        mode_loc = get_coord_TL_tsne(g = g,TL = igraph::V(g)$TL,beta = beta,
                                              TL_tsne.config = TL_tsne.config)
        rownames(mode_loc) = igraph::V(g)$name
      }
    }
  }
  if(mode == "group-TL-tsne"){
    attr_names = igraph::vertex_attr_names(g)
    if(length(grep(paste0("group_layout_x_beta",beta),attr_names)) == 0){
      stop("to use 'group-TL-tsne', you need to attach a layout for the desired beta value
           and resolution, see attach_layout function with 'group-TL-tsne' mode")
    }else{
      mode_loc = cbind(igraph::get.vertex.attribute(g,attr_names[
        match(paste0("group_layout_x_beta",beta),attr_names)]),
                       igraph::get.vertex.attribute(g,attr_names[
                         match(paste0("group_layout_y_beta",beta),attr_names)]))
      rownames(mode_loc) = igraph::V(g)$name
    }
  }
  }
  
  if(sum(class(mode) %in% c('matrix','data.frame'))>0){
    #custom layout
    if(!(prod(dim(mode) == c(igraph::vcount(g),2)))){
      stop("mode must be of dimension: node number of g * 2")
    } else{
      message("mode is custom")
      mode_loc = mode[V(g)$name,]
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
        net = GGally::ggnet2(g_Network,mode = mode_loc,color = "color",edge.color = edge_color_loc,
                             size = size_loc,edge.size = "weight",
                             label = ggnet.config$label, label.size = ggnet.config$label.size,
                             max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                             arrow.size = ggnet.config$arrow.size,
                             arrow.gap = ggnet.config$arrow.gap,
                             alpha = alpha,
                             edge.alpha = edge_alpha_vec*ifelse(edge_color_loc == 'black',ggnet.config$edge.alpha,ggnet.config$edge.alpha_diff),
                             legend.position = ggnet.config$legend.position,
                             palette = c("only pres in g1" = "#8c2a00", "more ab in g1" = "#d55e00","more ab in g2" = "#56b4e9",
                                         "only pres in g2" = "#1a2d4d","shared" = "grey75")) +
                             ggplot2::theme(legend.box = "vertical")
        return(net)
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
        net = GGally::ggnet2(g_Network,mode = mode_loc,color = "color",edge.color = edge_color_loc,
                             size = as.numeric(igraph::V(g)$ab/max(igraph::V(g)$ab)),
                             label = ggnet.config$label, 
                             label.size = ggnet.config$label.size*alpha,
                             max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                             edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                             arrow.gap = ggnet.config$arrow.gap,
                             alpha = alpha,
                             edge.alpha = edge_alpha_vec,
                             legend.position = ggnet.config$legend.position,
                             palette = c("only pres in g1" = "#8c2a00", "more ab in g1" = "#d55e00","more ab in g2" = "#56b4e9",
                                         "only pres in g2" = "#1a2d4d","shared" = "grey75")) +
          ggplot2::theme(legend.box = "vertical")
        return(net)
      }
    } else{
      #check level number of attribute "color", bug in ggnet2 with palette if single level
      if(length(unique( network::get.vertex.attribute(g_Network,"color"))) == 1){
        net = GGally::ggnet2(g_Network,mode = mode_loc,color = "color",edge.color = edge_color_loc,
                             size = size_loc,edge.size = "weight",
                             label = ggnet.config$label, label.size = ggnet.config$label.size,
                             max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                             arrow.size = ggnet.config$arrow.size,
                             arrow.gap = ggnet.config$arrow.gap,
                             alpha = alpha,
                             edge.alpha = ifelse(edge_color_loc == 'black',ggnet.config$edge.alpha,ggnet.config$edge.alpha_diff),
                             legend.position = ggnet.config$legend.position) +
          ggplot2::theme(legend.box = "vertical")
      }else{
        net = GGally::ggnet2(g_Network,mode = mode_loc,color = "color",edge.color = edge_color_loc,
                             size = size_loc,edge.size = "weight",
                             label = ggnet.config$label, label.size = ggnet.config$label.size,
                             max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                             arrow.size = ggnet.config$arrow.size,
                             arrow.gap = ggnet.config$arrow.gap,
                             alpha = alpha,
                             edge.alpha = ifelse(edge_color_loc == 'black',ggnet.config$edge.alpha,ggnet.config$edge.alpha_diff),
                             legend.position = ggnet.config$legend.position,
                             palette = c("only pres in g1" = "#8c2a00", "more ab in g1" = "#d55e00","more ab in g2" = "#56b4e9",
                                         "only pres in g2" = "#1a2d4d","shared" = "grey75")) +
          ggplot2::theme(legend.box = "vertical")
      }
      return(net)
    }
  }else{
    #assigning max_size if single abundance value
    size_loc = as.numeric(igraph::V(g)$ab/max(igraph::V(g)$ab))
    if(length(unique(size_loc)) == 1){
      size_loc = ggnet.config$max_size
    }
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
          #make label size proportional to alpha
          
          
          net = GGally::ggnet2(g_Network,mode = mode_loc,color = color_loc,
                         size = as.numeric(igraph::V(g)$ab/max(igraph::V(g)$ab)),
                         label = ggnet.config$label, 
                         label.size = ggnet.config$label.size*alpha,
                         max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                         edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                         arrow.gap = ggnet.config$arrow.gap,
                         alpha = alpha,
                         edge.alpha = edge_alpha_vec,
                         palette = mycolors,
                         legend.position = ggnet.config$legend.position) +
            ggplot2::theme(legend.box = "vertical")
          return(net)
        } else{ #use colors and shapes
          shapes_colors = assign_shapes_colors(color_loc,metanetwork,g,legend,
                                               beta,nrep_ly,mode_loc,flip_coords,ggnet.config)
          colors_loc_nodes = shapes_colors$colors
          shapes_loc_nodes = shapes_colors$shapes
          
          # shape_color_table = data.frame(shape = shapes_loc,color = shapes_loc)
          # rownames(shape_color_table) = groups_loc
          net = GGally::ggnet2(g_Network,mode = mode_loc,
                         color = colors_loc_nodes,
                         size = as.numeric(igraph::V(g)$ab/max(igraph::V(g)$ab)),
                         label = ggnet.config$label, 
                         label.size = ggnet.config$label.size*alpha,
                         max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                         edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                         arrow.gap = ggnet.config$arrow.gap,
                         alpha = alpha,
                         edge.alpha = edge_alpha_vec,
                         legend.position = ggnet.config$legend.position,
                         shape = shapes_loc_nodes)
          return(net)
        }
      }else{
        nodes_focal = alpha_per_node$nodes
        alpha_focal = alpha_per_node$alpha_focal
        alpha_hidden = alpha_per_node$alpha_hidden
        alpha = ifelse(V(g)$name %in% nodes_focal,
                       as.numeric(alpha_focal),as.numeric(alpha_hidden))
        
        edge_alpha = (alpha %*% t(alpha)) * igraph::get.adjacency(g) * ggnet.config$edge.alpha
        g_Network = network::set.edge.value(g_Network,"edge_alpha",edge_alpha)
        edge_alpha_vec = network::get.edge.attribute(g_Network,'edge_alpha')
        
        nb_cols = length(unique(color_loc))
        #use colors only if reasonible number of colors
        if(nb_cols < 15){
          mycolors = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, ggnet.config$palette))(nb_cols)
          names(mycolors) = unique(color_loc)
          net = GGally::ggnet2(g_Network,mode = mode_loc,color = color_loc,
                         size = size_loc,
                         label = ggnet.config$label, 
                         label.size = ggnet.config$label.size*alpha,
                         max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                         edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                         arrow.gap = ggnet.config$arrow.gap,
                         alpha = alpha,
                         edge.alpha = edge_alpha_vec,
                         palette = mycolors,
                         legend.position = ggnet.config$legend.position) +
            ggplot2::theme(legend.box = "vertical")
          return(net)
        } else{ #use colors and shapes
          shapes_colors = assign_shapes_colors(color_loc,metanetwork,g,legend,beta,
                                               nrep_ly,mode_loc,flip_coords,ggnet.config)
          colors_loc_nodes = shapes_colors$colors
          shapes_loc_nodes = shapes_colors$shapes
          
          net = GGally::ggnet2(g_Network,mode = mode_loc,
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
                               shape = shapes_loc_nodes)  +
            ggplot2::theme(legend.box = "vertical")
          # shape_color_table = data.frame(shape = shapes_loc,color = shapes_loc)
          # rownames(shape_color_table) = groups_loc
          return(net)
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
            net = GGally::ggnet2(g_Network,mode = mode_loc,
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
                                 palette = mycolors)  +
              ggplot2::theme(legend.box = "vertical")
            return(net)
          } else{ #use colors and shapes
            shapes_colors = assign_shapes_colors(color_loc,metanetwork,g,legend,beta,nrep_ly,
                                                 mode_loc,flip_coords,ggnet.config)
            colors_loc_nodes = shapes_colors$colors
            shapes_loc_nodes = shapes_colors$shapes
            # shape_color_table = data.frame(shape = shapes_loc,color = shapes_loc)
            # rownames(shape_color_table) = groups_loc
            net = GGally::ggnet2(g_Network,mode = mode_loc,
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
                                 shape = shapes_loc_nodes)  +
              ggplot2::theme(legend.box = "vertical")
            return(net)
          }
    } else {
        net = GGally::ggnet2(g_Network,mode = mode_loc,node.color = ggnet.config$default.color,
                             size = size_loc,
                             label = ggnet.config$label, label.size = ggnet.config$label.size,
                             max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                             edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                             arrow.gap = ggnet.config$arrow.gap, alpha = ggnet.config$alpha,
                             edge.alpha = ggnet.config$edge.alpha, 
                             legend.position = ggnet.config$legend.position)  +
          ggplot2::theme(legend.box = "vertical")
      return(net)
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
          net = GGally::ggnet2(g_Network,mode = mode_loc,
                               size = as.numeric(igraph::V(g)$ab/max(igraph::V(g)$ab)),
                               label = ggnet.config$label, 
                               label.size = ggnet.config$label.size*alpha,
                               max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                               edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                               arrow.gap = ggnet.config$arrow.gap,
                               alpha = alpha,
                               edge.alpha = edge_alpha_vec,
                               legend.position = ggnet.config$legend.position)  +
            ggplot2::theme(legend.box = "vertical")
          return(net)
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
          net = GGally::ggnet2(g_Network,mode = mode_loc,color = color_loc,
                               size = as.numeric(igraph::V(g)$ab/max(igraph::V(g)$ab)),
                               label = ggnet.config$label, 
                               label.size = ggnet.config$label.size,
                               max_size = ggnet.config$max_size,size.cut = ggnet.config$size.cut,
                               edge.size = "weight",arrow.size = ggnet.config$arrow.size,
                               arrow.gap = ggnet.config$arrow.gap,
                               alpha = alpha,
                               edge.alpha = edge_alpha_vec,
                               palette = mycolors,
                               legend.position = ggnet.config$legend.position)  +
            ggplot2::theme(legend.box = "vertical")
          return(net)
        }
      }
    }
  }
}    
     
#function to build a legend mixing shapes and colors    
assign_shapes_colors <- function(color_loc,metanetwork,g,legend,
                                 beta,nrep_ly,mode_loc,flip_coords,ggnet.config){
  groups_loc = unique(color_loc)
  #re-order with group trophic levels for shapes
  networks = extract_networks(metanetwork)
  agg_net_loc = networks[[
    which(sapply(networks,
                 function(x) c(x$res == legend,x$name == g$name)) %>%
            colSums() == 2)]]
  if(!(is.null(V(agg_net_loc)$TL))){
    groups_loc = groups_loc[order(V(agg_net_loc)$TL)]
  }
  n_groups_loc = length(groups_loc)
  nb_cols = length(unique(color_loc))
  
  #shapes: use of 15,16,17,18 and 25
  k = floor(n_groups_loc/5)
  shapes_loc = c(rep(15,k),rep(16,k),rep(17,k),
                 rep(18,k),rep(20,n_groups_loc-4*k)) 
  names(shapes_loc) = groups_loc
  shapes_pool = unique(shapes_loc)
  #generate colors
  mycolors = grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(8, ggnet.config$palette))(floor(nb_cols/5)+4)
  
  #assign colors according to TL-tsne coordinate (if existing) of the aggregated network
  ind_attr_loc = grep(paste0("layout_beta",beta), igraph::vertex_attr_names(agg_net_loc))
  if(length(ind_attr_loc) == 1){
    lay_agg_loc = igraph::get.vertex.attribute(agg_net_loc)[[ind_attr_loc]]
  } else if(length(ind_attr_loc) == 0){
    if(flip_coords){
      lay_agg_loc = mode_loc[,1]
    } else{
      lay_agg_loc = mode_loc[,2]
    }
    #repet layout case 
  } else if(length(ind_attr_loc>1)){
    ind_attr_loc = ind_attr_loc[nrep_ly]
    lay_agg_loc = igraph::get.vertex.attribute(agg_net_loc)[[ind_attr_loc]]
  }
  
  names(lay_agg_loc) = igraph::V(agg_net_loc)$name
  lay_agg_loc = lay_agg_loc[groups_loc]
  
  #filling colors
  comp = 0
  for(j in 1:length(shapes_pool)){
    # print(comp)
    # print(j)
    if(j == 1){
      order_colors = order(lay_agg_loc[which(shapes_loc == shapes_pool[j])])
      colors_loc = mycolors[1:length(order_colors)]
      comp = comp + length(order_colors)
    }else{
      order_colors_loc = order(lay_agg_loc[which(shapes_loc == shapes_pool[j])])
      colors_loc = c(colors_loc,
                     mycolors[1:length(order_colors_loc)])
      order_colors = c(order_colors,order_colors_loc + comp)
      comp = comp + length(order_colors_loc)
    }
  }
  names(colors_loc) = groups_loc[order_colors]
  colors_loc = colors_loc[groups_loc]
  
  # colors_loc = c(mycolors[1:k],mycolors[1:k],mycolors[1:k],mycolors[1:k],mycolors[1:(n_groups_loc-4*k)])
  # names(colors_loc) = groups_loc
  
  #legend identical to network resolution
  if(g$res == legend){
    shapes_loc_nodes = sapply(igraph::V(g)$name,
                              function(x) shapes_loc[x])
    names(shapes_loc_nodes) = igraph::V(g)$name
    #colors
    colors_loc_nodes = sapply(igraph::V(g)$name,
                              function(x) colors_loc[x])
    names(colors_loc_nodes) = igraph::V(g)$name
    
  }else{ #must be at the original resolution 
    shapes_loc_nodes = sapply(metanetwork$trophicTable[igraph::V(g)$name,legend],
                              function(x) shapes_loc[x])
    names(shapes_loc_nodes) = igraph::V(g)$name
    #colors
    colors_loc_nodes = sapply(metanetwork$trophicTable[igraph::V(g)$name,legend],
                              function(x) colors_loc[x])
    names(colors_loc_nodes) = igraph::V(g)$name
  }
  message("too many groups, assigning colors and shapes")
  # print(shapes_loc_nodes)
  # print(colors_loc_nodes)
  #plot legend if resolution of the current network is different than legend
  if(g$res != legend){
    if(ggnet.config$legend.big.nets){
      p  = plot_legend(shapes_loc,colors_loc,metanetwork,
                       legend,order_colors,ggnet.config)
      plot(p) 
    }
  }
  return(list(shapes = shapes_loc_nodes,
              colors = colors_loc_nodes))
}


plot_legend <- function(shapes_loc,colors_loc,metanetwork,
                        legend,order_colors,ggnet.config){
  message("plotting legend")
  shapes_colors_df = cbind(shapes_loc,colors_loc)
  colnames(shapes_colors_df) = c('shape','color')
  shapes_colors_df = as.data.frame(shapes_colors_df)
  shapes_colors_df =  cbind(shapes_colors_df,
                            group = rownames(shapes_colors_df),
                            x = rep(0,nrow(shapes_colors_df)),
                            y = rep(0,nrow(shapes_colors_df)))
  
  shape_pool = unique(shapes_colors_df$shape)
  comp = rep(0,5)
  #re-order with colors
  shapes_colors_df = shapes_colors_df[names(colors_loc)[order_colors],]
  
  for(k in 1:nrow(shapes_colors_df)){
    if(shapes_colors_df$shape[k] == 15){
      comp[1] = comp[1] + 1
      shapes_colors_df[k,"x"] = 1
      shapes_colors_df[k,"y"] = comp[1]
    } else if(shapes_colors_df$shape[k] == 16){
      comp[2] = comp[2] + 1
      shapes_colors_df[k,"x"] = 6
      shapes_colors_df[k,"y"] = comp[2]
    } else if(shapes_colors_df$shape[k] == 17){
      comp[3] = comp[3] + 1
      shapes_colors_df[k,"x"] = 11
      shapes_colors_df[k,"y"] = comp[3]
    } else if(shapes_colors_df$shape[k] == 18){
      comp[4] = comp[4] + 1
      shapes_colors_df[k,"x"] = 16
      shapes_colors_df[k,"y"] = comp[4]
    } else if(shapes_colors_df$shape[k] == 20){
      comp[5] = comp[5] + 1
      shapes_colors_df[k,"x"] = 20
      shapes_colors_df[k,"y"] = comp[5]
    }
  }
  
  if(is.null(ggnet.config$img_PATH)){
    legend_plot = ggplot2::ggplot(shapes_colors_df, aes(.data$x, .data$y, label = .data$group)) +
      ggplot2::geom_point(aes(shape = .data$shape, color = .data$color)
                          , size = 5) +
      ggplot2::scale_shape_manual(values = unique(as.numeric(shapes_colors_df$shape))) +
      ggplot2::scale_color_manual(breaks = shapes_colors_df[shapes_colors_df$shape == 20,]$color,
                                  values = shapes_colors_df[shapes_colors_df$shape == 20,]$color) + 
      ggplot2::geom_text(hjust=0, vjust=2) + 
      ggplot2::theme(legend.position = "none",
                     panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     # Modifier le trait des axes
                     axis.line = element_line(colour = "white"),
                     plot.background = element_rect(fill = "white"),
                     panel.grid = element_blank(),
                     panel.background = element_blank(),
                     panel.grid.major.x = element_blank(),
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank()) + 
      ggplot2::ggtitle("legend") + xlim(c(0,24)) + ylim(c(0,10)) 
  } else{
    message("building a legend with images")
    message(paste0("image PATH: "),ggnet.config$img_PATH)
    img_PATH = ggnet.config$img_PATH
    img = list.files(path = img_PATH,
                     pattern="png", full.names=TRUE)
    #reorder
    img = img[sapply(shapes_colors_df$group,function(x) grep(paste0("/",x,"_"),img))]
    shapes_colors_df = cbind(shapes_colors_df,image = img)
    
    legend_plot = ggplot2::ggplot(shapes_colors_df, aes(.data$x, .data$y, label = .data$group)) +
      ggplot2::geom_point(aes(shape = .data$shape, color = .data$color)
                          , size = 5) +
      ggplot2::scale_shape_manual(values = unique(as.numeric(shapes_colors_df$shape))) +
      ggplot2::scale_color_manual(breaks = shapes_colors_df[shapes_colors_df$shape == 20,]$color,
                                  values = shapes_colors_df[shapes_colors_df$shape == 20,]$color) + 
      ggplot2::geom_text(hjust=0, vjust=2) + 
      ggplot2::theme(legend.position = "none",
                     panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     # Modifier le trait des axes
                     axis.line = element_line(colour = "white"),
                     plot.background = element_rect(fill = "white"),
                     panel.grid = element_blank(),
                     panel.background = element_blank(),
                     panel.grid.major.x = element_blank(),
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank()) + 
      ggplot2::ggtitle("legend") + xlim(c(0,24)) + ylim(c(0,10)) + 
      ggimage::geom_image(aes(x = .data$x + 2,y=.data$y,image = .data$image))
  }
  return(legend_plot)
}


