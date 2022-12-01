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

#' Default configuration for visNetwork
#'
#' A list with parameters customizing visNetwork visualisation (see visNetwork documentations)
#'
#' @examples
#' # display all default settings
#' visNetwork.default
#'
#' # create a new settings
#' visNetwork.custom = visNetwork.default
#' visNetwork.custom$label.size = 10
#' visNetwork.custom
#' 
#' @export
visNetwork.default = list(
  label = TRUE,
  label.size = 10,
  visEvent = T,
  edge.size = 1)

  # arrow.size = 6, 
  # arrow.gap = 0.015,
  # alpha = 0.8,
  # edge.alpha = 0.5,
  # alpha_diff = 0.8,
  # edge.alpha_diff = 0.8,
  # size.cut = 5,
  # palette = "Set2",
  # default.color = "grey75",
  # legend.position = "bottom")

class(visNetwork.default) = 'metanetwork_config'

#custom layout function for visnetwork
customLayout <- function(graph,g,metanetwork,mode,beta,x_y_range,
                         nrep_ly, layout_metaweb, flip_coords, diff_plot_bool,
                         TL_tsne.config = TL_tsne.default){
  igraphlayout = list(type = "square")
  # g = graph_from_data_frame(graph$x$edges)
  n = igraph::vcount(g)
  TL_loc = igraph::V(g)$TL
  names(TL_loc) = igraph::V(g)$name
  if(mode == 'TL-tsne'){
    if(layout_metaweb == T){
      #check for a layout associated to the metaweb
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
        coords_metaweb = cbind(igraph::V(metaweb)$TL,
                                 igraph::get.vertex.attribute(
                                   metaweb,attr_names[grep(paste0("beta",beta),attr_names)[nrep_ly]]))
        rownames(coords_metaweb) = igraph::V(metaweb)$name
        coords = coords_metaweb[igraph::V(g)$name,]
      }
    }else{
      #check for existing saved layout
      attr_names = igraph::vertex_attr_names(g)
      #check if a layout is attached to the focal network
      if(length(grep(paste0("beta",beta),attr_names))>0){
        coords = cbind(igraph::V(g)$TL,
                       igraph::get.vertex.attribute(g,attr_names[grep(paste0("beta",beta),attr_names)[nrep_ly]]))
        rownames(coords) = igraph::V(g)$name
      }else{
        coords = get_coord_TL_tsne(g = g,TL = igraph::V(g)$TL,beta = beta,
                                            TL_tsne.config = TL_tsne.config)
        rownames(coords) = igraph::V(g)$name
      }
    }
  }
  rownames(coords) = names(TL_loc)
  coords = scale(coords,center = F,scale = T)
  #flip coordinates
  if(flip_coords){
    graph$x$nodes$x = coords[graph$x$nodes$id, 2]*4*x_y_range[2]
    graph$x$nodes$y = -coords[graph$x$nodes$id, 1]*4*x_y_range[1]
    graph$x$nodes$y = graph$x$nodes$y - min(graph$x$nodes$y)
  } else{
    graph$x$nodes$x = coords[graph$x$nodes$id, 1]*8*x_y_range[1]
    graph$x$nodes$y = coords[graph$x$nodes$id, 2]*2*x_y_range[2]
  }
  message(paste0('x_max = ',max(graph$x$nodes$x)))
  message(paste0('y_max = ',max(graph$x$nodes$y)))
  graph$x$options$layout = igraphlayout
  graph %>% visNetwork::visNodes(physics = FALSE) %>% visNetwork::visEdges(smooth = F,shadow = F)
}

#' vismetaNetwork
#'
#' Function that provides network dynamic representation  (using 'visNetwork') from a 
#' 'metanetwork' object with a layout based on a diffusion kernel
#'
#' @param metanetwork object of class metanetwork
#' @param g network (igraph object) to represent, default is metaweb
#' @param beta the diffusion parameter of the diffusion kernel, a positive scalar controlling the 
#' vertical squeezing of the network
#' @param legend resolution for the legend, legend resolution must be a coarser resolution than the resolution of g, default is NULL
#' @param mode mode used for layout, 'TL-tsne' for trophic level t-sne. Default is 'TL-tsne'
#' @param edge_thrs if non-null, a numeric (between 0 and 1) indicating an edge threshold for the representation
#' @param layout_metaweb a boolean indicating whether the layout of the metaweb should be used to represent the network
#' to use metaweb layout = T, you need first to compute metaweb layout for this beta value using \code{attach_layout()}
#' @param nrep_ly If several layouts for this beta value are attached to the metaweb (if \code{layout_metaweb = T}), index of the layout to use, see \code{attach_layout()}
#' @param flip_coords a boolean indicating whether coordinates should be flipped. 
#' In that case, y-axis is the trophic level and x-axis is the layout axis
#' @param x_y_range a two dimension numeric vector, indicating dilatation of x,y axis
#' @param visNetwork.config configuration list for visNetwork representation, default is visNetwork.default
#' @param TL_tsne.config configuration list for mode 'TL-tsne', default is TL_tsne.default
#' @param diff_plot_bool boolean, do not edit by hand
#' 
#' @return object of class 'visNetwork', dynamic representation of the current network
#'
#' @examples
#' library(metanetwork)
#' library(igraph)
#' data("meta_angola")
#' ## Return htmlwidget
#' # on angola dataset
#' meta_angola = attach_layout(meta_angola, beta = 0.05)
#' vismetaNetwork(meta_angola, beta = 0.05)
#' 
#'
#' @export
vismetaNetwork <- function(metanetwork,g = NULL,beta = 0.1,
                           legend = NULL,mode = 'TL-tsne',
                           edge_thrs = NULL,
                           layout_metaweb = FALSE,nrep_ly = 1,
                           flip_coords = FALSE,
                           diff_plot_bool = FALSE,
                           x_y_range = c(100,100),
                           visNetwork.config = visNetwork.default,
                           TL_tsne.config = TL_tsne.default){
  `%>%` <- magrittr::`%>%`
  if(is.null(g)){
    g = metanetwork$metaweb
  }
  #edge weight threshold
  if(!(is.null(edge_thrs))){
    g = igraph::delete.edges(graph = g,edges =  which(igraph::E(g)$weight < edge_thrs))
  }
  if(!(diff_plot_bool)){
    if(is.null(legend)){
      #get resources and consumers
      resources_list = sapply(igraph::V(g)$name,
                              function(x) igraph::neighbors(g,x,mode = 'in')$name)
      resources_list[which(lapply(resources_list,length)==0)] = ""
      consumers_list = sapply(igraph::V(g)$name,
                              function(x) igraph::neighbors(g,x,mode = 'out')$name)
      consumers_list[which(lapply(consumers_list,length)==0)] = "" 
      
      nodes = data.frame(id = igraph::V(g)$name, value = igraph::V(g)$ab,
                         font.size = visNetwork.config$label.size,
                         resources = paste("\n",unname(unlist(lapply(
                           resources_list,
                           paste,collapse = ', '))),"\n"),
                         consumers = paste("\n",unname(unlist(lapply(
                           consumers_list,
                           paste,collapse = ', '))),"\n"),
                         TL = paste('\n Trophic level',round(igraph::V(g)$TL,digits = 3)),
                         width = ifelse(is.null(igraph::E(g)$weight),1,igraph::E(g)$weight),
                         stringsAsFactors = F)
      #removing labels
      if(!(visNetwork.config$label)){
        nodes = dplyr::mutate(nodes,label = NA)
      }else{
        nodes = dplyr::mutate(nodes,label = igraph::V(g)$name)
      }
      edges = data.frame(igraph::get.edgelist(g),
                         arrows = c("to"),
                         width = igraph::E(g)$weight*visNetwork.config$edge.size)
      colnames(edges) = c('from','to','arrows','width')
      
      network_loc = visNetwork::visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
        #visNodes(shape = "square") %>%  # square for all nodes
        visNetwork::visNodes(shapeProperties = list(useBorderWithImage = TRUE)) %>% 
        visNetwork::visOptions(highlightNearest = TRUE) %>%
        customLayout(metanetwork = metanetwork,g = g,mode = mode,beta = beta,x_y_range = x_y_range,
                     nrep_ly = nrep_ly,layout_metaweb = layout_metaweb,flip_coords = flip_coords) %>%
        visNetwork::visLegend(width = 0.1, position = "right")
      if(!(visNetwork.config$visEvent)){
        return(network_loc)
      } else{
        network_loc = visNetwork::visEvents(network_loc,selectNode = "function(properties) {
      alert(this.body.data.nodes.get(properties.nodes[0]).id +
              ' consumes: ' + this.body.data.nodes.get(properties.nodes[0]).resources +
              this.body.data.nodes.get(properties.nodes[0]).id +
              ' is consumed by: ' + this.body.data.nodes.get(properties.nodes[0]).consumers +
              this.body.data.nodes.get(properties.nodes[0]).TL);}")
        return(network_loc)
      }
    } else{
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
      if(!(which(colnames(metanetwork$trophicTable) == res_local) <
           which(colnames(metanetwork$trophicTable) == legend))){
        stop("legend must be a coarser resolution than resolution of the current network")
      }else{
        #get resources and consumers
        resources_list = sapply(igraph::V(g)$name,
                                function(x) igraph::neighbors(g,x,mode = 'in')$name)
        resources_list[which(lapply(resources_list,length)==0)] = ""
        consumers_list = sapply(igraph::V(g)$name,
                                function(x) igraph::neighbors(g,x,mode = 'out')$name)
        consumers_list[which(lapply(consumers_list,length)==0)] = "" 
        
        nodes = data.frame(id = igraph::V(g)$name, value = igraph::V(g)$ab,
                           font.size = visNetwork.config$label.size,
                           group = metanetwork$trophicTable[igraph::V(g)$name,legend],
                           resources = paste("\n",unname(unlist(lapply(
                             resources_list,
                             paste,collapse = ', '))),"\n"),
                           consumers = paste("\n",unname(unlist(lapply(
                             consumers_list,
                             paste,collapse = ', '))),"\n"),
                           TL = paste('\n Trophic level',round(V(g)$TL,digits = 3)),
                           stringsAsFactors = F)
        
        #removing labels
        if(!(visNetwork.config$label)){
          nodes = dplyr::mutate(nodes,label = NA)
        }else{
          nodes = dplyr::mutate(nodes,label = igraph::V(g)$name)
        }
        
        edges = data.frame(igraph::get.edgelist(g),
                           # arrows
                           arrows = c("to"),
                           width = igraph::E(g)$weight*visNetwork.config$edge.size)
        colnames(edges) = c('from','to','arrows','width')
        
        network_loc = visNetwork::visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
          #visNodes(shape = "square") %>%  # square for all nodes
          visNetwork::visNodes(shapeProperties = list(useBorderWithImage = TRUE)) %>% # images
          visNetwork::visOptions(highlightNearest = TRUE) %>%
          customLayout(metanetwork = metanetwork,g = g,mode = mode,beta = beta,x_y_range = x_y_range,
                       nrep_ly = nrep_ly,layout_metaweb = layout_metaweb,flip_coords = flip_coords,
                       TL_tsne.config = TL_tsne.config) %>%
          visNetwork::visLegend(width = 0.1, position = "right") 
        if(!(visNetwork.config$visEvent)){
          return(network_loc)
        } else{
          network_loc = visNetwork::visEvents(network_loc,selectNode = "function(properties) {
      alert(this.body.data.nodes.get(properties.nodes[0]).id +
              ' consumes: ' + this.body.data.nodes.get(properties.nodes[0]).resources +
              this.body.data.nodes.get(properties.nodes[0]).id +
              ' is consumed by: ' + this.body.data.nodes.get(properties.nodes[0]).consumers +
              this.body.data.nodes.get(properties.nodes[0]).TL);}")
          return(network_loc)
        }
        return(network_loc)
      }
    }
  }else{ #diff_plot_bool
    #metaweb of the two networks (meta_diff)
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
    
    nodes = data.frame(id = igraph::V(g)$name, value = igraph::V(g)$ab,
                       group = color_loc,
                       font.size = visNetwork.config$label.size,
                       # resources = paste("\n",unname(unlist(lapply(
                       #   resources_list,
                       #   paste,collapse = ', '))),"\n"),
                       # consumers = paste("\n",unname(unlist(lapply(
                       #   consumers_list,
                       #   paste,collapse = ', '))),"\n"),
                       # TL = paste('\n Trophic level',round(igraph::V(g)$TL,digits = 3)),
                       width = ifelse(is.null(igraph::E(g)$weight),1,igraph::E(g)$weight),
                       stringsAsFactors = F)
    
    #removing labels
    if(!(visNetwork.config$label)){
      nodes = dplyr::mutate(nodes,label = NA)
    }else{
      nodes = dplyr::mutate(nodes,label = igraph::V(g)$name)
    }
    edges = data.frame(igraph::get.edgelist(g),
                       arrows = c("to"),
                       width = igraph::E(g)$weight*visNetwork.config$edge.size,color = edge_color_loc)
    colnames(edges) = c('from','to','arrows','width','color')
    
    network_loc = visNetwork::visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
      #visNodes(shape = "square") %>%  # square for all nodes
      visNetwork::visNodes(shapeProperties = list(useBorderWithImage = TRUE)) %>% # images
      visNetwork::visOptions(highlightNearest = TRUE) %>%
      customLayout(metanetwork = metanetwork,g = g,mode = mode,beta = beta,x_y_range = x_y_range,
                   nrep_ly = nrep_ly,layout_metaweb = layout_metaweb,flip_coords = flip_coords,
                   TL_tsne.config = TL_tsne.config) %>%
      visNetwork::visLegend(width = 0.1, position = "right",main = 
                              list(text = "Legend",
                                   style = "font-family:Comic Sans MS;color:#6d799c;font-size:12px;text-align:center;")) %>%
      visNetwork::visGroups(groupname = "shared", color = "grey50") %>%
      visNetwork::visGroups(groupname = "more ab in g1", color = "#a1d99b") %>%
      visNetwork::visGroups(groupname = "more ab in g2", color = "#fc9272") %>%
      visNetwork::visGroups(groupname = "only pres in g1", color = "#31a354") %>%
      visNetwork::visGroups(groupname = "only pres in g2", color = "#de2d26") 
      return(network_loc)
  }
}
