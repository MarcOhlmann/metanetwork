#function to build a legend mixing shapes and colors    
assign_shapes_colors <- function(color_loc,metanetwork,g,legend,ggnet.config,
                                 beta,nrep_ly,mode_loc,flip_coords){
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
    p  = plot_legend(shapes_loc,colors_loc,metanetwork,
                     legend,order_colors)
    plot(p)
  }
  return(list(shapes = shapes_loc_nodes,
              colors = colors_loc_nodes))
}


plot_legend <- function(shapes_loc,colors_loc,metanetwork,
                        legend,order_colors){
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
    legend_plot = ggplot2::ggplot(shapes_colors_df, aes(x, y, label = group)) +
      ggplot2::geom_point(aes(shape = shape, color = color)
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
    img = list.files(path = img_PATH,
                     pattern="png", full.names=TRUE)
    #reorder
    img = img[sapply(shapes_colors_df$group,function(x) grep(paste0("/",x,"_"),img))]
    shapes_colors_df = cbind(shapes_colors_df,image = img)
    
    legend_plot = ggplot2::ggplot(shapes_colors_df, aes(x, y, label = group)) +
      ggplot2::geom_point(aes(shape = shape, color = color)
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
      ggplot2::ggtitle("legend") + xlim(c(0,24)) + ylim(c(0,10)) + geom_image(aes(x = x + 2,y=y,image = image))
  }
  return(legend_plot)
}

