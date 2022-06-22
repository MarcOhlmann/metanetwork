# metanetwork = meta_norway
# beta = 0.1
# g = metanetwork$metaweb
# res = "taxa"
# metanetwork = append_agg_nets(metanetwork)
# metanetwork = compute_TL(metanetwork)
# metanetwork = attach_layout(metanetwork, g = metanetwork$metaweb_taxa)
# 
# get_coord_group_TL_tsne <- function(g,metanetwork,group_layout_config = NULL,res,beta){
#   #need for a trophic Table
#   if(is.null(metanetwork$trophicTable)){
#     stop("no trophicTable provided, group layout is impossible !")
#   }
#   networks = metanetwork[lapply(metanetwork,class) == "igraph"]
#   if(!(res %in% sapply(networks,function(g) g$res))){
#     stop("to use group-TL-tsne layout, you need to compute aggregated networks first, see append_agg_nets")
#   }
#   #group-TL-tsne layout is available only at the original resolution
#   if(!(g$res == colnames(metanetwork$trophicTable)[1])){
#     stop("group layout is only available at the orginal resolution")
#   }
#   
#   
#   #get the network at the desired res
#   name_res_tab = sapply(networks,
#                         function(g) list(res = g$res, name = g$name))
#   g_agg = networks[[which(colSums(name_res_tab == c(res,g$name)) == 2)]]
#   if(is.null(igraph::V(g_agg)$TL)){
#     stop("to use group-TL-tsne layout, you need to compute trophic levels first, see compute_TL")
#   }
#   #check if TL-tsne layout is already computed for g_agg
#   if(is.null(igraph::get.vertex.attribute(g_agg,paste0("layout_beta",beta)))){
#     message(paste0("computing TL-tsne layout for ",g$name,"at resolution: ",res,". See attach_layout to store it."))
#     g_agg = attach_layout_g(g_agg,metanetwork,mode = 'TL-tsne',beta,TL_tsne.config = TL_tsne.default)
#   }
#   
#   nbreaks_group = 1
#   group_height = 1
#   group_width = 2
#   group_height = 5
#   group_width = 5
#   
#   groups = V(g_agg)$name
#   if(nbreaks_group == 1){
#     coords_list = c()
#     for(k in 1:length(groups)){
#       sp_loc = metanetwork$trophicTable[
#         which(metanetwork$trophicTable[,res] == groups[k]),
#         1]
#       g_loc = induced_subgraph(metanetwork$metaweb,sp_loc)
#       layout_loc = layout_with_graphopt(g_loc)
#       rownames(layout_loc) = V(g_loc)$name
#       #centering and scaling
#       layout_loc[,1] = group_height*(layout_loc[,1] - mean(layout_loc[,1]))/(sd(layout_loc[,1]))
#       layout_loc[,2] = group_width*(layout_loc[,2] - mean(layout_loc[,2]))/(sd(layout_loc[,2]))
#       #adding the coordinate of the current group in the aggregated layout
#       layout_loc[,1] = layout_loc[,1] + 100*V(g_agg)$TL[k]/max(abs(V(g_agg)$TL))
#       layout_loc[,2] = layout_loc[,2] + 100*igraph::get.vertex.attribute(
#         g_agg,paste0("layout_beta",beta))[k]/
#         max(abs(igraph::get.vertex.attribute(
#           g_agg,paste0("layout_beta",beta))))
#       coords_list = c(coords_list,list(layout_loc))
#     }
#     names(coords_list) = groups
#     #rbinding coordinates
#     coords_mat = do.call(rbind,coords_list)
#   } else{
#     #get group sizes
#     group_sizes = table(metanetwork$trophicTable[,res])[V(g_agg)$name]
#     sapply(1:(nbreaks_group-1),function(k) quantile(group_sizes,k/nbreaks_group))
#   }
#   return(coords_mat)
# }
