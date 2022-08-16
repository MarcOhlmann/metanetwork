library(metanetwork)
library(edgebundle)
library(igraph)
library(ggplot2)
library(ggraph)

data("meta_norway")

meta_norway = append_agg_nets(meta_norway)
meta_norway = compute_TL(meta_norway)

beta = 0.2

meta_norway = attach_layout(g = meta_norway$metaweb_trophic_class,beta = beta,metanetwork = meta_norway)
ggmetanet(g = meta_norway$metaweb_trophic_class,metanetwork = meta_norway,beta = beta,legend = "trophic_class")

group_layout.custom = group_layout.default
group_layout.custom$nbreaks_group = 2
group_layout.custom$group_height = c(5,7)
group_layout.custom$group_width = c(5,7)
meta_norway = attach_layout(metanetwork = meta_norway,beta = beta,mode ="group-TL-tsne",
                            res = "trophic_class",group_layout.config = group_layout.custom)
ggmetanet(meta_norway,beta = beta,mode = "group-TL-tsne",legend = "trophic_class")

ggnet.custom = ggnet.default
ggnet.custom$label = F
ggnet.custom$arrow.size = 3
ggnet.custom$edge.alpha = 0.3
diff_plot(g1 = meta_norway$high,g2 = meta_norway$low,
          meta_norway,beta = beta,mode = "group-TL-tsne",
          layout_metaweb = T,ggnet.config = ggnet.custom)

g = meta_norway$metaweb

#threshold on the edges
g = delete.edges(g, which(E(g)$weight<0.15))
xy = cbind(get.vertex.attribute(g,"group_layout_x_beta0.2"),get.vertex.attribute(g,"group_layout_y_beta0.2"))

fbundle <- edge_bundle_force(g,xy,compatibility_threshold = 0.5,K = 1)
pbundle <- edge_bundle_path(g, xy, max_distortion = 50, weight_fac = 1000, segments = 3000)

# sbundle <- edge_bundle_stub(g,xy)
# hbundle <- edge_bundle_hammer(g,xy,bw = 0.7,decay = 0.5)

ggplot(fbundle)+
  geom_path(aes(x,y),
            size = 0.1,show.legend = FALSE,alpha = 0.5)+
  geom_point(data=as.data.frame(xy),aes(V1,V2,label = "re"),size=2)+
  theme_void()
