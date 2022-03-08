setwd("~/Bureau/metanetwork/")
devtools::install(dependencies = F)

library(metanetwork)
library(igraph)

setwd("~/Desktop/metanetwork_project/metanetwork/vignette/")

#the data
print(barents)
meta_barents = barents

# append aggregated networks and compute trophic levels
meta_barents = append_aggregated_networks(meta_barents)
print(meta_barents)
meta_barents = compute_trophic_levels(meta_barents)

#plot_trophicTable
plot_trophicTable(meta_barents)

#ggmetanet#
##########
g = meta_barents$R0.20
ggmetanet(g = g,metanetwork = meta_barents,legend = 'group',mode = 'TL-tsne',beta = 0.006)

#more
ggmetanet(metanetwork = meta_barents,legend = 'group',groups = c('Basal','Seabirds','Mammals'))

#diffplot#
g1 = meta_barents$R0.18
g2 = meta_barents$R0.12
diff_plot(g1,g2,meta_barents,mode = 'TL-tsne',beta = 0.006)

#vismetaNetwork
vismetaNetwork(g,metanetwork = meta_barents,legend = 'group',mode = 'TL-tsne',
               beta = 0.006,x_y_range = c(10,0.01))
