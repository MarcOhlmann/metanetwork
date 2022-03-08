# library(GGally)
# library(network)
# library(sna)
# library(ggplot2)
# library(intergraph)
# library(dplyr)
# library(igraph)
# library(Matrix)
# library(Matrix.utils)
# library(ade4)


#with metanetwork package

#installation
library(devtools)
setwd("~/Desktop/metanetwork_project/metanetwork/")
devtools::install()

library(metanetwork)
library(igraph)
search()
setwd("~/Desktop/globnets")


#import the globnets files

load("data/abTable.RData")
trophic_groups_class <- read.table('data/trophicTable.csv',sep = ",",header = T,stringsAsFactors = F)
load("data/metaweb.RData")
metaweb <- g.metaweb.t
metaweb = igraph::permute(metaweb,order(order(V(metaweb)$name)))
abTable <- t(abTable.t)
abTable = abTable[,V(metaweb)$name]
trophicTable = cbind(trophic_groups_class$trophic_group,trophic_groups_class$trophic_class,trophic_groups_class$taxa)
colnames(trophicTable) = c('trophic_group','trophic_class','taxa')

V(metaweb)$name = trophic_groups_class$trophic_group[order(trophic_groups_class$trophic_group)]
colnames(abTable) = V(metaweb)$name
abTable = abTable[1:5,]

covariable = c('high','low','low','low','high')
names(covariable) = rownames(abTable)

metanetwork = build_metanetwork(metaweb,abTable,trophicTable,compute_local_networks = T,covariable)
metanetwork::print(metanetwork)
metanetwork = append_aggregated_networks(metanetwork)
metanetwork::print(metanetwork)

metanetwork = compute_trophic_levels(metanetwork)
V(metanetwork$metaweb)$name
V(metanetwork$metaweb)$TL

metanetwork = append_mean_networks(metanetwork)
metanetwork::print(metanetwork)
names(metanetwork)

ggmetanet(metanetwork =  metanetwork)
ggmetanet(metanetwork =  metanetwork,legend = 'taxa',mode = 'TL-tsne',
          alpha_per_group = list(groups = c('Bacteria','Animals','Fungi'),alpha_focal = 0.8,alpha_hidden =0),beta = 0.006)

ggnet.custom = ggnet.default
ggnet.custom$label.size = 1.5*ggnet.default$label.size
ggmetanet(metanetwork =  metanetwork,legend = 'taxa',mode = "TL-tsne",alpha_interactive = T,
          beta = 0.006,ggnet.config = ggnet.custom)

#for small networks more regular
ggmetanet(g = 'metaweb_taxa',metanetwork =  metanetwork)

#diffplot
diff_plot(g1 = metanetwork$mean_net_high,g2 = metanetwork$mean_net_low,
          metanetwork)

diff_plot(g1 = metanetwork$mean_net_high,g2 = metanetwork$mean_net_low,
          metanetwork)

diff_plot(g1 = metanetwork$mean_net_high,g2 = metanetwork$mean_net_low,
          metanetwork,mode = 'TL-tsne')

ggnet.custom = ggnet.default
ggnet.custom$edge.alpha = 0.05
diff_plot(g1 = metanetwork$mean_net_high,g2 = metanetwork$mean_net_low,
          metanetwork,mode = 'TL-tsne',ggnet.config = ggnet.custom)

#some constrains
diff_plot(g1 = metanetwork$mean_net_high_Trophic_class,g2 = metanetwork$mean_net_low,
          metanetwork)

ggnet.custom2 = ggnet.default
ggnet.custom2$edge.size = 2
diff_plot(g1 = metanetwork$mean_net_high_trophic_class,g2 = metanetwork$mean_net_low_trophic_class,
          metanetwork,ggnet.config = ggnet.custom2)

#plot_trophicTable
plot_trophicTable(metanetwork)

#vismetaNetwork
vismetaNetwork(metanetwork = metanetwork)




