library(GGally)
library(network)
library(sna)
library(ggplot2)
library(intergraph)
library(dplyr)
library(igraph)
library(Matrix)
library(Matrix.utils)
library(ade4)
library(visNetwork)

remove.packages("metanetwork", lib="~/R/x86_64-pc-linux-gnu-library/3.6/")
setwd("~/Bureau/metanetwork/")
devtools::install(dependencies = F)
library(metanetwork)

setwd("~/Desktop/metanetwork_project/metanetwork/vignette/")

set.seed(2458)
###lattice no colors#######
###########################

# 
#generate a directed lattice
g = make_lattice(dim = 2,length = 10,directed = T)
n = vcount(g)
V(g)$name = as.character(1:n)

#sampling a presence table
presence = rbind(rbinom(n,size = 1,prob = 0.4),
                 rbinom(n,size = 1,prob = 0.7),
                 rep(1,n))

rownames(presence) = c('a','b','c')
colnames(presence) = V(g)$name

#build_metanetwork#
#build metanetwork object (of S3 class 'metanetwork')
meta0 = build_metanetwork(g,abTable = presence)
class(g)
class(meta0)
#print#
print(meta0)
#abundances
V(meta0$c)$ab
V(meta0$metaweb)$ab


#compute_trophic_levels#
#compute trophic levels for metaweb and local networks
meta0 = compute_trophic_levels(meta0)
#trophic levels
V(meta0$metaweb)$name
V(meta0$metaweb)$TL

#ggmetanet#
#represent the complete lattice
beta = 20
ggmetanet(meta0$c,beta = beta,metanetwork = meta0)
#represent the metaweb (notice the different nodes ab)
beta = 20
ggmetanet(meta0$metaweb,beta = beta,metanetwork = meta0)
#represent fragmented lattice
beta = 20
ggmetanet(meta0$b,beta = beta,metanetwork = meta0)
#more fragmented
beta = 0.1
ggmetanet(meta0$a,beta =  beta,metanetwork = meta0)

#vismetanetwork#
lattice_metaweb = vismetaNetwork(g = meta0$metaweb,
                                 metanetwork = meta0,beta = 30,x_y_range = c(1,50))
lattice_b = vismetaNetwork(g = meta0$b,
                           metanetwork = meta0,beta = 10,x_y_range = c(2,100))
setwd("~/Desktop/metanetwork_project/test_visnetwork/")
visSave(lattice_metaweb, file = "lattice_metaweb.html")
visSave(lattice_b, file = "lattice_b.html")

#diff_plot#
#b - c
beta = 20
diff_plot(g1 = meta0$b,g2 = meta0$c,beta = beta,metanetwork = meta0)
ggmetanet(meta0$b,beta = 0.9,metanetwork = meta0)

# a - c
diff_plot(g1 = meta0$a,g2 = meta0$c,beta = beta,metanetwork = meta0)
ggmetanet(meta0$a,beta = 10,metanetwork = meta0)


#example lattice color
#####################

#add a trophic table
trophicTable = cbind(V(g)$name,
                     c(rep("A",16),rep("B",16),rep("C",16),rep("D",16)))
colnames(trophicTable) = c('name','group')

#generate a new metanetwork object
meta1 = build_metanetwork(metaweb = g,abTable = presence,
                          trophicTable = trophicTable)

#computing trophic levels
meta1 = compute_trophic_levels(meta1)

#plot_trophic_table
plot_trophicTable(meta1)

#ggmetanet
ggmetanet(g = meta1$c,beta = 30,metanetwork = meta1,
          legend = 'group')
ggmetanet(g = meta1$b,mode = 'TL-kpco',beta = 10,metanetwork = meta1,
          legend = 'group')
ggmetanet(g = meta1$a,mode = 'TL-kpco',beta =  0.3,metanetwork = meta1,
          legend = 'group')

#vismetanetwork
lattice_b = vismetaNetwork(g = meta1$b,metanetwork = meta1
                                 ,beta = 30,x_y_range = c(3,300),legend = 'group')
visSave(lattice_b, file = "lattice_b_colors.html")


###GLOBNETS DATA####
####################

#data loading
setwd(dir = "~/Desktop/globnets/data/")
load('metaweb_globnets.Rdata')
load('abTable_globnets.Rdata')
trophicTable = read.csv('trophicTable.csv')
colnames(abTable_globnets) = trophicTable$trophic_group
V(metaweb_globnets)$name = colnames(abTable_globnets)
abTable_globnets = abTable_globnets

#build metanetwork object
meta_globnets = build_metanetwork(metaweb = metaweb_globnets,abTable = abTable_globnets,
                                  trophicTable = trophicTable)

print(meta_globnets)

#get aggreageted_networks#
meta_globnets = append_aggregated_networks(meta_globnets)
print(meta_globnets)
#attributes of each network
meta_globnets$S1L1_01$name
meta_globnets$S1L1_01$res
#compute trophic levels
meta_globnets = compute_trophic_levels(meta_globnets)

#plot_trophicTable
plot_trophicTable(meta_globnets)
plot_trophicTable(meta_globnets,resolutions = c('taxa','trophic_group'))
plot_trophicTable(meta_globnets,resolutions = c('taxa','trophic_class'))

#ggmetanet with TL-tsne
ggmetanet(metanetwork = meta_globnets,mode = 'TL-tsne',legend = 'taxa',beta = 0.006)
ggmetanet(g = meta_globnets$metaweb_trophic_class,metanetwork = meta_globnets)

#vismetaNetwork
vismetaNetwork(metanetwork = meta_globnets,mode = 'TL-tsne',
               beta = 0.006,x_y_range = c(10,0.025),legend = 'taxa')

#adding covariable#
##################

#disturbance
disturb = rep("low",nrow(meta_globnets$abTable))
names(disturb) = rownames(meta_globnets$abTable)
disturb[grep(x = names(disturb),pattern = 'S1|S3')] = "high"

meta_globnets = build_metanetwork(metaweb = metaweb_globnets,abTable = abTable_globnets,
                                  trophicTable = trophicTable,covariable = disturb)

meta_globnets = append_aggregated_networks(meta_globnets)
meta_globnets = append_mean_networks(meta_globnets)
#compute trophic levels (again, bug to fix)
meta_globnets = compute_trophic_levels(meta_globnets)
diff_plot(g2 = meta_globnets$mean_net_high,g1 = meta_globnets$mean_net_low,
          mode = 'TL-tsne',beta = 0.006,metanetwork = meta_globnets)

#
diff_plot(g1 = meta_globnets$S1L1_03,g2 = meta_globnets$S2L1_06,
mode = 'TL-tsne',beta = 0.006,metanetwork = meta_globnets)

#ggnet custom parameters
ggnet.custom = ggnet.default
ggnet.custom$edge.alpha_diff = 0.5
ggnet.custom$edge.alpha = 0.1
diff_plot(g1 = meta_globnets$S1L1_03,g2 = meta_globnets$S2L1_07,
          mode = 'TL-tsne',beta = 0.006,metanetwork = meta_globnets,ggnet.config = ggnet.custom)

#at another resolution
ggnet.custom$edge.alpha_diff = 1
ggnet.custom$edge.alpha = 1
diff_plot(g1 = meta_globnets$mean_net_high_trophic_class,g2 = meta_globnets$mean_net_low_trophic_class,
          beta = 0.1,metanetwork = meta_globnets,ggnet.config = ggnet.custom)

