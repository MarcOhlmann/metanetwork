library(metanetwork)
library(igraph)

n = 7
g = make_lattice(c(n,n),directed = T)
V(g)$name = as.character(1:n^2)
meta0 = build_metanetwork(g)
meta0 = compute_trophic_levels(metanetwork = meta0)

abTable =  matrix(0,nrow = 2, ncol = vcount(g))
colnames(abTable) = V(g)$name
abTable[,V(g)[which(V(meta0$metaweb)$TL > 3)]$name] = 1
abTable[2,c("19","20","26","27")] = 0
rownames(abTable)  = c('a','b')

meta_pyramid = build_metanetwork(metaweb = g,abTable = abTable)
meta_pyramid = compute_trophic_levels(meta_pyramid)
ggmetanet(metanetwork = meta_pyramid,g = meta_pyramid$b,beta = 3.5)
ggmetanet(metanetwork = meta_pyramid,g = meta_pyramid$b,beta = 150)
