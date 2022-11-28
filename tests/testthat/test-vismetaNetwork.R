data("meta_angola")

library(igraph)

metaweb = meta_angola$metaweb
abTable = meta_angola$abTable
trophicTable = meta_angola$trophicTable

#igraph example
n = 10
#directed ring
g = make_ring(n,directed = T)
V(g)$name = as.character(1:n) 
#disconnected local networks
abTable_g = matrix(rep(1,n*2),ncol = n)
abTable_g[1,c(1:2,7:9)] = 0
colnames(abTable_g) = V(g)$name 
rownames(abTable_g) = c("A","B")

test_that("test errors in input of vismetaNetwork",  {
  expect_error(vismetaNetwork(meta_angola,g = meta_angola$X2003,layout_metaweb = T))
})

test_that("output of vismetaNetwork",{
  expect_is(vismetaNetwork(meta_angola),"visNetwork")
  #layout metaweb
  meta_angola = attach_layout(meta_angola)
  expect_is(vismetaNetwork(meta_angola),"visNetwork")
  #legend
  expect_is(vismetaNetwork(meta_angola,legend = "Phylum"),"visNetwork")
})

test_that("diffplot with visnetwork",{
  expect_is(diff_plot(meta_angola,vis_tool = "visNetwork",g1 = meta_angola$X1986,
                      g2 = meta_angola$X2003),"visNetwork")
})
