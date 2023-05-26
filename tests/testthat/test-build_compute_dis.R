data("meta_angola")

library(igraph)

metaweb = meta_angola$metaweb
abTable = meta_angola$abTable
trophicTable = meta_angola$trophicTable

#igraph example
n = 10
g = make_ring(n,directed = T)
A = as.matrix(get.adjacency(g))
meta0 = build_metanet(g)


test_that("no local network",  {
  #no computation of trophic levels
  expect_error(compute_dis(meta0))
})

dis_angola = compute_dis(meta_angola)

test_that("dimension of dis table",{
  expect_equal(length(dis_angola),2)
  expect_equal(length(dis_angola$Species),2)
  expect_equal(dim(dis_angola$Species$nodes),c(2,2))
  dis_angola_Phylum = compute_dis(meta_angola,res = c("Phylum"))
  expect_equal(length(dis_angola_Phylum),1)
})


