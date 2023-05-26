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
  expect_error(compute_div(meta0))
})

div_angola = compute_div(meta_angola)

test_that("dimension of div table",{
  expect_equal(length(div_angola),2)
  expect_equal(dim(div_angola$nodes),c(5,2))
})

