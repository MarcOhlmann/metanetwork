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


test_that("TL computation",  {
  #no computation of trophic levels
  expect_error(compute_metrics(meta0))
})

meta0 = compute_TL(meta0)

test_that("dimension of metrics table",{
  expect_equal(length(compute_metrics(meta0)),4)
  expect_equal(class(compute_metrics(meta0)),"numeric")
  metrics_angola = compute_metrics(meta_angola)
  expect_equal(length(metrics_angola),2)
  expect_equal(dim(metrics_angola$Species),c(3,4))
})

