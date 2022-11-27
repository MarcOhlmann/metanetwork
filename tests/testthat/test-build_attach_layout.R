data("meta_angola")

library(igraph)

metaweb <- meta_angola$metaweb
abTable <- meta_angola$abTable
trophicTable <- meta_angola$trophicTable

#igraph example
n <- 10
g <- make_ring(n,directed = T)
A <- as.matrix(get.adjacency(g))
V(g)$name <- as.character(1:n)
g_u <- make_ring(n)
V(g_u)$name <- as.character(1:n) 
g_d <- erdos.renyi.game(n,n-2,type = "gnm",directed = T)
V(g_d)$name <- as.character(1:n) 

test_that("test if enough output are provided",{
  meta0 <- build_metanet(g)
 # are trophic levels computed ?
  expect_error(attach_layout(meta0))
  #correct resolution
})

test_that("test of the computation of TL-tsne axis",{
  meta0 <- build_metanet(g)
  meta0 = compute_TL(meta0)
  meta0 = attach_layout(meta0)
  expect_equal(length(vertex_attr_names(meta0$metaweb)),4)
  meta0 = attach_layout(meta0)
  expect_equal(length(vertex_attr_names(meta0$metaweb)),5)
})



