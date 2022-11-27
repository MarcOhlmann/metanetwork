data("meta_angola")

library(igraph)

metaweb = meta_angola$metaweb
abTable = meta_angola$abTable
trophicTable = meta_angola$trophicTable

#igraph example
n = 10
g = make_ring(n,directed = T)
A = as.matrix(get.adjacency(g))
g_u = make_ring(n)
V(g_u)$name = as.character(1:n) 
g_d = erdos.renyi.game(n,n-2,type = "gnm",directed = T)
V(g_d)$name = as.character(1:n) 

test_that("property of the metaweb",  {
  ## No names for metaweb (igraph)
  expect_warning(build_metanet(g))
  ## No names for metaweb (Matrix)
  expect_warning(build_metanet(A))
  #metaweb must be directed
  expect_error(build_metanet(g_u))
  #metaweb must be connected
  expect_error(build_metanet(g_d))
})

test_that("correspondance of node names between metaweb, abTable and trophic table",{
  abTable_r = abTable
  colnames(abTable_r)[1] = "aa"
  trophicTable_r = trophicTable
  trophicTable_r[1,1] = "aa"
  expect_error(build_metanet(metaweb,abTable_r,trophicTable))
  expect_error(build_metanet(metaweb,abTable,trophicTable_r))
  expect_error(build_metanet(metaweb,abTable_r,trophicTable_r))
})

test_that("extract networks",{
  # test on the number of local networks
  expect_equal(length(extract_networks(meta_angola)),6)
  expect_equal(length(extract_networks(build_metanet(metaweb,abTable))),3)
})