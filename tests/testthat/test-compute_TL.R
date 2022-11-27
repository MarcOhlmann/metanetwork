data("meta_angola")

library(igraph)

metaweb <- meta_angola$metaweb
abTable <- meta_angola$abTable
trophicTable <- meta_angola$trophicTable

#igraph example
n <- 10
#directed ring
g <- make_ring(n,directed = T)
V(g)$name <- as.character(1:n) 
#disconnected local networks
abTable_g <- matrix(rep(1,n*2),ncol = n)
abTable_g[1,c(1:2,7:9)] = 0
colnames(abTable_g) <- V(g)$name 
rownames(abTable_g) <- c("A","B")

test_that("Computation of trophic levels of norway dataset",  {
  #test of computation of angola metaweb trophic levels
  metaweb <- delete_vertex_attr(metaweb,"TL")
  meta0 <- build_metanet(metaweb)
  meta0 = compute_TL(meta0)
  expect_equal(V(meta0$metaweb)$TL,V(meta_angola$metaweb)$TL)
  #disconnected local network
  meta00 <- build_metanet(g,abTable_g)
  meta00 = compute_TL(meta00)
  expect_equal(V(meta00$metaweb)$TL[10],V(meta00$A)$TL[which(V(meta00$A)$name == 10)])
})

