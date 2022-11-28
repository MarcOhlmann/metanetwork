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

test_that("Test of metanet_pipe",  {
  meta0 <- build_metanet(g,abTable_g)
  meta0 = metanet_pipe(meta0,beta = 0.1)
  g0 <-meta0$metaweb
  expect_type(V(g0)$TL,"double")
  expect_type(V(g0)$layout_beta0.1,"double")
  
  #on angola data
  meta_ango_loc <- build_metanet(metaweb,abTable,trophicTable)
  meta_ango_loc = metanet_pipe(meta_ango_loc,beta = 0.1)
  g_loc <-meta_ango_loc$metaweb
  expect_type(V(g_loc)$TL,"double")
  expect_type(V(g_loc)$layout_beta0.1,"double")
  
  
})

test_that("Test of metanet_build_pipe",  {
  meta0 <- build_metanet(g,abTable_g)
  meta0 = metanet_build_pipe(g,abTable_g,beta = 0.1)
  g0 <-meta0$metaweb
  expect_type(V(g0)$TL,"double")
  expect_type(V(g0)$layout_beta0.1,"double")   
  
    #on angola data
  meta_ango_loc = metanet_build_pipe(metaweb,abTable,trophicTable,beta = 0.1)
  g_loc <-meta_ango_loc$metaweb
  expect_type(V(g_loc)$TL,"double")
  expect_type(V(g_loc)$layout_beta0.1,"double")
})
