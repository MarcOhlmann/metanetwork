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

test_that("test errors in input of plot_trophicTable",  {
  meta0 = build_metanet(g)
  expect_error(plot_trophicTable(meta0))
  #on angola dataset
  expect_error(plot_trophicTable(meta_angola,res = c('Species')))
  expect_error(plot_trophicTable(meta_angola,res = c('Species','aaa')))
})

test_that("test if output is ggplot2 object",{
  expect_is(plot_trophicTable(meta_angola),"ggplot")
  expect_is(plot_trophicTable(meta_angola,res = c('Species','Phylum')),"ggplot")
})

