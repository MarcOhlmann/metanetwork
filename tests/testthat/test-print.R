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

test_that("test print",{
  meta0 <- build_metanet(g)
 # with single resolution
  expect_output(print(meta0),
                paste("metaweb has 10 nodes and 10 edges ",
                      "single network ",
                      "single resolution available ",
                      sep = "\\n"))
  #on meta_angola
  expect_output(print(meta_angola),
                paste("metaweb has 28 nodes and 127 edges ",
                      "2 local networks ",
                      "available resolutions are: Species Phylum ",
                      sep = "\\n"))
})

