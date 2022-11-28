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
abTable_g[1,c(1:2,7:9)] <- 0
colnames(abTable_g) <- V(g)$name 
rownames(abTable_g) <- c("A","B")

test_that("test errors in input of diffplot",  {
  expect_error(diff_plot(meta_angola,meta_angola$X1986,g))
  expect_error(diff_plot(meta_angola,meta_angola$X1986,meta_angola$X2003_Phylum))
})

test_that("test of diffplot mode and layout metaweb",{
  expect_is(diff_plot(meta_angola,g1 = meta_angola$X1986,g2 = meta_angola$X2003),"ggplot")
  #using layout_metaweb
  expect_error(diff_plot(meta_angola,g1 = meta_angola$X1986,g2 = meta_angola$X2003,
                         layout_metaweb = T))
  meta_angola2 = attach_layout(meta_angola)
  expect_is(diff_plot(meta_angola2,g1 = meta_angola2$X1986,g2 = meta_angola2$X2003,
                      layout_metaweb = T),"ggplot")
})

test_that("group-TL-tsne layout",{
  expect_error(diff_plot(meta_angola2,g1 = meta_angola2$X1986,
                         g2 = meta_angola2$X2003,mode = "group-TL-tsne"))
  expect_error(diff_plot(meta_angola2,g1 = meta_angola2$X1986,
                         g2 = meta_angola2$X2003,mode = "group-TL-tsne",
                         layout_metaweb = T))
  meta_angola = attach_layout(meta_angola,mode = "group-TL-tsne",
                              res = "Phylum")
  expect_is(diff_plot(meta_angola,g1 = meta_angola$X1986,
                      g2 = meta_angola$X2003,mode = "group-TL-tsne",
                      layout_metaweb = T),"ggplot")
})

test_that("test of alpha options",{
  expect_is(diff_plot(meta_angola,g1 = meta_angola$X1986,g2 = meta_angola$X2003,
                      alpha_per_node = list(nodes = c("Detritus"),
                                                                  alpha_focal = 0.5,
                                                                  alpha_hidden = 0.2)),
            "ggplot")
  expect_is(diff_plot(meta_angola,g1 = meta_angola$X1986,g2 = meta_angola$X2003,
                      alpha_per_group = list(resolutions = "Phylum",
                                             groups = c("Mammals"),
                                            alpha_focal = 0.5,
                                            alpha_hidden = 0.2)),
            "ggplot")
})
