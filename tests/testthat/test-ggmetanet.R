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

test_that("test errors in input of ggmetanet",  {
  expect_error(ggmetanet(g))
  expect_error(ggmetanet(meta_angola,mode = "aa"))
  expect_error(ggmetanet(g,legend = "aa"))
  expect_error(ggmetanet(meta_angola,legend = "aa"))
})

test_that("test of ggmetanet mode and layout metaweb",{
  expect_is(ggmetanet(meta_angola),"ggplot")
  expect_is(ggmetanet(meta_angola,legend = "Phylum"),"ggplot")
  expect_is(ggmetanet(meta_angola,mode = "fr"),"ggplot")
  expect_is(ggmetanet(meta_angola,mode = "kk"),"ggplot")
  expect_is(ggmetanet(meta_angola,mode = "circle"),"ggplot")
  #using layout_metaweb
  expect_error(ggmetanet(meta_angola,g = meta_angola$X2003,layout_metaweb = T))
  meta_angola2 = attach_layout(meta_angola)
  expect_is(ggmetanet(meta_angola2,g = meta_angola2$X2003,layout_metaweb = T),"ggplot")
  expect_error(ggmetanet(meta_angola2,nrep_ly = 2))
  meta_angola2 = attach_layout(meta_angola2)
  expect_is(ggmetanet(meta_angola2,nrep_ly = 2),"ggplot")
  #custom layout
  layout_loc = matrix(runif(2*vcount(meta_angola$metaweb)),ncol = 2)
  rownames(layout_loc) = V(meta_angola$metaweb)$name
  layout_err_loc = matrix(runif(3*vcount(meta_angola$metaweb)),ncol = 3)
  rownames(layout_err_loc) = V(meta_angola$metaweb)$name
  expect_is(ggmetanet(meta_angola,mode = layout_loc),"ggplot")
  expect_error(ggmetanet(meta_angola,mode = layout_err_loc))
})

test_that("group-TL-tsne layout",{
  expect_error(ggmetanet(meta_angola2,mode = "group-TL-tsne"))
  meta_angola = attach_layout(meta_angola,mode = "group-TL-tsne",
                              res = "Phylum")
  expect_is(ggmetanet(meta_angola,mode = "group-TL-tsne"),"ggplot")
})

test_that("test of alpha options",{
  expect_is(ggmetanet(meta_angola,legend = "Phylum",
                      alpha_per_node = list(nodes = c("Detritus"),
                                                                  alpha_focal = 0.5,
                                                                  alpha_hidden = 0.2)),
            "ggplot")
  expect_is(ggmetanet(meta_angola,legend = "Phylum",
                      alpha_per_group = list(resolutions = "Phylum",
                                             groups = c("Mammals"),
                                            alpha_focal = 0.5,
                                            alpha_hidden = 0.2)),
            "ggplot")
})

data("meta_vrtb")

test_that("test legend for large networks",{
  #single warning for size cut
  expect_warning(ggmetanet(meta_vrtb,legend = "group",beta = 4e-06))
  expect_warning(ggmetanet(meta_vrtb,legend = "group",alpha_per_node = list(nodes = c("Alytes_cisternasii"),
                                                                       alpha_focal = 0.5,
                                                                       alpha_hidden = 0.2),beta = 4e-06))
  expect_warning(ggmetanet(meta_vrtb,legend = "group",alpha_per_node = list(resolutions = "group",
                                                                            groups = c("1"),
                                                                            alpha_focal = 0.5,
                                                                            alpha_hidden = 0.2),beta = 4e-06))
  ggnet.custom = ggnet.default
  ggnet.custom$img_PATH = "aa"
  expect_error(ggmetanet(meta_vrtb,beta = 4e-06,legend = "group",ggnet.config = ggnet.custom))
})
