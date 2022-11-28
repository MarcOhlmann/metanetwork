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

test_that("computation of TL-tsne axis",{
  meta0 <- build_metanet(g)
  meta0 = compute_TL(meta0)
  meta0 = attach_layout(meta0)
  expect_equal(length(vertex_attr_names(meta0$metaweb)),4)
  meta0 = attach_layout(meta0)
  expect_equal(length(vertex_attr_names(meta0$metaweb)),5)
})

test_that("computation of group-TL-tsne layout",{
  #on angola dataset
  #error when no resolution is provided
  expect_error(attach_layout(meta_angola,mode = "group-TL-tsne"))
  #attach group layout but not TL-tsne layout
  meta1 <- attach_layout(meta_angola,mode = "group-TL-tsne",res = 'Phylum')
  expect_equal(length(grep("beta",vertex_attr_names(meta1$metaweb))),2)
  #attach group layout and TL-tsne layout
  meta2 <- attach_layout(meta_angola)
  meta2 = attach_layout(meta2,mode = "group-TL-tsne",res = 'Phylum')
  expect_equal(length(grep("beta",vertex_attr_names(meta2$metaweb))),3)
  #add different sizes
  group_layout.custom = group_layout.default
  group_layout.custom$nbreaks_group = 2
  group_layout.custom$group_height = c(5,2)
  group_layout.custom$group_width = c(2,5)
  meta2 = attach_layout(meta2,mode = "group-TL-tsne",res = 'Phylum',
                        group_layout.config = group_layout.custom)
  expect_equal(length(grep("beta",vertex_attr_names(meta2$metaweb))),3)
})


