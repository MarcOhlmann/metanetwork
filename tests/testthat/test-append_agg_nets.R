data("meta_angola")

library(igraph)

metaweb <- meta_angola$metaweb
abTable <- meta_angola$abTable
trophicTable <- meta_angola$trophicTable

test_that("check computation of the aggregated networks",  {
  meta0 <- build_metanet(metaweb,abTable)
  meta1 <- build_metanet(metaweb,abTable,trophicTable)
  ## No possible aggregation
  expect_error(append_agg_nets(meta0))
  #test for computation of aggregated networks
  meta_1app <- append_agg_nets(meta1)
  res_net1app <- sapply(extract_networks(meta_1app),function(g) g$res) %>%
                unique()
  expect_equal(res_net1app,colnames(trophicTable))
  #no action when aggregated networks are already computed
  expect_equal(meta_angola,append_agg_nets(meta_angola))
})

