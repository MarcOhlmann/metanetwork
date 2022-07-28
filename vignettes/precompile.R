# Must manually move image files from metanetwork/ to metanetwork/vignettes/ after knit
setwd("~/Desktop/metanetwork_project/github/metanetwork/")
library(knitr)
knit("vignettes/vertebrates0.Rmd",output =  "vignettes/vertebrates.Rmd")
