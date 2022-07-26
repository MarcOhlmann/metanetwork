# Must manually move image files from metanetwork/ to metanetwork/vignettes/ after knit

library(knitr)
knit("vignettes/vertebrates.Rmd.orig",output =  "vignettes/vertebrates.Rmd")
