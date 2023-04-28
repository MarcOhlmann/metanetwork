# Must manually move image files from metanetwork/ to metanetwork/vignettes/ after knit
library(knitr)
#knit("vertebrates0.Rmd",output =  "vertebrates.Rmd")
knit("angola0.Rmd",output =  "angola.Rmd")
knit("norway0.Rmd",output =  "norway.Rmd")
knit("pyramid0.Rmd",output =  "pyramid.Rmd")
