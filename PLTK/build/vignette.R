#install.packages("devtools")
library(devtools)

test.package=file.path("~/git", "PLTK/PLTK")
#devtools::create(test.package)
## Modify all package files here

devtools::document(test.package)
#devtools::check(test.package)

#devtools::use_vignette("pltk-vignette", test.package)
#devtools::build_vignettes(test.package)
devtools::build(test.package)

install.packages("~/git/PLTK/PLTK_0.0.0.9000.tar.gz",
                 repos = NULL, type = "source")
#install.packages("roxygen2")

detach("package:PLTK", unload=TRUE)
library(PLTK)



screen(1)
par(mar=c(0.5, 4.1, 2, 2.1))
suppressWarnings(plotGrMetadata(gr.expr, plot.settings, data.type='expr', 
                                yrange=c(0,5), side=3,
                                add.axis=TRUE, axis.mark=4, add.y.axis=TRUE,
                                anno.track=0.4, add.annotations=TRUE,
                                cex=0.5, min.anno.gap=60, add.anno.lines=TRUE,
                                target.chr=t.chr))
close.screen(all.screens=TRUE)



