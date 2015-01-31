#!/usr/bin/Rscript

# Copy or soft link this to a directory with a bunch of .grd layers,
# then run it ( `./view-layers.R` )
# and it will produce an html doc with the layers in.

rmd <- '
```{r setup, include=FALSE}
    require(raster)
    fig.dim <- 8
    opts_chunk$set(fig.width=fig.dim,fig.height=fig.dim,fig.align="center")
    layer.files <- list.files(pattern="[.]((grd)|(tif)|(asc))$")
    layer.names <- gsub("[.]((grd)|(tif)|(asc))$","",layer.files)
    tort.loc.obj <- load(file.path(gsub("tortoisescape.*","tortoisescape",normalizePath(getwd())),"tort_272_info/geog_coords.RData"))
    assign("tort.locs",get(tort.loc.obj))
    dem.name <- list.files("../expanded/64x",pattern="dem_30.grd",full.names=TRUE)
    do.contour <- ( length(dem.name) > 0 )
    if (do.contour) { dem <- raster(dem.name[1]) }
```
```{r plot_layers, echo=FALSE, warning=FALSE}
    for (k in seq_along(layer.files)) {
        layer <- raster(layer.files[k])
        if (!compareRaster(layer,dem,stopiffalse=FALSE)) {
            dem <- projectRaster(dem,layer)
        }
        contour(dem,col=adjustcolor("black",0.2),main=layer.names[k])
        plot(layer,add=TRUE)
        points(tort.locs,pch=20,cex=0.5)
    }
```
'

require(knitr)
html <- knit2html(text=rmd,options="base64_images")
cat(html, file="layers.html")
