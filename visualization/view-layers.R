#!/usr/bin/Rscript

# Copy or soft link this to a directory with a bunch of .grd layers,
# then run it ( `./view-layers.R` )
# and it will produce an html doc with the layers in.

rmd <- '
```{r setup, include=FALSE}
    require(raster)
    fig.dim <- 8
    opts_chunk$set(fig.width=fig.dim,fig.height=fig.dim,fig.align="center")
    layer.names <- gsub(".grd$","",list.files(pattern="*.grd"))
    tort.loc.obj <- load(file.path(gsub("tortoisescape.*","tortoisescape",normalizePath(getwd())),"tort_272_info/geog_coords.RData"))
    assign("tort.locs",get(tort.loc.obj))
    dem.name <- grep("dem_30$", layer.names, value=TRUE)
    do.contour <- ( length(dem.name) > 0 )
    if (do.contour) { dem <- raster(dem.name) }
```
```{r plot_layers, echo=FALSE, warning=FALSE}
    for (ln in layer.names) {
        layer <- raster(ln)
        plot(layer,main=ln)
        points(tort.locs,pch=20,cex=0.5)
        if (do.contour) { contour(dem,add=TRUE,col=adjustcolor("black",0.2)) }
    }
```
'

require(knitr)
html <- knit2html(text=rmd,options="base64_images")
cat(html, file="layers.html")
