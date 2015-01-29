#!/usr/bin/Rscript

# Copy or soft link this to a directory with a bunch of .grd layers,
# then run it ( `./view-layers.R` )
# and it will produce an html doc with the layers in.

rmd <- '
```{r setup, include=FALSE}
    require(raster)
    require(rgdal)
    require(colorspace)
    fig.dim <- 8
    opts_chunk$set(fig.width=fig.dim,fig.height=fig.dim,fig.align="center")
    layer.files <- list.files(pattern="((grd)|(tif)|(asc))$")
    layer.names <- gsub(".((grd)|(tif)|(asc))$","",layer.files)
    sample.loc.obj <- load(file.path(gsub("tortoisescape.*","tortoisescape",normalizePath(getwd())),"tort_272_info/geog_coords.RData"))
    assign("sample.locs",get(sample.loc.obj))
    tort.ids <- row.names(sample.locs)
    pcs <- read.csv(file.path(gsub("tortoisescape.*","tortoisescape",normalizePath(getwd())),"tort_272_info/pcs.csv"),header=TRUE)
    # pc.cols <- adjustcolor( ifelse( pcs$PC1[match(tort.ids,pcs$etort)] > 0, "blue", "purple" ), .75 )
    pc.cols <- adjustcolor(diverge_hcl(64,h=c(240,300),c=100,l=c(100,50))[cut(pcs$PC1,breaks=64)],.75)
    dem.name <- grep("dem_30$", layer.names, value=TRUE)
    do.contour <- ( length(dem.name) > 0 )
    if (do.contour) { dem <- raster(dem.name) }
    layer <- raster(layer.files[1])
    sample.locs <- spTransform( sample.locs, CRSobj=CRS(proj4string(layer)))
```
```{r plot_layers, echo=FALSE, warning=FALSE}
    for (k in seq_along(layer.files)) {
        layer <- raster(layer.files[k])
        plot(layer,main=layer.names[k])
        if (do.contour) { contour(dem,add=TRUE,col=adjustcolor("black",0.2)) }
        points(sample.locs,pch=21,cex=1,lwd=0.5)
        points(sample.locs,pch=20,cex=1,col=pc.cols)
    }
```
'

require(knitr)
html <- knit2html(text=rmd,options="base64_images")
cat(html, file="layers.html")
