northdir <- "../tort_272_north_info"
southdir <- "../tort_272_south_info"
dir.create(northdir)
dir.create(southdir)

pcs <- read.csv("pcs.csv",header=TRUE,stringsAsFactors=FALSE)

pc.breaks <- c(-Inf,-.02,0,.02,Inf)
pc.fac <- cut( pcs$PC1, breaks=pc.breaks )
pc.cols <- c("blue","lightblue","mediumpurple","purple")[as.numeric(pc.fac)]

n.torts <- pcs$etort[as.numeric(pc.fac)==4]
s.torts <- pcs$etort[as.numeric(pc.fac)==1]

# save pcs and locations
single.csvs <- c("pcs.csv","272torts_metadata.csv")
for (csvf in single.csvs) {
    csv <- read.csv(csvf,header=TRUE,stringsAsFactors=FALSE)
    n.csv <- subset( csv, (csv[,1] %in% n.torts) )
    s.csv <- subset( csv, (csv[,1] %in% s.torts) )
    stopifnot(nrow(n.csv)==length(n.torts))
    stopifnot(nrow(s.csv)==length(s.torts))
    write.csv(n.csv,file=file.path(northdir,csvf),row.names=FALSE)
    write.csv(s.csv,file=file.path(southdir,csvf),row.names=FALSE)
}


require(raster)
require(rgdal)
dem <- raster("../geolayers/expanded/64x/dem_30.grd")
tort.loc.obj <- load("geog_coords.RData")
tort.locs <- spTransform(get(tort.loc.obj),CRSobj=CRS(proj4string(dem)))
stopifnot(all(row.names(tort.locs)==pcs$etort))

# save location objects
n.tort.locs <- tort.locs[row.names(tort.locs)%in%n.torts]
s.tort.locs <- tort.locs[row.names(tort.locs)%in%s.torts]
save(n.tort.locs,file=file.path(northdir,"geog_coords.RData"))
save(s.tort.locs,file=file.path(southdir,"geog_coords.RData"))

# plots
pdf(file=file.path(northdir,"north-torts.pdf"),width=10,height=5,pointsize=10)
layout(t(1:2))
plot(dem)
points(tort.locs,col=adjustcolor(pc.cols,0.75),cex=0.5,pch=20)
points(tort.locs[as.numeric(pc.fac)==4])
plot(dem)
text(tort.locs[as.numeric(pc.fac)==4],col=adjustcolor(pc.cols,0.75),cex=0.5,labels=gsub(" *(sheared)","",gsub("etort-","",n.torts)))
dev.off()

# plots
pdf(file=file.path(southdir,"south-torts.pdf"),width=10,height=5,pointsize=10)
layout(t(1:2))
plot(dem)
points(tort.locs,col=adjustcolor(pc.cols,0.75),cex=0.5,pch=20)
points(tort.locs[as.numeric(pc.fac)==1])
plot(dem)
text(tort.locs[as.numeric(pc.fac)==1],col=adjustcolor(pc.cols,0.75),cex=0.5,labels=gsub(" *(sheared)","",gsub("etort-","",s.torts)))
dev.off()

# save pairwise csvs
pairwise.csvs <- c("all_angsd_snps.pwp.csv","covariance.csv","geog_distance.csv")
for (csvf in pairwise.csvs) {
    csv <- read.csv(csvf,header=TRUE,stringsAsFactors=FALSE)
    n.csv <- subset( csv, (etort1 %in% n.torts) & (etort2 %in% n.torts) )
    s.csv <- subset( csv, (etort1 %in% s.torts) & (etort2 %in% s.torts) )
    write.csv(n.csv,file=file.path(northdir,csvf),row.names=FALSE)
    write.csv(s.csv,file=file.path(southdir,csvf),row.names=FALSE)
}

