# Watershed data

## Watershed Boundary Dataset (WBD)


Downloaded national data from [rockyftp.usgs](ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/WBD/FileGDB931/),
to file geodatabase (their shapefile was missing, grrr!!) [WBD_National_931.gdb/](WBD_National_931.gdb)
What we want out of this are the layers named
`WBDHUn`, where `n` is even and ranges from 2 (coarsest) to 12 (finest).
This *can* be read with `readOGR()`.

## Hydrography

Downloaded national data from the [National Hydrography Dataset](http://nhd.usgs.gov/data.html)
(of which WBD is actually a part)
for subbasins 1807, 1809, 1810, 1606, 1501, and 1503 to [NHD/](NHD/):
```
BASEURL="ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/Hydrography/NHD/HU4/HighResolution/Shape/NHD_H_"
for x in 1807 1809 1810 1606 1501 and 1503
do
    wget ${BASEURL}${x}_Shape.jpg
    wget ${BASEURL}${x}_Shape.xml
    wget ${BASEURL}${x}_Shape.zip
done
```


# Creating migration matrices from polygons

These would be possible, using tools in rGEOS.
For instance, as [mentioned here](http://stackoverflow.com/questions/26499010/finding-adjacent-polygons-in-r-neighbors)
a quick way to find adjacent polygons is
```
nb <- spdep::poly2nb( vegmap, foundInBox=rgeos::gUnarySTRtreeQuery(vegmap) )
```
This took 7 minutes on the western desert_veg dataset, with 47,514 polygons.


