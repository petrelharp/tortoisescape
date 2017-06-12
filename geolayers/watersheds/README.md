# Watershed data

## Watershed Boundary Dataset (WBD)


Downloaded national data from [rockyftp.usgs](ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/Hydrography/WBD/National/GDB/),
to file geodatabase (their shapefile was missing, grrr!!) [WBD_National_GDB.gdb/](WBD_National_GDB.gdb)
What we want out of this are the layers named
`WBDHUn`, where `n` is even and ranges from 2 (coarsest) to 12 (finest).
This *can* be read with `readOGR()`.
```
mkdir -p WBD; cd WBD
BASEURL="ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/Hydrography/WBD/National/GDB"
wget ${BASEURL}/WBD_National_GDB.zip
mkdir -p WBD_Metadata; cd WBD_Metadata
for x in HU2 HU4 HU6 HU8 HU10 HU12 HU14 HU16 Line _FeatureDataset
do
    wget ${BASEURL}/WBD_Metadata/WBD${x}_FGDC.xml
done
```

## Hydrography

Downloaded national data from the [National Hydrography Dataset](http://nhd.usgs.gov/data.html)
(of which WBD is actually a part)
for subbasins 1807, 1809, 1810, 1606, 1501, and 1503 to [NHD/](NHD/):
```
mkdir -p NHD; cd NHD
BASEURL="ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/Hydrography/NHD/HU4/HighResolution/Shape/NHD_H_"
for x in 1807 1809 1810 1606 1501 and 1503
do
    wget ${BASEURL}${x}_Shape.jpg
    wget ${BASEURL}${x}_Shape.xml
    wget ${BASEURL}${x}_Shape.zip
done
for x in *zip; do y=$(echo $x | sed -e "s/_Shape.zip//"); mkdir -p $y; unzip -d $y $x; done
```


# Creating migration matrices from polygons

These would be possible, using tools in rGEOS.
For instance, as [mentioned here](http://stackoverflow.com/questions/26499010/finding-adjacent-polygons-in-r-neighbors)
a quick way to find adjacent polygons is
```
nb <- spdep::poly2nb( vegmap, foundInBox=rgeos::gUnarySTRtreeQuery(vegmap) )
```
This took 7 minutes on the western desert_veg dataset, with 47,514 polygons.


