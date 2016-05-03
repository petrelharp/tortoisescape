In this directory:

- [map-utils.R](map-utils.R) : source this with `chdir=TRUE` to get

    * `get_counties(x)` : returns a county lines object with same projection as x
    * `get_elev(x)` : returns an elevation raster with same projection as x (from dem_30)
    * `get_contours(x)` : returns elevation contours with same projection as x (from dem_30)
    * `get_shading(x)` : returns an elevation *shading* with same projection as x (from Natural Earth)
    * `.raster.crs` : the projection used by our rasters; project here using `spTransform( , .raster.crs)`


# Notes on mapping:

## Coordinates and projections.

Geographic coordinate systems tells you where something is on the *surface of the earth*;
the way (most?) of these work is to: 

    1.  transfer the point-in-the-universe to the nearest point on an ellipsoid (the **datum**) 
        that approximates the surface of the earth
        (this transfer is a projection in the mathematical sense)
    2.  transfer the point-on-the-big-ellipse into a (Cartesian) coordinate system
        (this transfer is called the **projection**)

So, to say where a thing is in the world, you need to give the coordinates, the projection, and the datum.
All this stuff is specified by the *coordinate reference system* (**CRS**).

Maps have such info associated as well.
So, if you have *coordinates*, i.e., two numbers that say where the thing is on its ellipse,
to put it on a map, you need to

    1.  figure out where the thing is on its ellipsoid
    2.  transfer from its ellipsoid to the map's ellipsoid
    3.  project back onto the map's projection


## Reading maps

Most things can be read with `readOGR`.
Stay away from "file geodatabase"s if possible: some of these are *not* readable by anything non-ESRI.
To find out what layers are in a file, use `ogrinfo` on the command line, or `ogrListLayers()` in R.
