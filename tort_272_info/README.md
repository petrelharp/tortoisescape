Useful files in this directory:


- [pcs.csv](pcs.csv) with first two PCs for each sample.
- [covariance.csv](covariance.csv) matrix of covariances for each pair of samples.
- [geog_distance.csv](geog_distance.csv) matrix of geographic distances for each pair of samples

- [geog_coords.RData](geog_coords.RData) saved R file of the SpatialPoints object with sample locations


# Data processing notes

- The file `torts_272_coordinates.csv` had one tortoise location as being in UTM zone 12.  
    It is clearly in zone 11; and coordinates are also present in [272torts_metadata.csv](272torts_metadata.csv),
    so `torts_272_coordinates.csv` was deleted.  
    This error may be still be present in the upstream `.xslx` file, `DTfieldsamples.xls`.

- The tortoise `etort-231` is in the middle of the Colorado; it's sample ID in `Tortoise Locations 2010-2012.xslx` is `2012082`,
    which is grouped along with `2012083` and `2012004` in that file, which correspond to `etort-256` and `etort-278`.
    The closest relatives of `etort-231` are in the area of these two latter tortoises.
    `etort-231`'s Northing is listed as `3923493`; while the others are `3821051` and `3817182`, respectively;
    changing the Northing to `3823493` puts the tortoise right next to the others.
    Therefore, this was done.
    This may be still present in `Tortoise Locations 2010-2012.xslx`.
