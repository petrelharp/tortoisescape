1. `setup-real-G.R` will set up the generator matrix for the non-NA entries of a particular set of layers.  For instance, 
```
  Rscript setup-real-G.R ../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_
```
will produce the files 
- `crop_resampled_masked_aggregated_100x_G.RData` : the generator matrix and associated things:
  - `G` : generator matrix, indexed by **nonmissing** raster pixels
  - `layers` : matrix of **nonmissing* raster values (just two layers; replace this in scripts)
  - `update.G` : function returning `G@x`
  - `ndelta` : required for `update.G`
  - `ngamma` : required for `update.G`
  - `transfn` : transformation function, required for `update.G`
  - `valfn` : computes linear combination of layers, required for `update.G`
 
- `crop_resampled_masked_aggregated_100x_nonmissing.RData` : the necessary info to translate from nonmissing index back to position in the layer
  - `nonmissing` : indices of nonmissing items in layers


2. `setup-tort-locs.R` will figure out which raster cell each tortoise fell in, and save this to, e.g.
```
  Rscript setup-tort-locs.R ../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_
```
will save `locs` and `locs.ij` in `crop_resampled_masked_aggregated_100x_tortlocs.RData`.
  - This also saves `all.locs.ij` in `crop_resampled_masked_aggregated_100x_alllocs.RData` (for use in smoothing).
