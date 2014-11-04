Setup
====

The script `setup-everything.sh` will, for each of several spatial resolutions:

1. `make-overlap-na-layer.R` : figures out the common set of nonmissing locations
```
  Rscript make-overlap-na-layer.R ../geolayers/TIFF/500x/500x_ six-raster-list
```
will create the `six-raster-list-na.grd` layer in `../geolayers/TIFF/500x/`.

2. `setup-real-G.R` will set up the generator matrix for the non-NA entries of a particular set of layers.  For instance, 
```
  Rscript setup-real-G.R ../geolayers/TIFF/500x/500x_ 500x six-raster-list
```
will produce the files 
- `500x/500x__six-raster-list_G.RData` : the generator matrix and associated things:
  - `G` : generator matrix, indexed by **nonmissing** raster pixels
  - `layers` : matrix of **nonmissing** raster values, used to update G
  - `update.G` : function returning `G@x`
  - `ndelta` : required for `update.G`
  - `ngamma` : required for `update.G`
  - `transfn` : transformation function, required for `update.G`
  - `valfn` : computes linear combination of layers, required for `update.G`
 
- `500x/500x__six-raster-list_nonmissing.RData` : the necessary info to translate from nonmissing index back to position in the layer
  - `nonmissing` : indices of nonmissing items in layers


3. `setup-tort-locs.R` will figure out which raster cell each tortoise fell in, and their neighborhood, e.g.
```
  Rscript setup-tort-locs.R ../geolayers/TIFF/500x/500x_ 500x six-raster-list
```
will produce the following files, each in the same order as tortoise locations in `../tort.coords.rasterGCS.Robj`:

- `500x/500x_tortlocs.RData`: contains `locs`, which gives the indices of **nonmissing** raster cells that tortoises fall in.
- `500x/500x_six-raster-list_neighborhoods.RData`: contains `neighborhoods`, which is a list with indices of **nonmissing** raster cells within 15km of that tortoise; nearby raster cells not falling on the layer are recorded as `NA`.


4. `setup-inference.R` will load these things and `pimat` up in an `.RData`, e.g.
```
Rscript setup-inference.R ../geolayers/TIFF/500x/500x_ 500x six-raster-list ../pairwisePi/alleleCounts_1millionloci.pwp
```
will save everything plus the kitchen sink to:
- `500x/six-raster-list-500x_setup.RData`


Computing resistance distances
==============================

To find mean hitting times on a fine grid at a given set of parameters, as in `dem-hitting-time.sh`, after running `setup-dem-inference.sh`:
1. Use `make-resistance-distances.R` on a small grid to find initial guesses with analytic inference, e.g.
```
Rscript make-resistance-distances.R ../geolayers/TIFF/500x/500x_ 500x ../inference/six-raster-list simple-init-params-six-raster-list.tsv analytic 
```
2. Push these up to a finer grid, using `disaggregate-ht.R`, e.g.
```
Rscript disaggregate-ht.R ../geolayers/TIFF/500x/500x_ ../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_ 500x 100x six-raster-list 500x/six-raster-list-hitting-times.tsv 5
```
which produces `./100x/500x-aggregated-hitting-times.tsv`.
3. Use these as a starting point for inference on the finer grid, e.g.
```
Rscript make-resistance-distances.R ../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_ 100x ../inference/six-raster-list simple-init-params-six-raster-list.tsv CG 100x/500x-aggregated-hitting-times.tsv 120
```



Inference
=========

**Outline:**
0. Begin with reasonable guess at parameter values, by solving very-low-resolution problem exactly, and construct `G` matrix.
1. Interpolate observed mean pairwise divergence times using `G` to get estimate of full matrix of divergence times, as in `interp-inference.R`.
  1. use `initial-hitting-times.R` as e.g.
```
Rscript initial-hitting-times.R ../geolayers/TIFF/500x/500x_ 500x six-raster-list
```
  which produces `500x_six-raster-list-init-hts.RData`
2. Given full matrix of divergence times to infer parameter values, as in `exponential-transform.R`.
  1. use `fit-exponential-model.R` as e.g.
```
Rscript fit-exponential-model.R ../geolayers/TIFF/500x/500x_ 500x six-raster-list
```
3. Return to (1) if necessary.
Also, do this for sequentially finer grids, using previously inferred parameter values to start the next,
multiplied by the square of the ratio of the two grid sizes.



