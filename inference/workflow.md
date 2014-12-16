Setup
====

The script `setup-everything.sh` will, for each of several spatial resolutions:

1. `make-overlap-na-layer.R` : figures out the common set of nonmissing locations
```
  Rscript make-overlap-na-layer.R ../geolayers/multigrid/512x/crm_ six-raster-list
```
will create the `six-raster-list-na.grd` layer in `../geolayers/multigrid/512x/".

2. `setup-real-G.R` will set up the generator matrix for the non-NA entries of a particular set of layers.  For instance, 
```
  Rscript setup-real-G.R ../geolayers/multigrid/512x/crm_ 512x six-raster-list
```
will produce the files 
- `512x/crm__six-raster-list_G.RData` : the generator matrix and associated things:
  - `G` : generator matrix, indexed by **nonmissing** raster pixels
  - `layers` : matrix of **nonmissing** raster values, used to update G
  - `update.G` : function returning `G@x`
  - `ndelta` : required for `update.G`
  - `ngamma` : required for `update.G`
  - `transfn` : transformation function, required for `update.G`
  - `valfn` : computes linear combination of layers, required for `update.G`
 
- `512x/crm__six-raster-list_nonmissing.RData` : the necessary info to translate from nonmissing index back to position in the layer
  - `nonmissing` : indices of nonmissing items in layers


3. `setup-tort-locs.R` will figure out which raster cell each tortoise fell in, and their neighborhood, e.g.
```
  Rscript setup-tort-locs.R ../geolayers/multigrid/512x/crm_ 512x six-raster-list
```
will produce the following files, each in the same order as tortoise locations in `../tort.coords.rasterGCS.Robj`:

- `512x/crm_tortlocs.RData`: contains `locs`, which gives the indices of **nonmissing** raster cells that tortoises fall in.
- `512x/crm_six-raster-list_neighborhoods.RData`: contains `neighborhoods`, which is a list with indices of **nonmissing** raster cells within 15km of that tortoise; nearby raster cells not falling on the layer are recorded as `NA`.


4. `setup-inference.R` will load these things and `pimat` up in an `.RData`, e.g.
```
Rscript setup-inference.R ../geolayers/multigrid/512x/crm_ 512x six-raster-list ../pairwisePi/alleleCounts_1millionloci.pwp
```
will save everything plus the kitchen sink to:
- `512x/crm_six-raster-list-setup.RData`


Computing resistance distances
==============================

To find mean hitting times on a fine grid at a given set of parameters, as in `dem-hitting-time.sh`, after running `setup-dem-inference.sh`:
1. Use `make-resistance-distances.R` on a small grid to find initial guesses with analytic inference, e.g.
```
Rscript make-resistance-distances.R ../geolayers/multigrid/512x/crm_ 512x six-raster-list multigrid-six-raster-list.tsv analytic
```
2. Push these up to a finer grid, using `disaggregate-ht.R`, e.g.
```
Rscript disaggregate-ht.R ../geolayers/multigrid/512x/crm_ ../geolayers/multigrid/256x/crm_ 512x 256x six-raster-list 512x/six-raster-list-hitting-times.tsv 2
```
which produces `./100x/500x-six-raster-list-aggregated-hitting-times.tsv`.
3. Use these as a starting point for inference on the finer grid, e.g.
```
Rscript make-resistance-distances.R ../geolayers/multigrid/512x/crm_ 256x six-raster-list multigrid-six-raster-list.tsv numeric 256x/512x-six-raster-list-aggregated-hitting-times.tsv 120
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



