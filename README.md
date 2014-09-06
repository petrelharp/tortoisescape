tortoisescape
=============


Files
=====


`tortoise_gis_notes.R`: figuring out how `raster` works and how much time/memory it takes

`resistance-timing.R`: looks at time it takes to use and solve large sparse grid-based generator matrices by various methods

`resistance-optimization.tex`: calculations.

`interp-inference.R`: does the interpolation-estimation method; tested in `test-interp-inference.R`

`resistance-fns.R`: functions for estimating hitting time distributions, creating generator matrices for random walks on a grid, and computing mean hitting times analytically

`sim-resistances.R`: computes analytic mean hitting times on a grid, plots these, and some forays into estimation

`gradient-ascent.R`: implements the gradient ascent method (inferring a given full hitting times) from the latex document.

`interp-resistances.R`: given correct alpha, interpolate resistance maps from observed values.

`rw-testing.R`: tests underlying math relating inverse of submatrix of generator to mean hitting times, comparing empirically estimated mean to theory. (it works)

`rw-testing-fns.R`: functions for simulating random walks on a graph (not efficient)
