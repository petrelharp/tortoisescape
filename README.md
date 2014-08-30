tortoisescape
=============


Files
=====


`tortoise_gis_notes.R`: figuring out how `raster` works and how much time/memory it takes

`sim-resistances.R`: computes analytic mean hitting times on a grid, plots these, and tries to infer the underlying parameters from these.

`rw-testing-fns.R`: functions for simulating random walks on a graph (not efficient), estimating hitting time distributions, creating generator matrices for random walks on a grid, and computing mean hitting times analytically

`rw-testing.R`: tests underlying math relating inverse of submatrix of generator to mean hitting times, comparing empirically estimated mean to theory. (it works)

`resistance-timing.R`: looks at time it takes to use and solve large sparse grid-based generator matrices by various methods

`resistance-optimization.tex`: calculations.
