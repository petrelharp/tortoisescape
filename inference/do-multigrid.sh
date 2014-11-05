#!/bin/bash


#  Rscript make-resistance-distances.R (layer prefix) (subdir) (layer file) (parameter file) (method) [initial guess] [max running time] [output file]

Rscript make-resistance-distances.R ../geolayers/multigrid/512x/crm_ 512x six-raster-list multigrid-six-raster-list.tsv analytic
Rscript disaggregate-ht.R ../geolayers/multigrid/512x/crm_ ../geolayers/multigrid/256x/crm_ 512x 256x six-raster-list 512x/six-raster-list-hitting-times.tsv 2

Rscript make-resistance-distances.R ../geolayers/multigrid/256x/crm_ 256x six-raster-list multigrid-six-raster-list.tsv numeric 256x/512x-aggregated-hitting-times.tsv
Rscript disaggregate-ht.R ../geolayers/multigrid/256x/crm_ ../geolayers/multigrid/128x/crm_ 256x 128x six-raster-list 256x/six-raster-list-hitting-times.tsv 2

Rscript make-resistance-distances.R ../geolayers/multigrid/128x/crm_ 128x six-raster-list multigrid-six-raster-list.tsv numeric 128x/256x-aggregated-hitting-times.tsv
Rscript disaggregate-ht.R ../geolayers/multigrid/128x/crm_ ../geolayers/multigrid/64x/crm_ 128x 64x six-raster-list 128x/six-raster-list-hitting-times.tsv 2

Rscript make-resistance-distances.R ../geolayers/multigrid/64x/crm_ 64x six-raster-list multigrid-six-raster-list.tsv numeric 64x/128x-aggregated-hitting-times.tsv
