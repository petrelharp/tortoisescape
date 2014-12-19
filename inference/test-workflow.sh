#!/bin/bash

# Stuff to test inference with (in progress)
Rscript make-resistance-distances.R ../geolayers/multigrid/256x/crm_ 256x six-raster-list test01/six-params.tsv analytic test01/256x/six-raster-list-hitting-times-full.tsv 

# fit?
Rscript fit-logistic-model.R ../geolayers/multigrid/256x/crm_ 256x six-raster-list test01/256x/six-raster-list-hitting-times-full.tsv test01/six-params.tsv test01/inferred-six-params.tsv
