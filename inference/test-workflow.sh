#!/bin/bash

# Simulate hitting times then try to fit parameters

Rscript make-resistance-distances.R ../geolayers/multigrid/256x/crm_ 256x six-raster-list multigrid-six-raster-list.tsv analytic

Rscript sim-hitting-times.R ../geolayers/multigrid/256x/crm_ 256x  six-raster-list six-params.tsv 0.05 256x/six-raster-list-sim-hts.tsv

