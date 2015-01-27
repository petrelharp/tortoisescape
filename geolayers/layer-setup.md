We have a new set of layers, and need to set up inference on them.
These are in `geolayers/expanded/expanded-TIFF/` (not really, but symlinked in).
First, we need to make the scaled down layers: in the directory `geolayers/`, run:
```
./aggregate-batch.sh expanded/expanded-TIFF/ expanded/
```
This calls the script `aggregate-layers.R` on each layer in the directory `expanded/expanded-TIFF/`,
aggregating by a factor of 2, and repeating.

Then, we need to put together the consensus masking NA layer (this is in `make-masking-layer.R`):
```
Rscript make-masking-layer.R expanded/2x 32
```


