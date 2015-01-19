Configuration file syntax
=========================

Configuration files are in JSON.  Some gotchas:
- always put quotes around names
- non-integer numbers need digits before *and* after decimals (`0.01` *not* `.01`)

These can be read in with `read.json.config()`.

Standard parts:

- `layer_names` : names of the layers involved
- `setup_files` : these will be loaded at the start of inference
- `observed_ht_file` : *tsv* file with observed hitting times (or divergences)
- `fitted_ht_file` : write out full hitting times at the end to this file
- `alpha` : controls the tradeoff between interpolation and hitting-time-ness
- `params` : parameter values, as a named list
- `paramscale` : scale on which to optimize over parameters
- `paramtol` : tolerance; stop optimization doesn't change parameters by more than this much
- `maxit` : maximum number of steps in optimization for hitting time interpolation
- `maxstep` : maximum number of ( parameter inference ) <-> ( hitting time interpolation ) steps
