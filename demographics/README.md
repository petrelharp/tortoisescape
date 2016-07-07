# Parameterized models

## Parameters

1. Demographics:

    1. population density
    2. dispersal distance

2. Refugia:

    1. location
    2. radii

3. Contraction
    
    1. end time
    2. speed

4. Re-expansion

    1. begin time
    2. speed

## Implementation

Functions to do this are in [grid-refugia-fns.R](grid-refugia-fns.R):

- `model_setup()` : sets up refugia model, as in `model-setup-msarg.R` (takes about 15 seconds)
- `sim_data()` : simulates data and returns the mean pairwise distance matrix (takes longer)

Running `fit_refugia_params.R` produces many iterations with parameters fairly widely resampled,
and diagnostic plots for each.  Results can be examined with
`parse-runs.R`.


# Ad hoc models

1.  **refugia-demography** : expansion in the north (and south?) that moves into the north and recontacts the south

    - [refugia.ms/] : results

2.  **barrier-demography** : constant lightly leaky barrier between the two

    - [barrier.ms/] : results

## What we're aiming for:

Barrier:
    - 18Kya  (Afton Canyon)
    - leakiness: do calc from panmictic pops

Expansion:
    - 200 Kya ??


## Run the code:

1.  Models are set up in `model-setup-msarg.R`, producing files
    that contain for each scenario `barrier.demog` and `sample.config`:
    
    * [refugia-demography-msarg.RData]
    * [barrier-demography-msarg.RData]

2.  In `run-ms.R`, trees are simulated with `ms`, these are processed, and information is written to e.g. `barrier.ms`.

3.  PNGs of IBD plots can be made by running e.g. `Rscript ms-maps.R barrier.ms`.
