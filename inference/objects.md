Code structure
==============

These are the things we work with:

*Landscape* : geographical information for a set of layers, plus sample information:
    + *@baselayer* : a reference `Raster` object that can be used for plotting
    + *@n.layers* : number of layers carried along
    + *@layer.names* : character vector of the layer names
    + *@layers* : numeric matrix of nonmissing values of the layers, of dimensions `n.nomissing , n.layers`
    + *@nonmissing* : vector of length `n.nonmissing` of indices of nonmissing layer values in `values(@baselayer)`
    + *@n.samples* : number of samples
    + *@sample.coords* : `SpatialPoints` object of the locations of the samples
    + *@locs* : indices of the nonmissing cells that `@sample.coords` fall in
    + *@neighborhoods* : list of vectors of indices of nonmissing cells nearby to `@sample.coords`
    + *@neighborhood.radius* : maximal distance from `@sample.coords` to a cell in `@neighborhoods`
    + *@resolution* : numeric resolution, used for adjusting the overall migration rates to be comparable across different resolutions

A `Landscape` has methods:
    + *plot(lname)* : plots the layer whose name is `lname`
    + *n.nonmissing()* : number of nonmissing cell values
    + *plot.samples()* : plot sample IDs on some background

*MigrationModel* : inherits from *Landscape*, plus migration model based on it
    + *@beta* : an overall scaling rate to the migration (nonzero)
    + *@gamma* : a vector of length `n.layers` of coefficients that determine the stationary distribution
    + *@delta* : a vector of length `n.layers` of coefficients that determine the relative migration rates
    + *@transfn* : transformation function to take sums of landscape values to rates
    + *@G* : generator matrix for the random walk
    + *@Gjj* : 1-based indices of columns of `@G@x`

A `MigrationModel` has methods:
    + *coef()* : returns the vector `c(@beta,@gamma,@delta)`
    + *hitting.times([locs=@locs])* : returns a `Landscape` object of mean hitting times to `locs` (whic is a vector of nonmissing entries as above)
    + *update(beta,gamma,delta)* : updates the parameters and the nonzero entries of the generator matrix `@G`
    + *gamma.layer([gamma=@gamma])* : returns a `Raster` whose values are the linear combination produced by `gamma` that determines the stationary distribution
    + *delta.layer([delta=@delta])* : returns a `Raster` whose values are the linear combination produced by `delta` that determines the relative migration rates
    + *stationary([beta=@beta],[gamma=@gamma],[delta=@delta])* : returns a `Raster` whose values are proportional to the stationary distribution
    + *migration([beta=@beta],[gamma=@gamma],[delta=@delta])* : returns a `Raster` whose values are the total migration rate out of that cell

*MigrationModelFit* : inherits from *MigrationModel*, plus observed divergences and methods to fit the parameters
    + *@observed.hts* : `n.samples * n.samples` matrix of "observed" hitting times we are trying to fit
    + *@fitted.hts* : `n.nonmissing * n.layers` matrix of fitted hitting times, interpolated across the entire landscape
    + *@interp.objective* : list of information needed to interpolate `@observed.hts`, given the parameters (objective function, gradient, shared environment)
    + *@param.objective* : list of information needed to find parameters that best fit `@fitted.hts` (objective function, gradient, shared environment)
