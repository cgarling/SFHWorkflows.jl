# SFHWorkflows.jl
[StarFormationHistories.jl](https://github.com/cgarling/StarFormationHistories.jl) contains highly modular components for measuring resolved star formation histories from high-precision color magnitude diagrams. This package provides standardized workflows for making these measurements that can be configured and ran from simple YAML configuration files.

## Usage
The main entry point for measuring resolved star formation histories from catalogs of photometry and artificial star tests is `fit_sfh`. This function reads a YAML configuration file that defines all relevant parameters for the fit. An example YAML configuration file is given in `examples/config.yaml`. In this documentation, placeholders like `<output.path>` indicates the value given under the `output` section in the `path` variable of the YAML configuration file. All text files written will use the same file extension as provided in `<output.filename>` for consistency.

## Output
### Fit Parameters
`fit_sfh` will write a number of output files containing results of the fit. The main output file, written to `<output.path>/<output.filename>`, will be a whitespace-delimited table with column names given in the first row. An example of this file is given in `examples/output/results.txt`. Each row in the file specifies SFH parameters (e.g., cumulative SFH, SFR, metallicity, etc.) in a bin of logarithmic age (defined as `log10(age [yr])`) defined by the left and right bin edges in the first two columns. 

The next two columns give the lower and upper bounds on the systematic error in the cumulative SFH, defined as the interval that bounds the cumulative SFH for _all_ combinations of stellar track libraries and bolometric corrections grids used in the fit. The next two columns give the lower and upper bounds on the systematic error of the metallicity of stars forming in that logarithmic age bin, defined in the same way. 

The rest of the columns contain SFH parameters for individual combinations of different stellar tracks and bolometric correction grids that were defined in the YAML configuration file. The naming of the columns follows the convention `<quantity>_<stellar track>_<BC grid>` for best-fit quantities, and the lower and upper random uncertainty bounds are given by `<quantity>_lower_<stellar track>_<BC grid>` and `<quantity>_upper_<stellar track>_<BC grid>`, respectively. These are the actual estimates of the values of these parameters 1-σ below and above the best-fit value, _not_ errors on the best-fit value, so the random uncertainty range for the SFR would be from `sfr_lower_<track>_<bc>` to `sfr_upper_<track>_<bc>`, for example. SFRs are in solar masses / yr, MH is logarithmic metallicity (\[M/H\] is the same as \[Fe/H\] for scaled-solar abundance patterns).

### Total Stellar Mass Formed
A plain text file `<output.path>/<output.filename>_mass` will be written containing the total stellar mass formed and 1-σ upper and lower estimates. Values are given for each combination of stellar track library and bolometric correction grid. Masses are in solar masses.

### Hess Diagrams
`fit_sfh` will also write a plain text file containing the observed Hess diagram given the binning scheme defined in the YAML configuration file under the `data.binning` section. An example file is given in `examples/output/results_obshess.txt`. The first row gives the x-axis histogram edges and the second row gives the y-axis histogram edges. The rest of the file is the histogram -- each row will have length equal to the number of y-axis bins minus 1 (`length(ybins) - 1`) and there are a number of columns equal to the number of x-axis bins minus 1 (`length(xbins) - 1`). This layout follows Julia's column-major array layout; you may need to transpose this matrix if you read it into a Python NumPy array, which is row-major. These files can be parsed into instances of `StatsBase.Histogram` with `SFHWorkflows.SFHFitting.read_histogram(<filename>)`.

Plain text files with identical data layouts are also written containing the best-fit model Hess diagram for each combination of stellar track library and bolometric correction grid defined in the YAML configuration file. These files have naming convention `<output.filename>_modelhess_<stellar track>_<BC grid>` and can also be read with `SFHWorkflows.SFHFitting.read_histogram(<filename>)`.

### Plots
A PDF file `<output.path>/diagnostics.pdf` illustrating the photometric error and completeness models measured from the artificial star tests will be written if `plotting.diagnostics: true` in the YAML configuration file.

A convenience function for making a 4-panel Hess diagram (observed Hess, model Hess, observed - model, residual significance) is provided in `SFHWorkflows.SFHFitting.Plotting.plot_cmd_residuals`. An example of its usage is given in `examples/fit_sfh.jl`.

A convenience function for making a 2-panel cumulative SFH and AMR plot is provided in `SFHWorkflows.SFHFitting.Plotting.plot_cumsfh_sys`. An example of its usage is given in `examples/fit_sfh.jl`.
