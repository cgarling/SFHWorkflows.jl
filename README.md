# SFHWorkflows.jl
[StarFormationHistories.jl](https://github.com/cgarling/StarFormationHistories.jl) contains highly modular components for measuring resolved star formation histories from high-precision color magnitude diagrams. This package provides standardized workflows for making these measurements that can be configured and run from simple YAML configuration files:

 - [`fit_sfh`](#measuring-star-formation-histories-fit_sfh) measures the star formation history (SFH) of a resolved stellar population from a catalog of photometry and a model of the observational errors and completeness.
 - [`fit_ssp`](#fitting-single-stellar-populations-fit_ssp) fits the age, metallicity, distance, extinction, and binary fraction of a single stellar population, such as a star cluster, and its birth stellar mass, with uncertainties from posterior sampling.
 - [`simulate_catalog`](#simulating-catalogs-simulate_catalog) samples a mock photometric catalog from a known SFH, optionally applies observational errors and incompleteness, and optionally fits the result with `fit_sfh` so you can check how well the input SFH is recovered.

Users may also be interested in [BolometricCorrections.jl](https://github.com/cgarling/BolometricCorrections.jl) and [StellarTracks.jl](https://github.com/cgarling/StellarTracks.jl), which are used to interpolate isochrones on the fly, and [InitialMassFunctions.jl](https://github.com/cgarling/InitialMassFunctions.jl), which provides the initial mass function models.

## Contents
 - [Installation](#installation): [Julia](#julia), [Projects](#projects), [Package](#package), [Data files](#data-files), [Performance](#performance)
 - [Examples](#examples)
 - [Measuring star formation histories: `fit_sfh`](#measuring-star-formation-histories-fit_sfh)
   - [Fit Parameters](#fit-parameters)
   - [Total Stellar Mass Formed](#total-stellar-mass-formed)
   - [Hess Diagrams](#hess-diagrams)
   - [Plots](#plots)
 - [Fitting single stellar populations: `fit_ssp`](#fitting-single-stellar-populations-fit_ssp)
   - [`fit_ssp` outputs](#fit_ssp-outputs)
   - [`fit_ssp` plots](#fit_ssp-plots)
 - [Simulating catalogs: `simulate_catalog`](#simulating-catalogs-simulate_catalog)

## Installation

### Julia
If you need to install Julia, it is recommended to do so through [Juliaup](https://github.com/julialang/juliaup), which is Julia's version multiplexer (similar to how `uv` lets you manage multiple Python versions). On Mac, Linux, and FreeBSD, install it from the command line with

```
curl -fsSL https://install.julialang.org | sh
```

and on Windows with

```
winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore
```

See the [Juliaup repository](https://github.com/julialang/juliaup) for other installation options.

Once `juliaup` is installed it should install the most recent stable version. You can install other versions if needed with `juliaup add <X.X.X>` such as `juliaup add 1.13.0`. To set a default version that will be used when you launch `julia`, set `juliaup default <X.X.X>`. To launch a non-default version of julia, execute, for example, `julia +1.13.0`. SFHWorkflows.jl requires Julia 1.10 or newer.

### Projects
Julia has built-in "projects" which serve the same purpose as virtual environments in Python, allowing for self-contained dependency resolution on a per-project basis. When you start Julia as `julia`, it launches into the **global** environment. Launching with `julia --project=my_project/` will start Julia using the Project.toml dependency file located in the directory `my_project` (e.g., `julia --project=.` to use the Project.toml in the current directory, or make a new one if none exists). The Project.toml is the user-facing file where project dependencies are recorded and can have loose version specifiers for package dependencies, or no version specifiers at all. The Manifest.toml is computer-generated and records the **exact versions** of the dependencies that actually got installed. Only the Manifest.toml guarantees an **exactly** reproducible environment, but in practice Project.toml is what people typically host and share.

We recommend installing SFHWorkflows.jl into specific projects where you need it rather than into the global environment. Packages installed in a project are only available when Julia is started with that project, so you will always have to start Julia with `julia --project=<path to project>`.

### Package
Start Julia (preferably pointing it to a project via the `julia --project=<path to project>` syntax) and enter the package manager by hitting `]` at the REPL prompt. You will see the prompt change from `julia> ` to `(my_project) pkg>` or similar, indicating you are in the package manager. Typing `help` here lists the available operations; the most common are `add` to add a package, `rm` to remove a package, and `update` to update installed packages.

SFHWorkflows.jl is registered in the Julia General Registry. At the `pkg>` prompt, enter

```
add SFHWorkflows CairoMakie YAML
```

CairoMakie and YAML are not needed to call `fit_sfh`, `fit_ssp`, or `simulate_catalog`, but the example scripts use them to make figures and read the configuration files. All of SFHWorkflows.jl's dependencies will be resolved and installed alongside it.

Return to the main Julia REPL prompt by hitting backspace when at the root level of the package manager prompt. Type `using SFHWorkflows` at the REPL prompt and it should complete successfully. You are ready to use SFHWorkflows.jl in your project.

### Data files
The stellar track libraries and bolometric correction grids are not bundled with the packages. Each one is downloaded the first time it is used and cached in your Julia depot (usually `~/.julia/scratchspaces`), so later uses do not download it again. Most downloads ask for confirmation at the REPL; answer `y`. To accept these prompts automatically (e.g., when running scripts non-interactively), set the environment variable `DATADEPS_ALWAYS_ACCEPT=true`.

If you want to download data before running a workflow (e.g., before going offline), you can trigger the downloads directly. For example, the data used by the `simulate_catalog` example can be downloaded with

```julia
using StellarTracks, BolometricCorrections # add these packages first if they are not in your project
PARSECLibrary()      # PARSEC stellar tracks
YBCGrid("acs_wfc")   # YBC bolometric corrections for HST ACS/WFC
YBCGrid("Roman2024") # YBC bolometric corrections for Roman
```

The other libraries are available as `MISTv1Library`, `MISTv2Library`, and `BaSTIv2Library` (see the [StellarTracks.jl documentation](https://cgarling.github.io/StellarTracks.jl/stable/)) and `MISTv1BCGrid` and `MISTv2BCGrid` (see the [BolometricCorrections.jl documentation](https://cgarling.github.io/BolometricCorrections.jl/stable/)).

### Performance
The workflows use Julia's multithreading. Julia starts with the number of threads given by the `--threads` command line option or the `JULIA_NUM_THREADS` environment variable; start Julia with `julia --threads=auto` to use all available cores, or set `export JULIA_NUM_THREADS=8` (for example) in your `~/.bashrc` or equivalent. You can check the number of threads Julia is using by running `Threads.nthreads()` at the REPL.

Faster linear algebra (BLAS) libraries are available for some systems. On Apple silicon with macOS 13.4 or later, `add AppleAccelerate` and load it with `using AppleAccelerate`; on Intel and AMD processors, the same applies to `MKL`. To load one of these in every Julia session, add the `using` line to `~/.julia/config/startup.jl` (create the file if it does not exist). You can check which BLAS library is in use with `import LinearAlgebra: BLAS; BLAS.get_config()`.

## Examples
Runnable examples for each workflow, with configuration files and scripts that make the standard figures, are in [`examples/`](examples). The [`simulate_catalog` example](examples/simulate_catalog) and the [`fit_ssp` example](examples/fit_ssp) run without any external data; the `simulate_catalog` example is a good place to start.

## Measuring star formation histories: `fit_sfh`
`fit_sfh` measures the SFH from a catalog of photometry. It reads a YAML configuration file that defines all relevant parameters for the fit; an example with every option documented is given in [`examples/fit_sfh/config.yml`](examples/fit_sfh/config.yml).

```julia
using SFHWorkflows
result, h = fit_sfh("config.yml") # `h` is the observed Hess diagram as a StatsBase.Histogram
```

The configuration defines
 - `data`: the photometry file, the Hess diagram binning, and the observational model for each filter used in the fit. Optional polygons in color and magnitude (`data.binning.gates`) exclude regions of the Hess diagram, such as those with foreground contamination, from the fit. The observational model (photometric error, bias, and completeness) is measured either from artificial star tests (`data.ASTs`) or from a tabulated signal-to-noise ratio or magnitude error curve for each filter (`data.filters`).
 - `stellartracks` and `bolometriccorrections`: the stellar track libraries and bolometric correction grids, and the grid of ages and metallicities the SFH is measured on. The SFH is measured once for every combination of stellar track library and bolometric correction grid, and the spread between these results is reported as a systematic uncertainty.
 - `imf`, `binaries`, `properties` (distance, extinction, and an approximate stellar mass), and `metallicity` (the age-metallicity or mass-metallicity relation, with initial guesses for its parameters).
 - `plotting` and `output`.

In this documentation, placeholders like `<output.path>` indicate the value given under the `output` section in the `path` variable of the YAML configuration file. All text files written will use the same file extension as provided in `<output.filename>` for consistency.

### Fit Parameters
`fit_sfh` will write a number of output files containing results of the fit. The main output file, written to `<output.path>/<output.filename>`, will be a whitespace-delimited table with column names given in the first row. An example of this file is given in `examples/fit_sfh/output/results.txt`. Each row in the file specifies SFH parameters (e.g., cumulative SFH, SFR, metallicity, etc.) in a bin of logarithmic age (defined as `log10(age [yr])`) defined by the left and right bin edges in the first two columns. 

The next two columns give the lower and upper bounds on the systematic error in the cumulative SFH, defined as the interval that bounds the cumulative SFH for _all_ combinations of stellar track libraries and bolometric corrections grids used in the fit. The next two columns give the lower and upper bounds on the systematic error of the metallicity of stars forming in that logarithmic age bin, defined in the same way. 

The rest of the columns contain SFH parameters for individual combinations of different stellar tracks and bolometric correction grids that were defined in the YAML configuration file. The naming of the columns follows the convention `<quantity>_<stellar track>_<BC grid>` for best-fit quantities, and the lower and upper random uncertainty bounds are given by `<quantity>_lower_<stellar track>_<BC grid>` and `<quantity>_upper_<stellar track>_<BC grid>`, respectively. These are the actual estimates of the values of these parameters 1-σ below and above the best-fit value, _not_ errors on the best-fit value, so the random uncertainty range for the SFR would be from `sfr_lower_<track>_<bc>` to `sfr_upper_<track>_<bc>`, for example. SFRs are in solar masses / yr, MH is logarithmic metallicity (\[M/H\] is the same as \[Fe/H\] for scaled-solar abundance patterns).

### Total Stellar Mass Formed
A plain text file `<output.path>/<output.filename>_mass` will be written containing the total stellar mass formed and 1-σ upper and lower estimates. Values are given for each combination of stellar track library and bolometric correction grid. Masses are in solar masses. An example file is given in `examples/fit_sfh/output/results_mass.txt`.

### Hess Diagrams
`fit_sfh` will also write a plain text file containing the observed Hess diagram given the binning scheme defined in the YAML configuration file under the `data.binning` section. An example file is given in `examples/fit_sfh/output/results_obshess.txt`. The first row gives the x-axis histogram edges and the second row gives the y-axis histogram edges. The rest of the file is the histogram -- each row will have length equal to the number of y-axis bins minus 1 (`length(ybins) - 1`) and there are a number of columns equal to the number of x-axis bins minus 1 (`length(xbins) - 1`). This layout follows Julia's column-major array layout; you may need to transpose this matrix if you read it into a Python NumPy array or want to plot it with matplotlib, as these expect row-major matrices. These files can be parsed into instances of `StatsBase.Histogram` with `SFHWorkflows.SFHFitting.read_histogram(<filename>)`.

Plain text files with identical data layouts are also written containing the best-fit model Hess diagram for each combination of stellar track library and bolometric correction grid defined in the YAML configuration file. These files have naming convention `<output.filename>_modelhess_<stellar track>_<BC grid>` and can also be read with `SFHWorkflows.SFHFitting.read_histogram(<filename>)`.

### Plots
All plotting uses [Makie.jl](https://github.com/MakieOrg/Makie.jl).

If `plotting.diagnostics: true` in the YAML configuration file, a PDF file `<output.path>/diagnostics.pdf` is written illustrating the observational model used in the fit. When artificial star tests are given, it shows the photometric residuals, errors, and completeness measured from them in each filter. It also shows the completeness used in the fit over the Hess diagram, when the Hess diagram's magnitude is one of the two filters in its color.

A convenience function for making a 4-panel Hess diagram (observed Hess, model Hess, observed - model, residual significance) is provided in `SFHWorkflows.SFHFitting.Plotting.plot_cmd_residuals`. An example of its usage is given in `examples/fit_sfh/fit_sfh.jl`.

A convenience function for making a 2-panel cumulative SFH and AMR plot is provided in `SFHWorkflows.SFHFitting.Plotting.plot_cumsfh_sys`. An example of its usage is given in `examples/fit_sfh/fit_sfh.jl`.

## Fitting single stellar populations: `fit_ssp`
`fit_ssp` fits the color-magnitude diagram of a single stellar population (one age and one metallicity), such as a star cluster, plus a population of background stars. The parameters are the age, metallicity, distance modulus, V-band extinction, binary fraction, and the fraction of observed stars in the background; the stellar mass formed in the population is derived from them. Unlike `fit_sfh`, which fits a fixed grid of templates, `fit_ssp` computes the model Hess diagram for any values of the parameters, finds the best fit by optimization, and samples the posterior with Markov chain Monte Carlo. It reads a YAML configuration file; an example is given in [`examples/fit_ssp/config.yml`](examples/fit_ssp/config.yml), and the [example README](examples/fit_ssp/README.md) describes the fitting method, the prior distributions, and the degeneracies between parameters in more detail.

```julia
using SFHWorkflows
results, h = fit_ssp("config.yml") # `h` is the observed Hess diagram as a StatsBase.Histogram
```

The configuration defines
 - `data`: the photometry, Hess diagram binning, gates, and observational model, as for `fit_sfh`. The optional `data.background` gives the shape of the background over the Hess diagram from the photometry of a nearby field or a Hess diagram file; without it, the background is uniform.
 - `stellartracks`, `bolometriccorrections`, `imf`, and `binaries`, as for `fit_sfh`, except that no grid of ages and metallicities is needed and the binary fraction is given in `parameters`.
 - `parameters`: each of `logAge` (or `age` in Gyr), `MH`, `distance_modulus` (or `distance` in parsecs), `Av`, `binary_fraction`, and `background_fraction` is either a number, which fixes it, or a prior distribution written as in [Distributions.jl](https://juliastats.org/Distributions.jl/stable/univariate/), e.g., `Normal(24.95, 0.1)`, which makes it a free parameter.
 - `fit` (the starting grid in age and metallicity for the optimization) and `sampling` (the number of walkers and steps of the sampler and the random seed, or `run: false` to stop after the best fit).
 - `plotting` and `output`.

The fit is run for every combination of stellar track library and bolometric correction grid, in parallel over threads. `results` is a matrix indexed by stellar track library first and bolometric correction grid second, like `result.results` from `fit_sfh`; each entry has fields `model`, `fit` (the best fit), and `chain` (the posterior samples as an `MCMCChains.Chains`, or `nothing` if sampling was not run).

### `fit_ssp` outputs
`fit_ssp` writes the following files to `<output.path>`, where `<base>` and `<ext>` are the base name and extension of `<output.filename>` and `<tracks>_<bcs>` names a combination of stellar track library and bolometric correction grid (e.g., `PARSEC_YBC`):
 - `<output.filename>`: a whitespace-delimited summary table with one row per combination. It gives the best fit (`<parameter>_best`, the highest posterior density found by the optimizer or the sampler), its log posterior density (`lp_best`), and the stellar mass formed at the best fit (`mass_best`), followed by the 16th, 50th, and 84th percentiles of the posterior (`<quantity>_lower`, `<quantity>_median`, and `<quantity>_upper`) of each free parameter, the age in Gyr (`age_Gyr`), the stellar mass formed in solar masses (`mass`), and the expected number of background stars (`background_stars`).
 - `<base>_chain_<tracks>_<bcs><ext>`: the posterior samples, with one row per step of each walker after burn-in, giving the free parameters, `mass`, `background_stars`, and the log posterior density (`lp`).
 - `<base>_obshess<ext>` and `<base>_modelhess_<tracks>_<bcs><ext>`: the observed and best-fit model Hess diagrams, in the same format as the [`fit_sfh` Hess diagram files](#hess-diagrams).
 - `input.yml`: a copy of the configuration, including the random seed used for sampling.
 - `diagnostics.pdf`, if `plotting.diagnostics: true`, as for `fit_sfh`.

### `fit_ssp` plots
Two plotting functions are provided in `SFHWorkflows.SSPFitting`. Each takes an entry of `results`, saves the figure, and returns the Makie figure and axes for further changes. Their optional `truth` keyword marks known values, e.g., for a simulated cluster.
 - `plot_ssp_corner(result, output_file; truth=nothing)` makes a corner plot of the posterior samples. It can also be made later from a chain file, with the best fit read from the summary table by `read_best`.
 - `plot_ssp_cmd(result, color, mag, output_file; truth=nothing)` plots the observed color-magnitude diagram with the best-fit isochrone.

The 4-panel Hess diagram plot of `fit_sfh`, `SFHWorkflows.SFHFitting.Plotting.plot_cmd_residuals`, also accepts a model Hess diagram matrix. Examples of all three are given in [`examples/fit_ssp/fit_ssp.jl`](examples/fit_ssp/fit_ssp.jl).

## Simulating catalogs: `simulate_catalog`
`simulate_catalog` samples a mock photometric catalog from a model stellar population with a known SFH. It reads a YAML configuration file; an example with every option documented is given in [`examples/simulate_catalog/config.yml`](examples/simulate_catalog/config.yml). Sections shared with the `fit_sfh` configuration (`imf`, `binaries`, `bolometriccorrections`, `metallicity`, `output`) use the same keys where possible.

```julia
using SFHWorkflows
catalog, truth, fit = simulate_catalog("config.yml")
```

The configuration defines
 - `stellartrack` and `bolometriccorrections`: one stellar track library, the grid of ages and metallicities used to discretize the SFH, and the bolometric correction grids. The catalog contains magnitudes in every filter of every listed bolometric correction grid (or a selected subset).
 - `properties`: the total birth stellar mass (or, alternatively, the integrated absolute magnitude in one filter), the distance, and the extinction.
 - `sfh` and `metallicity`: the input SFH, either a constant star formation rate or a tabulated cumulative SFH, and the age-metallicity or mass-metallicity relation.
 - `imf`, `binaries`, and `sampling` (random seed and an optional magnitude limit).
 - `mockobservations` (optional): observational models, from artificial star tests or tabulated signal-to-noise ratio curves, used to add photometric errors and bias and to draw which stars are detected.
 - `fit` (optional): settings for fitting the simulated catalog with `fit_sfh`.

`simulate_catalog` writes the following files to `<output.path>`:
 - `<output.filename>`: the catalog. Each row is one star or unresolved binary system that is still alive (and brighter than the magnitude limit, if one is given), with its initial mass(es), age, metallicity, and true apparent magnitudes. Mock-observed filters have an extra `<filter>_obs` column with the observed magnitude, or `NaN` if the star was not detected.
 - `<base>_truth<ext>`, where `<base>` and `<ext>` are the base name and extension of `<output.filename>`: the input SFH binned on the age grid, with the same columns as the `fit_sfh` results (`sfr`, `cum_sfh`, `MH`) so the two can be compared directly.
 - `input.yml`: a copy of the configuration, including the random seed, so that the catalog can be regenerated exactly.
 - `fit/` (if `fit.run: true`): the generated `fit_sfh` configuration (`input.yml`), the photometry of the detected stars (`phot.dat`), and all `fit_sfh` outputs. The generated configuration can be edited and rerun with `fit_sfh` directly.

It returns the catalog and the truth table as `TypedTables.Table`s, and `fit`, which is `nothing` unless `fit.run: true`, in which case it is `(result, h)` as returned by `fit_sfh`.
