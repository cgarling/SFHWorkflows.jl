# Examples

Each example is a directory containing a YAML configuration file (`config.yml`) and a script that runs the workflow on it and makes the standard figures.

 - [`simulate_catalog/`](simulate_catalog): simulates the photometric catalog of a dwarf irregular galaxy with a known star formation history, fits it with `fit_sfh`, and compares the fit to the input. This example needs no external data, so it is the best place to start.
 - [`fit_sfh/`](fit_sfh): measures the star formation history of a galaxy from its photometry and artificial star tests. You need to supply your own photometry and artificial star test files.

## Getting the examples
The installed copy of the package is read-only, so download the examples from this repository, either by cloning it with `git clone https://github.com/cgarling/SFHWorkflows.jl` or by downloading the example directory you want from GitHub. The examples on the `main` branch may use features newer than the latest release; if they do not run with your installed version, use the examples from the tagged release matching your version. Copy the example directory wherever you want to work; results are written to a `results/` subdirectory of the directory you run it from.

## Running an example
First install SFHWorkflows.jl, CairoMakie, and YAML into a Julia project, as described in the [installation instructions](../README.md#installation). Then change into the example directory and start Julia with that project and multiple threads,

```
cd simulate_catalog
julia --project=<path to project> --threads=auto
```

and run the script at the REPL:

```julia
include("simulate_catalog.jl")
```

You can also run the script directly from the command line with `julia --project=<path to project> --threads=auto simulate_catalog.jl`, but running it from the REPL keeps the results (e.g., the simulated catalog and the fit result) available for further exploration, and running it again in the same session is much faster because the code is already compiled.

The first run downloads the stellar track libraries and bolometric correction grids used by the configuration and asks you to confirm some of the downloads; see [Data files](../README.md#data-files).

The scripts assume that they are run from their own directory, because the configuration and output paths in `config.yml` are relative paths.
