# `simulate_catalog` example

This example simulates HST ACS/WFC and Roman photometry of a dwarf irregular galaxy similar to Aquarius or Leo A, at a distance modulus of 24.95, with 2 million solar masses of stars formed. Star formation is low before 8 Gyr ago, rises between 8 and 4 Gyr ago, and continues until 100 Myr ago, the youngest age in the grid. The catalog is mock-observed in F475W, F606W, and F814W using tabulated signal-to-noise ratio curves, and then fit with `fit_sfh` in F814W vs. F475W - F814W to check how well the input SFH is recovered. It needs no external data.

See [the examples README](../README.md) for how to download and run the examples.

## Files
 - `config.yml`: the configuration, with every option documented in comments. The `mockobservations.ASTs` section shows how to use artificial star tests instead of signal-to-noise ratio curves for two of the filters.
 - `simulate_catalog.jl`: runs `simulate_catalog` on `config.yml` and makes the figures below.

## Output
`simulate_catalog` writes the catalog (`results/catalog.txt`), the input SFH binned on the age grid (`results/catalog_truth.txt`), a copy of the configuration including the random seed (`results/input.yml`), and the fit to `results/fit/`. See the [main README](../../README.md#simulating-catalogs-simulate_catalog) for details.

The script makes these figures:
 - `results/catalog_cmd.pdf`: color-magnitude diagrams of the true magnitudes and of the mock-observed magnitudes of the detected stars, limited to the Hess diagram used in the fit. Sparse regions are shown as individual stars and dense regions as a density map. Any gates in `fit.binning` are outlined; they exclude regions from the fit only, so the catalog contains all stars.
 - `results/fit/results_hess.pdf`: the simulated and best-fit model Hess diagrams with their residuals.
 - `results/fit/results_cumsfh.pdf`: the fitted cumulative SFH and mean metallicity with the input SFH overlaid.
 - `results/fit/results_sfr.pdf`: the fitted star formation rate in each age bin, with its 16th to 84th percentile range, and the input star formation rate.
 - `results/fit/diagnostics.pdf`: the completeness used in the fit over the Hess diagram.

If `fit.run` is set to `false`, the script instead plots the input cumulative SFH and mean metallicity to `results/truth_cumsfh.pdf`.

## Things to try
 - Change `sfh` to simulate a different star formation history, or `properties.stellar_mass` to change the number of stars.
 - Change `sampling.seed` to draw a different random realization of the same population.
 - Change `fit.binning` to fit a different region of the color-magnitude diagram, or give `fit.stellartracks` or `fit.bolometriccorrections` to fit with different models than were used to simulate the catalog and see the effect of model systematics.
 - Edit the generated `results/fit/input.yml` and rerun the fit alone with `fit_sfh("results/fit/input.yml")`.
