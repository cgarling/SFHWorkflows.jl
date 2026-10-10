# Bimodal `fit_ssp` example

This example fits a single stellar population to a catalog that contains two, so that the posterior has two separate modes, and shows how `fit_ssp` reports them. It needs no external data. See [the `fit_ssp` example](../README.md) for how the fit and the posterior sampling work, the configuration options, and the output files, and [the examples README](../../README.md) for how to run the examples.

The catalog combines two simulated populations, both at a distance modulus of 18.5 with Av = 0.1 and without binaries, mock-observed in F475W and F814W with the signal-to-noise ratio curves of the `fit_ssp` example:
 - a younger, metal-rich population, 3.2 Gyr old (logAge 9.5) with [M/H] = -0.3 and 3000 solar masses of stars formed;
 - an older, metal-poor population, 12.6 Gyr old (logAge 10.1) with [M/H] = -1.6. It is simulated with 12000 solar masses of stars formed, and a random subset of 1764 of its detected stars is kept, a number chosen so that the two modes have similar probabilities.

Only stars brighter than F814W = 23.5, where the photometric errors are about 0.01 mag, are fit. A single population can describe either the young or the old stars, and the stars of the other population then go into the uniform background, so the posterior has one mode near each population. The background is needed for this: without it, every star that the population cannot produce makes the fit far worse, and the posterior is set by those stars rather than by the population that fits.

## Files
 - `simulate_young.yml` and `simulate_old.yml`: the `simulate_catalog` configurations of the two populations.
 - `config.yml`: the `fit_ssp` configuration.
 - `bimodal.jl`: simulates both populations, writes the combined catalog of detected stars to `results/simulation/phot.dat`, runs `fit_ssp` on `config.yml`, and makes the figures below.

## Results
In `results/fit/`, besides the output files described in the `fit_ssp` example:
 - `results_modes_PARSEC_YBC.txt` lists the two modes with their probabilities (`weight`). The probabilities are those of the two single-population descriptions of the catalog; neither is the true answer, since the catalog has two populations.
 - `results_cmd_mode1.pdf` and `results_cmd_mode2.pdf` show the combined catalog with the isochrone of each mode (blue) and the true isochrone of the population closest in age (red).
 - `results_corner.pdf` shows the posterior samples of both modes, with the best fit in blue and the location of each mode as an orange circle. The `mode` column of `results_samples_PARSEC_YBC.txt` says which mode each sample was drawn from.

The percentiles in `results.txt` mix both modes, so a median can fall between them where the posterior density is low. For a multimodal posterior, summarize each mode separately, e.g., from the samples with each value of `mode`.
