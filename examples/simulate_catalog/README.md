# `simulate_catalog` example

This example simulates HST ACS/WFC and Roman photometry of a dwarf irregular galaxy similar to Aquarius or Leo A, at a distance modulus of 24.95, with 2 million solar masses of stars formed. Star formation is low before 8 Gyr ago, rises between 8 and 4 Gyr ago, and continues until 100 Myr ago, the youngest age in the grid. The catalog is mock-observed in F475W, F606W, and F814W using tabulated signal-to-noise ratio curves, and then fit with `fit_sfh` in F814W vs. F475W - F814W to check how well the input SFH is recovered. It needs no external data.

Note that with `dAv` or `dAvy`, `properties.absolute_magnitude` (used instead of `stellar_mass`) is the integrated magnitude of the population at extinction `Av`, before the spread is applied, so the simulated catalog is fainter than `absolute_magnitude`.

See [the examples README](../README.md) for how to download and run the examples.

## Files
 - `config.yml`: the configuration, with every option documented in comments. The `mockobservations.ASTs` section shows how to use artificial star tests instead of signal-to-noise ratio curves for two of the filters.
 - `simulate_catalog.jl`: runs `simulate_catalog` on `config.yml` and makes the figures below.
 - `m31.yml` and `m31_roman.yml`: configurations for an M31 disk region observed with HST and with Roman, described [below](#m31-disk-region-with-differential-extinction). To run one, change `config` in `simulate_catalog.jl` to `"./m31.yml"` or `"./m31_roman.yml"`; their outputs are written to `results_m31/` and `results_m31_roman/`.

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

## M31 disk region with differential extinction
`m31.yml` simulates a ~100 pc region of the M31 disk like those whose recent SFHs were measured from PHAST photometry by [Wainer et al. 2026](https://arxiv.org/abs/2605.12721), and fits it self-consistently. Each star's extinction is drawn uniformly from `Av = 0.4` to `Av + dAv = 1.4`. The SFR declines over the last 500 Myr, following their combined PHAT and PHAST SFH, and is constant from 13.7 Gyr to 1 Gyr ago. This old population gives ~0.8 stars per arcsec² with 21.5 < F814W < 23, the stellar density of their regions with this depth. The mock photometry has 50% completeness at F475W ≈ 27.2 and F814W ≈ 25.9. The fit uses F475W vs. F475W - F814W down to the F475W 50% completeness limit, with a gate excluding F475W - F814W > 1.25 and F475W > 21, which removes the red giant branch and red clump. About 2000 stars fall in the fit region, close to the 1789 in their region 835. The catalog has ~15 times as many stars; most are fainter than the fit's limit or inside the gate. Simulating all ages is cheap, because stars fainter than `sampling.mag_lim` are not stored.

Because the fit region contains mainly upper main sequence and helium burning stars, the SFH is well constrained only over the last ~500 Myr. For this seed, the mass formed in the last 100 Myr is recovered within 6% and in the last 500 Myr within 22%. The oldest bin, from 1.26 Gyr to 13.7 Gyr, is constrained only by the few turnoff stars near the faint limit, so the total stellar mass and the cumulative SFH are not reliable. A fit of this region with 0.1 dex age bins out to 3.2 Gyr did not converge, so `fit.stellartracks` combines all ages older than 1.26 Gyr into one bin.

`m31_roman.yml` observes the same stars (same population and random seed) in Roman F087 and F158, with the same signal-to-noise ratio as a function of apparent magnitude as F475W and F814W in `m31.yml`. Stars are brighter at these longer wavelengths, so the same depth reaches the main sequence turnoffs of populations a few Gyr old, and the red giant branch and red clump are bright and well measured. The full F087 vs. F087 - F158 CMD is fit, with all age bins. For this seed, the fit recovers the SFH and the age-metallicity relation over all ages, and the total stellar mass within 2%. A gate analogous to the HST exclusion region is included in the configuration, commented out; in these filters it does not separate the faint red giant branch from the main sequence, and without the red giant branch the metallicity is poorly constrained.
