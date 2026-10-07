# `fit_sfh` example

This example measures the star formation history of the Aquarius dwarf irregular galaxy from HST ACS/WFC photometry in F475W and F814W, using artificial star tests to model the photometric errors and completeness. The fit is run for every combination of four stellar track libraries (PARSEC, MIST v1.2, MIST v2.5, and BaSTI v2) and three bolometric correction grids (YBC, MIST v1.2, and MIST v2.5), and the spread between these results is reported as a systematic uncertainty.

The photometry and artificial star test files are not included, so to run this example you need to supply your own and edit `config.yml` to match (see below). If you do not have data at hand, start with the [`simulate_catalog` example](../simulate_catalog), which generates a mock catalog and fits it with `fit_sfh`.

See [the examples README](../README.md) for how to download and run the examples.

## Files
 - `config.yml`: the configuration, with every option documented in comments.
 - `fit_sfh.jl`: runs `fit_sfh` on `config.yml` and makes the figures below.
 - `output/`: example output files, described in the [main README](../../README.md#measuring-star-formation-histories-fit_sfh).

## Using your own data
Edit the `data` section of `config.yml`:
 - `path`: the directory containing your photometry and artificial star test files.
 - `photometry`: the photometry file, a whitespace-delimited text file with one column of apparent magnitudes per filter, and the filter name of each column.
 - `ASTs`: the artificial star test file, with columns (input magnitude 1, input magnitude 2, output - input magnitude 1, output - input magnitude 2), and the value of the output - input columns that indicates a non-detection. Alternatively, comment out `ASTs` and give a tabulated signal-to-noise ratio or magnitude error curve for each filter under `filters`.
 - `binning`: the filters and the range and bin sizes of the Hess diagram to fit, and optionally `gates`, polygons in color and magnitude whose bins are excluded from the fit.
 - `background` (optional): the shape of the background of contaminating stars, from the photometry of a nearby control field or a Hess diagram file. The default is a uniform background, which handles contamination spread over the Hess diagram; give a control field if the contamination is concentrated, e.g., in faint background galaxies (see the [main README](../../README.md#background)).

Then update `properties` (distance modulus, extinction, and an approximate stellar mass) for your galaxy, and `bolometriccorrections` for your filter system.

## Differential extinction
By default every star has extinction `properties.Av`. Set `properties.dAv` to spread the extinction of all stars uniformly from `Av` to `Av + dAv`, as the `-dAv` option of MATCH does, and `properties.dAvy` for an additional, independent uniform spread from 0 to `dAvy` for young stars, as `-dAvy` does (MATCH's default is 0.5). The young-star spread is full for ages below `properties.dAvy_t1` and tapers linearly to zero at `properties.dAvy_t2` (in Gyr; defaults 0.04 and 0.1, as in MATCH). `dAv` and `dAvy` are fixed, not fit; to choose them, compare the best fit over a few values. The templates follow each star along its reddening vector, so the extinctions `Av + dAv + dAvy` must be within the range of every bolometric correction grid. Templates with differential extinction take several times longer to build than those without.

Fitting every combination of stellar track library and bolometric correction grid takes a while. For a quicker first run, comment out all but one entry under `stellartracks` and `bolometriccorrections`.

## Output
`fit_sfh` writes its results to `results/` (see the [main README](../../README.md#measuring-star-formation-histories-fit_sfh)), including `results/diagnostics.pdf`, which shows the photometric errors and completeness measured from the artificial star tests and the completeness used in the fit over the Hess diagram.

The script additionally makes these figures:
 - `results/results_hess.pdf`: the observed and best-fit model Hess diagrams with their residuals, for the first stellar track library and bolometric correction grid. Change `idx` in the script to plot a different combination.
 - `results/results_cumsfh.pdf`: the cumulative SFH and mean metallicity, with random uncertainties for the same combination and systematic uncertainties across all combinations.

The script sets `galaxy_name = "Aquarius"` for the figure labels; change it for your galaxy.
