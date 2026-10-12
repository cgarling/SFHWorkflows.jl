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

## Binaries
`binaries.model` sets how unresolved binary stars are modeled, and `binaries.binary_fraction` is the fraction of stellar systems that are binaries. Like the extinction, the binary model is fixed, not fit.
 - `NoBinaries`: every system is a single star.
 - `RandomBinaryPairs`: a binary is two stars drawn independently from the IMF.
 - `BinaryMassRatio`: the IMF gives the total mass of each system, so a system IMF such as `Chabrier2003System` is the natural choice. A binary splits its mass between two stars with a secondary-to-primary mass ratio drawn from `binaries.mass_ratio`, a distribution written as in [Distributions.jl](https://juliastats.org/Distributions.jl/stable/univariate/), e.g., `Beta(2, 2)` (default `Uniform(0.1, 1)`). Distributions extending outside [0, 1], e.g., `Normal(0.5, 0.25)`, are truncated to [0, 1].

## Differential extinction
By default every star has extinction `properties.Av`, e.g., from the Milky Way foreground. The optional `differential_extinction` section adds a layer of dust within the galaxy, so that each star also receives the extinction of the part of the layer in front of it:

| Key | Default | Meaning |
|---|---|---|
| `dAv` | 0 | Total extinction through the dust layer. With the other defaults, extinctions are uniform from `Av` to `Av + dAv`, as with the `-dAv` option of MATCH. |
| `scale_height_ratio` | 1 | Scale height of the stars divided by that of the dust, both with exponential vertical profiles. At 1, stars and dust are perfectly mixed, giving the uniform distribution above. Below 1, the stars are embedded in the dust and most have about half of `dAv`. Above 1, most stars are either in front of or behind the dust, approaching a thin dust layer, with extinctions near `Av` and `Av + dAv`. |
| `f_red` | 0.5 | Fraction of the stars behind the midplane of the dust. |
| `dAv_sigma` | 0 | Width of a log-normal spread of the total extinction of the layer across the observed field, in natural log; `dAv` is then its median. It is truncated 3σ above the median, so extinctions reach `Av + dAv exp(3 dAv_sigma)`. |
| `dAvy` | 0 | An additional, independent uniform spread from 0 to `dAvy` for young stars, as the `-dAvy` option of MATCH (MATCH's default is 0.5). It is full for ages below `dAvy_t1` and tapers linearly to zero at `dAvy_t2`. |
| `dAvy_t1`, `dAvy_t2` | 0.04, 0.1 | Ages [Gyr] bounding the taper of `dAvy`, as in MATCH. |

`dAv` and `scale_height_ratio` can also depend on age, given as a list of `[age_Gyr, value]` pairs, interpolated linearly in log age and constant before the first and after the last age. For example, young stars embedded in the dust and older stars in a thicker disk than the dust:
```yaml
differential_extinction:
  dAv: 1.0
  scale_height_ratio: [[0.09, 1.0], [0.11, 20.0]]
```
A table of `dAv` with `scale_height_ratio: 1` instead gives uniform distributions whose width changes with age, e.g., for young stars that lie in dustier parts of the field than old stars. A varying scale height ratio changes the shape of the distribution of extinctions but not its range or (for `f_red = 0.5`) its mean, while a varying `dAv` changes its range and mean.

These are the parameters of the `-diskav` model of MATCH as described by [Williams et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJ...846..145W): `Av` is its foreground extinction, `dAv`, `dAv_sigma`, and `f_red` are the median and width of the log-normal extinction and the reddened fraction that [Dalcanton et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015ApJ...814....3D) measured in M31, and a two-entry table of `scale_height_ratio` gives its scale height ratios of young and old stars and their transition age (Williams et al. used 20 for old stars and a transition at 0.1 Gyr). MATCH's exact functional form is not published, so this is a physically consistent realization of that model, not necessarily identical to it.

The parameters are fixed, not fit; to choose them, compare the best fit over a few values. The templates follow each star along its reddening vector, so the largest extinction must be within the range of every bolometric correction grid. Templates with differential extinction take several times longer to build than those without, and those with `dAv_sigma` (which also cover a wider range of extinction) about twice as long again.

Fitting every combination of stellar track library and bolometric correction grid takes a while. For a quicker first run, comment out all but one entry under `stellartracks` and `bolometriccorrections`.

## Output
`fit_sfh` writes its results to `results/` (see the [main README](../../README.md#measuring-star-formation-histories-fit_sfh)), including `results/diagnostics.pdf`, which shows the photometric errors and completeness measured from the artificial star tests and the completeness used in the fit over the Hess diagram.

The script additionally makes these figures:
 - `results/results_hess.pdf`: the observed and best-fit model Hess diagrams with their residuals, for the first stellar track library and bolometric correction grid. Change `idx` in the script to plot a different combination.
 - `results/results_cumsfh.pdf`: the cumulative SFH and mean metallicity, with random uncertainties for the same combination and systematic uncertainties across all combinations.

The script sets `galaxy_name = "Aquarius"` for the figure labels; change it for your galaxy.
