# Changelog

## Version v2.0.0 -- 2026-10-11
Requires StarFormationHistories.jl v2, StellarTracks.jl v1.2, BolometricCorrections.jl v0.2.3, InitialMassFunctions.jl v0.1.7

### Upgrade Instructions
 - The on-disk layout of the BaSTI v1 and v2 tracks has changed. You will get the new versions when you upgrade. If you wish to remove the old versions, you will need to do so manually (the old versions take up about 150 MB total). Unless you have defined `DATADEPS_LOAD_PATH` to overwrite the default, running the following Julia code should delete the old tracks:

```julia
for name in ("BaSTIv1", "BaSTIv2")
    rm(joinpath(DEPOT_PATH[1], "scratchspaces", "124859b0-ceae-595e-8997-d05f6a7a8dfe", "datadeps", name); recursive=true, force=true)
end
```

### New features
 - `fit_ssp` fits the color-magnitude diagram of a single stellar population, such as a star cluster, for its age, metallicity, distance, extinction, binary fraction, and background fraction, with priors on each, and derives its birth stellar mass, with uncertainties from posterior sampling that includes every mode of the posterior and reports the probability of each. `plot_ssp_corner` and `plot_ssp_cmd` plot the results. See `examples/fit_ssp`.
 - `simulate_catalog` samples a mock photometric catalog from a known star formation history, optionally applies observational errors and incompleteness, and optionally fits the catalog with `fit_sfh` to check how well the input is recovered. See `examples/simulate_catalog`.
 - Gates: polygons in color and magnitude (`data.binning.gates`) exclude regions of the Hess diagram, such as those with foreground contamination, from fits.
 - Differential extinction, in the new `differential_extinction` section: a dust layer of total extinction `dAv` in addition to the foreground `properties.Av`, with the scale height of the stars relative to the dust (`scale_height_ratio`), the fraction of stars behind the dust (`f_red`), a log-normal spread of `dAv` across the field (`dAv_sigma`), and MATCH's additional spread for young stars (`dAvy`, tapering between ages `dAvy_t1` and `dAvy_t2`). By default extinctions are uniform from `Av` to `Av + dAv`, as the `-dAv` and `-dAvy` options of MATCH; together the parameters describe the `-diskav` model of MATCH. `dAv` and `scale_height_ratio` can depend on age. `fit_sfh` builds its templates with the extinction distribution of each age. `simulate_catalog` draws each star's extinction from the age-dependent extinction distribution, records it in the catalog, and passes the section to its follow-up fit. `fit_ssp` does not yet support differential extinction.
 - `fit_sfh` fits a background of contaminating stars (e.g., background galaxies or foreground stars) alongside the stellar populations, uniform over the Hess diagram by default or shaped like a control field given in `data.background`; `data.background: none` turns it off (also accepted by `fit_ssp`). The expected number of background stars is reported in the `_mass` output file, and the model Hess diagrams include the background.
 - `fit_sfh` writes a `_convergence` output file reporting, for each combination of stellar track library and bolometric correction grid, whether its MAP and MLE optimizations converged, why they stopped, their numbers of iterations and objective evaluations, the largest final gradient component, the MLE log likelihood, and the time taken; it warns about optimizations that did not converge, including those that stalled with a largest gradient component above 1.
 - AST processing: error, bias, and completeness models are now built from provided AST files in a 2-D format that better captures filter covariances.
 - Observational models can be given for each filter from a tabulated signal-to-noise ratio or magnitude error curve (`data.filters`) instead of from artificial star tests, so the `ASTs` section is optional.
 - The `BinaryMassRatio` binary model is now fully supported, with a configurable distribution of mass ratios (`binaries.mass_ratio`).
 - The age-metallicity relations can be specified by their [M/H] at two lookback times (`metallicity.constraints`).
 - `properties.T_max` sets the lookback time at which star formation begins, i.e., the oldest edge of the age grid, which was previously fixed at 13.7 Gyr.
 - The `_mass` output file of `fit_sfh` gives the present-day stellar mass (`mstar_present`, with `mstar_present_lower` and `mstar_present_upper`), the initial mass of the stars that are still alive, alongside the total stellar mass formed.

### Changes affecting results
 - The total stellar mass formed in the `_mass` output file of `fit_sfh` (`mstar`, `mstar_lower`, `mstar_upper`) is now the 16th, 50th, and 84th percentiles of the total masses of sampled star formation histories, rather than the integrals of the percentiles of the star formation rate in each age bin, which ignored the correlations between age bins. The uncertainties are typically narrower and the median can change slightly.
 - Completeness, bias, and error measured from artificial star tests are now joint functions of the magnitudes in both AST filters, so selection that depends on both filters, such as a signal-to-noise cut in each, is modeled correctly. Both AST filters must be used in the Hess diagram.
 - The photometric errors measured from artificial star tests are now the standard deviation of the residuals for Gaussian errors, about 1.48 times larger than the median absolute deviation used before, which underestimated the errors.
 - The SFH templates change slightly with StarFormationHistories.jl v2, which treats isochrones whose initial masses decrease between evolutionary phases differently and samples isochrones more finely relative to the photometric errors.
 - `fit_sfh` fits a uniform background by default. Results for catalogs without contamination are essentially unchanged, since the fitted background then goes to zero, while contaminated catalogs no longer have their contaminants fit as extra (mostly old) stars.
 - The diagnostics plot shows the completeness used in the fit over the Hess diagram.
 - The `fit_sfh` example moved to `examples/fit_sfh`.

## Version v1.1.0 -- 2026-05-20
Features of the released package:
 - `fit_sfh` measures the star formation history and the age-metallicity or mass-metallicity relation of a resolved stellar population from its photometry and artificial star tests, configured by a YAML file. Templates are computed from isochrones interpolated on the fly from the PARSEC, MIST v1, MIST v2 (including v2.5), BaSTI v1, and BaSTI v2 stellar track libraries and the YBC, MIST v1, and MIST v2 bolometric correction grids. The fit is repeated for every combination of the listed stellar track libraries and bolometric correction grids, and the spread between the results is reported as a systematic uncertainty alongside the random uncertainties.
 - Metallicity models `LinearAMR`, `LogarithmicAMR`, and `PowerLawMZR` with Gaussian metallicity dispersion; the IMFs of InitialMassFunctions.jl; and the `NoBinaries` and `RandomBinaryPairs` binary models.
 - Outputs: a table of the cumulative SFH, star formation rate, and metallicity in each age bin with random and systematic uncertainties, the total stellar mass formed, the observed and best-fit model Hess diagrams, a copy of the configuration, and an optional diagnostics plot of the artificial star tests.
 - `plot_cmd_residuals` and `plot_cumsfh_sys` plot the Hess diagram residuals and the cumulative SFH and metallicity evolution.
