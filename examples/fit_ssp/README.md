# `fit_ssp` example

This example simulates a star cluster with `simulate_catalog` and fits its color-magnitude diagram with `fit_ssp` to check how well the true properties are recovered. The cluster is a single stellar population 5 Gyr old with [M/H] = -1.2, 5 × 10<sup>5</sup> solar masses of stars formed, a binary fraction of 0.3, Av = 0.1, and a distance modulus of 24.95. It is mock-observed in F475W and F814W with the same signal-to-noise ratio curves as the `simulate_catalog` example and fit in F814W vs. F475W - F814W. It needs no external data.

See [the examples README](../README.md) for how to download and run the examples.

## Files
 - `simulate.yml`: the `simulate_catalog` configuration for the cluster. See the [`simulate_catalog` example](../simulate_catalog) for every option.
 - `config.yml`: the `fit_ssp` configuration.
 - `fit_ssp.jl`: simulates the cluster, writes the detected stars to `results/simulation/phot.dat`, runs `fit_ssp` on `config.yml`, and makes the figures below.
 - `output/`: example output files of the fit (`results.txt` and the mode and sample files), described under [Output](#output).
 - [`bimodal/`](bimodal): a second example, which fits a catalog of two populations whose posterior has two modes.

## How the fit works
`fit_ssp` fits the observed Hess diagram with the model Hess diagram of a single stellar population plus a background, for every combination of the stellar track libraries and bolometric correction grids in the configuration. The parameters are the age (`logAge`), metallicity (`MH`), distance modulus (`distance_modulus`), V-band extinction (`Av`), the binary fraction (`binary_fraction`), and the fraction of the observed stars in the background (`background_fraction`). The number of stars in the model is not a fit parameter: for the best-fit shape it is set by the number of observed stars, and the birth stellar mass of the cluster is computed from it in closed form for each posterior sample. The background stars are spread uniformly over the Hess diagram unless a field is given in `data.background`.

The fit has two steps:
 1. The best fit (the maximum of the posterior density) is found in two stages. First, the posterior density is evaluated at `fit.ncandidates` candidate starting points spread evenly over the free age, metallicity, and distance modulus within the bulk of their priors, with the extinction at its prior median and the binary and background fractions optimized at each. Then all free parameters are optimized from the best candidate on each separate peak of the posterior density, at most `fit.nstarts` of them (a warning says if more peaks were found). Every one of these local solutions, not only the best, is passed to the sampler, so it can sample each mode of the posterior. Fits for different combinations of stellar tracks and bolometric corrections run in parallel on threads. Because the extinction is not searched, give it a prior no wider than about 0.2 mag (e.g., from foreground dust maps) if you can.
 2. The posterior is sampled with `sample_ssp` from StarFormationHistories.jl. Set `sampling.run` to `false` to stop after the best fit.

### Sampling the posterior
The posterior of a cluster fit can have more than one mode (peak), e.g., a younger, more metal-rich solution and an older, more metal-poor one. Step 1 finds a local maximum of the posterior density from each of several starting points, and the sampler starts from every local maximum whose log posterior density is within 20 of the best one. A Markov chain Monte Carlo (MCMC) sampler started in several modes rarely moves between them, so the fraction of its samples in each mode does not measure the mode's probability. `sample_ssp` instead combines two methods:
 1. **Pilot runs.** From each local maximum, the affine-invariant ensemble sampler of [KissMCMC.jl](https://github.com/mauro3/KissMCMC.jl), which is similar to Python's [emcee](https://emcee.readthedocs.io/en/stable/), runs `sampling.pilot_walkers` walkers for `sampling.pilot_steps` steps. The second half of each run gives the location and width of its mode, including the correlations between parameters. Local maxima that lie inside each other's modes are merged into one mode.
 2. **Importance sampling.** Draws are taken from a mixture of broad multivariate Student's t distributions, one per mode, and each draw is weighted by the posterior density divided by the density of the mixture. Because the weights come from the posterior density itself, the total weight in each mode measures its probability correctly. The sampler takes `sampling.is_rounds` rounds (at least 3) of `sampling.is_samples` draws, refitting the location, width, and weight of each mode to the weighted draws after each round. The draws of all rounds are then combined and resampled into `sampling.is_rounds × sampling.is_samples` equally weighted posterior samples, in which a sample can appear more than once.

The posterior density of the draws is evaluated in parallel on threads, so start Julia with several threads, e.g., `julia --threads=auto`. The sampler also estimates the evidence (marginal likelihood) of the fit and of each mode, from which the probability of each mode follows.

The defaults (`pilot_steps: 300`, `is_rounds: 3`, `is_samples: 2000`) were enough in tests on simulated clusters, including one with two modes. After sampling, the weights tell whether the settings were adequate, and `fit_ssp` warns, and records the warning in `results.txt`, in two cases:
 - **Not converged**: the last round of importance sampling still moved the location or width of a mode substantially, or a few draws carry most of the weight (the Pareto k̂ diagnostic is above 0.7). The mixture does not yet describe the posterior well, so the samples can be biased. Increase `is_rounds` or `pilot_steps`.
 - **Imprecise**: the weighted draws are worth fewer than 1000 independent samples. The samples are unbiased but noisy. Increase `is_samples`.

Each round of importance sampling costs `is_samples` evaluations of the posterior density and each pilot run `pilot_steps × pilot_walkers`.

## Parameters and priors
Each entry in `parameters` is either a number, which fixes that parameter, or a prior distribution, which makes it a free parameter. Priors are written as in [Distributions.jl](https://juliastats.org/Distributions.jl/stable/univariate/), with numeric arguments only. These distributions are available:
 - `Uniform(a, b)`, `Normal(μ, σ)`, `LogNormal(μ, σ)`, `Beta(α, β)`, `Gamma(α, θ)`, and `Exponential(θ)`;
 - `truncated(d, lo, hi)`, which truncates the distribution `d` to the range from `lo` to `hi`; use `Inf` or `-Inf` for an open end, as in `truncated(Normal(0.1, 0.05), 0, Inf)`.

Age can be given either as `logAge` (log10 of the age in years) or as `age` in Gyr, and distance either as `distance_modulus` or as `distance` in parsecs, as in `simulate_catalog`. A prior on `age` or `distance` is a prior on that quantity, not on `logAge` or `distance_modulus`; e.g., `age: Uniform(1, 13.5)` gives every age between 1 and 13.5 Gyr the same prior probability, while `logAge: Uniform(9, 10.13)` favors younger ages. The binary fraction is required unless `binaries.model` is `NoBinaries`, and the background fraction defaults to `Uniform(0, 1)`. To fit without a background, set `data.background: none` (or fix `background_fraction: 0`), which is an error if you *also* set a prior or a nonzero value for `background_fraction`.

Priors on the age, metallicity, and extinction are truncated to the ranges that every listed stellar track library and bolometric correction grid can evaluate, with a message saying so; a fixed value outside those ranges is an error.

### [M/H] versus [Fe/H]
`MH` is the total metallicity [M/H] in the chemical composition of the stellar track library. For libraries with solar-scaled compositions, such as PARSEC and the default options of the other libraries, [M/H] equals [Fe/H]. For α-enhanced compositions (MISTv2, BaSTIv1, and BaSTIv2 with nonzero `alpha_fe`), [M/H] is larger than [Fe/H]. A prior from a spectroscopic [Fe/H] should be converted to [M/H] for these libraries, e.g., with the approximation of [Salaris, Chieffi, & Straniero (1993)](https://ui.adsabs.harvard.edu/abs/1993ApJ...414..580S), [M/H] ≈ [Fe/H] + log10(0.638 × 10<sup>[α/Fe]</sup> + 0.362), which gives [M/H] ≈ [Fe/H] + 0.29 for [α/Fe] = 0.4.

## Degeneracies
Some combinations of parameters change the model Hess diagram in similar ways, so the data constrain their combination better than each one alone. Their posterior samples are correlated, which shows as tilted ellipses in `results_corner.pdf`, and a fit can land away from the true values along these directions. The reported uncertainties include these correlations; narrowing them takes independent information, given as priors.

### Metallicity and extinction
Lowering the metallicity makes the stars of a population bluer, and more extinction makes them redder and fainter, so a lower [M/H] combined with a higher Av gives a similar color-magnitude diagram. The distance modulus and age are also involved, since extinction dims the stars as a larger distance would and age and metallicity both change the color and brightness of the main-sequence turnoff. In the fit of this example, the posterior correlation between [M/H] and Av is -0.82, and between Av and the distance modulus -0.46. The fitted values are offset from the true ones along this direction, by 0.5 standard deviations lower in [M/H] and 1.1 higher in Av, which is within what is expected for one random realization of the cluster. If you have an independent estimate, use it as a prior, e.g., Av from foreground dust maps or [M/H] from spectroscopy (see [M/H] versus [Fe/H] above for the conversion).

### Binary fraction and cluster mass
The fitted birth stellar mass is correlated with the binary fraction (with a correlation of 0.31 in this example, where the binary fraction is well constrained). Unresolved binaries appear as single stars, so each one hides the mass of its companion, and a higher binary fraction implies more mass for the same observed stars. The data constrain the binary fraction mainly through the sequence of unresolved binaries that lies brighter and redder than the main sequence. This constraint is weak when there are few stars on the main sequence, e.g., for old or sparse clusters, or when the photometric errors are comparable to the separation of the binary sequence from the main sequence.

When the data say little about the binary fraction, its posterior follows its prior, and the mass follows along. The reported uncertainties include this, but they assume the prior is trustworthy. With a prior of `Uniform(0, 1)` the posterior median of a poorly constrained binary fraction is near 0.5, which biases the mass high if the true fraction is lower. In tests on simulated clusters like the one in this example, with a true binary fraction of 0.3, the posterior median binary fraction ranged from 0.18 to 0.78 and the median mass ranged from 3% below to 10% above the true mass, closely following the binary fraction. If you have prior knowledge of the binary fraction, e.g., from other clusters of similar age and mass, use an informative prior such as `Beta(3, 7)` (mean 0.3, standard deviation 0.14) or fix it to a value, and report how the mass depends on that choice.

## Output
`simulate_catalog` writes the simulated catalog to `results/simulation/`. `fit_ssp` writes to `results/fit/`, where `<tracks>_<bcs>` names each combination of stellar track library and bolometric correction grid (e.g., `PARSEC_YBC`):
 - `results.txt`: one row per combination, with the best fit (`_best`, including the birth stellar mass at the best fit, `mass_best`) and the 16th, 50th, and 84th percentiles of the posterior (`_lower`, `_median`, `_upper`) of each free parameter, and percentiles of the age in Gyr, the birth stellar mass of the cluster, and the number of background stars. The best fit is the highest posterior density found by either the optimizer or the sampler, and `lp_best` is its log posterior density. `logZ` is the log evidence, `ess` the number of independent samples the posterior samples are worth, and `modes` the number of modes; warnings from the sampler are listed in the header.
 - `results_samples_<tracks>_<bcs>.txt`: the equally weighted posterior samples, with the free parameters, the birth stellar mass (`mass`), the expected number of background stars (`background_stars`), the log posterior density (`lp`), and the mode each sample was drawn from (`mode`).
 - `results_modes_<tracks>_<bcs>.txt`: one row per mode of the posterior, with its probability (`weight`), log evidence (`logZ`), effective number of samples (`ess`), and the highest local maximum found in it by the optimizer (`lp` and the free parameters).
 - `results_obshess.txt` and `results_modelhess_<tracks>_<bcs>.txt`: the observed Hess diagram and the best-fit model Hess diagram.
 - `input.yml`: a copy of the configuration, including the random seed used for the fit and the sampler.
 - `diagnostics.pdf`: the completeness used in the fit over the Hess diagram.

The script makes these figures in `results/fit/` for the first stellar track library and bolometric correction grid:
 - `results_cmd.pdf`: the color-magnitude diagram of the detected stars with the best-fit and true isochrones.
 - `results_hess.pdf`: the observed and best-fit model Hess diagrams with their residuals.
 - `results_corner.pdf`: the posterior distribution of each free parameter and of the birth stellar mass, and of each pair of them, with the true values (red), the best fit (blue), and the local maximum of each mode (orange circle). The best-fit mass is that of the best-fit parameters with the expected number of stars equal to the number observed.

The CMD and corner plots are made by `plot_ssp_cmd` and `plot_ssp_corner` from `SFHWorkflows.SSPFitting`, whose optional `truth` keyword marks the true values; omit it for real data. `plot_ssp_corner` can also be run later on a saved samples file, with the best fit read from the summary table by `read_best`, as shown in comments in `fit_ssp.jl`.

## Things to try
 - Change `properties.stellar_mass` in `simulate.yml` to see how the uncertainties shrink as the number of stars grows, or `stellartrack.logAge` and `stellartrack.MH` to simulate a different cluster.
 - Fix a parameter in `config.yml`, e.g., `distance_modulus: 24.95`, or give it a tighter prior, and see how the uncertainties of the others change; age, metallicity, distance, and extinction are correlated.
 - Change the prior on `binary_fraction`, e.g., to `Beta(3, 7)` or a fixed `0.3`, and compare the fitted mass in `results.txt` and the binary fraction and mass panels of `results_corner.pdf`.
 - Use a smaller distance modulus, so that the simulated CMD has a stronger main sequence, and observe the increase in precision in estimating the binary fraction and other parameters.
 - Add `track2` with a different stellar track library to `stellartracks` in `config.yml` to see the effect of model systematics on the fitted age and metallicity.
