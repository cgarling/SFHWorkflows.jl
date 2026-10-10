# This example simulates a catalog that combines two single stellar populations of different ages and metallicities
# (configured by simulate_young.yml and simulate_old.yml), fits it with a single stellar population with `fit_ssp`
# (configured by config.yml), and plots both modes of the posterior that result.
# The plotting code also requires CairoMakie and YAML in the active environment: `]add CairoMakie YAML`
using SFHWorkflows
using Random: Xoshiro, randperm

# Simulate both populations and keep the stars detected in both filters
function detected_stars(file)
    catalog, _, _ = simulate_catalog(file)
    detected = isfinite.(catalog.F475W_obs) .& isfinite.(catalog.F814W_obs)
    return catalog.F475W_obs[detected], catalog.F814W_obs[detected]
end
young, old = detected_stars("simulate_young.yml"), detected_stars("simulate_old.yml")
# Keep a random subset of the old stars, sized so that the two modes have similar probabilities.
keep = randperm(Xoshiro(53), length(old[1]))[1:1764]
F475W, F814W = vcat(young[1], old[1][keep]), vcat(young[2], old[2][keep])
open(joinpath("results", "simulation", "phot.dat"), "w") do io
    foreach((b, i) -> println(io, b, " ", i), F475W, F814W)
end
@info "Combined $(length(young[1])) young and $(length(keep)) old stars"

# Fit the combined catalog. `results[1, 1].samples.modes` holds the probability and evidence of each mode.
results, h = fit_ssp("config.yml")

# Below is code to make figures of the fit. If you'd rather make your own from the output files, you can disregard the
# code below.
###########################################
using SFHWorkflows.SSPFitting: plot_ssp_cmd, plot_ssp_corner, SFH
import YAML
r = results[1, 1]
# True properties of the two populations; without binaries, the binary fraction is not a parameter
truths = map(("simulate_young.yml", "simulate_old.yml")) do file
    sim = YAML.load_file(file)
    (logAge = only(eval(Meta.parse(sim["stellartrack"]["logAge"]))), MH = only(eval(Meta.parse(sim["stellartrack"]["MH"]))),
     dmod = sim["properties"]["distance_modulus"], Av = sim["properties"]["Av"])
end

# Observed CMD of the combined catalog with the isochrone of the highest local maximum in each mode, and the true
# isochrone of the population with the closest age
for (k, m) in enumerate(r.samples.modes)
    p = SFH.ssp_params(r.model, m.θ)
    plot_ssp_cmd(r.model, p, F475W .- F814W, F814W, joinpath("results", "fit", "results_cmd_mode$k.pdf"); truth = argmin(t -> abs(t.logAge - p.logAge), truths), xlabel = "F475W - F814W",
                 ylabel = "F814W", title = "Mode $k, probability $(round(m.weight; digits=2))")
end

# Corner plot of the posterior samples, with the best fit in blue and the modes as orange circles
plot_ssp_corner(r, joinpath("results", "fit", "results_corner.pdf"))
