# This example simulates a star cluster with `simulate_catalog` (configured by simulate.yml), fits it with `fit_ssp`
# (configured by config.yml), and plots the fit compared to the true properties of the simulated cluster.
# The plotting code also requires CairoMakie and YAML in the active environment: `]add CairoMakie YAML`
using SFHWorkflows
import YAML

# Simulate the cluster. `catalog` contains every simulated star; the mock-observed magnitudes of stars that were not
# detected are NaN.
catalog, _, _ = simulate_catalog("simulate.yml")
# Write the detected stars to the photometry file read by config.yml, with one column per filter
sim = YAML.load_file("simulate.yml")
detected = isfinite.(catalog.F475W_obs) .& isfinite.(catalog.F814W_obs)
open(joinpath(sim["output"]["path"], "phot.dat"), "w") do io
    foreach((b, i) -> println(io, b, " ", i), catalog.F475W_obs[detected], catalog.F814W_obs[detected])
end
@info "Detected $(count(detected)) of $(length(detected)) simulated stars"

# Fit the cluster. `results` is a matrix indexed by stellar track library first and bolometric correction grid second,
# like `result.results` from fit_sfh; each entry has fields `model`, `fit` (the best fit), and `chain` (the posterior
# samples). `h` is the observed Hess diagram.
results, h = fit_ssp("config.yml")

# Below is code to make figures comparing the fit to the simulated cluster. If you'd rather make your own from the
# output files, you can disregard the code below.
###########################################
using SFHWorkflows.SFHFitting.Parsing: strip_whitespace, parse_gates
using SFHWorkflows.SSPFitting: best_fit, plot_ssp_cmd, plot_ssp_corner, SFH
using CairoMakie
cluster_name = "Simulated cluster"

dict = YAML.load_file("config.yml")
output_path = dict["output"]["path"]
xcolor = split(strip_whitespace(dict["data"]["binning"]["xcolor"]), ",")
yfilter = dict["data"]["binning"]["yfilter"]
gates = parse_gates(dict["data"]["binning"])
# True properties of the simulated cluster; it has no background stars. With real data, omit the `truth` keyword below.
truth = (logAge = only(eval(Meta.parse(sim["stellartrack"]["logAge"]))), MH = only(eval(Meta.parse(sim["stellartrack"]["MH"]))),
         dmod = sim["properties"]["distance_modulus"], Av = sim["properties"]["Av"], binary_fraction = sim["binaries"]["binary_fraction"],
         background_fraction = 0.0, mass = sim["properties"]["stellar_mass"])

# Plot the result for the first stellar track library and first bolometric correction grid in config.yml
idx = [1, 1]
r = results[idx...]
label = dict["stellartracks"]["track"*string(idx[1])]["name"] * " + " * dict["bolometriccorrections"]["bc"*string(idx[2])]["name"]

# Observed CMD of the detected stars with the best-fit and true isochrones
plot_ssp_cmd(r, catalog.F475W_obs[detected] .- catalog.F814W_obs[detected], catalog.F814W_obs[detected], joinpath(output_path, "results_cmd.pdf");
             truth, xlabel = join(xcolor, " - "), ylabel = yfilter, title = "$cluster_name, $label", gates)

# Observed and best-fit model Hess diagrams with their residuals. The best fit is the highest posterior density found by
# the optimizer or the sampler, as in the summary table.
best, _ = best_fit(r.fit, r.chain)
plot_path = joinpath(output_path, "results_hess.pdf")
fig, axs = SFHWorkflows.SFHFitting.plot_cmd_residuals(h, SFH.ssp_hess(r.model, best), join(xcolor, " - "), yfilter, cluster_name, plot_path; gates)
fig[0, :] = Label(fig, label, fontsize = 22, halign = :center)
save(plot_path, fig)

# Corner plot of the posterior samples, with the true values in red and the best fit in blue. It can also be made later
# from the output files, e.g., for the chain file of PARSEC + YBC,
# plot_ssp_corner(joinpath(output_path, "results_chain_PARSEC_YBC.txt"), "corner.pdf";
#                 best = SFHWorkflows.SSPFitting.read_best(joinpath(output_path, "results.txt"), "PARSEC_YBC"))
isnothing(r.chain) || plot_ssp_corner(r, joinpath(output_path, "results_corner.pdf"); truth)
