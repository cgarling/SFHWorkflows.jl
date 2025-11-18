# This example runs the SFH fitting workflow `fit_sfh` and shows available plotting utilities.
using SFHWorkflows

config = "/home/cgarling/Development/julia/SFHWorkflows.jl/examples/config.yml"
# Perform the SFH fit
result, h = fit_sfh(config); # `result` contains results, h is a StatsBase.Histogram containing the observed Hess diagram

# Below is code to make common figures (Hess diagram / residual plots and cumulative SFH / AMR). If you'd rather make your own from the SFH result files, you can disregard the code below.
using SFHWorkflows.SFHFitting.Parsing: strip_whitespace
using CairoMakie
import YAML
# using LaTeXStrings: @L_str
galaxy_name = "Aquarius"

# Here we'll parse the x-axis color label and y-axis magnitude label from the provided config.
# You can also just hard-code it if you want.
dict = YAML.load_file(config) # Load config into dictionary so we can access filter names
xcolor = split(dict["data"]["binning"]["xcolor"], ",")
xcolor = join(strip_whitespace.(xcolor), " - ") # creates a string like "F475W - F814W" 
yfilter = dict["data"]["binning"]["yfilter"]    # creates a string like "F814W"
plot_path = joinpath("results", "results_hess.pdf") # file path to save figure to
# `result.results` contains the SFH fit for every combination of stellar track and bolometric correction listed in `config`.
# The results are stored in a 2-D Matrix indexed by stellar track first, bolometric correction grid second, so
# result.results[1,2] gives you the result for the first listed stellar track and the second listed bolometric correction grid.
# You can specify which result you want to plot via the idx keyword argument, which follows this convention;
# result.results[idx...] will be the result used to make the plot.
idx = [1,1] # Use the result for the first stellar track and first bolometric correction grid listed in `config`
fig, axs = SFHWorkflows.SFHFitting.plot_cmd_residuals(h, result, xcolor, yfilter, galaxy_name, plot_path; idx=idx, 
    c_clim=(-30,30), d_clim=(-12,12)); # Set limits on residual and significance panels
# This function saves the figure to the provided `plot_path`, but also 
# returns the figure `fig` and axes `axs`, which we could further modify if we wanted to.
# These are Makie.jl objects, see their documentation to further modify.
# We'll add a figure title labeling which stellar track and bolometric correction library we used
fig[0, :] = Label(fig, dict["stellartracks"]["track"*string(idx[1])]["name"] * " + " * dict["bolometriccorrections"]["bc"*string(idx[2])]["name"], fontsize=22, halign=:center)
# Save the updated figure
save(plot_path, fig)


# Now we'll make the cumulative SFH and AMR plot. 
plot_path2 = joinpath("results", "results_cumsfh.pdf")
fig, axs = SFHWorkflows.SFHFitting.plot_cumsfh_sys(result, plot_path2; idx=[1,1]);
# This function also returns the plot figure and axes so you can further modify them.
# We'll add a galaxy name label to the first axis.
tl_kws = Dict(:space => :relative, :text_align => (:left, :top), :background_color => (:white, 1.0), :text_color => :black, :strokewidth => 0)
textlabel!(axs[1], 0.1, 0.95; text = "$galaxy_name", tl_kws...)
# Save the updated figure
save(plot_path2, fig)
