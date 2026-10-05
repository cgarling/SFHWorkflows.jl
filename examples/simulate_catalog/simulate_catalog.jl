using SFHWorkflows

config = "./config.yml"
# Simulate the catalog. `catalog` and `truth` are tables that are also written to the output path; `fit` is `nothing`
# unless `fit.run` is true in `config`, in which case it is `(result, h)` as returned by `fit_sfh`
catalog, truth, fit = simulate_catalog(config)

# Below is code to make common figures (CMDs of the catalog and the input SFH, compared to the fit if one was run).
# If you'd rather make your own from the output files, you can disregard the code below.
###########################################
using SFHWorkflows.SFHFitting.Parsing: strip_whitespace
using CairoMakie
import YAML
galaxy_name = "Simulated"

# Return a boolean mask indicating which points are in sparsely populated bins.
function sparse_mask(x, y; bins = 80, threshold = 25)
    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)

    nx = ny = bins

    # Map each point to a rectangular density bin.
    # Rectangular bins are fine here; they're just being used to decide
    # which points are worth drawing individually.
    ix = clamp.(floor.(Int, (x .- xmin) ./ (xmax - xmin) .* nx) .+ 1, 1, nx)
    iy = clamp.(floor.(Int, (y .- ymin) ./ (ymax - ymin) .* ny) .+ 1, 1, ny)

    counts = zeros(Int, nx, ny)

    @inbounds for i in eachindex(ix)
        counts[ix[i], iy[i]] += 1
    end

    keep = BitVector(undef, length(x))

    @inbounds for i in eachindex(keep)
        keep[i] = counts[ix[i], iy[i]] < threshold
    end

    return keep
end

dict = YAML.load_file(config)
output_path = dict["output"]["path"]
# CMD filters and limits come from fit.binning when present; otherwise set them here
binning = get(get(dict, "fit", Dict()), "binning", nothing)
if isnothing(binning)
    xcolor, yfilter = ["F475W", "F814W"], "F814W"
    limits = (nothing, nothing) # Full range of the data
else
    xcolor, yfilter = split(strip_whitespace(binning["xcolor"]), ","), binning["yfilter"]
    limits = (extrema(eval(Meta.parse(binning["xbins"]))), extrema(eval(Meta.parse(binning["ybins"]))))
end
xlabel = join(xcolor, " - ")

# CMD panel; sparse regions are drawn as individual points and dense regions as a hexbin density
function cmd_panel!(fig, col, suffix, title)
    x = getproperty(catalog, Symbol(xcolor[1], suffix)) .- getproperty(catalog, Symbol(xcolor[2], suffix))
    y = getproperty(catalog, Symbol(yfilter, suffix))
    sel = isfinite.(x) .& isfinite.(y) # Undetected stars have NaN observed magnitudes
    isnothing(limits[1]) || (sel .&= (limits[1][1] .<= x .<= limits[1][2]) .& (limits[2][1] .<= y .<= limits[2][2]))
    xs, ys = x[sel], y[sel]
    ax = Axis(fig[1, 2col-1]; xlabel, ylabel = yfilter, title, limits, yreversed = true)
    smask = sparse_mask(xs, ys; bins = 80)
    scatter!(ax, xs[smask], ys[smask]; color = :black, markersize = 2, rasterize = true)
    # The hexbin counts only the dense points (weight 1) so each star is drawn once, while binning all points keeps the
    # hexagon grid spanning the full data range; sparse CMDs have no dense points and are drawn only as scatter
    if !all(smask)
        h = hexbin!(ax, xs, ys; weights = .!smask, threshold = 1, bins = 80, colorscale = log10)
        Colorbar(fig[1, 2col], h; label = "Stars")
    end
    return ax
end

# CMDs of the intrinsic magnitudes and, if the CMD filters were mock-observed, the observed magnitudes
observed = all(f -> Symbol(f, "_obs") in propertynames(catalog), [xcolor; yfilter])
fig = Figure(size = (observed ? 1200 : 600, 600))
cmd_panel!(fig, 1, "", "Intrinsic")
observed && cmd_panel!(fig, 2, "_obs", "Mock-observed")
save(joinpath(output_path, "catalog_cmd.pdf"), fig)

# Input SFH (truth) as cumulative SFH and mean metallicity vs. lookback time, with the same conventions as plot_cumsfh_sys
function plot_truth!(ax1, ax2; kws...)
    x = exp10.(vcat(truth.logAge_lower[1], truth.logAge_upper) .- 9) # Lookback time in Gyr
    lines!(ax1, x, vcat(truth.cum_sfh, 0.0); kws...)
    # Mean [M/H] is only defined where stars formed; NaN breaks the line in bins without star formation
    lines!(ax2, x[begin:end-1], ifelse.(truth.sfr .> 0, truth.MH, NaN); kws...)
end

if isnothing(fit)
    fig = Figure(size = (750, 450))
    axs = (Axis(fig[1, 1], xlabel = "Lookback Time [Gyr]", ylabel = "Cumulative SF", xreversed = true),
           Axis(fig[1, 2], xlabel = "Lookback Time [Gyr]", ylabel = "⟨[M/H]⟩", xreversed = true))
    plot_truth!(axs...; color = :black)
    save(joinpath(output_path, "truth_cumsfh.pdf"), fig)
else
    # Fit figures as in examples/fit_sfh/fit_sfh.jl, written to the fit's output path with its configuration
    result, h = fit
    # simulate_catalog writes the fit to <output.path>/fit unless the config gives fit.output
    fit_path = haskey(dict["fit"], "output") ? dict["fit"]["output"]["path"] : joinpath(output_path, "fit")
    fitdict = YAML.load_file(joinpath(fit_path, "input.yml"))
    # result.results is indexed by stellar track first, bolometric correction grid second
    idx = [1, 1]
    plot_path = joinpath(fit_path, "results_hess.pdf")
    fig, axs = SFHWorkflows.SFHFitting.plot_cmd_residuals(h, result, xlabel, yfilter, galaxy_name, plot_path; idx = idx)
    fig[0, :] = Label(fig, fitdict["stellartracks"]["track"*string(idx[1])]["name"] * " + " * fitdict["bolometriccorrections"]["bc"*string(idx[2])]["name"], fontsize = 22, halign = :center)
    save(plot_path, fig)

    # Fitted cumulative SFH and AMR with the input SFH overlaid
    plot_path2 = joinpath(fit_path, "results_cumsfh.pdf")
    fig, axs = SFHWorkflows.SFHFitting.plot_cumsfh_sys(result, plot_path2; idx = idx)
    plot_truth!(axs...; color = :red, linestyle = :dash, label = "Truth")
    axislegend(axs[1]; position = :rb)
    save(plot_path2, fig)
end
