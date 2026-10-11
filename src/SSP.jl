"""Module to fit the photometry of a single stellar population of one age and metallicity, such as a star cluster. Main function is `fit_ssp`, which takes in a properly formatted YAML configuration file."""
module SSPFitting

export fit_ssp, plot_ssp_corner, plot_ssp_cmd, read_best

import StarFormationHistories as SFH
import CSV
import StellarTracks
import YAML
using BolometricCorrections: gridname, select_filters
using OrderedCollections: OrderedDict
using Printf: Format, format
using StableRNGs: StableRNG
using StatsBase: Histogram, quantile
using TypedTables: Table
using CairoMakie: Figure, Axis, Legend, LineElement, MarkerElement, Point2f, WilkinsonTicks, cgrad, save, axislegend, @L_str,
    hist!, scatter!, hexbin!, lines!, vlines!, hlines!, ylims!, hidedecorations!, hidexdecorations!, hideydecorations!,
    linkxaxes!, rowgap!
using ..SFHFitting: observational_models, write_diagnostics, gate_mask, obs_hess, write_histogram, background_hess
using ..SFHFitting.Parsing: parse_ssp_config, restrict
using ..SFHFitting.Systematics: iso_filters
using ..SFHFitting.Plotting: sparse_mask

# Whether `isochrone` can be evaluated; parameters outside a library's range throw a DomainError or ArgumentError
function evaluable(tl, bc, logAge, mh, Av)
    try
        StellarTracks.isochrone(tl, bc, logAge, mh, Av)
        return true
    catch e
        e isa Union{DomainError, ArgumentError} || rethrow()
        return false
    end
end

# The logAge closest to `bound`, between the valid `la0` and `bound`, at which `isochrone` can be evaluated, by bisection;
# track libraries do not report their age range directly
function age_limit(tl, bc, la0, bound, mh, Av)
    evaluable(tl, bc, bound, mh, Av) && return bound
    good, bad = la0, bound
    while abs(bad - good) > 1e-3
        mid = (good + bad) / 2
        evaluable(tl, bc, mid, mh, Av) ? (good = mid) : (bad = mid)
    end
    return good
end

# Center and bounds of a fixed parameter or prior
center(p::Real) = p
center(p) = quantile(p, 0.5)
bounds(p::Real) = (p, p)
bounds(p) = (minimum(p), maximum(p))

"""
    restrict_parameters(parameters, stellar_tracks, bcs)
Truncates the priors on `logAge`, `MH`, and `Av` (or checks their fixed values) to the ranges that every stellar track
library and bolometric correction grid can evaluate: [M/H] to the range of the track libraries, Av to the range of the
BC grids, and logAge to the ages at which every track library returns an isochrone.
"""
function restrict_parameters(parameters, stellar_tracks, bcs)
    mh_lo, mh_hi = maximum(first(extrema(StellarTracks.MH(tl))) for tl in stellar_tracks), minimum(last(extrema(StellarTracks.MH(tl))) for tl in stellar_tracks)
    av_lo, av_hi = maximum(Float64(first(extrema(bc).Av)) for bc in bcs), minimum(Float64(last(extrema(bc).Av)) for bc in bcs)
    MH = restrict(parameters.MH, mh_lo, mh_hi, "parameters.MH")
    Av = restrict(parameters.Av, av_lo, av_hi, "parameters.Av")
    mh, av, la0 = center(MH), center(Av), center(parameters.logAge)
    # Search only as far as ages that are plausible for stellar populations
    lo_bound, hi_bound = max(bounds(parameters.logAge)[1], 5.0), min(bounds(parameters.logAge)[2], 10.5)
    la_lo = maximum(age_limit(tl, first(bcs), la0, lo_bound, mh, av) for tl in stellar_tracks)
    la_hi = minimum(age_limit(tl, first(bcs), la0, hi_bound, mh, av) for tl in stellar_tracks)
    return merge(parameters, (logAge=restrict(parameters.logAge, la_lo, la_hi, "parameters.logAge"), MH, Av))
end

# Fits one combination of stellar track library and BC grid, and samples the posterior if requested, with random numbers
# from `seed`
function fit_combination(tl, bc, h, config, parameters, err, completeness, bias, background, mask, seed)
    iso_symb, yidx, xidxs = iso_filters(bc, config.ystring, config.xstrings)
    grid = select_filters(bc, iso_symb) # BC tables interpolate only the filters that are used
    isofunc(logAge, mh, Av) = (iso = StellarTracks.isochrone(tl, grid, logAge, mh, Av); (iso.m_ini, [getproperty(iso, s) for s in iso_symb]))
    model = SFH.SSPModel(isofunc, h.weights, h.edges, err, yidx, xidxs, config.imf, completeness, bias;
                         parameters..., binary_model=config.binary_model, background, background_floor=config.background.floor, mask)
    rng = StableRNG(seed)
    fit = SFH.fit_ssp(model; config.ncandidates, config.nstarts, config.restarts, rng)
    config.run_sampling || return (; model, fit, samples=nothing)
    samples = SFH.sample_ssp(fit; config.pilot_steps, config.pilot_walkers, config.is_rounds, config.is_samples, rng)
    return (; model, fit, samples)
end

# Number formats for the output tables; fractions, masses, and star counts can span orders of magnitude (e.g., a background
# fraction near 0), so they are written in exponential format
column_format(name) = name in (:mode, :modes) ? "%d" : name == :ess ? "%.0f" : name in (:lp, :lp_best, :logZ) ? "%.3f" :
    any(s -> occursin(s, string(name)), ("fraction", "mass", "background_stars", "weight")) ? "%.4e" : "%.5f"
function write_table(fname, table, comments)
    names = collect(propertynames(table))
    open(fname, "w") do io
        foreach(c -> println(io, "# ", c), comments)
        # append keeps CSV.write from seeking back to the start of `io` and overwriting the comments
        CSV.write(io, table; append=true, writeheader=true, delim=' ', transform=(col, val) -> val isa Number ? format(Format(column_format(names[col])), val) : val)
    end
end

# Table of the modes of the posterior samples `s` (StarFormationHistories.SSPSamples): the posterior probability, log
# evidence, and effective sample size of each mode, and the local solution of fit_ssp it was found from
modes_table(s) = Table([merge((mode=k, weight=m.weight, logZ=m.logZ, ess=m.ess, lp=m.logdensity), NamedTuple{Tuple(s.fit.model.free)}(m.θ))
                        for (k, m) in enumerate(s.modes)])

# Highest-density parameters found and their log density: the MAP from the optimizer, or the posterior sample with the
# highest log density if it is higher
function best_fit(fit, samples)
    isnothing(samples) && return fit.map, fit.logdensity
    lp = samples.samples.lp
    i = argmax(lp)
    lp[i] <= fit.logdensity && return fit.map, fit.logdensity
    return SFH.ssp_params(fit.model, [samples.samples[k][i] for k in fit.model.free]), lp[i]
end

# Summary row: best fit and 16th, 50th, and 84th percentiles of each free parameter and derived quantity, the stellar mass
# formed at the best fit, and the log evidence, effective sample size, and number of modes of the posterior samples
function summary_row(name, fit, s)
    best, lp = best_fit(fit, s)
    row = OrderedDict{Symbol, Any}(:name => name, :lp_best => lp)
    for k in fit.model.free
        row[Symbol(k, :_best)] = best[k]
        samples = isnothing(s) ? Float64[] : s.samples[k]
        for (suffix, q) in zip((:_lower, :_median, :_upper), (0.16, 0.5, 0.84))
            row[Symbol(k, suffix)] = isempty(samples) ? NaN : quantile(samples, q)
        end
    end
    row[:mass_best] = SFH.ssp_mass(fit.model, best)
    if !isnothing(s)
        derived = Dict(:mass => s.samples.mass, :background_stars => s.samples.background_stars)
        :logAge in fit.model.free && (derived[:age_Gyr] = exp10.(s.samples.logAge .- 9))
        for k in (:age_Gyr, :mass, :background_stars), (suffix, q) in zip((:_lower, :_median, :_upper), (0.16, 0.5, 0.84))
            haskey(derived, k) && (row[Symbol(k, suffix)] = quantile(derived[k], q))
        end
        row[:logZ], row[:ess], row[:modes] = s.logZ, s.ess, length(s.modes)
    end
    return NamedTuple(row)
end

"""
    (results, h) = fit_ssp(config_file::AbstractString)
    (results, h) = fit_ssp(config::NamedTuple)

Fits the Hess diagram of the photometry described by the configuration with a single stellar population of one age and
metallicity plus a background, for every combination of the stellar track libraries and bolometric correction grids in
the configuration (see `examples/fit_ssp/config.yml`). The combinations are fit in parallel over threads. For each, the
maximum a posteriori (MAP) parameters are found with `StarFormationHistories.fit_ssp` and the posterior is sampled with
`StarFormationHistories.sample_ssp` unless `sampling.run` is `false`.

Writes to `output.path` a summary table (`output.filename`) with one row per combination, the posterior samples and the
modes of the posterior of each combination, the observed Hess diagram and the best-fit model Hess diagram of each
combination, and a copy of the configuration including the random seed; warnings from the sampler are also recorded in
the summary table. Returns the matrix `results` indexed by stellar track library and BC grid, whose entries are
`NamedTuple`s with fields `model` (`StarFormationHistories.SSPModel`), `fit` (`StarFormationHistories.SSPFit`), and
`samples` (`StarFormationHistories.SSPSamples`, or `nothing`), and the observed Hess diagram `h`.
"""
function fit_ssp(config::NamedTuple)
    output_path = config.output_path
    mkpath(output_path)
    # The seed is always recorded in the saved input.yml so the fit and samples can be reproduced
    seed = isnothing(config.seed) ? rand(0:typemax(Int)) : Int(config.seed)
    cfg = deepcopy(config.config)
    cfg["seed"] = seed
    YAML.write_file(joinpath(output_path, "input.yml"), cfg)

    edges = (config.xbins, config.ybins)
    mask = gate_mask(edges, config.gates)
    completeness, bias, err = observational_models(config.xstrings, config.ystring, config.ast_file, config.ast_filters, config.filter_models;
                                                   config.badval, config.minerr, config.maxerr, config.plot_diagnostics, output_path)
    config.plot_diagnostics && write_diagnostics(completeness, config.ast_file, edges, config.xstrings, config.ystring, output_path)
    filters = string.(config.filters)
    h = obs_hess(config.phot_file, filters, config.xstrings, config.ystring, edges)
    background = background_hess(config.background, filters, config.xstrings, config.ystring, edges)
    parameters = restrict_parameters(config.parameters, config.stellar_tracks, config.bcs)

    # One task per combination; each fit's starting grid and sampling are also threaded, and Julia schedules the nested
    # threaded loops across all threads. Seeds are drawn serially so the results do not depend on scheduling.
    combos = collect(Iterators.product(eachindex(config.stellar_tracks), eachindex(config.bcs)))
    seeds = rand(StableRNG(seed), UInt64, length(combos))
    @info "Fitting $(length(combos)) combination(s) of stellar tracks and bolometric corrections"
    tasks = map(zip(combos, seeds)) do ((i, j), s)
        Threads.@spawn fit_combination(config.stellar_tracks[i], config.bcs[j], h, config, parameters, err, completeness, bias, background, mask, s)
    end
    results = reshape(fetch.(tasks), size(combos))

    base, ext = splitext(config.output_filename)
    write_histogram(h, joinpath(output_path, base * "_obshess" * ext))
    fixed = join(("$k = $v" for (k, v) in pairs(parameters) if v isa Real), ", ")
    warnings = String[]
    rows = map(CartesianIndices(results)) do I
        r = results[I]
        name = gridname(config.stellar_tracks[I[1]]) * "_" * gridname(config.bcs[I[2]])
        write_histogram(Histogram(edges, SFH.ssp_hess(r.model, first(best_fit(r.fit, r.samples)))), joinpath(output_path, base * "_modelhess_" * name * ext))
        if !isnothing(r.samples)
            write_table(joinpath(output_path, base * "_samples_" * name * ext), Table(r.samples.samples),
                        ["$name equal-weight posterior samples, worth $(round(Int, r.samples.ess)) independent samples", "Fixed parameters: $fixed", "Seed: $seed"])
            write_table(joinpath(output_path, base * "_modes_" * name * ext), modes_table(r.samples),
                        ["$name modes of the posterior: probability (weight), log evidence (logZ), effective sample size (ess), and the best local solution of the optimizer in each (lp and the free parameters)", "Seed: $seed"])
            append!(warnings, ("Warning for $name: $w" for w in r.samples.warnings))
        end
        summary_row(name, r.fit, r.samples)
    end
    write_table(joinpath(output_path, config.output_filename), Table(vec(rows)),
                ["Best fit (highest log posterior density found by the optimizer or the sampler) and 16th, 50th, and 84th percentiles",
                 "Fixed parameters: $fixed", "Seed: $seed", warnings...])
    return results, h
end
fit_ssp(config_file::AbstractString) = fit_ssp(parse_ssp_config(config_file))

"""
    read_best(summary_file::AbstractString, name::AbstractString)

Returns the best fit of the combination of stellar tracks and bolometric corrections `name` (e.g., `"PARSEC_YBC"`) from the
summary table written by `fit_ssp`, as a `NamedTuple` of the free parameters and `mass`, for use with
[`plot_ssp_corner`](@ref).
"""
function read_best(summary_file::AbstractString, name::AbstractString)
    row = only(r for r in CSV.File(summary_file; delim=' ', comment="#") if r.name == name)
    return (; (Symbol(chopsuffix(string(k), "_best")) => getproperty(row, k) for k in propertynames(row) if endswith(string(k), "_best") && k != :lp_best)...)
end

# Best fit of an entry of the `results` returned by fit_ssp, with the stellar mass formed at the best fit
function best_with_mass(result)
    best, _ = best_fit(result.fit, result.samples)
    return merge(best, (mass = SFH.ssp_mass(result.model, best),))
end

# Free parameters and stellar mass formed at the local solution of fit_ssp from which each mode of `s` was found
mode_points(s) = [merge(NamedTuple{Tuple(s.fit.model.free)}(m.θ), (mass = SFH.ssp_mass(s.fit.model, SFH.ssp_params(s.fit.model, m.θ)),))
                  for m in s.modes]

# Columns of posterior samples by name, from the StarFormationHistories.SSPSamples returned by sample_ssp or a samples file
# written by fit_ssp
sample_columns(s::SFH.SSPSamples) = Dict(pairs(s.samples))
function sample_columns(file::AbstractString)
    f = CSV.File(file; delim=' ', comment="#")
    return Dict(k => collect(Float64, getproperty(f, k)) for k in propertynames(f))
end

const CORNER_LABELS = Dict(:logAge => L"\log_{10}(\mathrm{Age} / \mathrm{yr})", :MH => L"[\mathrm{M/H}]", :dmod => L"\mu\ [\mathrm{mag}]",
                           :Av => L"A_V\ [\mathrm{mag}]", :binary_fraction => L"f_\mathrm{binary}",
                           :background_fraction => L"f_\mathrm{background}", :mass => L"\log_{10}(M_\star / M_\odot)")

"""
    (fig, axs) = plot_ssp_corner(result::NamedTuple, output_file::AbstractString; truth=nothing)
    (fig, axs) = plot_ssp_corner(samples, output_file::AbstractString; best=nothing, truth=nothing, modes)

Makes a corner plot of the posterior samples of an SSP fit and saves it to `output_file`: the distribution of each free
parameter and of the stellar mass formed (as log10) on the diagonal, and of each pair of them below it. `result` is an
entry of the `results` returned by `fit_ssp`, whose best fit (including its mass) is marked in blue and the modes of its
posterior with orange circles. Alternatively, `samples` is the `StarFormationHistories.SSPSamples` returned by
`StarFormationHistories.sample_ssp` or the path to a samples file written by `fit_ssp`, with an optional best fit `best`
given as a `NamedTuple` (e.g., from [`read_best`](@ref)) and the modes `modes` as a vector of `NamedTuple`s (by default
those of an `SSPSamples`, and none for a file). If given, the values in the `NamedTuple` `truth` (e.g., the true
parameters of a simulated cluster, including `mass`) are marked in red; parameters missing from `truth` are not marked.
Returns the figure and the matrix of axes, indexed by row and column.
"""
function plot_ssp_corner(result::NamedTuple, output_file::AbstractString; truth=nothing)
    isnothing(result.samples) && throw(ArgumentError("`result` has no posterior samples; run fit_ssp with sampling.run set to true."))
    return plot_ssp_corner(result.samples, output_file; best=best_with_mass(result), truth)
end
function plot_ssp_corner(samples, output_file::AbstractString; best=nothing, truth=nothing,
                         modes=samples isa SFH.SSPSamples ? mode_points(samples) : nothing)
    cols = sample_columns(samples)
    pars = [[k for k in SFH.SSP_PARAMETERS if haskey(cols, k)]; :mass]
    # The mass is plotted as log10(mass) so its tick labels stay short
    tf(k) = k == :mass ? log10 : identity
    samples = [tf(k).(cols[k]) for k in pars]
    value(p, k) = isnothing(p) || !haskey(p, k) ? nothing : tf(k)(p[k])
    n = length(pars)
    fig = Figure(size = (170n, 170n))
    axs = Matrix{Any}(nothing, n, n)
    for i in 1:n, j in 1:i
        ax = axs[i, j] = Axis(fig[i, j]; xlabel = CORNER_LABELS[pars[j]], ylabel = CORNER_LABELS[pars[i]], xticks = WilkinsonTicks(3),
                              yticks = WilkinsonTicks(3), xticklabelrotation = π / 4)
        tj, ti, bj, bi = value(truth, pars[j]), value(truth, pars[i]), value(best, pars[j]), value(best, pars[i])
        if i == j
            hist!(ax, samples[i]; bins = 30, color = :gray)
            ylims!(ax, 0, nothing)
            isnothing(ti) || vlines!(ax, ti; color = :red)
            isnothing(bi) || vlines!(ax, bi; color = :dodgerblue, linestyle = :dash)
            # The bottom histogram keeps its x-axis, which labels the last column
            i == n ? hideydecorations!(ax) : hidedecorations!(ax)
        else
            # Sparse regions are drawn as individual samples and dense regions as a hexbin density
            smask = sparse_mask(samples[j], samples[i]; bins = 40, threshold = 10)
            scatter!(ax, samples[j][smask], samples[i][smask]; color = :gray50, markersize = 3, strokewidth = 0)
            all(smask) || hexbin!(ax, samples[j], samples[i]; weights = .!smask, threshold = 1, bins = 30, colormap = cgrad([:gray70, :black]), colorscale = log10)
            isnothing(tj) || vlines!(ax, tj; color = :red)
            isnothing(ti) || hlines!(ax, ti; color = :red)
            # Modes are open circles, so both show when the best fit is also the location of a mode
            if !isnothing(modes)
                mj, mi = [value(m, pars[j]) for m in modes], [value(m, pars[i]) for m in modes]
                ok = .!isnothing.(mj) .& .!isnothing.(mi)
                any(ok) && scatter!(ax, Float64.(mj[ok]), Float64.(mi[ok]); color = :transparent, strokecolor = :darkorange, strokewidth = 2, markersize = 14)
            end
            isnothing(bi) || isnothing(bj) || scatter!(ax, [bj], [bi]; color = :dodgerblue, marker = :xcross, markersize = 12)
            j > 1 && hideydecorations!(ax; grid = false)
        end
        i < n && hidexdecorations!(ax; grid = false)
    end
    foreach(j -> linkxaxes!(filter(!isnothing, axs[:, j])...), 1:n)
    entries = Tuple{Any, String}[]
    isnothing(truth) || push!(entries, (LineElement(color = :red), "Truth"))
    isnothing(best) || push!(entries, (MarkerElement(color = :dodgerblue, marker = :xcross, markersize = 12), "Best fit"))
    isnothing(modes) || isempty(modes) || push!(entries, (MarkerElement(marker = :circle, color = :transparent, strokecolor = :darkorange, strokewidth = 2, markersize = 14), length(modes) == 1 ? "Mode" : "Modes"))
    isempty(entries) || Legend(fig[1, max(2, n - 1):n], first.(entries), last.(entries); framevisible = false)
    rowgap!(fig.layout, 0)
    save(output_file, fig)
    return fig, axs
end

# Color and magnitude of the isochrone of `model` for the parameters `p`; `isofunc` returns the initial masses and the
# absolute magnitudes in the filters of the model, which are indexed by `y_index` and `color_indices`
function isochrone_cmd(model, p)
    _, mags = model.isofunc(p.logAge, p.MH, p.Av)
    return mags[model.color_indices[1]] .- mags[model.color_indices[2]], mags[model.y_index] .+ p.dmod
end

"""
    (fig, ax) = plot_ssp_cmd(result::NamedTuple, color, mag, output_file::AbstractString; truth=nothing, xlabel="",
                             ylabel="", title="", gates=[])
    (fig, ax) = plot_ssp_cmd(model::StarFormationHistories.SSPModel, best::NamedTuple, color, mag,
                             output_file::AbstractString; kws...)

Plots the color-magnitude diagram of the observed stars, with colors `color` and magnitudes `mag`, and the isochrone of
the best fit of `result`, an entry of the `results` returned by `fit_ssp`, over the range of its Hess diagram, and saves
it to `output_file`. Alternatively, give the `model` and the best-fit parameters `best` (e.g., from [`read_best`](@ref)). If given, the isochrone for the parameters in the `NamedTuple` `truth` (e.g., the true parameters of
a simulated cluster) is also plotted; parameters missing from `truth` are taken from the best fit. The polygons `gates`
(e.g., from `data.binning.gates`) are outlined. Returns the figure and axis.
"""
plot_ssp_cmd(result::NamedTuple, color, mag, output_file::AbstractString; kws...) =
    plot_ssp_cmd(result.model, first(best_fit(result.fit, result.samples)), color, mag, output_file; kws...)
function plot_ssp_cmd(model::SFH.SSPModel, best::NamedTuple, color, mag, output_file::AbstractString; truth=nothing, xlabel="", ylabel="", title="", gates=[])
    # Fixed parameters are not in `best` when it is read from the summary table
    best = merge(model.params, best)
    fig = Figure(size = (550, 650))
    ax = Axis(fig[1, 1]; xlabel, ylabel, title, limits = (extrema(model.edges[1]), extrema(model.edges[2])), yreversed = true)
    scatter!(ax, color, mag; color = :black, markersize = 2, rasterize = true)
    lines!(ax, isochrone_cmd(model, best)...; color = :dodgerblue, linewidth = 3, label = "Best fit")
    # Drawn after the best fit and dashed so it shows where the two overlap
    isnothing(truth) || lines!(ax, isochrone_cmd(model, merge(best, truth))...; color = :red, linestyle = :dash, label = "Truth")
    for g in gates
        lines!(ax, [Point2f(v[1], v[2]) for v in vcat(g, g[1:1])]; color = :darkorange, linewidth = 1.5)
    end
    axislegend(ax; position = :rt)
    save(output_file, fig)
    return fig, ax
end

end # module
