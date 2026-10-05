"""Module with code to sample simulated photometric catalogs from model stellar populations. Main function is `simulate_catalog` which takes in a properly formatted YAML configuration file."""
module Simulate

export simulate_catalog

import StarFormationHistories as SFH
import CSV
import YAML
using ArgCheck: @argcheck
using BolometricCorrections: filternames
using DelimitedFiles: writedlm
using OrderedCollections: OrderedDict
using StableRNGs: StableRNG
using StellarTracks: isochrone
using TypedTables: Table
using ..SFHFitting: fit_sfh
using ..SFHFitting.Parsing: parse_config, parse_imf, parse_binaries, parse_tracks, parse_bcs, parse_metallicity, parse_filter_models, parse_to_vector
using ..SFHFitting.ASTs: snr_model, process_ast_file
using ..SFHFitting.Systematics: mag_select

"""
    sfh_mass_fractions(sfh, logAge)
Returns the fraction of the total birth stellar mass formed in each bin `[logAge[i], logAge[i+1])` of the sorted grid
`logAge`, where the oldest bin ends at the lookback time `sfh["T_max"]` (Gyr). `sfh["model"]` is either `constant`, for a
constant star formation rate from `T_max` to the youngest grid age, or `cumulative`, in which case the cumulative fraction
of mass formed `sfh["cum_sfh"]` at each `sfh["logAge"]` is interpolated linearly in lookback time, so the star formation
rate is constant between the given points (and between `T_max`, where the cumulative fraction is 0, and the oldest point).
"""
function sfh_mass_fractions(sfh, logAge)
    @argcheck issorted(logAge) && allunique(logAge)
    T_max = Float64(sfh["T_max"])
    edges = vcat(exp10.(logAge .- 9), T_max) # Bin edges in lookback time [Gyr], youngest first
    if edges[end] <= edges[end-1]
        error("Invalid configuration: sfh.T_max must be older than the oldest stellartrack.logAge.")
    end
    model = lowercase(sfh["model"])
    cumf = if model == "constant"
        t -> (T_max - t) / (T_max - edges[1])
    elseif model == "cumulative"
        t, c = exp10.(Float64.(sfh["logAge"]) .- 9), Float64.(sfh["cum_sfh"])
        if length(t) != length(c) || isempty(t)
            error("Invalid configuration: sfh.logAge and sfh.cum_sfh must be nonempty and have equal lengths.")
        end
        p = sortperm(t; rev=true) # Oldest first
        tn, cn = vcat(T_max, t[p]), vcat(0.0, c[p])
        if !allunique(tn) || !issorted(tn; rev=true)
            error("Invalid configuration: sfh.logAge values must be distinct and younger than sfh.T_max.")
        elseif !issorted(cn) || !isapprox(cn[end], 1)
            error("Invalid configuration: sfh.cum_sfh must increase toward the present and end at 1.")
        elseif tn[end] < edges[1]
            error("Invalid configuration: the youngest sfh.logAge is younger than the youngest stellartrack.logAge, so stars formed after it would fall outside the grid.")
        end
        function (e)
            e <= tn[end] && return cn[end]
            k = findlast(>=(e), tn) # tn[k] >= e > tn[k+1]
            return cn[k] + (cn[k+1] - cn[k]) * (tn[k] - e) / (tn[k] - tn[k+1])
        end
    else
        error("Invalid configuration: sfh.model $(sfh["model"]) invalid; valid options are constant and cumulative.")
    end
    return [cumf(edges[i]) - cumf(edges[i+1]) for i in eachindex(logAge)]
end

# Mock-observes the pure magnitudes in `cols` for `filters`, with one detection draw per star for the group, and stores
# columns <filter>_obs that are NaN for undetected stars
function observe!(cols, filters, model, rng)
    mags = [[cols[Symbol(f)][i] for f in filters] for i in eachindex(cols[:m_ini])]
    obs, good = SFH.model_cmd(mags, model.err, model.completeness, model.bias; rng, ret_idxs=true)
    for (j, f) in enumerate(filters)
        c = fill(NaN, length(mags))
        c[good] .= getindex.(obs, j)
        cols[Symbol(f, "_obs")] = c
    end
    return cols
end

simulate_catalog(file::AbstractString) = simulate_catalog(YAML.load_file(file; dicttype=OrderedDict{String, Any}))

"""
    simulate_catalog(config_file::AbstractString)
    simulate_catalog(config::AbstractDict)
Samples a simulated photometric catalog as described by a YAML configuration file (see `examples/simulate_catalog/config.yml`), or
by the equivalent dictionary, and writes it to `config["output"]["path"]`. Returns a `NamedTuple` with the `catalog` and
`truth` tables and the `fit_sfh` result `fit` (`nothing` unless `fit.run` is true).
"""
function simulate_catalog(config::AbstractDict)
    config = deepcopy(config)
    output_path = get(config["output"], "path", ".")
    sampling = get!(config, "sampling", OrderedDict{String, Any}())
    # The seed is always recorded in the saved input.yml so the catalog can be regenerated
    seed = get!(sampling, "seed", rand(0:typemax(Int)))
    rng = StableRNG(seed)

    # Stellar models and output filters
    st = config["stellartrack"]
    tracklib = parse_tracks(st)
    logAge = sort(collect(Float64, eval(Meta.parse(st["logAge"]))))
    MH_grid = collect(Float64, eval(Meta.parse(st["MH"])))
    bcdicts = config["bolometriccorrections"]
    @info "Loading bolometric corrections"
    bcs = parse_bcs(config)
    bcfilters = map(zip(bcs, values(bcdicts))) do (bc, d)
        available = String[string(f) for f in filternames(bc)]
        haskey(d, "filters") || return available
        f = parse_to_vector(d["filters"])
        missing_filters = setdiff(f, available)
        if !isempty(missing_filters)
            error("Invalid configuration: filters $(join(missing_filters, ", ")) not found in bolometric correction grid $(d["name"]) $(d["filterset"]); available filters are $(join(available, ", ")).")
        end
        f
    end
    mag_names = reduce(vcat, bcfilters)
    if !allunique(mag_names)
        dupes = unique(f for f in mag_names if count(==(f), mag_names) > 1)
        error("Invalid configuration: output filters $(join(dupes, ", ")) appear in more than one bolometric correction grid; use `filters` to select each from one grid.")
    end
    imf = parse_imf(config)
    binary_model = parse_binaries(config)

    # Population properties
    props = config["properties"]
    if haskey(props, "stellar_mass") == haskey(props, "absolute_magnitude")
        error("Invalid configuration: give exactly one of properties.stellar_mass or properties.absolute_magnitude.")
    elseif haskey(props, "distance_modulus") == haskey(props, "distance")
        error("Invalid configuration: give exactly one of properties.distance_modulus or properties.distance.")
    end
    dmod = haskey(props, "distance_modulus") ? Float64(props["distance_modulus"]) : SFH.distance_modulus(Float64(props["distance"]))
    Av = Float64(props["Av"])
    absmag_name = get(props, "absolute_magnitude_filter", nothing)
    if haskey(props, "absolute_magnitude") && absmag_name ∉ mag_names
        error("Invalid configuration: properties.absolute_magnitude_filter $absmag_name not found in the output filters $(join(mag_names, ", ")).")
    end
    T_max = Float64(config["sfh"]["T_max"])
    age_fracs = sfh_mass_fractions(config["sfh"], logAge)
    MH_model, disp_model = parse_metallicity(config; T_max)
    # One simple stellar population (SSP) per (logAge, MH) pair; unique(ssp_logAge) == logAge as calculate_coeffs requires
    ssp_logAge = repeat(logAge; inner=length(MH_grid))
    ssp_MH = repeat(MH_grid; outer=length(logAge))
    ssp_masses(M) = SFH.calculate_coeffs(MH_model, disp_model, M .* age_fracs, ssp_logAge, ssp_MH)

    mag_lim = Float64(get(sampling, "mag_lim", Inf))
    mag_lim_name = get(sampling, "mag_lim_filter", nothing)
    if isfinite(mag_lim) && mag_lim_name ∉ mag_names
        error("Invalid configuration: sampling.mag_lim_filter $mag_lim_name not found in the output filters $(join(mag_names, ", ")).")
    end
    min_ssp_mass = Float64(get(sampling, "min_ssp_mass", 0.01))

    # Observational models
    mock = get(config, "mockobservations", nothing)
    asts = isnothing(mock) ? nothing : get(mock, "ASTs", nothing)
    ast_filters = isnothing(asts) ? String[] : parse_to_vector(asts["filters"])
    if !isnothing(asts) && (length(ast_filters) != 2 || !allunique(ast_filters))
        error("Invalid configuration: mockobservations.ASTs.filters must list exactly two different filters.")
    end
    filter_models = isnothing(mock) ? OrderedDict{String, NamedTuple}() : parse_filter_models(mock)
    for f in Iterators.flatten((ast_filters, keys(filter_models)))
        f ∈ mag_names || error("Invalid configuration: mock-observed filter $f not found in the output filters $(join(mag_names, ", ")).")
    end
    for f in intersect(ast_filters, keys(filter_models))
        error("Invalid configuration: filter $f has both an AST model and a per-filter `filters` model in mockobservations; give only one.")
    end
    models = Any[(filters=[f], model=(completeness=(m.completeness,), err=(m.err,), bias=(m.bias,)))
              for (f, m) in ((f, snr_model(s.mag, s.snr; s.bias, s.minerr, s.snr50, s.width)) for (f, s) in filter_models)]
    if !isnothing(asts)
        model = process_ast_file(asts["ast_file"], ast_filters, get(asts, "badval", 99.999), Float64(get(asts, "minerr", 0.0)), Float64(get(asts, "maxerr", Inf)), false, output_path)
        pushfirst!(models, (filters=ast_filters, model=model))
    end
    obs_filters = vcat(ast_filters, collect(keys(filter_models)))

    # fit_sfh configuration, built before sampling so configuration errors surface early
    fit = get(config, "fit", OrderedDict{String, Any}())
    run_fit = Bool(get(fit, "run", false))
    if run_fit
        binning = fit["binning"]
        fit_filters = unique(vcat(binning["yfilter"], parse_to_vector(binning["xcolor"])))
        for f in fit_filters
            f ∈ obs_filters || error("Invalid configuration: fit filter $f is not mock-observed; add it to mockobservations.")
        end
        if count(in(fit_filters), ast_filters) == 1
            error("Invalid configuration: the AST models are joint functions of both AST filters, so both $(join(ast_filters, " and ")) must be used in fit.binning, or neither.")
        end
        fitdir = joinpath(output_path, "fit")
        data = OrderedDict{String, Any}("path" => fitdir, "photometry" => OrderedDict{String, Any}("photometry_file" => "phot.dat", "filters" => join(fit_filters, ", ")), "binning" => binning)
        if any(in(fit_filters), ast_filters)
            data["ASTs"] = merge(OrderedDict{String, Any}(asts), OrderedDict{String, Any}("ast_file" => abspath(asts["ast_file"]), "badval" => get(asts, "badval", 99.999)))
        end
        fit_models = OrderedDict{String, Any}(f => mock["filters"][f] for f in fit_filters if haskey(filter_models, f))
        isempty(fit_models) || (data["filters"] = fit_models)
        metallicity = OrderedDict{String, Any}(config["metallicity"])
        tracks = OrderedDict{String, Any}("track1" => OrderedDict{String, Any}(k => v for (k, v) in st if k ∉ ("logAge", "MH")), "logAge" => st["logAge"], "MH" => st["MH"])
        # Default to every BC grid containing the fit filters; mag_select throws an ArgumentError when a filter is missing
        fit_bcs = OrderedDict{String, Any}()
        for (bc, (key, d)) in zip(bcs, bcdicts)
            has_filters = try
                mag_select(bc, binning["yfilter"], parse_to_vector(binning["xcolor"])); true
            catch e
                e isa ArgumentError || rethrow()
                false
            end
            has_filters && (fit_bcs[key] = OrderedDict{String, Any}(k => v for (k, v) in d if k != "filters"))
        end
        fitcfg = OrderedDict{String, Any}("data" => data, "imf" => config["imf"], "binaries" => config["binaries"],
                                          "stellartracks" => tracks, "bolometriccorrections" => fit_bcs,
                                          "properties" => OrderedDict{String, Any}("Av" => Av, "distance_modulus" => dmod, "T_max" => T_max),
                                          "metallicity" => metallicity, "plotting" => OrderedDict{String, Any}("diagnostics" => true),
                                          "output" => OrderedDict{String, Any}("path" => fitdir, "filename" => "results" * splitext(config["output"]["filename"])[2]))
        # Any fit_sfh section given under `fit` replaces the generated default
        for (k, v) in fit
            k ∈ ("run", "stellar_mass_error", "binning") || (fitcfg[k] = deepcopy(v))
        end
        get!(fitcfg["metallicity"], "T_max", T_max)
        if isempty(fitcfg["bolometriccorrections"])
            error("Invalid configuration: no bolometric correction grid contains all fit filters; give fit.bolometriccorrections.")
        end
        # A finite mag_lim removes stars that the fit_sfh templates expect, unless they are essentially undetectable
        if isfinite(mag_lim)
            i = findfirst(m -> mag_lim_name ∈ m.filters, models)
            if !isnothing(i)
                m = models[i]
                cmax = if length(m.filters) == 1
                    only(m.model.completeness)(mag_lim)
                else
                    maximum(x -> mag_lim_name == m.filters[1] ? m.model.completeness(mag_lim, x) : m.model.completeness(x, mag_lim), range(-10, 50; step=0.01))
                end
                if cmax > 0.01
                    @warn "sampling.mag_lim = $mag_lim in $mag_lim_name removes stars that are detected with probability up to $(round(cmax; sigdigits=3)); the fit_sfh templates include these stars, so the fit will be biased. Use a fainter mag_lim."
                end
            end
        end
    end

    mkpath(output_path)
    YAML.write_file(joinpath(output_path, "input.yml"), config)

    # Isochrones of the SSPs to sample; SSPs below min_ssp_mass are skipped because sampling draws at least one star per SSP
    # Names local to these closures must not match locals of simulate_catalog, which closures would capture and share
    function ssp_isochrone(k)
        grids = [isochrone(tracklib, bc, ssp_logAge[k], ssp_MH[k], Av) for bc in bcs]
        m_ini = collect(Float64, first(grids).m_ini)
        all(g -> g.m_ini == m_ini, grids) || error("Isochrones from different bolometric correction grids have different initial masses.")
        return m_ini, [collect(Float64, getproperty(g, Symbol(f))) for (g, fs) in zip(grids, bcfilters) for f in fs]
    end
    function ssp_isochrones(idxs)
        out = Vector{Tuple{Vector{Float64}, Vector{Vector{Float64}}}}(undef, length(idxs))
        Threads.@threads for i in eachindex(idxs)
            out[i] = ssp_isochrone(idxs[i])
        end
        return out
    end
    @info "Computing isochrones"
    if haskey(props, "stellar_mass")
        M = Float64(props["stellar_mass"])
        masses = ssp_masses(M)
        kept = findall(>=(min_ssp_mass), masses)
        isos = ssp_isochrones(kept)
    else
        # Total birth stellar mass whose expected integrated luminosity matches absolute_magnitude; the PowerLawMZR SSP
        # masses depend on the total mass, so iterate to a fixed point (the AMRs converge in one step)
        candidates = findall(>(0), repeat(age_fracs; inner=length(MH_grid)))
        isos = ssp_isochrones(candidates)
        j = findfirst(==(absmag_name), mag_names)
        lpm = [SFH.luminosity_per_mass(m_ini, mags, imf)[j] for (m_ini, mags) in isos]
        L = SFH.mag2flux(Float64(props["absolute_magnitude"]))
        M = 1e6
        for iter in 1:100
            M_new = M * L / sum(ssp_masses(M)[candidates] .* lpm)
            converged = abs(M_new / M - 1) < 1e-10
            M = M_new
            converged && break
            iter == 100 && error("Total stellar mass for properties.absolute_magnitude did not converge.")
        end
        masses = ssp_masses(M)
        keep = masses[candidates] .>= min_ssp_mass
        kept, isos = candidates[keep], isos[keep]
    end

    @info "Sampling stars"
    kws = isfinite(mag_lim) ? (; mag_lim, mag_lim_name) : (;)
    massvec, magvec = if haskey(props, "stellar_mass")
        SFH.generate_stars_mass_composite(first.(isos), last.(isos), mag_names, sum(masses[kept]), masses[kept], imf; binary_model, rng, dist_mod=dmod, kws...)
    else
        SFH.generate_stars_mag_composite(first.(isos), last.(isos), mag_names, Float64(props["absolute_magnitude"]), absmag_name, masses[kept], imf; frac_type=:mass, binary_model, rng, dist_mod=dmod, kws...)
    end
    systems = reduce(vcat, massvec)
    allmags = reduce(vcat, magvec)
    nstars = length.(massvec)
    cols = OrderedDict{Symbol, Vector{Float64}}(
        :m_ini => first.(systems),
        :m_ini2 => [length(s) > 1 ? s[2] : 0.0 for s in systems], # 0 for single stars
        :logAge => reduce(vcat, fill.(ssp_logAge[kept], nstars)),
        :MH => reduce(vcat, fill.(ssp_MH[kept], nstars)))
    for (j, f) in enumerate(mag_names)
        cols[Symbol(f)] = getindex.(allmags, j)
    end
    for m in models
        observe!(cols, m.filters, m.model, rng)
    end

    @info "Writing catalog"
    catalog = Table(NamedTuple{Tuple(keys(cols))}(Tuple(values(cols))))
    CSV.write(joinpath(output_path, config["output"]["filename"]), catalog; delim=' ')
    _, cum_sfh, sfr, mean_MH = SFH.calculate_cum_sfr(masses, ssp_logAge, ssp_MH, T_max; sorted=true)
    truth = Table(logAge_lower=logAge, logAge_upper=vcat(logAge[begin+1:end], log10(T_max) + 9), sfr=sfr, cum_sfh=cum_sfh, MH=mean_MH)
    base, ext = splitext(config["output"]["filename"])
    CSV.write(joinpath(output_path, base * "_truth" * ext), truth; delim=' ')

    result = nothing
    if run_fit
        @info "Fitting the simulated catalog"
        mkpath(fitdir)
        phot = reduce(hcat, cols[Symbol(f, "_obs")] for f in fit_filters)
        writedlm(joinpath(fitdir, "phot.dat"), phot[vec(all(!isnan, phot; dims=2)), :])
        fitcfg["properties"]["Mstar"] = M * (1 + Float64(get(fit, "stellar_mass_error", 0.0)) * randn(rng))
        result = fit_sfh(parse_config(fitcfg))
    end
    return (catalog=catalog, truth=truth, fit=result)
end

end # module
