"""Module with code to measure resolved star formation histories from provided photometry. Main function is `fit_sfh` which takes in a properly formatted YAML configuration file."""
module SFHFitting

export fit_sfh

using BolometricCorrections: gridname
import StarFormationHistories as SFH
using ArgCheck: @argcheck, @check
using DelimitedFiles: readdlm
using PDFmerger: merge_pdfs
using StatsBase: Histogram

# Includes
include("Parsing.jl") # Code to parse configuration file
using .Parsing
include("ASTs.jl") # Code to analyze artificial star tests
using .ASTs
include("Systematics.jl")
using .Systematics
include("Plotting.jl")
using .Plotting

function write_histogram(hist::Histogram, filename::String)
    @argcheck length(hist.edges) == 2
    @argcheck ndims(hist.weights) == 2
    open(filename, "w") do io
        # Write bin edges: each dimension on its own line
        println(io, join(collect(hist.edges[1]), " "))
        println(io, join(collect(hist.edges[2]), " "))

        # Write weights: each row of the 2D array on its own line
        for row in eachrow(hist.weights)
            println(io, join(row, " "))
        end
    end
end

function read_histogram(filename::String)
    lines = readlines(filename)

    # Parse bin edges from first two lines
    edges1 = parse.(Float64, split(lines[1]))
    edges2 = parse.(Float64, split(lines[2]))

    # Parse weights from remaining lines
    weights_data = lines[3:end]
    weights = [parse.(Float64, split(line)) for line in weights_data]
    weights_array = reduce(vcat, [reshape(row, 1, :) for row in weights])

    # Construct histogram
    return Histogram((edges1, edges2), weights_array)
end

# Hess diagram bins excluded by `gates`, with a warning for each gate that contains no bin centers
function gate_mask(edges, gates)
    for (i, g) in enumerate(gates)
        any(SFH.hess_mask(edges, [g])) || @warn "Gate $i contains no Hess diagram bin centers; its vertices must be (color, magnitude) pairs inside the binning range."
    end
    mask = SFH.hess_mask(edges, gates)
    all(mask) && error("The gates cover every Hess diagram bin, leaving nothing to fit.")
    return mask
end

"""
    (completeness, bias, err) = observational_models(xstrings, ystring, astfile, ast_filters, filter_models;
                                                     badval, minerr, maxerr, plot_diagnostics, output_path)
Returns the completeness, bias, and error models of the magnitudes used for the Hess diagram with x-axis color `xstrings`
and y-axis filter `ystring`, in the order of the isochrone magnitudes from `Systematics.iso_filters`. The two filters
`ast_filters` are modeled jointly from the artificial star tests in `astfile`, and every other filter by its per-filter
model in `filter_models`; `astfile` may be `nothing` when every filter has a per-filter model.
"""
function observational_models(xstrings, ystring, astfile, ast_filters, filter_models; badval, minerr, maxerr, plot_diagnostics, output_path)
    ast_filters = isnothing(astfile) ? String[] : string.(ast_filters)
    # Isochrone magnitudes are ordered this way by Systematics.iso_filters; models must match
    mag_order = ystring in xstrings ? xstrings : [ystring; xstrings]
    Parsing.check_filter_models(mag_order, ast_filters, keys(filter_models))
    models = Dict(name => snr_model(s.mag, s.snr; s.bias, s.minerr, s.snr50, s.width) for (name, s) in filter_models)
    isnothing(astfile) && return Tuple([getproperty(models[f], k) for f in mag_order] for k in (:completeness, :bias, :err))
    @argcheck length(ast_filters) == 2
    asts = process_ast_file(astfile, ast_filters, badval, minerr, maxerr, plot_diagnostics, output_path)
    # Joint models of the magnitudes in `mag_order`: the AST models of the two AST filters, combined with the
    # per-filter models of any other filter
    ia, ib = (findfirst(==(f), mag_order) for f in ast_filters)
    other = [(k, models[mag_order[k]]) for k in eachindex(mag_order) if k ∉ (ia, ib)]
    completeness = (m...) -> asts.completeness(m[ia], m[ib]) * prod(o.completeness(m[k]) for (k, o) in other; init=1.0)
    joint(f, g) = (m...) -> (a = f(m[ia], m[ib]); ntuple(k -> k == ia ? a[1] : k == ib ? a[2] : getproperty(models[mag_order[k]], g)(m[k]), length(m)))
    return completeness, joint(asts.bias, :bias), joint(asts.err, :err)
end

# Merges the artificial star test diagnostics written by `observational_models` and the completeness over the Hess
# diagram into `output_path/diagnostics.pdf`
function write_diagnostics(completeness, astfile, edges, xstrings, ystring, output_path)
    pdfs = isnothing(astfile) ? String[] : ["residuals1.pdf", "residuals2.pdf", "error.pdf", "completeness.pdf"]
    # With a third filter, a Hess diagram bin does not determine every magnitude the completeness depends on
    if ystring in xstrings
        plot_completeness_hess(completeness, edges, xstrings, ystring, joinpath(output_path, "completeness_hess.pdf"))
        push!(pdfs, "completeness_hess.pdf")
    end
    isempty(pdfs) || merge_pdfs(joinpath.(output_path, pdfs), joinpath(output_path, "diagnostics.pdf"); cleanup=true)
end

# Observed Hess diagram from the photometry file, whose columns are the magnitudes in `filters`
function obs_hess(obsfile, filters, xstrings, ystring, edges)
    data = readdlm(obsfile, Float64)
    yidx = findfirst(==(ystring), filters)
    xidxs = [findfirst(==(x), filters) for x in xstrings]
    return SFH.bin_cmd(view(data, :, xidxs[1]) .- view(data, :, xidxs[2]), view(data, :, yidx); edges=edges)
end

# Top-level functions
# `astfile` may be `nothing` when every filter used in the fit has a model in `filter_models`
function fit_sfh(obsfile::AbstractString, astfile::Union{AbstractString, Nothing}, filters, xstrings, ystring, edges,
                 MH_model0::SFH.AbstractMetallicityModel, disp_model0::SFH.AbstractDispersionModel, Mstar::Number, stellar_tracks, bcs,
                 dmod::Number, Av::Number, imf, MH, logAge, binary_model::SFH.AbstractBinaryModel, output_filename::AbstractString;
                 badval::Number=99.999, minerr::Number=0.0, maxerr::Number=0.2, plot_diagnostics::Bool=true, output_path::AbstractString=".",
                 ast_filters=first(string.(filters), 2), filter_models=Dict{String, NamedTuple}(), T_max::Number=13.7,
                 gates=Vector{NTuple{2, Float64}}[]) # filters=("mag1", "mag2")
    @argcheck length(xstrings) == 2
    @argcheck Mstar > 0
    mask = gate_mask(edges, gates)
    completeness, bias, err = observational_models(xstrings, ystring, astfile, ast_filters, filter_models; badval, minerr, maxerr,
                                                   plot_diagnostics, output_path)
    plot_diagnostics && write_diagnostics(completeness, astfile, edges, xstrings, ystring, output_path)
    h = obs_hess(obsfile, string.(filters), xstrings, ystring, edges)
    out_file = joinpath(output_path, output_filename)
    result = systematics(MH_model0, disp_model0, Mstar, vec(h.weights), stellar_tracks, bcs, xstrings, ystring, dmod, Av, err, completeness, bias, imf, MH, logAge, edges; binary_model=binary_model, output=out_file, T_max, mask)
    # Write histograms to files
    ext = splitext(output_filename)[2]
    write_histogram(h, joinpath(output_path, splitext(output_filename)[1]*"_obshess"*ext))
    for i in eachindex(stellar_tracks)
        for j in eachindex(bcs)
            # Construct best-fit model histogram
            coeffs = SFH.calculate_coeffs(result.results[i,j], result.logAge[i,j], result.MH[i,j])
            model_hess = sum(coeffs .* result.templates[i,j]) #  ./ normalize_value)
            fname = joinpath(output_path, splitext(output_filename)[1]*"_modelhess_"*gridname(stellar_tracks[i])*"_"*gridname(bcs[j])*ext)
            write_histogram(Histogram(edges, model_hess), fname)
        end
    end
    return result, h
end

fit_sfh(@nospecialize(config::NamedTuple)) = fit_sfh(config.phot_file, config.ast_file, config.filters, config.xstrings, config.ystring, (config.xbins, config.ybins), config.MH_model0, config.disp_model0, config.Mstar, config.stellar_tracks, config.bcs, config.dmod, config.Av, config.imf, config.MH, config.logAge, config.binary_model, config.output_filename; badval=config.badval, minerr=config.minerr, maxerr=config.maxerr, plot_diagnostics=config.plot_diagnostics, output_path=config.output_path, config.ast_filters, config.filter_models, config.T_max, config.gates)
fit_sfh(config_file::AbstractString) = fit_sfh(parse_config(config_file))


end # module