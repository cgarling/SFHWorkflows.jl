"""Module with code to measure resolved star formation histories from provided photometry. Main function is `fit_sfh` which takes in a properly formatted YAML configuration file."""
module SFHFitting

export fit_sfh

using BolometricCorrections: gridname
import StarFormationHistories as SFH
using ArgCheck: @argcheck, @check
using DelimitedFiles: readdlm
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

# Top-level functions
function fit_sfh(obsfile::AbstractString, astfile::AbstractString, filters, xstrings, ystring, edges,
                 MH_model0::SFH.AbstractMetallicityModel, disp_model0::SFH.AbstractDispersionModel, Mstar::Number, stellar_tracks, bcs,
                 dmod::Number, Av::Number, imf, MH, logAge, binary_model::SFH.AbstractBinaryModel, output_filename::AbstractString; 
                 badval::Number=99.999, minerr::Number=0.0, maxerr::Number=0.2, plot_diagnostics::Bool=true, output_path::AbstractString=".",
                 dtype::Type{<:AbstractFloat}=Float64) # filters=("mag1", "mag2")
    @argcheck length(xstrings) == 2
    @argcheck Mstar > 0
    filters = string.(filters)
    completeness, bias, err = process_ast_file(astfile, filters, badval, minerr, maxerr, plot_diagnostics, output_path; dtype=dtype)
    data = readdlm(obsfile, ' ', dtype)
    yidx = findfirst(==(ystring), filters)
    xidxs = [findfirst(==(x), filters) for x in xstrings]
    h = SFH.bin_cmd(view(data, :, xidxs[1]) .- view(data, :, xidxs[2]), view(data, :, yidx); edges=edges)
    out_file = joinpath(output_path, output_filename)
    result = systematics(MH_model0, disp_model0, Mstar, vec(h.weights), stellar_tracks, bcs, xstrings, ystring, dmod, Av, err, completeness, bias, imf, MH, logAge, edges; binary_model=binary_model, output=out_file, dtype=dtype)
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

fit_sfh(@nospecialize(config::NamedTuple)) = fit_sfh(config.phot_file, config.ast_file, config.filters, config.xstrings, config.ystring, (config.xbins, config.ybins), config.MH_model0, config.disp_model0, config.Mstar, config.stellar_tracks, config.bcs, config.dmod, config.Av, config.imf, config.MH, config.logAge, config.binary_model, config.output_filename; badval=config.badval, minerr=config.minerr, maxerr=config.maxerr, plot_diagnostics=config.plot_diagnostics, output_path=config.output_path, dtype=config.dtype)
fit_sfh(config_file::AbstractString) = fit_sfh(parse_config(config_file))


end # module